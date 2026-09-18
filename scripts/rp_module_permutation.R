library(tidyverse)
library(here)
source(here("scripts", "heatmap_helpers.R"))

# Was the ribosomal protein module simply the most visible slice of a global
# effect? RP genes are the most abundant class in the transcriptome, so any
# broad change would surface there first. Two controls:
#
#   1. compare RP genes with non-RP genes drawn from the same expression bins
#   2. compare the module with random gene sets matched for size AND abundance
#
# The module was chosen after inspecting genome-wide GSEA output, so these test
# whether that choice is doing the work. They do not remove the selection bias
# in the effect sizes themselves, which remain upper estimates.

set.seed(1)
n_permutations <- 2000

counts <- read_corr(here("output", "pseudobulk_counts_by_phase.csv"))
meta <- tibble(sample_id = colnames(counts)) |>
  separate_wider_delim(
    sample_id,
    delim = "|", names = c("cell_line", "guide", "phase"), cols_remove = FALSE
  ) |>
  filter(phase == "G1")

cpm <- prop.table(counts[, meta[["sample_id"]]], margin = 2) * 1e6
colnames(cpm) <- str_c(meta[["cell_line"]], meta[["guide"]], sep = "|")

mean_logfc_per_line <- function(line) {
  samples <- colnames(cpm) |> str_subset(str_glue("^{line}\\|"))
  controls <- samples |> str_subset(regex("mruby", ignore_case = TRUE))
  targeting <- samples |>
    discard(\(s) str_detect(str_to_lower(s), "mruby")) |>
    discard(is_ineffective)

  baseline_all <- rowMeans(cpm[, controls])
  expressed <- names(baseline_all)[baseline_all >= 10]

  logfc <- targeting |>
    map(\(sample) {
      reference <- setdiff(controls, sample)
      log2(cpm[expressed, sample] + 1) -
        log2(rowMeans(cpm[expressed, reference, drop = FALSE]) + 1)
    }) |>
    reduce(`+`) / length(targeting)

  list(logfc = logfc, abundance = baseline_all[expressed])
}

results <- c("A549", "H1299") |>
  set_names() |>
  map(\(line) {
    fit <- mean_logfc_per_line(line)
    logfc <- fit[["logfc"]]
    abundance <- fit[["abundance"]]
    rp_genes <- names(logfc) |> str_subset("^RP[LS]\\d")

    # 20 abundance strata; each permutation draws a non-RP gene from the same
    # stratum as each RP gene, so size and expression profile are matched
    bin <- cut(
      log10(abundance),
      breaks = quantile(log10(abundance), seq(0, 1, 0.05)),
      include.lowest = TRUE, labels = FALSE
    )
    names(bin) <- names(abundance)
    pool <- split(setdiff(names(logfc), rp_genes), bin[setdiff(names(logfc), rp_genes)])
    rp_bins <- bin[rp_genes]

    null <- map_dbl(seq_len(n_permutations), \(i) {
      drawn <- map_chr(rp_bins, \(b) sample(pool[[as.character(b)]], 1))
      mean(logfc[drawn])
    })
    observed <- mean(logfc[rp_genes])

    tibble(
      cell_line = line,
      n_rp_genes = length(rp_genes),
      observed = observed,
      null_mean = mean(null),
      null_sd = sd(null),
      z = (observed - mean(null)) / sd(null),
      p_empirical = (sum(null <= observed) + 1) / (n_permutations + 1),
      sets_beaten = sum(null > observed)
    )
  }) |>
  list_rbind()

cat("=== ribosomal module vs abundance-matched random sets ===\n")
print(results |> mutate(across(where(is.numeric), \(x) signif(x, 3))))
write_csv(results, here("output", "rp_module_permutation.csv"))
