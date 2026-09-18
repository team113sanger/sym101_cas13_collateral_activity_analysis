library(tidyverse)
library(here)
library(edgeR)
library(limma)
library(fgsea)
library(msigdbr)

# Differential expression and GSEA for the shared Cas13d collateral response.
# Analysis only: writes the DE and GSEA tables that the figure scripts read.
# Guide-level pseudobulks are the replicates: 8 targeting vs 3 non-targeting
# (mRuby) per cell line, tested within line and as a cross-line average.

set.seed(1)

fdr_cutoff <- 0.05
min_cells <- 50

counts <- read_csv(here("output", "pseudobulk_counts.csv"), show_col_types = FALSE) |>
  column_to_rownames(var = "...1") |>
  as.matrix()

# CD9_3 achieves no knockdown of a well-measured target, so it is dropped from
# the targeting group rather than diluting it
ineffective_guides <- c("CD9_3")

metadata <- read_csv(here("output", "pseudobulk_metadata.csv"), show_col_types = FALSE) |>
  filter(n_cells >= min_cells, !str_to_upper(guide) %in% str_to_upper(ineffective_guides)) |>
  mutate(
    condition = if_else(
      str_starts(str_to_lower(guide), "mruby"), "non_targeting", "targeting"
    ),
    target = guide |>
      str_to_upper() |>
      str_remove("_\\d+$") |>
      recode(CD95 = "FAS", MRUBY = "non-targeting"),
    group = str_c(cell_line, condition, sep = "_")
  ) |>
  # the lines sit at very different baseline phase mixes (A549 ~91% G1, H1299
  # ~33%), so centre within line: the covariate then carries only the
  # guide-to-guide shift, not the between-line difference
  mutate(
    across(c(frac_S, frac_G2M), \(x) x - mean(x), .names = "{.col}_centred"),
    .by = cell_line
  ) |>
  arrange(cell_line, condition, guide)

counts <- counts[, metadata[["sample_id"]]]

dge <- DGEList(counts = counts, samples = metadata, group = metadata[["group"]])
keep <- filterByExpr(dge, group = dge[["samples"]][["group"]])
dge <- dge[keep, , keep.lib.sizes = FALSE] |> calcNormFactors()

# targeting guides shift the cell cycle mix (notably in H1299), so S and G2M
# fractions are covariates: without them a composition shift reads as an
# E2F/G2M expression signature
design <- model.matrix(~ 0 + group + frac_S_centred + frac_G2M_centred, data = metadata)
colnames(design) <- str_remove(colnames(design), "^group")

# pseudobulks differ several-fold in cell number, so weight samples by quality
fit <- voomWithQualityWeights(dge, design, plot = FALSE) |>
  lmFit(design)

# the two lines respond in opposite directions, so "shared" is reported only
# alongside the interaction that tests whether the responses actually differ
contrasts <- makeContrasts(
  A549 = A549_targeting - A549_non_targeting,
  H1299 = H1299_targeting - H1299_non_targeting,
  shared = (A549_targeting - A549_non_targeting +
    H1299_targeting - H1299_non_targeting) / 2,
  interaction = (H1299_targeting - H1299_non_targeting) -
    (A549_targeting - A549_non_targeting),
  levels = design
)

fit <- contrasts.fit(fit, contrasts) |> eBayes()

results <- colnames(contrasts) |>
  set_names() |>
  map(\(contrast) {
    topTable(fit, coef = contrast, number = Inf, sort.by = "P") |>
      rownames_to_column(var = "gene") |>
      as_tibble() |>
      mutate(contrast = contrast)
  })

results |>
  list_rbind() |>
  write_csv(here("output", "de_targeting_vs_control.csv"))

walk2(results, names(results), \(res, contrast) {
  n_sig <- sum(res[["adj.P.Val"]] < fdr_cutoff)
  cat(str_glue("{contrast}: {n_sig} genes at FDR < {fdr_cutoff}"), "\n")
})

# Positive control: the targeting group pools guides against three different
# genes, so on-target knockdown is diluted there. Test each guide family against
# the controls separately, where the on-target gene must move.
on_target <- tibble(family = c("CD9", "CD95", "CD274"), gene = c("CD9", "FAS", "CD274")) |>
  mutate(
    result = map(family, \(fam) {
      samples <- metadata |>
        filter(str_starts(str_to_upper(guide), fam) | condition == "non_targeting")
      fam_dge <- dge[, samples[["sample_id"]]]
      fam_design <- model.matrix(~ cell_line + condition, data = samples)
      fam_fit <- voomWithQualityWeights(fam_dge, fam_design, plot = FALSE) |>
        lmFit(fam_design) |>
        eBayes()
      topTable(fam_fit, coef = "conditiontargeting", number = Inf) |>
        rownames_to_column(var = "gene") |>
        as_tibble()
    })
  ) |>
  mutate(
    on_target_stats = map2(result, gene, \(res, g) filter(res, gene == g) |> select(logFC, adj.P.Val)),
    n_sig = map_int(result, \(res) sum(res[["adj.P.Val"]] < fdr_cutoff))
  ) |>
  unnest(on_target_stats, keep_empty = TRUE) |>
  select(family, gene, logFC, adj.P.Val, n_sig)

cat("\n=== on-target knockdown per guide family ===\n")
print(on_target)
write_csv(on_target, here("output", "de_on_target_control.csv"))

# Sensitivity analysis: repeat the per-line contrasts on G1 cells only, so the
# phase shift cannot contribute at all
g1_counts <- read_csv(here("output", "pseudobulk_counts_g1.csv"), show_col_types = FALSE) |>
  column_to_rownames(var = "...1") |>
  as.matrix()
g1_metadata <- metadata |> filter(sample_id %in% colnames(g1_counts))
g1_counts <- g1_counts[, g1_metadata[["sample_id"]]]

g1_dge <- DGEList(counts = g1_counts, samples = g1_metadata)
g1_keep <- filterByExpr(g1_dge, group = g1_metadata[["group"]])
g1_dge <- g1_dge[g1_keep, , keep.lib.sizes = FALSE] |> calcNormFactors()

g1_design <- model.matrix(~ 0 + group, data = g1_metadata)
colnames(g1_design) <- str_remove(colnames(g1_design), "^group")
g1_contrasts <- makeContrasts(
  A549 = A549_targeting - A549_non_targeting,
  H1299 = H1299_targeting - H1299_non_targeting,
  levels = g1_design
)
g1_fit <- voomWithQualityWeights(g1_dge, g1_design, plot = FALSE) |>
  lmFit(g1_design) |>
  contrasts.fit(g1_contrasts) |>
  eBayes()

g1_results <- colnames(g1_contrasts) |>
  set_names() |>
  map(\(contrast) {
    topTable(g1_fit, coef = contrast, number = Inf) |>
      rownames_to_column(var = "gene") |>
      as_tibble() |>
      mutate(contrast = contrast)
  })

g1_results |>
  list_rbind() |>
  write_csv(here("output", "de_targeting_vs_control_g1only.csv"))

cat("\n=== G1-only sensitivity: concordance with all-cell logFC ===\n")
c("A549", "H1299") |>
  walk(\(contrast) {
    joined <- inner_join(
      results[[contrast]] |> select(gene, logFC_all = logFC),
      g1_results[[contrast]] |> select(gene, logFC_g1 = logFC),
      by = "gene"
    )
    n_sig_g1 <- sum(g1_results[[contrast]][["adj.P.Val"]] < fdr_cutoff)
    cat(str_glue(
      "{contrast}: r = {round(cor(joined[['logFC_all']], joined[['logFC_g1']]), 3)}, ",
      "{n_sig_g1} genes at FDR < {fdr_cutoff} in G1 only"
    ), "\n")
  })

# GSEA on the moderated t statistic, so ranking uses effect size and its variance
gene_sets <- list(
  hallmark = msigdbr(species = "Homo sapiens", collection = "H"),
  reactome = msigdbr(
    species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"
  ),
  go_bp = msigdbr(
    species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"
  )
) |>
  map(\(sets) split(sets[["gene_symbol"]], sets[["gs_name"]]))

run_fgsea <- function(res, pathways) {
  ranks <- res |>
    filter(!is.na(t)) |>
    arrange(desc(t)) |>
    pull(t, name = gene)

  fgsea(pathways, ranks, minSize = 15, maxSize = 500, nPermSimple = 10000) |>
    as_tibble() |>
    arrange(padj)
}

# H1299's phase shift is collinear with condition, so its all-cell contrast is
# unreliable; GSEA is run on the G1-restricted fit as well and both are reported
all_fits <- c(
  results |> set_names(str_c(names(results), "_allcells")),
  g1_results |> set_names(str_c(names(g1_results), "_g1only"))
)

gsea <- expand_grid(
  contrast = names(all_fits),
  collection = names(gene_sets)
) |>
  mutate(
    result = map2(contrast, collection, \(x, y) run_fgsea(all_fits[[x]], gene_sets[[y]]))
  ) |>
  unnest(result) |>
  mutate(leadingEdge = map_chr(leadingEdge, \(x) str_c(x, collapse = ";")))

write_csv(gsea, here("output", "gsea_targeting_vs_control.csv"))

gsea |>
  filter(padj < fdr_cutoff) |>
  count(contrast, collection, direction = if_else(NES > 0, "up", "down")) |>
  print(n = Inf)
