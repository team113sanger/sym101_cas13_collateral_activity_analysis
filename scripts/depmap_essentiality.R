library(tidyverse)
library(here)

# Are the targeted genes ones the cells depend on? If CD9 loss carried a fitness
# cost, the shared transcriptional change could be a stress response to losing it
# rather than collateral activity.
#
# Computed from DepMap 24Q4 Avana raw readcounts: log2 fold change of each gene's
# guides at 21 days against the matched pDNA pool, for the two lines used here.
# Scores are placed on an essential-to-non-essential scale so they can be read
# without reference to DepMap's own normalisation.

depmap <- "/home/ubuntu/projects/avana_depmap/DEPMAP24Q4_BAGEL_analysis/data"

screens <- list(
  A549 = c(
    "A549-311Cas9 Rep A p5_Avana-3", "A549-311Cas9 Rep B p5_Avana-3",
    "A549-311Cas9 Rep C p5_Avana-3", "A549-311Cas9 Rep D p5_Avana-3"
  ),
  H1299 = c("NCI-H1299-311Cas9_RepA_p5_Avana-3", "NCI-H1299-311Cas9_RepB_p5_Avana-3")
)
reference <- "pDNA_batch_Avana-3"
wanted <- c(unlist(screens, use.names = FALSE), reference)

counts <- read_csv(
  file.path(depmap, "AvanaRawReadcounts.csv"),
  col_select = all_of(c("...1", wanted)), show_col_types = FALSE
) |>
  rename(sgRNA = `...1`)

guide_map <- read_csv(
  file.path(depmap, "AvanaGuideMap.csv"),
  col_select = c("sgRNA", "Gene", "UsedByChronos"), show_col_types = FALSE
) |>
  filter(UsedByChronos) |>
  distinct(sgRNA, .keep_all = TRUE) |>
  mutate(symbol = str_remove(Gene, "\\s*\\(\\d+\\)$"))

joined <- inner_join(counts, guide_map |> select(sgRNA, symbol), by = "sgRNA")

# a small number of guides carry missing counts; they are dropped outright so
# that library totals are computed on a consistent guide set
mat <- joined |> select(all_of(wanted)) |> as.matrix()
complete <- stats::complete.cases(mat)
cat(str_glue("guides dropped for missing counts: {sum(!complete)}"), "\n")
mat <- mat[complete, , drop = FALSE]
joined <- joined[complete, ]

# CPM-normalise each library, then log2 ratio against the pDNA pool
cpm <- sweep(mat, 2, colSums(mat), "/") * 1e6
logfc <- log2(cpm + 1) - log2(cpm[, reference] + 1)
logfc <- logfc[, setdiff(wanted, reference), drop = FALSE]

guide_level <- imap(screens, \(cols, line) rowMeans(logfc[, cols, drop = FALSE]))

gene_level <- tibble(symbol = joined[["symbol"]]) |>
  bind_cols(as_tibble(guide_level)) |>
  summarise(across(c(A549, H1299), mean), .by = symbol)

# reference sets anchor the scale: 0 = non-essential, 1 = core essential
core_essential <- c(
  "RPL3", "RPL5", "RPS3", "RPS6", "POLR2A", "PSMA1", "PSMB2", "EIF3B", "SF3B1", "RAN"
)
non_essential <- c(
  "OR2T1", "OR5A1", "GPR15", "CD19", "MS4A1", "KRT4", "ACTL9", "DEFB1", "OR1A1", "TAS2R16"
)

anchors <- c("A549", "H1299") |>
  set_names() |>
  map(\(line) {
    list(
      essential = mean(gene_level[[line]][gene_level[["symbol"]] %in% core_essential]),
      non_essential = mean(gene_level[[line]][gene_level[["symbol"]] %in% non_essential])
    )
  })

targets <- gene_level |>
  filter(symbol %in% c("CD9", "CD274", "FAS")) |>
  pivot_longer(c(A549, H1299), names_to = "cell_line", values_to = "log2fc") |>
  mutate(
    scaled_essentiality = map2_dbl(cell_line, log2fc, \(line, value) {
      (value - anchors[[line]][["non_essential"]]) /
        (anchors[[line]][["essential"]] - anchors[[line]][["non_essential"]])
    }),
    percentile = map2_dbl(cell_line, log2fc, \(line, value) 100 * mean(gene_level[[line]] < value))
  )

cat("=== reference anchors (mean log2 FC vs pDNA, 21 days) ===\n")
print(
  imap(anchors, \(a, line) tibble(
    cell_line = line, core_essential = a[["essential"]], non_essential = a[["non_essential"]]
  )) |>
    list_rbind() |>
    mutate(across(where(is.numeric), \(x) round(x, 2)))
)

cat("\n=== targeted genes ===\n")
print(targets |> mutate(across(where(is.numeric), \(x) round(x, 3))))
write_csv(targets, here("output", "depmap_target_essentiality.csv"))
write_csv(gene_level, here("output", "depmap_lfc_a549_h1299.csv"))
