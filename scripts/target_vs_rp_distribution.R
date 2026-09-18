library(tidyverse)
library(here)
library(patchwork)
source(here("scripts", "heatmap_helpers.R"))

# Each guide's own target, placed on the same log2 FC axis as the ribosomal
# protein genes it depletes. The target is a single large change; the ribosomal
# genes are a tight cloud of small ones. G1 cells, vs mRuby controls.

theme_set(theme_bw(base_size = 9))

counts_by_phase <- read_corr(here("output", "pseudobulk_counts_by_phase.csv"))
phase_meta <- tibble(sample_id = colnames(counts_by_phase)) |>
  separate_wider_delim(
    sample_id,
    delim = "|", names = c("cell_line", "guide", "phase"), cols_remove = FALSE
  ) |>
  filter(phase == "G1")

cpm <- prop.table(counts_by_phase[, phase_meta[["sample_id"]]], margin = 2) * 1e6
colnames(cpm) <- str_c(phase_meta[["cell_line"]], phase_meta[["guide"]], sep = "|")

rp_genes <- rownames(cpm) |>
  str_subset("^RP[LS]\\d") |>
  keep(\(g) mean(cpm[g, ]) >= 10)

target_cpm_in_line <- function(gene, line) {
  columns <- colnames(cpm) |> str_subset(str_glue("^{line}\\|"))
  if (gene %in% rownames(cpm)) mean(cpm[gene, columns]) else 0
}

changes <- phase_meta |>
  distinct(cell_line) |>
  mutate(
    rows = map(cell_line, \(line) {
      samples <- colnames(cpm) |> str_subset(str_glue("^{line}\\|"))
      controls <- samples |> str_subset(regex("mruby", ignore_case = TRUE))

      samples |>
        discard(is_ineffective) |>
        discard(\(s) str_detect(str_to_lower(s), "mruby")) |>
        map(\(sample) {
          reference <- setdiff(controls, sample)
          log_change <- log2(cpm[, sample] + 1) -
            log2(rowMeans(cpm[, reference, drop = FALSE]) + 1)

          guide <- str_split_i(sample, "\\|", 2)
          target <- guide |>
            str_to_upper() |>
            str_remove("_\\d+$") |>
            recode(CD95 = "FAS")
          level <- target_cpm_in_line(target, line)

          bind_rows(
            tibble(
              gene = rp_genes,
              log2fc = log_change[rp_genes],
              class = "ribosomal protein genes"
            ),
            tibble(
              gene = target,
              log2fc = if (target %in% names(log_change)) log_change[[target]] else NA_real_,
              class = "the guide's own target"
            )
          ) |>
            mutate(
              guide = guide,
              target = target,
              target_cpm = level,
              quantifiable = level >= 10
            )
        }) |>
        list_rbind()
    })
  ) |>
  unnest(rows)

summary_table <- changes |>
  summarise(
    value = if (first(class) == "the guide's own target") first(log2fc) else mean(log2fc),
    .by = c(cell_line, guide, target, target_cpm, quantifiable, class)
  ) |>
  pivot_wider(names_from = class, values_from = value) |>
  mutate(
    ratio = `the guide's own target` / `ribosomal protein genes`
  )

cat("=== target log2 FC vs mean ribosomal log2 FC ===\n")
print(summary_table |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf)
write_csv(summary_table, here("output", "target_vs_rp_log2fc.csv"))

plot_changes <- function(data, title, subtitle, limits) {
  prepared <- data |>
    mutate(
      guide_label = if_else(quantifiable, guide, str_c(guide, " *")),
      guide_label = fct_reorder(guide_label, log2fc, .fun = \(x) mean(x, na.rm = TRUE))
    )

  labels <- prepared |>
    filter(
      class == "the guide's own target", !is.na(log2fc),
      log2fc >= limits[[1]], log2fc <= limits[[2]]
    )

  prepared |>
    ggplot(aes(log2fc, guide_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_boxplot(
      data = \(d) filter(d, class == "ribosomal protein genes"),
      aes(colour = class), outlier.size = 0.3, width = 0.6, linewidth = 0.35
    ) +
    geom_point(
      data = \(d) filter(d, class == "the guide's own target"),
      aes(colour = class), size = 2.6
    ) +
    geom_text(
      data = labels, aes(label = round(log2fc, 2)),
      vjust = -1.1, size = 2.3, colour = "#4292C6"
    ) +
    facet_wrap(vars(cell_line), scales = "free_y") +
    coord_cartesian(xlim = limits) +
    scale_colour_manual(
      values = c(
        `ribosomal protein genes` = "#B2182B",
        `the guide's own target` = "#4292C6"
      ), name = NULL
    ) +
    labs(title = title, subtitle = subtitle, x = "log2 FC vs mRuby controls", y = NULL)
}

full_view <- plot_changes(
  changes,
  "Target versus ribosomal genes on one scale",
  "box = 104 ribosomal protein genes; point = the guide's own target (* = target under 10 CPM, unreliable)",
  c(-3.1, 1.2)
)

zoom_view <- plot_changes(
  changes,
  "Same data, zoomed to the ribosomal range",
  "the targets sit far off this scale; the ribosomal shift is small but consistently negative",
  c(-0.35, 0.35)
)

ggsave(
  here("output", "target_vs_rp_distribution.png"),
  full_view / zoom_view + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom"),
  width = 11, height = 9, dpi = 300
)
