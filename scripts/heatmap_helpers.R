library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(here)

# Shared plotting helpers for the SYM10X0 correlation heatmaps.
# Sourced by plot_correlation_heatmaps.R and correlation_by_phase.R.

read_corr <- function(path) {
  read_csv(path, show_col_types = FALSE) |>
    column_to_rownames(var = "...1") |>
    as.matrix()
}

# CPM + log1p, then leave-one-out logFC against the mRuby controls of the same
# cell line, matching the contrast used throughout the project
logfc_vs_controls <- function(counts) {
  logcpm <- log1p(prop.table(counts, margin = 2) * 1e6)
  logcpm <- logcpm[rowMeans(exp(logcpm) - 1) >= 1, ]

  samples <- colnames(logcpm)
  line_of <- str_split_i(samples, "\\|", 1)
  is_control <- str_detect(str_to_lower(samples), "mruby")

  samples |>
    set_names() |>
    map(\(sample) {
      controls <- samples[line_of == line_of[samples == sample] & is_control & samples != sample]
      logcpm[, sample] - rowMeans(logcpm[, controls, drop = FALSE])
    }) |>
    bind_cols() |>
    as.matrix() |>
    (\(x) {
      rownames(x) <- rownames(logcpm)
      x
    })()
}

# CD9_3 reaches only -0.13 log2 against CD9, a gene measured at >130 CPM in both
# lines, so it is a non-functional guide rather than an unmeasurable one. It is
# excluded from every group-level statistic but kept, labelled, in per-guide
# displays where it acts as an internal inert control. CD95_1/cd95_2 are NOT
# excluded: they look flat only in H1299, where FAS sits at 3.7 CPM, and they
# work in A549.
ineffective_guides <- c("CD9_3")

is_ineffective <- function(sample_ids) {
  str_split_i(sample_ids, "\\|", 2) |>
    str_to_upper() |>
    (\(g) g %in% str_to_upper(ineffective_guides))()
}

# correlation between functional targeting guides; mRuby rows are fixed by
# construction once they define the baseline
correlate_targeting <- function(logfc, drop_ineffective = TRUE) {
  targeting <- colnames(logfc) |> str_subset("mruby", negate = TRUE)
  if (drop_ineffective) {
    targeting <- targeting[!is_ineffective(targeting)]
  }
  cor(logfc[, targeting])
}

# guides are named inconsistently (CD95_1 vs cd95_2), so fold case before parsing
build_annotation_table <- function(sample_ids) {
  tibble(sample_id = sample_ids) |>
    separate_wider_delim(
      sample_id,
      delim = "|",
      names = c("cell_line", "guide"),
      cols_remove = FALSE
    ) |>
    mutate(
      guide = str_to_upper(guide),
      target = str_remove(guide, "_\\d+$"),
      target = case_match(
        target,
        "CD95" ~ "FAS (CD95)",
        "MRUBY" ~ "mRuby (non-targeting)",
        .default = target
      ),
      target = fct_relevel(target, "mRuby (non-targeting)", after = Inf)
    )
}

make_palette <- function(values, palette) {
  levels <- levels(as.factor(values))
  set_names(colorRampPalette(RColorBrewer::brewer.pal(8, palette))(length(levels)), levels)
}

build_annotation <- function(annotation_table, which, colours,
                             fields = c("Cell line", "Target", "Guide")) {
  values <- list(
    `Cell line` = annotation_table[["cell_line"]],
    Target = annotation_table[["target"]],
    Guide = annotation_table[["guide"]]
  )[fields]

  rlang::exec(
    HeatmapAnnotation,
    !!!values,
    which = which,
    col = colours[fields],
    show_annotation_name = which == "column",
    annotation_name_side = if (which == "column") "left" else "top",
    annotation_name_gp = gpar(fontsize = 9),
    show_legend = which == "column",
    gap = unit(0.8, "mm"),
    simple_anno_size = unit(4, "mm")
  )
}

plot_correlation_heatmap <- function(corr, annotation_table, title, col_fun, legend_title,
                                     annotation_colours, legend_breaks = NULL,
                                     annotation_fields = c("Cell line", "Target", "Guide")) {
  diag(corr) <- NA
  Heatmap(
    corr,
    name = legend_title,
    col = col_fun,
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    top_annotation = build_annotation(
      annotation_table, "column", annotation_colours, annotation_fields
    ),
    left_annotation = build_annotation(
      annotation_table, "row", annotation_colours, annotation_fields
    ),
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    row_dend_width = unit(12, "mm"),
    column_dend_height = unit(12, "mm"),
    rect_gp = gpar(col = "white", lwd = 0.5),
    na_col = "grey95",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(corr[i, j])) {
        grid.text(sprintf("%.2f", corr[i, j]), x, y, gp = gpar(fontsize = 5.5))
      }
    },
    heatmap_legend_param = list(
      at = legend_breaks,
      title_gp = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    )
  )
}

save_heatmap <- function(heatmap, path, width = 11, height = 9) {
  png(path, width = width, height = height, units = "in", res = 300)
  draw(heatmap,
    merge_legend = TRUE, heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  dev.off()
}

split_by_cell_line <- function(corr) {
  rownames(corr) |>
    str_split_i("\\|", 1) |>
    unique() |>
    sort() |>
    set_names() |>
    map(\(line) {
      members <- str_subset(rownames(corr), str_glue("^{line}\\|"))
      corr[members, members]
    })
}

# guide is the only thing that varies within a panel, so drop the cell line bar
within_line_heatmap <- function(corr, line, palette, legend_title, diverging,
                                show_stats = TRUE) {
  annotation <- build_annotation_table(rownames(corr))
  values <- corr[upper.tri(corr)]

  col_fun <- if (diverging) {
    limit <- max(abs(values))
    colorRamp2(c(-limit, 0, limit), c("#2166AC", "white", "#B2182B"))
  } else {
    colorRamp2(range(values), c("#2D004B", "#FDE725"))
  }

  breaks <- if (diverging) {
    seq(-max(abs(values)), max(abs(values)), length.out = 5) |> round(2)
  } else {
    seq(min(values), max(values), length.out = 5) |> round(3)
  }

  panel_title <- if (show_stats) {
    str_glue("{line}  (n = {nrow(corr)} guides, r = {round(min(values), 2)}-{round(max(values), 2)})")
  } else {
    line
  }

  plot_correlation_heatmap(
    corr,
    annotation,
    panel_title,
    col_fun,
    legend_title,
    palette,
    legend_breaks = breaks,
    annotation_fields = c("Target", "Guide")
  )
}

save_per_line_panels <- function(corr, file_stem, legend_title, diverging,
                                annotation_colours, show_titles = TRUE) {
  panels <- split_by_cell_line(corr)

  iwalk(panels, \(matrix, line) {
    save_heatmap(
      within_line_heatmap(matrix, line, annotation_colours, legend_title, diverging),
      here("output", str_glue("{file_stem}_{line}.png")),
      width = 7.5, height = 6.5
    )
  })

  # ComplexHeatmap draws to a device, so the vector version is a second pass
  # through the same code rather than a conversion of the raster
  draw_panels <- function(open_device) {
    open_device()
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(
      2, length(panels),
      heights = unit.c(unit(if (show_titles) 8 else 1, "mm"), unit(1, "null"))
    )))
    if (show_titles) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:length(panels)))
      grid.text(
        str_glue("Within-cell-line correlation - {legend_title}"),
        gp = gpar(fontsize = 13, fontface = "bold")
      )
      popViewport()
    }
    iwalk(panels, \(matrix, line) {
      pushViewport(viewport(
        layout.pos.row = 2,
        layout.pos.col = which(names(panels) == line)
      ))
      draw(
        within_line_heatmap(
          matrix, line, annotation_colours, str_glue("r ({line})"), diverging,
          show_stats = show_titles
        ),
        newpage = FALSE, merge_legend = TRUE,
        heatmap_legend_side = "right", annotation_legend_side = "right"
      )
      popViewport()
    })
    popViewport()
    dev.off()
  }

  # Side by side, each panel drawn into its own viewport. They must NOT be
  # combined with `+`: a HeatmapList forces one shared row order across
  # heatmaps, and these two panels have entirely different samples.
  draw_panels(\() png(
    here("output", str_glue("{file_stem}_by_cell_line.png")),
    width = 16, height = 7, units = "in", res = 300
  ))
  draw_panels(\() cairo_pdf(
    here("output", str_glue("{file_stem}_by_cell_line.pdf")),
    width = 16, height = 7
  ))
}
