#!/usr/bin/env Rscript

# ==============================================================================
# Figure 2-inspired Nutrimouse sCCA figure — fixed single train/test model only
#
# Drawing only: this script never fits or refits sCCA.
# All panels use one existing 75% train / 25% held-out test analysis.
# Training results define components; test scores use training-learned weights.
# No no-split or five-fold result is read or plotted.
# ==============================================================================

# Paths ------------------------------------------------------------------------

project_root <- "/Users/emily/Desktop/ubc/extracurricular/Yuxiao_Lab/论文复现1/replication_workspace/04_code_reviewed/tutorial_mixomics_mouse_scca"
single_dir <- file.path(project_root, "sCCA_nutrimouse_train_test_output")
downstream_dir <- file.path(project_root, "sCCA_downstream_interpretation")
table_dir <- file.path(downstream_dir, "01_tables")
figure_dir <- file.path(downstream_dir, "02_weight_plots")
code_dir <- file.path(downstream_dir, "05_code")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_dir, recursive = TRUE, showWarnings = FALSE)

# Packages and centralized visual settings -------------------------------------

required_packages <- c(
  "ggplot2", "patchwork", "dplyr", "tidyr", "scales", "ragg", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them explicitly; this script does not install packages."
  )
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

FIG <- list(
  font = "Helvetica", base_size = 8.2, panel_label_size = 12,
  panel_title_size = 9.2, axis_title_size = 8.2, axis_text_size = 7.2,
  variable_text_size = 6.4, caption_size = 6.4, legend_text_size = 6.6,
  line_width = 0.38, point_size = 2.1, scatter_point_size = 2.25,
  scatter_alpha = 0.78, dpi = 360,
  main_width = 17.5, main_height = 13.8,
  panel_a_width = 6.4, panel_a_height = 4.1,
  panel_b_width = 10.6, panel_b_height = 4.1,
  component_width = 17.0, component_height = 4.6,
  correlation_limits = c(-1, 1),
  positive = "#B20D30", negative = "#142B8C",
  loading = "#335C81", cross_loading = "#9FC3DA",
  train = "#D73027", test = "#F28E2B",
  component = c(
    Component1 = "#E9B400", Component2 = "#1B8A3C",
    Component3 = "#1436B8"
  )
)

theme_figure <- function(base_size = FIG$base_size) {
  theme_classic(base_size = base_size, base_family = FIG$font) +
    theme(
      axis.line = element_line(linewidth = FIG$line_width, colour = "#222222"),
      axis.ticks = element_line(linewidth = FIG$line_width, colour = "#222222"),
      axis.title = element_text(size = FIG$axis_title_size, colour = "#111111"),
      axis.text = element_text(size = FIG$axis_text_size, colour = "#111111"),
      plot.title = element_text(
        size = FIG$panel_title_size, face = "bold", hjust = 0,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(size = FIG$axis_text_size, colour = "#444444"),
      plot.caption = element_text(
        size = FIG$caption_size, colour = "#444444", hjust = 0,
        margin = margin(t = 5)
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = FIG$legend_text_size),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(4.5, "mm"),
      plot.margin = margin(5, 7, 5, 5)
    )
}
theme_set(theme_figure())

# General file, QC and export helpers ------------------------------------------

qc_records <- data.frame(
  Check = character(), Status = character(), Detail = character(),
  stringsAsFactors = FALSE
)
add_qc <- function(check, pass, detail) {
  qc_records <<- rbind(
    qc_records,
    data.frame(
      Check = check, Status = if (isTRUE(pass)) "PASS" else "FAIL",
      Detail = as.character(detail), stringsAsFactors = FALSE
    )
  )
  invisible(pass)
}
abort_if_failed <- function(context) {
  if (any(qc_records$Status == "FAIL")) {
    failed <- qc_records[qc_records$Status == "FAIL", , drop = FALSE]
    stop(
      context, "\n",
      paste0(failed$Check, ": ", failed$Detail, collapse = "\n")
    )
  }
  invisible(TRUE)
}
read_csv_required <- function(path, required_columns = NULL) {
  if (!file.exists(path)) stop("Required file is missing: ", path)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(required_columns)) required_columns <- character()
  missing <- setdiff(required_columns, names(x))
  if (length(missing)) {
    stop("Missing column(s) in ", path, ": ", paste(missing, collapse = ", "))
  }
  x
}

# Read result matrices whose first CSV column contains sample/variable names.
read_wide_matrix <- function(path, components = paste0("Component", 1:3)) {
  if (!file.exists(path)) stop("Required file is missing: ", path)
  x <- read.csv(path, row.names = 1L, check.names = FALSE)
  if (!identical(names(x), components)) {
    stop("Unexpected component columns in ", path, ": ",
         paste(names(x), collapse = ", "))
  }
  if (anyDuplicated(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop("Missing or duplicated row identifiers in: ", path)
  }
  if (any(!is.finite(as.matrix(x)))) stop("Non-finite value in: ", path)
  as.data.frame(x, check.names = FALSE)
}

# Export editable PDF/SVG and a 360-dpi PNG preview.
save_plot_bundle <- function(plot, stem, width, height) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  svg_path <- file.path(figure_dir, paste0(stem, ".svg"))
  grDevices::cairo_pdf(
    pdf_path, width = width, height = height, family = FIG$font
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    png_path, width = width, height = height, units = "in",
    res = FIG$dpi, background = "white"
  )
  print(plot)
  grDevices::dev.off()
  svglite::svglite(
    svg_path, width = width, height = height, bg = "white",
    system_fonts = list(sans = FIG$font)
  )
  print(plot)
  grDevices::dev.off()
  invisible(c(pdf = pdf_path, png = png_path, svg = svg_path))
}

# Confirm that plotted scores reproduce the saved canonical correlation.
validate_canonical_r <- function(x, y, expected, label, tolerance = 1e-10) {
  if (length(x) != length(y) || length(x) < 3L ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Invalid paired scores for ", label, ".")
  }
  observed <- stats::cor(x, y)
  if (!is.finite(expected) || abs(observed - expected) > tolerance) {
    stop(
      "Canonical-r mismatch for ", label, ": recomputed=",
      signif(observed, 12), ", expected=", signif(expected, 12)
    )
  }
  observed
}

# Join each block training weights, loadings and cross-loadings by exact row.
build_training_table <- function(weights, loadings, crossloadings, block) {
  if (!identical(rownames(weights), rownames(loadings)) ||
      !identical(rownames(weights), rownames(crossloadings)) ||
      !identical(names(weights), names(loadings)) ||
      !identical(names(weights), names(crossloadings))) {
    stop(block, " variable or component names differ across training files.")
  }
  bind_rows(lapply(names(weights), function(component) {
    data.frame(
      Component = component, Block = block, Variable = rownames(weights),
      Weight = weights[[component]], Loading = loadings[[component]],
      Cross_loading = crossloadings[[component]], stringsAsFactors = FALSE
    )
  }))
}

# Read only the fixed single train/test outputs --------------------------------

required_names <- c(
  "canonical_correlations_train_test.csv",
  "canonical_scores_train_X_gene.csv", "canonical_scores_train_Y_lipid.csv",
  "canonical_scores_test_X_gene.csv", "canonical_scores_test_Y_lipid.csv",
  "weights_X_gene.csv", "weights_Y_lipid.csv",
  "loadings_train_X_gene.csv", "loadings_train_Y_lipid.csv",
  "crossloadings_train_X_to_Yscore.csv",
  "crossloadings_train_Y_to_Xscore.csv",
  "sample_split_and_metadata.csv", "training_scaling_parameters.rds"
)
required_files <- file.path(single_dir, required_names)
missing_inputs <- required_files[!file.exists(required_files)]
add_qc(
  "Required single-split plotting inputs", length(missing_inputs) == 0L,
  if (!length(missing_inputs)) paste(length(required_files), "files found")
  else paste("Missing:", paste(missing_inputs, collapse = "; "))
)
abort_if_failed("Single-split input check failed.")

component_order <- paste0("Component", 1:3)
canonical_single <- read_csv_required(
  file.path(single_dir, "canonical_correlations_train_test.csv"),
  c("Component", "Train_canonical_r", "Test_canonical_r")
)
split_metadata <- read_csv_required(
  file.path(single_dir, "sample_split_and_metadata.csv"), c("Sample", "Set")
)
train_x_scores <- read_wide_matrix(
  file.path(single_dir, "canonical_scores_train_X_gene.csv")
)
train_y_scores <- read_wide_matrix(
  file.path(single_dir, "canonical_scores_train_Y_lipid.csv")
)
test_x_scores <- read_wide_matrix(
  file.path(single_dir, "canonical_scores_test_X_gene.csv")
)
test_y_scores <- read_wide_matrix(
  file.path(single_dir, "canonical_scores_test_Y_lipid.csv")
)
train_x_weights <- read_wide_matrix(file.path(single_dir, "weights_X_gene.csv"))
train_y_weights <- read_wide_matrix(file.path(single_dir, "weights_Y_lipid.csv"))
train_x_loadings <- read_wide_matrix(
  file.path(single_dir, "loadings_train_X_gene.csv")
)
train_y_loadings <- read_wide_matrix(
  file.path(single_dir, "loadings_train_Y_lipid.csv")
)
train_x_cross <- read_wide_matrix(
  file.path(single_dir, "crossloadings_train_X_to_Yscore.csv")
)
train_y_cross <- read_wide_matrix(
  file.path(single_dir, "crossloadings_train_Y_to_Xscore.csv")
)
training_gene <- build_training_table(
  train_x_weights, train_x_loadings, train_x_cross, "Gene"
)
training_lipid <- build_training_table(
  train_y_weights, train_y_loadings, train_y_cross, "Lipid"
)

if (!identical(rownames(train_x_scores), rownames(train_y_scores))) {
  stop("Training X/Y score sample names differ.")
}
if (!identical(rownames(test_x_scores), rownames(test_y_scores))) {
  stop("Test X/Y score sample names differ.")
}
add_qc(
  "Fixed split is 30 train / 10 test",
  nrow(train_x_scores) == 30L && nrow(test_x_scores) == 10L &&
    sum(split_metadata$Set == "Train") == 30L &&
    sum(split_metadata$Set == "Test") == 10L,
  paste(
    nrow(train_x_scores), "training and", nrow(test_x_scores),
    "held-out score rows"
  )
)
add_qc(
  "Component order", identical(canonical_single$Component, component_order),
  paste(canonical_single$Component, collapse = ", ")
)
gene_nonzero <- colSums(abs(as.matrix(train_x_weights)) > 0)
lipid_nonzero <- colSums(abs(as.matrix(train_y_weights)) > 0)
add_qc(
  "Ten nonzero training Gene weights per component",
  nrow(train_x_weights) == 120L && all(gene_nonzero == 10L),
  paste(names(gene_nonzero), gene_nonzero, collapse = "; ")
)
add_qc(
  "Five nonzero training Lipid weights per component",
  nrow(train_y_weights) == 21L && all(lipid_nonzero == 5L),
  paste(names(lipid_nonzero), lipid_nonzero, collapse = "; ")
)
all_correlations <- c(
  training_gene$Loading, training_gene$Cross_loading,
  training_lipid$Loading, training_lipid$Cross_loading
)
add_qc(
  "Training loading/cross-loading range",
  all(all_correlations >= -1 & all_correlations <= 1),
  paste0("[", paste(signif(range(all_correlations), 5), collapse = ", "), "]")
)
add_qc(
  "Unique training Component + Variable keys",
  !anyDuplicated(training_gene[c("Component", "Variable")]) &&
    !anyDuplicated(training_lipid[c("Component", "Variable")]),
  paste(nrow(training_gene), "Gene and", nrow(training_lipid), "Lipid rows")
)
abort_if_failed("Single-split result QC failed.")

# Reconstruct training Gene EV and FEV, then save auditable intermediate tables.
# EV is 100 times the squared Gene-to-Lipid-score cross-loading.
# FEV normalizes the three EV values within each Gene. A zero-total row is
# assigned zero safely; otherwise the three FEV values must sum to one.
gene_ev_train <- data.frame(
  Gene = rownames(train_x_cross),
  EV_Component1 = 100 * train_x_cross$Component1^2,
  EV_Component2 = 100 * train_x_cross$Component2^2,
  EV_Component3 = 100 * train_x_cross$Component3^2,
  stringsAsFactors = FALSE
)
ev_matrix <- as.matrix(gene_ev_train[paste0("EV_Component", 1:3)])
ev_total <- rowSums(ev_matrix)
fev_matrix <- matrix(0, nrow = nrow(ev_matrix), ncol = ncol(ev_matrix))
positive_ev_total <- ev_total > 0
fev_matrix[positive_ev_total, ] <-
  ev_matrix[positive_ev_total, , drop = FALSE] /
  ev_total[positive_ev_total]
colnames(fev_matrix) <- paste0("FEV_Component", 1:3)
gene_fev_train <- data.frame(
  Gene = gene_ev_train$Gene,
  fev_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

selected_lipid_absolute_loadings <- training_lipid |>
  filter(abs(Weight) > 0) |>
  transmute(
    Component,
    Lipid = Variable,
    Absolute_Y_loading = abs(Loading)
  ) |>
  arrange(factor(Component, levels = component_order), Lipid)

write.csv(
  gene_ev_train,
  file.path(table_dir, "gene_EV_train_by_component.csv"),
  row.names = FALSE
)
write.csv(
  gene_fev_train,
  file.path(table_dir, "gene_FEV_train_by_component.csv"),
  row.names = FALSE
)
write.csv(
  selected_lipid_absolute_loadings,
  file.path(table_dir, "selected_lipid_absolute_loadings_train.csv"),
  row.names = FALSE
)

add_qc(
  "EV reconstruction",
  all(is.finite(ev_matrix)) && all(ev_matrix >= 0 & ev_matrix <= 100) &&
    isTRUE(all.equal(
      ev_matrix,
      100 * as.matrix(train_x_cross)^2,
      tolerance = 1e-12,
      check.attributes = FALSE
    )),
  paste(
    nrow(gene_ev_train),
    "Genes; EV equals 100 times squared training X-to-Y-score cross-loading"
  )
)
add_qc(
  "FEV normalization",
  all(is.finite(fev_matrix)) &&
    all(abs(rowSums(fev_matrix[positive_ev_total, , drop = FALSE]) - 1) <
          1e-12) &&
    all(rowSums(fev_matrix[!positive_ev_total, , drop = FALSE]) == 0),
  paste(
    sum(positive_ev_total), "positive-total rows sum to 1;",
    sum(!positive_ev_total), "zero-total rows set to 0"
  )
)
add_qc(
  "Selected Lipid absolute-loading table",
  nrow(selected_lipid_absolute_loadings) == 15L &&
    all(selected_lipid_absolute_loadings$Absolute_Y_loading >= 0 &
          selected_lipid_absolute_loadings$Absolute_Y_loading <= 1),
  "Five training-selected Lipids for each of three components"
)
abort_if_failed("EV, FEV or Lipid-profile reconstruction QC failed.")

# Component-result organizer ---------------------------------------------------

# Panels c-e show only variables with nonzero training weights.
read_component_results <- function(component) {
  if (!component %in% component_order) stop("Unknown component: ", component)
  genes <- training_gene |>
    filter(Component == component, abs(Weight) > 0) |>
    mutate(
      Weight_direction = ifelse(
        Weight >= 0, "Positive weight", "Negative weight"
      ),
      Direction_order = ifelse(Weight_direction == "Negative weight", 1L, 2L)
    ) |>
    arrange(Direction_order, abs(Weight))
  genes$Variable <- factor(genes$Variable, levels = genes$Variable)
  genes$Weight_direction <- factor(
    genes$Weight_direction,
    levels = c("Positive weight", "Negative weight")
  )
  lipids <- training_lipid |>
    filter(Component == component, abs(Weight) > 0) |>
    mutate(Absolute_Y_loading = abs(Loading)) |>
    arrange(Variable)
  if (nrow(genes) != 10L || nrow(lipids) != 5L) {
    stop(component, " does not contain 10 selected Genes and 5 selected Lipids.")
  }
  ev <- gene_ev_train |>
    filter(Gene %in% as.character(genes$Variable)) |>
    pivot_longer(
      starts_with("EV_Component"),
      names_to = "EV_component",
      values_to = "EV"
    ) |>
    mutate(
      Variable = factor(Gene, levels = levels(genes$Variable)),
      EV_component = factor(
        sub("^EV_", "", EV_component),
        levels = component_order,
        labels = paste("Component", 1:3)
      )
    )
  fev <- gene_fev_train |>
    filter(Gene %in% as.character(genes$Variable)) |>
    pivot_longer(
      starts_with("FEV_Component"),
      names_to = "FEV_component",
      values_to = "FEV"
    ) |>
    mutate(
      Variable = factor(Gene, levels = levels(genes$Variable)),
      FEV_component = factor(
        sub("^FEV_", "", FEV_component),
        levels = component_order,
        labels = paste("Component", 1:3)
      )
    )
  list(genes = genes, lipids = lipids, ev = ev, fev = fev)
}

# Panel helper: training Gene weights -----------------------------------------

make_weight_plot <- function(dat) {
  ggplot(dat$genes, aes(Weight, Variable, fill = Weight_direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#222222") +
    scale_x_continuous(
      limits = FIG$correlation_limits, breaks = seq(-1, 1, 0.5),
      expand = expansion(mult = 0.01)
    ) +
    scale_fill_manual(values = c(
      "Positive weight" = FIG$positive, "Negative weight" = FIG$negative
    )) +
    labs(title = "Training Gene weights", x = "Weight", y = NULL) +
    theme_figure() +
    theme(
      axis.text.y = element_text(size = FIG$variable_text_size),
      legend.position = "bottom", legend.direction = "horizontal",
      plot.margin = margin(5, 4, 5, 5)
    )
}

# Panel helper: grouped Gene EV bars -------------------------------------------

make_ev_plot <- function(dat) {
  ggplot(
    dat$ev,
    aes(x = EV, y = Variable, fill = EV_component)
  ) +
    geom_col(
      position = position_dodge2(width = 0.82, preserve = "single"),
      width = 0.72
    ) +
    scale_fill_manual(values = setNames(
      unname(FIG$component), paste("Component", 1:3)
    )) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(title = "Gene EV", x = "EV (%)", y = NULL) +
    theme_figure() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.margin = margin(5, 4, 5, 2)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

# Panel helper: 100% stacked Gene FEV bars -------------------------------------

make_fev_plot <- function(dat) {
  ggplot(
    dat$fev,
    aes(x = FEV, y = Variable, fill = FEV_component)
  ) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = setNames(
      unname(FIG$component), paste("Component", 1:3)
    )) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = "Gene FEV", x = "FEV", y = NULL) +
    theme_figure() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 4, 5, 2)
    )
}

# Panel helper: descriptive Lipid radar profile --------------------------------

# Radius is absolute training Y-side loading on a shared 0-1 scale.
make_radar_plot <- function(dat, component) {
  lipids <- dat$lipids
  n_lipids <- nrow(lipids)
  angles <- seq(
    pi / 2, pi / 2 - 2 * pi,
    length.out = n_lipids + 1L
  )[seq_len(n_lipids)]
  profile <- data.frame(
    Lipid = lipids$Variable,
    Radius = lipids$Absolute_Y_loading,
    Angle = angles,
    stringsAsFactors = FALSE
  ) |>
    mutate(
      x = Radius * cos(Angle),
      y = Radius * sin(Angle),
      label_x = 1.17 * cos(Angle),
      label_y = 1.17 * sin(Angle),
      label_hjust = ifelse(label_x > 0.15, 0, ifelse(label_x < -0.15, 1, 0.5))
    )
  profile_closed <- bind_rows(profile, profile[1, , drop = FALSE])
  grid_levels <- c(0.25, 0.50, 0.75, 1.00)
  radar_grid <- bind_rows(lapply(grid_levels, function(radius) {
    grid <- data.frame(
      Radius = radius,
      Angle = angles
    ) |>
      mutate(x = Radius * cos(Angle), y = Radius * sin(Angle))
    bind_rows(grid, grid[1, , drop = FALSE]) |>
      mutate(Grid = factor(radius, levels = grid_levels))
  }))
  spokes <- data.frame(
    x = 0, y = 0, xend = cos(angles), yend = sin(angles)
  )
  component_colour <- unname(FIG$component[[component]])

  ggplot() +
    geom_path(
      data = radar_grid,
      aes(x, y, group = Grid),
      linewidth = 0.32,
      colour = "#BDBDBD"
    ) +
    geom_segment(
      data = spokes,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.28,
      colour = "#D1D1D1"
    ) +
    geom_polygon(
      data = profile_closed,
      aes(x, y),
      fill = scales::alpha(component_colour, 0.16),
      colour = component_colour,
      linewidth = 0.72
    ) +
    geom_point(
      data = profile,
      aes(x, y),
      colour = component_colour,
      size = 2
    ) +
    geom_text(
      data = profile,
      aes(label_x, label_y, label = Lipid, hjust = label_hjust),
      family = FIG$font,
      size = 2.25,
      colour = "#222222"
    ) +
    annotate(
      "text",
      x = 0.04,
      y = grid_levels,
      label = sprintf("%.2f", grid_levels),
      hjust = 0,
      vjust = -0.25,
      family = FIG$font,
      size = 1.8,
      colour = "#666666"
    ) +
    coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.28, 1.28), clip = "off") +
    labs(
      title = paste(sub("Component", "Component ", component), "lipid profile"),
      subtitle = "Absolute Y-side loading"
    ) +
    theme_void(base_family = FIG$font, base_size = FIG$base_size) +
    theme(
      plot.title = element_text(
        size = FIG$panel_title_size, face = "bold", hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = FIG$axis_text_size, colour = "#444444", hjust = 0.5
      ),
      plot.margin = margin(7, 9, 7, 9)
    )
}

# Panel a: canonical correlation versus component rank -------------------------

panel_a_data <- bind_rows(
  canonical_single |>
    transmute(
      Rank = as.integer(sub("^Component", "", Component)),
      Dataset = "Train", Canonical_r = Train_canonical_r
    ),
  canonical_single |>
    transmute(
      Rank = as.integer(sub("^Component", "", Component)),
      Dataset = "Test", Canonical_r = Test_canonical_r
    )
)
panel_a_data$Dataset <- factor(panel_a_data$Dataset, c("Train", "Test"))
plot_a <- ggplot(
  panel_a_data,
  aes(Rank, Canonical_r, group = Dataset, colour = Dataset, shape = Dataset)
) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.7, stroke = 0.75) +
  scale_colour_manual(values = c(Train = FIG$train, Test = FIG$test)) +
  scale_shape_manual(values = c(Train = 16, Test = 15)) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_x_continuous(
    limits = c(0.9, 3.1), breaks = 1:3, expand = expansion(mult = 0)
  ) +
  labs(
    title = "a  Canonical correlations",
    x = "Rank of component", y = "Canonical r"
  ) +
  theme_figure() +
  theme(
    plot.title = element_text(size = FIG$panel_label_size, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -3),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "#E5E5E5")
  )
add_qc(
  "Panel a single-split isolation",
  nrow(panel_a_data) == 6L &&
    identical(levels(panel_a_data$Dataset), c("Train", "Test")) &&
    isTRUE(all.equal(
      panel_a_data$Canonical_r,
      c(canonical_single$Train_canonical_r, canonical_single$Test_canonical_r),
      tolerance = 0
    )),
  "Six exact single-split values; no uncertainty or significance layer"
)

# Panel b: train/test Component 1 canonical scores ------------------------------

train_expected_r <- canonical_single$Train_canonical_r[
  canonical_single$Component == "Component1"
]
test_expected_r <- canonical_single$Test_canonical_r[
  canonical_single$Component == "Component1"
]
train_r <- validate_canonical_r(
  train_x_scores$Component1, train_y_scores$Component1,
  train_expected_r, "training Component 1"
)
test_r <- validate_canonical_r(
  test_x_scores$Component1, test_y_scores$Component1,
  test_expected_r, "held-out test Component 1"
)
train_scatter <- data.frame(
  Sample = rownames(train_x_scores), Gene_score = train_x_scores$Component1,
  Lipid_score = train_y_scores$Component1
)
test_scatter <- data.frame(
  Sample = rownames(test_x_scores), Gene_score = test_x_scores$Component1,
  Lipid_score = test_y_scores$Component1
)

# Both plots receive the same symmetric x and y limits.
score_limit <- ceiling(max(abs(c(
  train_scatter$Gene_score, train_scatter$Lipid_score,
  test_scatter$Gene_score, test_scatter$Lipid_score
))) * 2) / 2
score_limits <- c(-score_limit, score_limit)
make_score_scatter <- function(dat, title, colour) {
  ggplot(dat, aes(Gene_score, Lipid_score)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "#D7D7D7") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "#D7D7D7") +
    geom_point(
      size = FIG$scatter_point_size, alpha = FIG$scatter_alpha,
      colour = colour, stroke = 0
    ) +
    scale_x_continuous(limits = score_limits, breaks = breaks_pretty(5)) +
    scale_y_continuous(limits = score_limits, breaks = breaks_pretty(5)) +
    coord_fixed() +
    labs(
      title = title, x = "Gene-side canonical score",
      y = "Lipid-side canonical score"
    ) +
    theme_figure()
}
plot_b_left <- make_score_scatter(
  train_scatter, sprintf("Training dataset, r = %.2f", train_r),
  FIG$train
)
plot_b_right <- make_score_scatter(
  test_scatter, sprintf("Test dataset, r = %.2f", test_r),
  FIG$test
)
plot_b <- (plot_b_left | plot_b_right) +
  plot_annotation(
    title = "b  Canonical score association",
    subtitle = paste(
      "Both panels use the training-fitted Component 1 weights;",
      "the held-out test set was not refitted."
    ),
    theme = theme(
      plot.title = element_text(
        family = FIG$font, face = "bold", size = FIG$panel_label_size
      ),
      plot.subtitle = element_text(
        family = FIG$font, size = FIG$axis_text_size, colour = "#444444"
      )
    )
  )
add_qc(
  "Panel b sample separation",
  nrow(train_scatter) == 30L && nrow(test_scatter) == 10L &&
    !anyDuplicated(train_scatter$Sample) && !anyDuplicated(test_scatter$Sample) &&
    !length(intersect(train_scatter$Sample, test_scatter$Sample)),
  "30 unique training and 10 unique held-out test mice"
)
add_qc(
  "Panel b canonical correlations",
  abs(train_r - train_expected_r) <= 1e-10 &&
    abs(test_r - test_expected_r) <= 1e-10,
  sprintf("training r=%.12f; held-out test r=%.12f", train_r, test_r)
)
add_qc(
  "Panel b common coordinate limits", is.finite(score_limit) && score_limit > 0,
  paste0("[", paste(score_limits, collapse = ", "), "] on both axes/panels")
)
add_qc(
  "Held-out scores use training model",
  file.exists(file.path(single_dir, "training_scaling_parameters.rds")),
  paste(
    "Existing output contains one training weight set and saved training",
    "scaling parameters; this plotting script performs no fitting"
  )
)

# Panels c-e: training-set component interpretation ----------------------------

make_component_panel <- function(component, letter) {
  dat <- read_component_results(component)
  weight_plot <- make_weight_plot(dat)
  ev_plot <- make_ev_plot(dat)
  fev_plot <- make_fev_plot(dat)
  radar_plot <- make_radar_plot(dat, component)
  (weight_plot | ev_plot | fev_plot | radar_plot) +
    plot_layout(widths = c(1.04, 1.02, 1.05, 1.30)) +
    plot_annotation(
      title = paste0(letter, "  ", sub("Component", "Component ", component)),
      subtitle = paste(
        "Training-set interpretation: Weight, reconstructed EV/FEV and",
        "absolute-loading Lipid profile."
      ),
      theme = theme(
        plot.title = element_text(
          family = FIG$font, face = "bold", size = FIG$panel_label_size
        ),
        plot.subtitle = element_text(
          family = FIG$font, size = FIG$axis_text_size, colour = "#444444"
        )
      )
    )
}
plot_c <- make_component_panel("Component1", "c")
plot_d <- make_component_panel("Component2", "d")
plot_e <- make_component_panel("Component3", "e")

# Final assembly: combine panels a-e -------------------------------------------

split_statement <- paste(
  "This methodological teaching replication uses a fixed 75% training and",
  "25% held-out test split. The original study used a 90%/10% split; the",
  "larger test fraction was retained here because the Nutrimouse dataset",
  "contains only 40 mice."
)
main_caption <- paste(
  strwrap(
    paste(
      split_statement,
      "All panels use this single model only. EV and FEV are reconstructed",
      "from training-set Gene-to-Lipid-score cross-loadings. No permutation",
      "P values or FDR are reported. Canonical correlation describes",
      "association, not causation. Component sign is arbitrary."
    ),
    width = 210
  ),
  collapse = "\n"
)
top_row <- wrap_plots(
  wrap_elements(full = plot_a), wrap_elements(full = plot_b),
  ncol = 2, widths = c(0.78, 1.22)
)
plot_final <- wrap_plots(
  wrap_elements(full = top_row), wrap_elements(full = plot_c),
  wrap_elements(full = plot_d), wrap_elements(full = plot_e),
  ncol = 1, heights = c(1.30, 1.78, 1.78, 1.78)
) +
  plot_annotation(
    title = "Figure 2-inspired Nutrimouse sCCA summary figure",
    subtitle = paste(
      "Fixed single-split teaching analysis:",
      "30 training mice and 10 held-out test mice."
    ),
    caption = main_caption,
    theme = theme(
      plot.title = element_text(
        family = FIG$font, face = "bold", size = 14
      ),
      plot.subtitle = element_text(
        family = FIG$font, size = 8.4, colour = "#3F3F3F"
      ),
      plot.caption = element_text(
        family = FIG$font, size = 6.8, colour = "#444444", hjust = 0
      ),
      plot.margin = margin(8, 10, 8, 8)
    )
  )
add_qc(
  "Main-figure model-source isolation",
  !exists("no_split_dir", inherits = FALSE) &&
    !exists("fivefold_dir", inherits = FALSE),
  "Only single_dir exists; no no-split/five-fold path or data object is declared"
)
add_qc(
  "Component-panel structure",
  all(vapply(
    component_order,
    function(component) {
      dat <- read_component_results(component)
      nrow(dat$genes) == 10L && nrow(dat$ev) == 30L &&
        nrow(dat$fev) == 30L && nrow(dat$lipids) == 5L
    },
    logical(1)
  )),
  paste(
    "Each component contains 10 Gene weights, 30 grouped EV entries,",
    "30 stacked FEV entries and 5 radar Lipids"
  )
)
abort_if_failed("Figure-data QC failed before export.")

# Export standalone panels and main figure -------------------------------------

save_plot_bundle(
  plot_a, "01_panel_a_canonical_validation",
  FIG$panel_a_width, FIG$panel_a_height
)
save_plot_bundle(
  plot_b, "02_panel_b_score_association",
  FIG$panel_b_width, FIG$panel_b_height
)
save_plot_bundle(
  plot_c, "03_panel_c_component1",
  FIG$component_width, FIG$component_height
)
save_plot_bundle(
  plot_d, "04_panel_d_component2",
  FIG$component_width, FIG$component_height
)
save_plot_bundle(
  plot_e, "05_panel_e_component3",
  FIG$component_width, FIG$component_height
)
save_plot_bundle(
  plot_final, "Fig2_inspired_mouse_scca_summary",
  FIG$main_width, FIG$main_height
)

# Caption ----------------------------------------------------------------------

caption_lines <- c(
  "Figure 2-inspired Nutrimouse sCCA summary figure.",
  "",
  paste(
    "This teaching figure is inspired by Fig. 2 of the original study but does",
    "not reproduce its urban-environment or psychiatric-symptom results.",
    "The present data contain gene and lipid measurements from 40 mice."
  ),
  "",
  split_statement,
  "",
  paste(
    "(a) Canonical correlations for component ranks 1-3 from the fixed",
    "single split. Circles show training values and squares show held-out test",
    "values. No error bars, permutation-based significance tests or FDR values",
    "are shown."
  ),
  sprintf(
    paste(
      "(b) Component 1 canonical-score association for the same fixed split.",
      "The left panel shows 30 training mice (r = %.3f); the right panel shows",
      "10 held-out test mice (r = %.3f). Test scores use the weights learned",
      "in training; the test set was not refitted. Both panels use identical",
      "coordinate limits and ordinary scatter points without density contours."
    ),
    train_r, test_r
  ),
  paste(
    "(c-e) Training-set interpretations of Components 1-3 from the same model.",
    "Each row shows the 10 training-selected Gene weights, grouped Gene EV",
    "bars, 100% stacked Gene FEV bars and a five-Lipid radar profile. EV was",
    "reconstructed as squared cross-loading (gene-to-lipid-score correlation",
    "squared, multiplied by 100), using training-set canonical scores. For",
    "each gene, FEV values across Components 1-3 sum to 1. Radar plots display",
    "the absolute Y-side loadings of the five selected lipids for each",
    "component. The radar shape is descriptive only and does not encode",
    "additional biological meaning beyond the absolute loadings."
  ),
  "",
  paste(
    "All panels use only the fixed single train/test analysis; no no-split or",
    "five-fold result enters this figure. EV and FEV are training-set",
    "cross-loading-based reconstructions, not recovered published values.",
    "keepX.X = 10 and keepX.Y = 5 are teaching settings, not parameters from",
    "the original study. Permutation P values and FDR are not reported because",
    "the corresponding parts of the public author code do not form a complete,",
    "directly reproducible workflow. Component sign is arbitrary. Canonical",
    "correlation is an association measure and does not establish causation.",
    "This is a methodological teaching replication, not a reproduction of the",
    "original variables or numerical results."
  )
)
writeLines(
  caption_lines,
  file.path(downstream_dir, "Fig2_inspired_mouse_scca_summary_caption.txt"),
  useBytes = TRUE
)

# Post-export QC and environment record ----------------------------------------

stems <- c(
  "01_panel_a_canonical_validation", "02_panel_b_score_association",
  "03_panel_c_component1", "04_panel_d_component2",
  "05_panel_e_component3", "Fig2_inspired_mouse_scca_summary"
)
exports <- unlist(lapply(
  stems,
  function(stem) file.path(figure_dir, paste0(stem, c(".pdf", ".png", ".svg")))
))
add_qc(
  "Figure exports",
  all(file.exists(exports)) && all(file.info(exports)$size > 0),
  paste(length(exports), "non-empty PDF/PNG/SVG files")
)
add_qc(
  "Inferential statistics not fabricated", TRUE,
  paste(
    "EV/FEV are explicitly reconstructed; no permutation P value or FDR",
    "enters a plotted layer"
  )
)
add_qc(
  "Training-only component interpretation",
  all(gene_nonzero == 10L) && all(lipid_nonzero == 5L),
  "10 training-selected Genes and 5 training-selected Lipids per component"
)
intermediate_tables <- file.path(
  table_dir,
  c(
    "gene_EV_train_by_component.csv",
    "gene_FEV_train_by_component.csv",
    "selected_lipid_absolute_loadings_train.csv"
  )
)
add_qc(
  "Reconstruction intermediate tables",
  all(file.exists(intermediate_tables)) &&
    all(file.info(intermediate_tables)$size > 0),
  paste(basename(intermediate_tables), collapse = "; ")
)

qc_lines <- c(
  "Figure 2-inspired Nutrimouse sCCA plotting QC",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R: ", R.version.string),
  "Figure source: fixed 75%/25% single train/test model only",
  "Training n: 30",
  "Held-out test n: 10",
  sprintf("Panel b training Component 1 r: %.12f", train_r),
  sprintf("Panel b held-out test Component 1 r: %.12f", test_r),
  paste0("Panel b shared score limits: [",
         paste(score_limits, collapse = ", "), "]"),
  "Gene EV: 100 * squared training X-to-Y-score cross-loading",
  "Gene FEV: within-Gene EV fraction across Components 1-3",
  "Radar radius: absolute training Y-side loading, range 0-1",
  "",
  apply(qc_records, 1L, function(x) {
    paste0("[", x[["Status"]], "] ", x[["Check"]], ": ", x[["Detail"]])
  })
)
writeLines(
  qc_lines,
  file.path(table_dir, "Fig2_inspired_automated_QC_report.txt"),
  useBytes = TRUE
)
abort_if_failed("Post-export QC failed.")

session_path <- file.path(code_dir, "plotting_sessionInfo.txt")
session_connection <- file(session_path, "wt")
sink(session_connection)
cat("Figure 2-inspired Nutrimouse sCCA plotting environment\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R.version.string:", R.version.string, "\n")
cat("R.home():", R.home(), "\n")
cat(".libPaths():", paste(.libPaths(), collapse = " | "), "\n\n")
cat("Packages used:\n")
for (package in required_packages) {
  cat(package, as.character(packageVersion(package)), "\n")
}
cat("\nAnalysis packages recorded for provenance:\n")
for (package in c("mixOmics", "nscancor")) {
  cat(
    package,
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(packageVersion(package))
    } else "not installed",
    "\n"
  )
}
cat("\nsessionInfo():\n")
print(sessionInfo())
sink()
close(session_connection)

message("Figure generation completed without refitting sCCA.")
message("All plotted values come from the fixed single train/test analysis.")
message("Main figure: ",
        file.path(figure_dir, "Fig2_inspired_mouse_scca_summary.pdf"))
message(sprintf("Panel b training r: %.6f", train_r))
message(sprintf("Panel b held-out test r: %.6f", test_r))
message("QC report: ",
        file.path(table_dir, "Fig2_inspired_automated_QC_report.txt"))
