#!/usr/bin/env Rscript

# ==============================================================================
# Nutrimouse sCCA downstream interpretation
#
# Scope:
#   - Read and summarize existing sCCA outputs.
#   - Do NOT fit or refit an sCCA model.
#   - No-split defines final components.
#   - Five-fold results describe selection stability and canonical-r variation.
#   - The single train/test split is supplementary validation.
#
# Important:
#   Component signs are arbitrary. Positive/negative weights are coefficient
#   directions under the current orientation, not risk/protective effects.
#   Components 2 and 3 may be mismatched across folds; no component matching or
#   sign alignment is performed, and fold weights/loadings are never averaged.
# ==============================================================================

project_root <- "/Users/emily/Desktop/ubc/extracurricular/Yuxiao_Lab/论文复现1/replication_workspace/04_code_reviewed/tutorial_mixomics_mouse_scca"

no_split_dir <- file.path(project_root, "sCCA_nutrimouse_no_split_output")
single_dir <- file.path(project_root, "sCCA_nutrimouse_train_test_output")
fivefold_dir <- file.path(project_root, "sCCA_nutrimouse_5fold_output")
downstream_dir <- file.path(project_root, "sCCA_downstream_interpretation")
table_dir <- file.path(downstream_dir, "01_tables")
weight_plot_dir <- file.path(downstream_dir, "02_weight_plots")
loading_plot_dir <- file.path(downstream_dir, "03_loading_crossloading_plots")
validation_plot_dir <- file.path(downstream_dir, "04_validation_plots")
code_dir <- file.path(downstream_dir, "05_code")

for (d in c(downstream_dir, table_dir, weight_plot_dir, loading_plot_dir,
            validation_plot_dir, code_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

qc_path <- file.path(table_dir, "automated_QC_report.txt")
session_path <- file.path(code_dir, "sessionInfo.txt")

required_packages <- c(
  "ggplot2", "patchwork", "dplyr", "tidyr", "scales", "ragg", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       ". No figures were generated.")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

qc_records <- data.frame(
  Check = character(),
  Status = character(),
  Detail = character(),
  stringsAsFactors = FALSE
)

add_qc <- function(check, pass, detail) {
  qc_records <<- rbind(
    qc_records,
    data.frame(
      Check = check,
      Status = if (isTRUE(pass)) "PASS" else "FAIL",
      Detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  )
  invisible(pass)
}

write_qc <- function() {
  header <- c(
    "Nutrimouse sCCA downstream automated QC report",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R: ", R.version.string),
    paste0("R.home(): ", R.home()),
    paste0(".libPaths(): ", paste(.libPaths(), collapse = " | ")),
    paste0("mixOmics: ",
           if (requireNamespace("mixOmics", quietly = TRUE))
             as.character(packageVersion("mixOmics")) else "not installed"),
    paste0("nscancor: ",
           if (requireNamespace("nscancor", quietly = TRUE))
             as.character(packageVersion("nscancor")) else "not installed"),
    "",
    "This QC validates existing outputs only. No sCCA model was fitted.",
    "Components 2 and 3 may be misaligned across folds; matching/sign alignment was not performed.",
    ""
  )
  body <- apply(qc_records, 1L, function(x) {
    paste0("[", x[["Status"]], "] ", x[["Check"]], ": ", x[["Detail"]])
  })
  writeLines(c(header, body), qc_path, useBytes = TRUE)
}

abort_if_failed <- function(context) {
  if (any(qc_records$Status == "FAIL")) {
    write_qc()
    failed <- qc_records[qc_records$Status == "FAIL", , drop = FALSE]
    stop(
      context, "\n",
      paste0(failed$Check, ": ", failed$Detail, collapse = "\n"),
      "\nSee: ", qc_path
    )
  }
  invisible(TRUE)
}

required_files <- c(
  file.path(no_split_dir, "weights_X_gene.csv"),
  file.path(no_split_dir, "weights_Y_lipid.csv"),
  file.path(no_split_dir, "loadings_X_gene.csv"),
  file.path(no_split_dir, "loadings_Y_lipid.csv"),
  file.path(no_split_dir, "crossloadings_X_to_Yscore.csv"),
  file.path(no_split_dir, "crossloadings_Y_to_Xscore.csv"),
  file.path(no_split_dir, "selected_variables_X_gene.csv"),
  file.path(no_split_dir, "selected_variables_Y_lipid.csv"),
  file.path(no_split_dir, "canonical_correlations.csv"),
  file.path(fivefold_dir, "selection_frequency_X_gene_by_component.csv"),
  file.path(fivefold_dir, "selection_frequency_Y_lipid_by_component.csv"),
  file.path(fivefold_dir, "canonical_correlations_all_folds.csv"),
  file.path(fivefold_dir, "canonical_correlations_5fold_summary.csv"),
  file.path(single_dir, "canonical_correlations_train_test.csv")
)

missing_files <- required_files[!file.exists(required_files)]
add_qc(
  "Required input files",
  length(missing_files) == 0L,
  if (length(missing_files) == 0L) {
    paste(length(required_files), "required files found")
  } else {
    paste("Missing:", paste(missing_files, collapse = "; "))
  }
)
abort_if_failed("Required input-file check failed.")

read_csv_checked <- function(path, required_columns = NULL) {
  x <- read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "Inf", "-Inf")
  )
  if (!is.null(required_columns)) {
    missing <- setdiff(required_columns, names(x))
    add_qc(
      paste0("Columns: ", basename(path)),
      length(missing) == 0L,
      if (length(missing) == 0L) {
        paste("columns =", paste(names(x), collapse = ", "))
      } else {
        paste("missing columns =", paste(missing, collapse = ", "))
      }
    )
  }
  x
}

read_wide_component <- function(path, expected_n) {
  x <- read.csv(
    path, row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE
  )
  expected_components <- paste0("Component", 1:3)
  add_qc(
    paste0("Dimensions/components: ", basename(path)),
    nrow(x) == expected_n && identical(names(x), expected_components),
    paste0(
      "observed ", nrow(x), "x", ncol(x),
      "; columns=", paste(names(x), collapse = ","),
      "; expected ", expected_n, "x3"
    )
  )
  add_qc(
    paste0("Unique variable names: ", basename(path)),
    !anyDuplicated(rownames(x)) && all(nzchar(rownames(x))),
    paste("unique variables =", length(unique(rownames(x))))
  )
  x
}

wide_to_long <- function(x, value_name) {
  out <- data.frame(Variable = rownames(x), x, check.names = FALSE)
  out <- tidyr::pivot_longer(
    out,
    cols = starts_with("Component"),
    names_to = "Component",
    values_to = value_name
  )
  out$Component_number <- as.integer(sub("^Component", "", out$Component))
  out
}

finite_numeric <- function(x) {
  is.numeric(x) && all(is.finite(x))
}

in_correlation_range <- function(x, tolerance = 1e-10) {
  finite_numeric(x) && all(x >= -1 - tolerance & x <= 1 + tolerance)
}

# ------------------------------------------------------------------------------
# Read no-split objects
# ------------------------------------------------------------------------------

w_x_wide <- read_wide_component(file.path(no_split_dir, "weights_X_gene.csv"), 120L)
w_y_wide <- read_wide_component(file.path(no_split_dir, "weights_Y_lipid.csv"), 21L)
l_x_wide <- read_wide_component(file.path(no_split_dir, "loadings_X_gene.csv"), 120L)
l_y_wide <- read_wide_component(file.path(no_split_dir, "loadings_Y_lipid.csv"), 21L)
cl_x_wide <- read_wide_component(
  file.path(no_split_dir, "crossloadings_X_to_Yscore.csv"), 120L
)
cl_y_wide <- read_wide_component(
  file.path(no_split_dir, "crossloadings_Y_to_Xscore.csv"), 21L
)

selected_x <- read_csv_checked(
  file.path(no_split_dir, "selected_variables_X_gene.csv"),
  c("Component", "Variable", "Weight")
)
selected_y <- read_csv_checked(
  file.path(no_split_dir, "selected_variables_Y_lipid.csv"),
  c("Component", "Variable", "Weight")
)
no_split_r <- read_csv_checked(
  file.path(no_split_dir, "canonical_correlations.csv"),
  c("Component", "Canonical_r")
)

# ------------------------------------------------------------------------------
# Read five-fold and single-split objects
# ------------------------------------------------------------------------------

freq_x <- read_csv_checked(
  file.path(fivefold_dir, "selection_frequency_X_gene_by_component.csv"),
  c("Component", "Variable", "Selected_in_n_folds", "Selection_frequency")
)
freq_y <- read_csv_checked(
  file.path(fivefold_dir, "selection_frequency_Y_lipid_by_component.csv"),
  c("Component", "Variable", "Selected_in_n_folds", "Selection_frequency")
)
fivefold_all <- read_csv_checked(
  file.path(fivefold_dir, "canonical_correlations_all_folds.csv"),
  c("Fold", "Component", "N_train", "N_test",
    "Train_canonical_r", "Test_canonical_r")
)
fivefold_summary_source <- read_csv_checked(
  file.path(fivefold_dir, "canonical_correlations_5fold_summary.csv"),
  c("Component", "Mean_train_r", "SD_train_r",
    "Mean_test_r", "SD_test_r")
)
single_r <- read_csv_checked(
  file.path(single_dir, "canonical_correlations_train_test.csv"),
  c("Component", "Train_canonical_r", "Test_canonical_r")
)

abort_if_failed("Input column or dimension check failed.")

# ------------------------------------------------------------------------------
# Structural and value QC
# ------------------------------------------------------------------------------

expected_components <- paste0("Component", 1:3)
add_qc(
  "No-split has three components",
  setequal(no_split_r$Component, expected_components) && nrow(no_split_r) == 3L,
  paste("components =", paste(no_split_r$Component, collapse = ", "))
)
add_qc(
  "X block has 120 genes",
  nrow(w_x_wide) == 120L,
  paste("observed =", nrow(w_x_wide))
)
add_qc(
  "Y block has 21 lipids",
  nrow(w_y_wide) == 21L,
  paste("observed =", nrow(w_y_wide))
)

nonzero_x <- colSums(w_x_wide != 0)
nonzero_y <- colSums(w_y_wide != 0)
add_qc(
  "No-split nonzero X weights per component",
  all(nonzero_x == 10L),
  paste(names(nonzero_x), nonzero_x, collapse = "; ")
)
add_qc(
  "No-split nonzero Y weights per component",
  all(nonzero_y == 5L),
  paste(names(nonzero_y), nonzero_y, collapse = "; ")
)

add_qc(
  "Five-fold contains folds 1-5",
  identical(sort(unique(fivefold_all$Fold)), 1:5) &&
    all(table(fivefold_all$Fold) == 3L),
  paste("fold rows =", paste(names(table(fivefold_all$Fold)),
                             as.integer(table(fivefold_all$Fold)),
                             sep = ":", collapse = "; "))
)
add_qc(
  "Five-fold sample sizes are 32 train and 8 test",
  all(fivefold_all$N_train == 32L) && all(fivefold_all$N_test == 8L),
  paste0(
    "N_train=", paste(sort(unique(fivefold_all$N_train)), collapse = ","),
    "; N_test=", paste(sort(unique(fivefold_all$N_test)), collapse = ",")
  )
)

freq_valid <- function(x) {
  finite_numeric(x$Selection_frequency) &&
    all(x$Selection_frequency >= 0 & x$Selection_frequency <= 1) &&
    all(x$Selected_in_n_folds %in% 0:5) &&
    all(abs(x$Selection_frequency - x$Selected_in_n_folds / 5) < 1e-12)
}
add_qc(
  "X selection frequencies are valid",
  freq_valid(freq_x),
  paste("range =", paste(range(freq_x$Selection_frequency), collapse = " to "))
)
add_qc(
  "Y selection frequencies are valid",
  freq_valid(freq_y),
  paste("range =", paste(range(freq_y$Selection_frequency), collapse = " to "))
)

numeric_objects <- list(
  weights_X_gene = as.matrix(w_x_wide),
  weights_Y_lipid = as.matrix(w_y_wide),
  loadings_X_gene = as.matrix(l_x_wide),
  loadings_Y_lipid = as.matrix(l_y_wide),
  crossloadings_X_to_Yscore = as.matrix(cl_x_wide),
  crossloadings_Y_to_Xscore = as.matrix(cl_y_wide)
)
for (nm in names(numeric_objects)) {
  add_qc(
    paste0("Finite values: ", nm),
    finite_numeric(numeric_objects[[nm]]),
    paste("n =", length(numeric_objects[[nm]]))
  )
}

for (nm in c("loadings_X_gene", "loadings_Y_lipid",
             "crossloadings_X_to_Yscore", "crossloadings_Y_to_Xscore")) {
  add_qc(
    paste0("Correlation range [-1,1]: ", nm),
    in_correlation_range(numeric_objects[[nm]]),
    paste("range =", paste(range(numeric_objects[[nm]]), collapse = " to "))
  )
}

check_variable_sets <- function(reference, other, label) {
  same <- setequal(rownames(reference), rownames(other))
  add_qc(
    paste0("Variable-name match: ", label),
    same,
    if (same) {
      paste("matched", nrow(reference), "variables")
    } else {
      paste0(
        "missing from second=",
        paste(setdiff(rownames(reference), rownames(other)), collapse = ","),
        "; extra in second=",
        paste(setdiff(rownames(other), rownames(reference)), collapse = ",")
      )
    }
  )
}
check_variable_sets(w_x_wide, l_x_wide, "X weights vs loadings")
check_variable_sets(w_x_wide, cl_x_wide, "X weights vs cross-loadings")
check_variable_sets(w_y_wide, l_y_wide, "Y weights vs loadings")
check_variable_sets(w_y_wide, cl_y_wide, "Y weights vs cross-loadings")

expected_selected_x <- wide_to_long(w_x_wide, "Weight") |>
  filter(Weight != 0) |>
  select(Component, Variable, Weight)
expected_selected_y <- wide_to_long(w_y_wide, "Weight") |>
  filter(Weight != 0) |>
  select(Component, Variable, Weight)

selected_matches <- function(expected, observed) {
  e <- arrange(expected, Component, Variable)
  o <- arrange(observed[, c("Component", "Variable", "Weight")],
               Component, Variable)
  nrow(e) == nrow(o) &&
    identical(e$Component, o$Component) &&
    identical(e$Variable, o$Variable) &&
    all(abs(e$Weight - o$Weight) < 1e-12)
}
add_qc(
  "Selected X file matches nonzero no-split weights",
  selected_matches(expected_selected_x, selected_x),
  paste("nonzero weights =", nrow(expected_selected_x),
        "; selected rows =", nrow(selected_x))
)
add_qc(
  "Selected Y file matches nonzero no-split weights",
  selected_matches(expected_selected_y, selected_y),
  paste("nonzero weights =", nrow(expected_selected_y),
        "; selected rows =", nrow(selected_y))
)

freq_key_valid <- function(freq, variables) {
  !anyDuplicated(freq[c("Component", "Variable")]) &&
    all(freq$Component %in% expected_components) &&
    all(freq$Variable %in% variables)
}
add_qc(
  "X five-fold frequency keys match no-split variable universe",
  freq_key_valid(freq_x, rownames(w_x_wide)),
  paste("frequency rows =", nrow(freq_x))
)
add_qc(
  "Y five-fold frequency keys match no-split variable universe",
  freq_key_valid(freq_y, rownames(w_y_wide)),
  paste("frequency rows =", nrow(freq_y))
)

add_qc(
  "Canonical correlations are finite and in [-1,1]",
  in_correlation_range(no_split_r$Canonical_r) &&
    in_correlation_range(single_r$Train_canonical_r) &&
    in_correlation_range(single_r$Test_canonical_r) &&
    in_correlation_range(fivefold_all$Train_canonical_r) &&
    in_correlation_range(fivefold_all$Test_canonical_r),
  "no-split, single split, and all five-fold r values checked"
)

abort_if_failed("Structural or numeric QC failed.")

# ------------------------------------------------------------------------------
# Build complete component summary tables
# ------------------------------------------------------------------------------

build_summary <- function(weights, loadings, crossloadings, frequencies, block) {
  w <- wide_to_long(weights, "Weight") |>
    select(Component, Component_number, Variable, Weight)
  l <- wide_to_long(loadings, "Loading") |>
    select(Component, Variable, Loading)
  cl <- wide_to_long(crossloadings, "Cross_loading") |>
    select(Component, Variable, Cross_loading)

  out <- w |>
    left_join(l, by = c("Component", "Variable")) |>
    left_join(cl, by = c("Component", "Variable")) |>
    left_join(
      frequencies[, c("Component", "Variable",
                       "Selected_in_n_folds", "Selection_frequency")],
      by = c("Component", "Variable")
    ) |>
    mutate(
      Block = block,
      Abs_weight = abs(Weight),
      Abs_loading = abs(Loading),
      Abs_cross_loading = abs(Cross_loading),
      Selected_in_no_split = Weight != 0,
      Selected_in_n_folds = ifelse(
        is.na(Selected_in_n_folds), 0L, as.integer(Selected_in_n_folds)
      ),
      Selection_frequency = ifelse(
        is.na(Selection_frequency), 0, Selection_frequency
      )
    ) |>
    select(
      Component, Component_number, Block, Variable,
      Weight, Abs_weight, Loading, Abs_loading,
      Cross_loading, Abs_cross_loading,
      Selected_in_no_split, Selected_in_n_folds, Selection_frequency
    ) |>
    arrange(Component_number, Variable)

  if (nrow(out) != nrow(weights) * 3L ||
      anyNA(out) || anyDuplicated(out[c("Component", "Variable")])) {
    add_qc(
      paste0("Merged summary integrity: ", block),
      FALSE,
      paste0("rows=", nrow(out), "; anyNA=", anyNA(out),
             "; duplicate keys=",
             anyDuplicated(out[c("Component", "Variable")]) > 0)
    )
    abort_if_failed("Summary-table merge failed.")
  }
  add_qc(
    paste0("Merged summary integrity: ", block),
    TRUE,
    paste("rows =", nrow(out), "; complete unique Component+Variable keys")
  )
  out
}

summary_x <- build_summary(w_x_wide, l_x_wide, cl_x_wide, freq_x, "Gene")
summary_y <- build_summary(w_y_wide, l_y_wide, cl_y_wide, freq_y, "Lipid")

selected_summary_x <- summary_x |> filter(Selected_in_no_split)
selected_summary_y <- summary_y |> filter(Selected_in_no_split)

selected_xy <- bind_rows(selected_summary_x, selected_summary_y) |>
  mutate(
    Block_order = ifelse(Block == "Gene", 1L, 2L)
  ) |>
  arrange(
    Component_number, Block_order,
    desc(Selection_frequency), desc(Abs_loading), desc(Abs_weight)
  ) |>
  select(-Block_order)

write.csv(
  summary_x |> select(-Component_number),
  file.path(table_dir, "component_summary_X_gene_all_variables.csv"),
  row.names = FALSE
)
write.csv(
  summary_y |> select(-Component_number),
  file.path(table_dir, "component_summary_Y_lipid_all_variables.csv"),
  row.names = FALSE
)
write.csv(
  selected_summary_x |> select(-Component_number),
  file.path(table_dir, "component_summary_X_gene_selected_only.csv"),
  row.names = FALSE
)
write.csv(
  selected_summary_y |> select(-Component_number),
  file.path(table_dir, "component_summary_Y_lipid_selected_only.csv"),
  row.names = FALSE
)
write.csv(
  selected_xy |> select(-Component_number),
  file.path(table_dir, "component_summary_selected_XY.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Stability tables
# ------------------------------------------------------------------------------

stability_x <- summary_x |>
  transmute(
    Component, Variable, Selected_in_n_folds, Selection_frequency,
    No_split_weight = Weight,
    No_split_loading = Loading,
    No_split_cross_loading = Cross_loading
  ) |>
  arrange(as.integer(sub("^Component", "", Component)),
          desc(Selection_frequency), Variable)

stability_y <- summary_y |>
  transmute(
    Component, Variable, Selected_in_n_folds, Selection_frequency,
    No_split_weight = Weight,
    No_split_loading = Loading,
    No_split_cross_loading = Cross_loading
  ) |>
  arrange(as.integer(sub("^Component", "", Component)),
          desc(Selection_frequency), Variable)

write.csv(stability_x, file.path(table_dir, "stability_summary_X_gene.csv"),
          row.names = FALSE)
write.csv(stability_y, file.path(table_dir, "stability_summary_Y_lipid.csv"),
          row.names = FALSE)

high_stability <- bind_rows(
  summary_x |> filter(Selection_frequency >= 0.8),
  summary_y |> filter(Selection_frequency >= 0.8)
) |>
  transmute(
    Component, Component_number, Block, Variable, Selection_frequency,
    Selected_in_n_folds, Weight, Loading, Cross_loading
  ) |>
  arrange(Component_number, Block, desc(Selection_frequency), Variable) |>
  select(-Component_number)

write.csv(
  high_stability,
  file.path(table_dir, "high_stability_variables_frequency_0.8.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Canonical-r comparison table
# ------------------------------------------------------------------------------

fivefold_recalculated <- fivefold_all |>
  group_by(Component) |>
  summarise(
    Fivefold_mean_train_r = mean(Train_canonical_r),
    Fivefold_sd_train_r = sd(Train_canonical_r),
    Fivefold_mean_test_r = mean(Test_canonical_r),
    Fivefold_sd_test_r = sd(Test_canonical_r),
    .groups = "drop"
  )

summary_agreement <- fivefold_summary_source |>
  inner_join(fivefold_recalculated, by = "Component") |>
  transmute(
    ok = abs(Mean_train_r - Fivefold_mean_train_r) < 1e-12 &
      abs(SD_train_r - Fivefold_sd_train_r) < 1e-12 &
      abs(Mean_test_r - Fivefold_mean_test_r) < 1e-12 &
      abs(SD_test_r - Fivefold_sd_test_r) < 1e-12
  )
add_qc(
  "Five-fold summary matches fold-level recalculation",
  all(summary_agreement$ok) && nrow(summary_agreement) == 3L,
  paste("components checked =", nrow(summary_agreement))
)
abort_if_failed("Canonical-r summary reconciliation failed.")

canonical_comparison <- no_split_r |>
  rename(No_split_r = Canonical_r) |>
  left_join(
    single_r |>
      rename(
        Single_split_train_r = Train_canonical_r,
        Single_split_test_r = Test_canonical_r
      ),
    by = "Component"
  ) |>
  left_join(fivefold_recalculated, by = "Component") |>
  mutate(Component_number = as.integer(sub("^Component", "", Component))) |>
  arrange(Component_number) |>
  select(-Component_number)

if (anyNA(canonical_comparison) || nrow(canonical_comparison) != 3L) {
  add_qc(
    "Canonical-r model comparison merge",
    FALSE,
    paste("rows =", nrow(canonical_comparison),
          "; anyNA =", anyNA(canonical_comparison))
  )
  abort_if_failed("Canonical-r comparison merge failed.")
}
add_qc(
  "Canonical-r model comparison merge",
  TRUE,
  "three complete components"
)
write.csv(
  canonical_comparison,
  file.path(table_dir, "canonical_r_model_comparison.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Plot contract and helpers
# ------------------------------------------------------------------------------

palette <- c(
  "Positive weight" = "#2C7FB8",
  "Negative weight" = "#D95F0E",
  "Loading" = "#2C7FB8",
  "Cross-loading" = "#D95F0E",
  "Train" = "#2C7FB8",
  "Test" = "#D95F0E",
  "No-split" = "#2B2B2B",
  "Single split train" = "#4C78A8",
  "Single split test" = "#F28E2B",
  "Five-fold mean train" = "#59A14F",
  "Five-fold mean test" = "#B07AA1"
)

theme_publication <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.title = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 1, hjust = 0),
      panel.grid = element_blank()
    )
}

save_plot <- function(plot, path_without_extension, width, height,
                      save_png = TRUE) {
  svglite::svglite(
    paste0(path_without_extension, ".svg"),
    width = width, height = height, bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(
    filename = paste0(path_without_extension, ".pdf"),
    width = width, height = height, family = "sans", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  if (save_png) {
    ggsave(
      paste0(path_without_extension, ".png"), plot,
      device = ragg::agg_png,
      width = width, height = height, units = "in",
      dpi = 300, bg = "white"
    )
  }
}

component_caption <- paste(
  "Weights are coefficients used to construct canonical scores;",
  "the global sign of a component is arbitrary."
)

plot_weight_panel <- function(data, block_label, show_legend = TRUE) {
  d <- data |>
    filter(Weight != 0) |>
    arrange(Weight) |>
    mutate(
      Variable_plot = factor(Variable, levels = Variable),
      Weight_direction = ifelse(
        Weight >= 0, "Positive weight", "Negative weight"
      )
    )

  ggplot(d, aes(x = Weight, y = Variable_plot, fill = Weight_direction)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#666666") +
    geom_col(width = 0.72) +
    scale_fill_manual(values = palette[c("Positive weight", "Negative weight")]) +
    labs(x = "Weight", y = NULL, title = block_label) +
    theme_publication(9) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      axis.text.y = element_text(size = 8)
    )
}

component_plots_weight <- list()
for (comp in expected_components) {
  comp_n <- sub("^Component", "", comp)
  px <- plot_weight_panel(filter(summary_x, Component == comp), "Gene weights")
  py <- plot_weight_panel(filter(summary_y, Component == comp), "Lipid weights")
  combined <- (px | py) +
    plot_layout(guides = "collect", widths = c(1.35, 1)) +
    plot_annotation(
      title = paste0("No-split final component — Component ", comp_n),
      caption = component_caption,
      theme = theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.caption = element_text(size = 8, hjust = 0)
      )
    ) &
    theme(legend.position = "bottom")
  component_plots_weight[[comp]] <- combined
  save_plot(
    combined,
    file.path(weight_plot_dir, paste0("Component", comp_n, "_weights_XY")),
    width = 12.5, height = 7.2, save_png = TRUE
  )
}

weight_overview_panels <- list()
for (comp in expected_components) {
  comp_n <- sub("^Component", "", comp)
  weight_overview_panels[[length(weight_overview_panels) + 1L]] <-
    plot_weight_panel(
      filter(summary_x, Component == comp),
      paste0("Component ", comp_n, " — Gene"), show_legend = FALSE
    )
  weight_overview_panels[[length(weight_overview_panels) + 1L]] <-
    plot_weight_panel(
      filter(summary_y, Component == comp),
      paste0("Component ", comp_n, " — Lipid"), show_legend = FALSE
    )
}
weight_overview <- wrap_plots(weight_overview_panels, ncol = 2) +
  plot_annotation(
    title = "No-split final components: Gene and Lipid weights",
    caption = component_caption,
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )
save_plot(
  weight_overview,
  file.path(weight_plot_dir, "all_components_weights_XY_overview"),
  width = 13.5, height = 17, save_png = TRUE
)

plot_loading_panel <- function(data, block_label, show_legend = TRUE) {
  d <- data |>
    filter(Selected_in_no_split) |>
    arrange(Abs_loading) |>
    mutate(Variable_plot = factor(Variable, levels = Variable))
  points <- bind_rows(
    d |> transmute(Variable_plot, Metric = "Loading", Value = Loading),
    d |> transmute(
      Variable_plot, Metric = "Cross-loading", Value = Cross_loading
    )
  )

  ggplot(d, aes(y = Variable_plot)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
    geom_segment(
      aes(x = Loading, xend = Cross_loading,
          yend = Variable_plot),
      linewidth = 0.55, colour = "#BDBDBD"
    ) +
    geom_point(
      data = points,
      aes(x = Value, shape = Metric, colour = Metric),
      size = 2.1, stroke = 0.7
    ) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    scale_colour_manual(values = palette[c("Loading", "Cross-loading")]) +
    scale_shape_manual(values = c("Loading" = 16, "Cross-loading" = 17)) +
    labs(x = "Correlation", y = NULL, title = block_label) +
    theme_publication(9) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      axis.text.y = element_text(size = 8),
      plot.margin = margin(7, 14, 7, 14)
    )
}

loading_caption <- paste0(
  "Loading = correlation with the same-block canonical score; ",
  "Cross-loading = correlation with the opposite-block score.\n",
  "Only variables with nonzero no-split weights are shown. ",
  "The global sign of a component is arbitrary."
)

loading_overview_panels <- list()
for (comp in expected_components) {
  comp_n <- sub("^Component", "", comp)
  px <- plot_loading_panel(
    filter(summary_x, Component == comp), "Gene block"
  )
  py <- plot_loading_panel(
    filter(summary_y, Component == comp), "Lipid block"
  )
  combined <- (px | py) +
    plot_layout(guides = "collect", widths = c(1.25, 1)) +
    plot_annotation(
      title = paste0("Component ", comp_n,
                     " — no-split loading and cross-loading"),
      caption = loading_caption,
      theme = theme(
        plot.title = element_text(
          face = "bold", size = 11, margin = margin(b = 8)
        ),
        plot.caption = element_text(
          size = 8, hjust = 0, lineheight = 1.15, margin = margin(t = 8)
        )
      )
    ) &
    theme(legend.position = "bottom")
  save_plot(
    combined,
    file.path(
      loading_plot_dir,
      paste0("Component", comp_n, "_loading_crossloading_XY")
    ),
    width = 14.0, height = 7.5, save_png = TRUE
  )
  loading_overview_panels[[length(loading_overview_panels) + 1L]] <-
    plot_loading_panel(
      filter(summary_x, Component == comp),
      paste0("Component ", comp_n, " — Gene block"), show_legend = FALSE
    )
  loading_overview_panels[[length(loading_overview_panels) + 1L]] <-
    plot_loading_panel(
      filter(summary_y, Component == comp),
      paste0("Component ", comp_n, " — Lipid block"), show_legend = FALSE
    )
}

loading_overview <- wrap_plots(loading_overview_panels, ncol = 2) +
  plot_annotation(
    title = "No-split selected variables: loadings and cross-loadings",
    caption = loading_caption,
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )
save_plot(
  loading_overview,
  file.path(
    loading_plot_dir, "all_components_loading_crossloading_overview"
  ),
  width = 14.0, height = 17, save_png = TRUE
)

# ------------------------------------------------------------------------------
# Selection-frequency plots
# ------------------------------------------------------------------------------

plot_stability <- function(data, block_label) {
  d <- data |>
    filter(Selection_frequency > 0) |>
    mutate(
      Component = factor(Component, levels = expected_components),
      Variable_plot = reorder(
        paste(Component, Variable, sep = "___"),
        Selection_frequency
      ),
      Frequency_label = paste0(Selected_in_n_folds, "/5"),
      Label_hjust = ifelse(Selection_frequency >= 0.9, 1.15, -0.15)
    )

  ggplot(d, aes(x = Selection_frequency, y = Variable_plot)) +
    geom_vline(
      xintercept = c(0.6, 0.8), linetype = "dashed",
      linewidth = 0.4, colour = "#666666"
    ) +
    geom_segment(
      aes(x = 0, xend = Selection_frequency, yend = Variable_plot),
      linewidth = 0.6, colour = "#BDBDBD"
    ) +
    geom_point(size = 2.1, colour = "#2C7FB8") +
    geom_text(
      aes(label = Frequency_label, hjust = Label_hjust),
      size = 2.6, colour = "#333333"
    ) +
    facet_wrap(
      ~ Component, scales = "free_y", ncol = 1,
      labeller = as_labeller(function(x) sub("Component", "Component ", x))
    ) +
    scale_y_discrete(labels = function(x) sub("^Component[123]___", "", x)) +
    scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      labels = c("0.0", "0.2\n(1/5)", "0.4\n(2/5)",
                 "0.6\n(3/5)", "0.8\n(4/5)", "1.0\n(5/5)")
    ) +
    labs(
      x = "Selection frequency across five folds",
      y = NULL,
      title = paste0(block_label, " variable-selection stability"),
      caption = paste(
        "Dashed lines at 0.6 and 0.8 are descriptive references,",
        "not significance thresholds.",
        "Components are grouped by the same numeric label;",
        "Components 2 and 3 may be mismatched across folds."
      )
    ) +
    theme_publication(8.5) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 7),
      panel.spacing.y = grid::unit(0.8, "lines")
    )
}

stability_plot_x <- plot_stability(stability_x, "Gene")
stability_plot_y <- plot_stability(stability_y, "Lipid")
save_plot(
  stability_plot_x,
  file.path(validation_plot_dir, "selection_frequency_X_gene"),
  width = 9.5, height = 16, save_png = TRUE
)
save_plot(
  stability_plot_y,
  file.path(validation_plot_dir, "selection_frequency_Y_lipid"),
  width = 9.0, height = 10.5, save_png = TRUE
)

# ------------------------------------------------------------------------------
# Canonical-correlation validation plots
# ------------------------------------------------------------------------------

fold_r_long <- fivefold_all |>
  select(Fold, Component, Train_canonical_r, Test_canonical_r) |>
  pivot_longer(
    cols = c(Train_canonical_r, Test_canonical_r),
    names_to = "Set", values_to = "Canonical_r"
  ) |>
  mutate(
    Set = recode(
      Set,
      Train_canonical_r = "Train",
      Test_canonical_r = "Test"
    ),
    Component = factor(Component, levels = expected_components)
  )

p_fold <- ggplot(
  fold_r_long,
  aes(x = Fold, y = Canonical_r, colour = Set, linetype = Set, group = Set)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#888888") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~ Component, nrow = 1,
    labeller = as_labeller(function(x) sub("Component", "Component ", x))
  ) +
  scale_colour_manual(values = palette[c("Train", "Test")]) +
  scale_linetype_manual(values = c("Train" = "solid", "Test" = "dashed")) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(
    x = "Fold", y = "Canonical r",
    title = "Canonical correlation by five-fold split",
    caption = paste(
      "Train and held-out test r are shown for every fold.",
      "No fold was selected as best.",
      "Components 2 and 3 may be mismatched across folds."
    )
  ) +
  theme_publication(9)

save_plot(
  p_fold,
  file.path(validation_plot_dir, "canonical_r_by_fold"),
  width = 11.5, height = 4.5, save_png = TRUE
)

summary_long <- fivefold_recalculated |>
  pivot_longer(
    cols = -Component,
    names_to = c(".value", "Set"),
    names_pattern = "Fivefold_(mean|sd)_(train|test)_r"
  ) |>
  rename(Mean = mean, SD = sd) |>
  mutate(
    Set = recode(Set, train = "Train", test = "Test"),
    Component = factor(Component, levels = expected_components)
  )

p_summary <- ggplot(
  summary_long,
  aes(x = Component, y = Mean, colour = Set, shape = Set)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#888888") +
  geom_errorbar(
    aes(ymin = pmax(-1, Mean - SD), ymax = pmin(1, Mean + SD)),
    width = 0.15, linewidth = 0.6,
    position = position_dodge(width = 0.38)
  ) +
  geom_point(size = 2.8, position = position_dodge(width = 0.38)) +
  scale_colour_manual(values = palette[c("Train", "Test")]) +
  scale_shape_manual(values = c("Train" = 16, "Test" = 17)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_x_discrete(labels = function(x) sub("Component", "Component ", x)) +
  labs(
    x = NULL, y = "Mean canonical r ± SD",
    title = "Five-fold canonical-correlation summary",
    caption = paste(
      "Points are fold means; error bars are ±1 SD across five folds.",
      "Components 2 and 3 may be mismatched across folds."
    )
  ) +
  theme_publication(9)

save_plot(
  p_summary,
  file.path(validation_plot_dir, "canonical_r_5fold_summary"),
  width = 7.2, height = 5.2, save_png = TRUE
)

comparison_long <- bind_rows(
  canonical_comparison |>
    transmute(
      Component, Model = "No-split",
      Canonical_r = No_split_r, SD = NA_real_
    ),
  canonical_comparison |>
    transmute(
      Component, Model = "Single split train",
      Canonical_r = Single_split_train_r, SD = NA_real_
    ),
  canonical_comparison |>
    transmute(
      Component, Model = "Single split test",
      Canonical_r = Single_split_test_r, SD = NA_real_
    ),
  canonical_comparison |>
    transmute(
      Component, Model = "Five-fold mean train",
      Canonical_r = Fivefold_mean_train_r, SD = Fivefold_sd_train_r
    ),
  canonical_comparison |>
    transmute(
      Component, Model = "Five-fold mean test",
      Canonical_r = Fivefold_mean_test_r, SD = Fivefold_sd_test_r
    )
) |>
  mutate(
    Component = factor(Component, levels = expected_components),
    Model = factor(
      Model,
      levels = c(
        "No-split", "Single split train", "Single split test",
        "Five-fold mean train", "Five-fold mean test"
      )
    )
  )

p_comparison <- ggplot(
  comparison_long,
  aes(x = Model, y = Canonical_r, colour = Model, shape = Model)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#888888") +
  geom_errorbar(
    data = filter(comparison_long, !is.na(SD)),
    aes(ymin = pmax(-1, Canonical_r - SD),
        ymax = pmin(1, Canonical_r + SD)),
    width = 0.12, linewidth = 0.55
  ) +
  geom_point(size = 2.7) +
  facet_wrap(
    ~ Component, nrow = 1,
    labeller = as_labeller(function(x) sub("Component", "Component ", x))
  ) +
  scale_colour_manual(values = palette[levels(comparison_long$Model)]) +
  scale_shape_manual(values = c(16, 15, 17, 18, 8)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(
    x = NULL, y = "Canonical r",
    title = "Canonical-r comparison across validation designs",
    caption = paste(
      "No-split r is in-sample; single train/test is one random split;",
      "five-fold points are means ± SD and reflect split-to-split variation.",
      "The final model must not be selected according to the largest r."
    )
  ) +
  theme_publication(8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

save_plot(
  p_comparison,
  file.path(validation_plot_dir, "canonical_r_model_comparison"),
  width = 12.0, height = 5.8, save_png = TRUE
)

# ------------------------------------------------------------------------------
# Final output QC and session information
# ------------------------------------------------------------------------------

expected_output_files <- c(
  file.path(table_dir, c(
    "component_summary_X_gene_all_variables.csv",
    "component_summary_Y_lipid_all_variables.csv",
    "component_summary_X_gene_selected_only.csv",
    "component_summary_Y_lipid_selected_only.csv",
    "component_summary_selected_XY.csv",
    "stability_summary_X_gene.csv",
    "stability_summary_Y_lipid.csv",
    "high_stability_variables_frequency_0.8.csv",
    "canonical_r_model_comparison.csv"
  )),
  file.path(weight_plot_dir, c(
    paste0("Component", 1:3, "_weights_XY.pdf"),
    paste0("Component", 1:3, "_weights_XY.png"),
    "all_components_weights_XY_overview.pdf",
    "all_components_weights_XY_overview.png"
  )),
  file.path(loading_plot_dir, c(
    paste0("Component", 1:3, "_loading_crossloading_XY.pdf"),
    paste0("Component", 1:3, "_loading_crossloading_XY.png"),
    "all_components_loading_crossloading_overview.pdf",
    "all_components_loading_crossloading_overview.png"
  )),
  file.path(validation_plot_dir, c(
    "selection_frequency_X_gene.pdf",
    "selection_frequency_X_gene.png",
    "selection_frequency_Y_lipid.pdf",
    "selection_frequency_Y_lipid.png",
    "canonical_r_by_fold.pdf",
    "canonical_r_by_fold.png",
    "canonical_r_5fold_summary.pdf",
    "canonical_r_5fold_summary.png",
    "canonical_r_model_comparison.pdf",
    "canonical_r_model_comparison.png"
  ))
)

missing_outputs <- expected_output_files[
  !file.exists(expected_output_files) | file.info(expected_output_files)$size <= 0
]
add_qc(
  "Expected downstream outputs generated",
  length(missing_outputs) == 0L,
  if (length(missing_outputs) == 0L) {
    paste(length(expected_output_files), "non-empty expected outputs found")
  } else {
    paste("missing/empty =", paste(missing_outputs, collapse = "; "))
  }
)

session_lines <- c(
  "Nutrimouse sCCA downstream session information",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R.version.string: ", R.version.string),
  paste0("R.home(): ", R.home()),
  paste0(".libPaths(): ", paste(.libPaths(), collapse = " | ")),
  paste0("mixOmics: ",
         if (requireNamespace("mixOmics", quietly = TRUE))
           as.character(packageVersion("mixOmics")) else "not installed"),
  paste0("nscancor: ",
         if (requireNamespace("nscancor", quietly = TRUE))
           as.character(packageVersion("nscancor")) else "not installed"),
  "",
  capture.output(sessionInfo())
)
writeLines(session_lines, session_path, useBytes = TRUE)

write_qc()
abort_if_failed("Final downstream output QC failed.")

message("Nutrimouse downstream interpretation completed without fitting sCCA.")
message("QC report: ", qc_path)
