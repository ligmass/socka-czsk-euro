#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

export_academic_table <- function(raw_df, output_filepath) {
  dt <- as.data.table(raw_df)
  required_cols <- c("outcome", "bin", "beta", "se", "p", "n", "r2")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols)) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  dt <- copy(dt)
  dt[, bin_key := sub("^bin", "", as.character(bin))]

  dt[, stars := fifelse(p < 0.01, "***",
                 fifelse(p < 0.05, "**",
                 fifelse(p < 0.1, "*", "")))]
  dt[, cell_value := sprintf("%.2f%s (%.2f)", beta, stars, se)]

  bin_labels <- c(
    "pre12_7" = "Pre-trend (-12 až -7)",
    "t0" = "Měsíc události (t=0)",
    "post1_6" = "1 až 6 měsíců po",
    "post7_12" = "7 až 12 měsíců po",
    "post13_18" = "13 až 18 měsíců po",
    "Pre_Tail" = "Pre-Tail",
    "Post_Tail" = "Post-Tail"
  )

  outcome_labels <- c(
    "exports_world_meur" = "Export (mil. EUR)",
    "imports_world_meur" = "Import (mil. EUR)",
    "log_exp" = "Export (log)",
    "log_imp" = "Import (log)",
    "hicp" = "HICP (index)",
    "hicp_yoy" = "Inflace YoY (p.b.)",
    "unemp_rate" = "Nezaměstnanost (p.b.)"
  )

  bin_order <- c("pre12_7", "t0", "post1_6", "post7_12", "post13_18", "Pre_Tail", "Post_Tail")
  dt[, bin_cz := fcoalesce(unname(bin_labels[bin_key]), bin_key)]
  dt[, bin_cz := factor(bin_cz, levels = unname(bin_labels[bin_order]), ordered = TRUE)]
  dt[, outcome_cz := fcoalesce(unname(outcome_labels[outcome]), outcome)]

  wide <- dcast(
    dt,
    bin_cz ~ outcome_cz,
    value.var = "cell_value",
    fun.aggregate = function(x) if (length(x)) x[1] else NA_character_
  )
  setnames(wide, "bin_cz", "Bin")

  n_stats <- dt[, .(n_val = as.integer(max(n, na.rm = TRUE))), by = .(outcome_cz)]
  r2_stats <- dt[, .(r2_val = max(r2, na.rm = TRUE)), by = .(outcome_cz)]

  all_outcome_cols <- setdiff(names(wide), "Bin")
  n_row <- data.table(Bin = "N")
  r2_row <- data.table(Bin = "R-squared")

  for (oc in all_outcome_cols) {
    n_value <- n_stats[outcome_cz == oc, n_val]
    r2_value <- r2_stats[outcome_cz == oc, r2_val]
    n_row[, (oc) := if (length(n_value) && is.finite(n_value)) as.character(n_value[1]) else NA_character_]
    r2_row[, (oc) := if (length(r2_value) && is.finite(r2_value)) sprintf("%.2f", r2_value[1]) else NA_character_]
  }

  out <- rbindlist(list(wide, n_row, r2_row), use.names = TRUE, fill = TRUE)

  out_dir <- dirname(output_filepath)
  if (!dir.exists(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(out, output_filepath)

  invisible(out)
}
