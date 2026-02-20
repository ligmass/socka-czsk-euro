#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
})

set.seed(42)

# ---------- config ----------
DATA_PATH <- Sys.getenv("DID_DATA", "data/monthly_panel_clean.csv")
OUT_DIR   <- file.path("result tables baseline", "adoption", "dynamic_2009_binned_tails"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTCOMES  <- c("hicp_yoy","unemp_rate","hicp","imports_world_meur","exports_world_meur","log_imp","log_exp")
TREAT_DATE <- as.IDate("2009-01-01")

# Window settings
TAU_MIN <- -12L; TAU_MAX <- 18L
REF_TAU <- -1L

# Define the "Inner Window" where we want granular monthly detail
INNER_MIN <- -12L
INNER_MAX <- 12L

# Czech translations for variable names
czech_labels <- c(
  "exports_world_meur" = "Celkový export (mil. EUR)",
  "imports_world_meur" = "Celkový import (mil. EUR)",
  "log_exp" = "Logaritmus exportu",
  "log_imp" = "Logaritmus importu",
  "hicp" = "Cenová hladina (HICP)",
  "hicp_yoy" = "Meziroční inflace (HICP YoY)",
  "unemp_rate" = "Míra nezaměstnanosti (%)"
)

# Czech y-axis labels with units and confidence intervals
czech_y_labels <- c(
  "exports_world_meur" = "Rozdíl [mil. EUR]",
  "imports_world_meur" = "Rozdíl [mil. EUR]",
  "log_exp" = "Efekt [log. diference]",
  "log_imp" = "Efekt [log. diference]",
  "hicp" = "Rozdíl [index. body]",
  "hicp_yoy" = "Rozdíl [p.b.]",
  "unemp_rate" = "Rozdíl [p.b.]"
)

# ---------- helpers ----------
auto_lag <- function(T){
  L <- floor(4 * (T/100)^(2/9))
  max(0L, min(L, 12L))
}

vcov_failsoft <- function(fit, L){
  V <- tryCatch(NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE), error = function(e) e)
  if (inherits(V, "error") || any(!is.finite(V))) V <- tryCatch(vcov(fit), error = function(e) e)
  if (inherits(V, "error") || any(!is.finite(V))) {
    X <- model.matrix(fit)
    XtX <- crossprod(X)
    lam <- 1e-8
    V <- tryCatch({
      sigma2 <- sum(resid(fit)^2) / max(1, nrow(X) - qr(X)$rank)
      sigma2 * solve(XtX + diag(lam, ncol(XtX)))
    }, error = function(e) matrix(NA_real_, ncol(XtX), ncol(XtX)))
    rownames(V) <- colnames(V) <- colnames(X)
  }
  V
}

run_one_series_dynamic_binned <- function(dat, y, break_date){
  # 1. Prepare Data
  need <- c("date","country", y)
  if (!all(need %in% names(dat))) return(NULL)
  d0 <- dat[, ..need]
  d0 <- dcast(d0, date ~ country, value.var = y)
  setnames(d0, c("CZ","SK"), c("cz","sk"), skip_absent = TRUE)
  d0[, diff := cz - sk]

  # 2. Create Event Time
  d0 <- d0[order(date)]
  d0[, time_id := as.integer(factor(date, levels = sort(unique(date))))]
  t0 <- unique(d0[date == break_date, time_id])
  if (length(t0) != 1L) return(NULL)
  d0[, tau := time_id - t0]
  
  # 3. Filter Window
  d0 <- d0[tau >= TAU_MIN & tau <= TAU_MAX]
  if (!nrow(d0)) return(NULL)

  # 4. Create "Hybrid" Bins
  # - Specific months for the inner window (-6 to +6)
  # - "Pre_Tail" for everything before -6
  # - "Post_Tail" for everything after +6
  d0[, bin_label := as.character(tau)]
  d0[tau < INNER_MIN, bin_label := "Pre_Tail"]
  d0[tau > INNER_MAX, bin_label := "Post_Tail"]
  
  # Create Factor with Reference
  # Reference is -1, which is inside our inner window
  d0[, tau_factor := factor(bin_label, levels = c("Pre_Tail", as.character(INNER_MIN:INNER_MAX), "Post_Tail"))]
  d0[, tau_factor := stats::relevel(tau_factor, ref = as.character(REF_TAU))]

  # 5. Run Model (Now we have enough degrees of freedom!)
  fit <- lm(diff ~ tau_factor, data = d0)
  
  # 6. Extract Diagnostics
  s <- summary(fit)
  r2 <- s$r.squared
  ar2 <- s$adj.r.squared
  
  bp <- bptest(fit)
  bp_p <- bp$p.value

  # 7. Extract Coefficients
  cf <- coef(fit)
  V  <- vcov_failsoft(fit, auto_lag(nrow(d0)))
  
  # Grab all tau_factor coefficients
  tau_names <- grep("^tau_factor", names(cf), value = TRUE)
  
  se <- sqrt(diag(V))[tau_names]
  beta <- cf[tau_names]
  z <- beta / se
  p <- 2 * pnorm(-abs(z))
  
  # Helper to map names back to tau/label
  clean_names <- sub("^tau_factor", "", tau_names)
  
  res <- data.table(
    outcome = y,
    label = clean_names,
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  )
  
  # Add reference row
  ref_row <- data.table(outcome = y, label = as.character(REF_TAU), beta = 0, se = 0, p = 1)
  res <- rbind(res, ref_row)
  
  # Add sorting index for nice display
  res[, sort_idx := as.numeric(label)]
  # Handle "Pre_Tail" and "Post_Tail" sorting
  res[label == "Pre_Tail", sort_idx := -999]
  res[label == "Post_Tail", sort_idx := 999]
  res <- res[order(sort_idx)]

  return(list(res = res, r2 = r2, ar2 = ar2, bp_p = bp_p))
}

# ---------- main ----------
dat <- fread(DATA_PATH)
if (!inherits(dat$date, "IDate")) dat[, date := as.IDate(date)]

es_list <- list()
cat("\n--- STARTING DYNAMIC RUN (BINNED TAILS) ---\n")

for (y in OUTCOMES[OUTCOMES %in% names(dat)]) {
  ret <- run_one_series_dynamic_binned(dat, y, TREAT_DATE)
  
  if (!is.null(ret)) {
    res <- ret$res
    # attach GOF diagnostics to the per-outcome table
    res[, r2 := ret$r2]
    res[, ar2 := ret$ar2]
    res[, bp_p := ret$bp_p]
    es_list[[y]] <- res
    
    cat(sprintf("\n========================\nOutcome: %s\n", y))
    cat(sprintf("R-squared: %.4f, Adj. R-squared: %.4f\n", ret$r2, ret$ar2))
    cat(sprintf("Heteroskedasticity (BP Test) p-value: %.4g\n", ret$bp_p))
    print(res[, .(label, beta = round(beta, 4), se = round(se, 4), p = format.pval(p, digits=3))], 
          row.names = FALSE) 
  }
}

es_dynamic <- rbindlist(es_list)
fwrite(es_dynamic, file.path(OUT_DIR, "es_dynamic_2009_binned.csv"))

# Create academic event-study plots per outcome
file_prefix <- "adoption"
event_name <- "Přijetí eura na Slovensku"
try({
  for (y in unique(es_dynamic$outcome)){
    sub <- es_dynamic[outcome == y]
    if (!nrow(sub)) next
    
    # Filter out Pre_Tail and Post_Tail
    sub <- sub[label != "Pre_Tail" & label != "Post_Tail"]
    if (!nrow(sub)) next
    
    # Convert sort_idx to character for proper numeric conversion
    sub[, x_pos := as.numeric(as.character(label))]
    sub <- sub[!is.na(x_pos)]
    
    # Calculate confidence intervals
    sub[, ymin := beta - 1.96 * se]
    sub[, ymax := beta + 1.96 * se]
    
    # Create plot with academic styling
    p <- ggplot(sub, aes(x = x_pos, y = beta)) +
      # Reference lines
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.5) +
      geom_vline(xintercept = -1, linetype = "dotted", color = "red", linewidth = 0.6) +
      # Geoms for point estimates and CIs
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.3, linewidth = 0.5) +
      geom_point(size = 2.5, color = "black") +
      geom_line(aes(group = 1), linewidth = 0.5, color = "black") +
      # Academic theme
      theme_classic(base_size = 14) +
      labs(
        title = event_name,
        subtitle = ifelse(y %in% names(czech_labels), czech_labels[y], y),
        x = "Měsíce relativně k události",
        y = ifelse(y %in% names(czech_y_labels), czech_y_labels[y], "Odhadovaný efekt")
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.grid = element_blank()
      )
    
    ggsave(filename = file.path(OUT_DIR, paste0("es_dynamic_", file_prefix, "_", y, ".png")), 
           plot = p, width = 8, height = 5, dpi = 300)
  }
}, silent = TRUE)

cat("\nDone. Event study plots saved in ", OUT_DIR, "\n", sep = "")