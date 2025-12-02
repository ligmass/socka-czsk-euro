#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
})

set.seed(42)

# ---------- config ----------
DATA_PATH <- Sys.getenv("DID_DATA", "data/monthly_panel_clean.csv")
OUT_DIR   <- "tables5_dynamic_ukraine"; dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTCOMES  <- c("hicp_yoy","unemp_rate","hicp","imports_world_meur","exports_world_meur","log_imp","log_exp")

# UKRAINE WAR DATE
TREAT_DATE <- as.IDate("2022-02-01")

# Window settings
TAU_MIN <- -12L; TAU_MAX <- 18L
REF_TAU <- -1L

# Inner Window for granular detail (-6 to +6 months)
INNER_MIN <- -6L
INNER_MAX <- 6L

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
  d0[, bin_label := as.character(tau)]
  d0[tau < INNER_MIN, bin_label := "Pre_Tail"]
  d0[tau > INNER_MAX, bin_label := "Post_Tail"]
  
  # Create Factor (ensure all levels exist or handle missing)
  # We explicitly define levels to ensure order
  d0[, tau_factor := factor(bin_label, levels = c("Pre_Tail", as.character(INNER_MIN:INNER_MAX), "Post_Tail"))]
  d0[, tau_factor := stats::relevel(tau_factor, ref = as.character(REF_TAU))]

  # 5. Run Model (Simple Spec: No controls to avoid saturation on short window)
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
  tau_names <- grep("^tau_factor", names(cf), value = TRUE)
  
  se <- sqrt(diag(V))[tau_names]
  beta <- cf[tau_names]
  z <- beta / se
  p <- 2 * pnorm(-abs(z))
  
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
  
  # Sort
  res[, sort_idx := as.numeric(label)]
  res[label == "Pre_Tail", sort_idx := -999]
  res[label == "Post_Tail", sort_idx := 999]
  res <- res[order(sort_idx)]

  return(list(res = res, r2 = r2, ar2 = ar2, bp_p = bp_p))
}

# ---------- main ----------
dat <- fread(DATA_PATH)
if (!inherits(dat$date, "IDate")) dat[, date := as.IDate(date)]

es_list <- list()
cat("\n--- STARTING DYNAMIC RUN (UKRAINE BINNED) ---\n")

for (y in OUTCOMES[OUTCOMES %in% names(dat)]) {
  ret <- run_one_series_dynamic_binned(dat, y, TREAT_DATE)
  
  if (!is.null(ret)) {
    res <- ret$res
    es_list[[y]] <- res
    
    cat(sprintf("\n========================\nOutcome: %s\n", y))
    cat(sprintf("R-squared: %.4f, Adj. R-squared: %.4f\n", ret$r2, ret$ar2))
    print(res[, .(label, beta = round(beta, 4), se = round(se, 4), p = format.pval(p, digits=3))], 
          row.names = FALSE) 
  }
}

es_dynamic <- rbindlist(es_list)
fwrite(es_dynamic, file.path(OUT_DIR, "es_dynamic_ukraine.csv"))
cat("\nDone.\n")