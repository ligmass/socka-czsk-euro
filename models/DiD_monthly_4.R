#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  wcb_ok <- requireNamespace("fwildclusterboot", quietly = TRUE)
  if (!wcb_ok) message("NOTE: fwildclusterboot not installed -> wild cluster p-values will be NA.")
  if (identical(Sys.getenv("USE_WCB"), "0")) { wcb_ok <- FALSE; message("NOTE: USE_WCB=0 -> skipping wild cluster bootstrap.") }
})

set.seed(42)
if (requireNamespace("dqrng", quietly = TRUE)) dqrng::dqset.seed(42L)

# -------------------- config --------------------
VERBOSE <- identical(Sys.getenv("VERBOSE", "1"), "1")
args <- commandArgs(trailingOnly = TRUE)
DATA_PATH <- if (length(args) >= 1) args[1] else Sys.getenv("DID_DATA", "data/monthly_panel_clean.csv")

OUTCOMES <- c("hicp_yoy","unemp_rate","hicp","imports_world_meur","exports_world_meur","mro","log_imp","log_exp")
CAND_CONTROLS <- c("fx_eur","fx_vol","repo","mro")

TREAT_DATE <- as.IDate("2009-01-01")
OUT_DIR <- "tables4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------- helpers --------------------
balance_by_time <- function(dt) {
  if (!nrow(dt)) return(dt)
  full_t <- dt[, .N, by = .(time_id, country)][, .N, by = time_id][N == 2L, time_id]
  dt[time_id %in% full_t]
}

miss_audit <- function(d, vars) {
  rbindlist(lapply(vars, function(v) {
    data.table(
      var           = v,
      pNA_all       = mean(is.na(d[[v]])),
      pNA_treatpost = mean(is.na(d[treated == 1 & post == 1][[v]])),
      pNA_treatpre  = mean(is.na(d[treated == 1 & post == 0][[v]])),
      pNA_control   = mean(is.na(d[treated == 0][[v]]))
    )
  }))
}

select_controls <- function(dat, y, cand_controls, min_n = 100L) {
  kept <- character(0)
  base_cols <- c("country","time_id","treated","post","D", y)
  for (cvar in cand_controls) {
    try_cols <- unique(c(base_cols, kept, cvar))
    if (!all(try_cols %in% names(dat))) next
    df_try <- na.omit(dat[, ..try_cols])
    df_try <- balance_by_time(df_try)
    if (nrow(df_try) >= min_n &&
        uniqueN(df_try$country) == 2L &&
        uniqueN(df_try$time_id)  >  1L &&
        uniqueN(df_try$D)        >  1L &&
        nrow(df_try[treated == 1 & post == 0]) > 0L) {
      kept <- c(kept, cvar)
    }
  }
  if (length(kept) == 0L) {
    for (k in length(cand_controls):0) {
      try_cols <- unique(c(base_cols, cand_controls[seq_len(k)]))
      df_try <- na.omit(dat[, ..try_cols])
      df_try <- balance_by_time(df_try)
      if (nrow(df_try) >= min_n &&
          uniqueN(df_try$country) == 2L &&
          uniqueN(df_try$time_id)  >  1L &&
          uniqueN(df_try$D)        >  1L &&
          nrow(df_try[treated == 1 & post == 0]) > 0L) {
        kept <- cand_controls[seq_len(k)]
        break
      }
    }
  }
  kept
}

mk_twfe_formula <- function(y, controls) {
  rhs <- if (length(controls)) paste("D +", paste(controls, collapse = " + ")) else "D"
  as.formula(paste0(y, " ~ ", rhs))
}

mk_sunab_formula <- function(y, controls, refp = -1) {
  rhs <- paste0("sunab(cohort, time_id, ref.p = ", refp, ")",
                if (length(controls)) paste0(" + ", paste(controls, collapse = " + ")) else "")
  as.formula(paste0(y, " ~ ", rhs))
}

wcb_p <- function(m, param) {
  if (!isTRUE(wcb_ok)) return(NA_real_)
  err <- NULL
  for (type in c("rademacher","mammen")) {
    bt <- tryCatch(
      fwildclusterboot::boottest(m, param = param, clustid = ~ country,
                                 bootcluster = "country", B = 999,
                                 impose_null = TRUE, type = type),
      error = function(e) { err <- conditionMessage(e); NULL }
    )
    if (!is.null(bt)) return(as.numeric(bt$p_val))
  }
  if (VERBOSE) message(sprintf("boottest() failed for '%s': %s", param, err))
  NA_real_
}

print_snapshot <- function(y, twfe_row, es_info = NULL, pre = NULL) {
  cat("\n========================\n")
  cat(sprintf("Outcome: %s\n", y))
  if (!is.null(twfe_row)) {
    cat("-- TWFE (hetero SEs; wild-cluster p by country) --\n")
    print(data.table(Estimate = twfe_row$coef, `Std. Error` = twfe_row$se,
                     `t value` = twfe_row$t, `Pr(>|t|)` = twfe_row$p_naive,
                     `p_wild` = twfe_row$p_wcb, n = twfe_row$n))
  }
  if (!is.null(es_info)) {
    if (!is.null(es_info$path) && nrow(es_info$path)) {
      cat("\n-- Event-study (tau in [-12,24]) --\n")
      print(head(es_info$path[tau >= -12 & tau <= 24, .(tau, beta, se, p_wcb)], 12))
      if (nrow(es_info$path) > 12) cat("... (full path in tables4/did_eventstudy_paths.csv)\n")
    } else {
      cat("\n-- Event-study --\n(no ES terms estimated — fell back to GAP-ITS if available)\n")
    }
  }
  if (!is.null(pre)) {
    cat("\n-- Pretrend joint test (all leads = 0) --\n")
    print(pre)
  }
  cat("\n")
}

# ---------- ES try: sunab with country FE only ----------
fit_es_sunab_countryFE <- function(dd, y, controls, refp = -1) {
  f_es  <- mk_sunab_formula(y, controls, refp)
  m_es  <- feols(f_es, data = dd, fixef = "country", warn = FALSE, notes = FALSE)
  cn    <- names(coef(m_es))
  keep  <- grepl("sunab\\(", cn) & grepl("::t=\\-?\\d+", cn)
  if (!any(keep)) return(NULL)
  es    <- data.table(term = cn[keep], beta = coef(m_es)[keep])
  es[, tau := as.integer(sub(".*::t=\\s*([-]?[0-9]+).*", "\\1", term))]
  V     <- vcov(m_es, vcov = "hetero")
  es[, se := sqrt(diag(V)[match(term, rownames(V))])]
  for (h in c(0,12)) if (h %in% es$tau) {
    par <- es[tau == h, term][1]
    es[tau == h, p_wcb := wcb_p(m_es, par)]
  }
  setorder(es, tau)
  list(model = m_es, path = es)
}

# ---------- GAP-ITS fallback ----------
# Build CZ–SK gap and run ITS with manual event time
fit_gap_its <- function(d, y, controls, treat_date, nw_lag = 12L) {
  # wide by country for y and controls
  vals <- unique(c(y, controls))
  w <- dcast(d, date + time_id ~ country, value.var = vals)
  # assume exactly "CZ" and "SK" columns exist for each var
  y_gap <- w[[paste0(y, "_CZ")]] - w[[paste0(y, "_SK")]]
  X <- NULL
  if (length(controls)) {
    X <- lapply(controls, function(cv) w[[paste0(cv, "_CZ")]] - w[[paste0(cv, "_SK")]])
    X <- as.data.table(setNames(X, paste0(controls, "_gap")))
  } else {
    X <- data.table()
  }
  # event time
  t0 <- unique(d[date == treat_date, time_id])
  if (length(t0) != 1L) stop("treatment date not found in time_id grid")
  et <- w$time_id - t0
  dt <- data.table(y_gap = y_gap, et = et, time_id = w$time_id, date = w$date, X)
  dt <- dt[!is.na(y_gap)]
  # ITS with manual event dummies; ref = -1
  f_rhs <- paste0("i(et, ref = -1)", if (ncol(X)) paste0(" + ", paste(names(X), collapse = " + ")) else "")
  fml <- as.formula(paste0("y_gap ~ ", f_rhs))
  m <- feols(fml, data = dt, warn = FALSE, notes = FALSE)
  s <- summary(m, vcov = "NW", nw = nw_lag)
  # pull ES-like path
  cn <- names(coef(m))
  keep <- grepl("^i\\(et, ref = -1\\)::", cn)
  es <- data.table(term = cn[keep], beta = coef(m)[keep])
  es[, tau := as.integer(sub("^i\\(et, ref = -1\\)::\\s*([-]?[0-9]+)$", "\\1", term))]
  V <- vcov(m, vcov = "NW", nw = nw_lag)
  es[, se := sqrt(diag(V)[match(term, rownames(V))])]
  setorder(es, tau)
  # joint pretrend (all tau<0)
  lead_idx <- which(es$tau < 0)
  pre <- if (length(lead_idx)) {
    LHS <- which(keep)[es$tau < 0]
    W <- wald(m, LHS ~ 0, vcov = "NW", nw = nw_lag)
    data.table(K = length(lead_idx), stat = unname(W$stat), p = unname(W$p), df = paste0(W$df1, ",", W$df2), n = nobs(m))
  } else data.table(K = 0, stat = NA_real_, p = NA_real_, df = NA, n = nobs(m))
  list(model = m, path = es, pre = pre, data = dt)
}

# ---------- main per-outcome block ----------
fit_block <- function(dat, y, controls, treat_date) {
  dd <- copy(dat)
  dd[, post := fifelse(date >= treat_date, 1L, 0L)]
  dd[, D    := treated * post]
  t0 <- unique(dd[date == treat_date, time_id])
  if (length(t0) != 1L) stop("treatment date not found in time_id grid")
  dd[, cohort := fifelse(treated == 1L, t0, 0L)]
  dd <- balance_by_time(na.omit(dd[, c("country","time_id","date","treated","post","D","cohort", y, controls), with = FALSE]))
  if (!nrow(dd) || uniqueN(dd$D) == 1L) {
    return(list(twfe_row = NULL, es_info = NULL, es_pre = NULL, diag = data.table(outcome = y, n = nrow(dd), ok = FALSE)))
  }
  # TWFE snapshot
  f_twfe <- mk_twfe_formula(y, controls)
  m_twfe <- feols(f_twfe, data = dd, fixef = c("country","time_id"), warn = FALSE, notes = FALSE)
  s_twfe <- summary(m_twfe, vcov = "hetero")
  d_row <- tryCatch(s_twfe$coeftable["D", , drop = FALSE], error = function(e) NULL)
  twfe_row <- if (is.null(d_row)) NULL else data.table(
    outcome = y, coef = unname(d_row[1]), se = unname(d_row[2]), t = unname(d_row[3]),
    p_naive = unname(d_row[4]), p_wcb = wcb_p(m_twfe, "D"), n = nobs(m_twfe), spec = "twfe"
  )
  # Try ES via sunab + country FE
  es_info <- tryCatch(fit_es_sunab_countryFE(dd, y, controls, refp = -1), error = function(e) NULL)
  if (is.null(es_info) || !nrow(es_info$path)) {
    # Fallback: GAP–ITS
    gap <- fit_gap_its(dd, y, controls, treat_date = treat_date, nw_lag = 12L)
    es_info <- list(model = gap$model, path = gap$path)
    es_pre  <- gap$pre
  } else {
    # Pretrend under sunab (cluster-robust HC)
    cn <- names(coef(es_info$model))
    lead_idx <- which(grepl("::t=\\-", cn))
    es_pre <- if (length(lead_idx)) {
      W <- wald(es_info$model, lead_idx ~ 0, vcov = "hetero")
      data.table(K = length(lead_idx), stat = unname(W$stat), p = unname(W$p), df = paste0(W$df1,",",W$df2), n = nobs(es_info$model))
    } else data.table(K = 0, stat = NA_real_, p = NA_real_, df = NA, n = nobs(es_info$model))
  }
  list(twfe_row = twfe_row, es_info = es_info, es_pre = es_pre, diag = data.table(outcome = y, n = nobs(m_twfe), ok = TRUE))
}

# -------------------- load & identifiers --------------------
dat <- fread(DATA_PATH)
if (!("date" %in% names(dat))) stop("expected 'date' column in data")
if (!inherits(dat$date, "IDate")) dat[, date := as.IDate(date)]
cat("\n-- columns present --\n"); print(names(dat))

if (!("time_id" %in% names(dat))) dat[, time_id := as.integer(factor(date))]
dat[, treated := fifelse(country == "CZ", 1L, 0L)]
dat[, post    := fifelse(date >= TREAT_DATE, 1L, 0L)]
dat[, D       := treated * post]

san <- dat[, .(n = .N, post_ones = sum(post == 1L), D_ones = sum(D == 1L)), by = country]
cat("\n-- treat/post sanity --\n"); print(san)

# -------------------- missingness policy --------------------
aud_pre <- miss_audit(dat, CAND_CONTROLS); cat("\n-- missingness audit (pre) --\n"); print(aud_pre)
aud_post <- miss_audit(dat, CAND_CONTROLS); cat("\n-- missingness audit (post) --\n"); print(aud_post)
bad <- aud_post[pNA_all > 0.40 | pNA_treatpost > 0.05 | pNA_treatpre > 0.05 | pNA_control > 0.05, var]
good_controls <- setdiff(CAND_CONTROLS, bad)
cat("\n-- controls dropped (NA policy): ", ifelse(length(bad), paste(bad, collapse = ", "), "(none)"),
    "\n-- controls kept: ", ifelse(length(good_controls), paste(good_controls, collapse = ", "), "(none)"), "\n", sep = "")

# -------------------- run per-outcome --------------------
outcomes_present <- OUTCOMES[OUTCOMES %in% names(dat)]
cat("\noutcomes used: ", paste(outcomes_present, collapse = ", "), "\n", sep = "")
cat("controls (candidates): ", paste(CAND_CONTROLS, collapse = ", "), "\n\n", sep = "")

rows_twfe <- list(); rows_es_pre <- list(); paths_es <- list(); diags <- list()

for (y in outcomes_present) {
  kept <- select_controls(dat, y, good_controls, min_n = 100L)
  dropped <- setdiff(good_controls, kept)
  cat(sprintf("[diagnostics: %s ] kept: %s | dropped: %s\n",
              y, ifelse(length(kept), paste(kept, collapse = ", "), "(none)"),
              ifelse(length(dropped), paste(dropped, collapse = ", "), "(none)")))

  res <- fit_block(dat, y, kept, TREAT_DATE)

  if (!is.null(res$twfe_row)) rows_twfe[[y]] <- res$twfe_row
  if (!is.null(res$es_pre))    rows_es_pre[[y]] <- cbind(outcome = y, res$es_pre)
  if (!is.null(res$es_info$path)) paths_es[[y]] <- cbind(outcome = y, res$es_info$path)
  if (!is.null(res$diag))      diags[[y]] <- res$diag

  print_snapshot(y, res$twfe_row, res$es_info, if (!is.null(res$es_pre)) res$es_pre else NULL)
}

main_twfe <- rbindlist(rows_twfe, use.names = TRUE, fill = TRUE)
es_pre    <- rbindlist(rows_es_pre, use.names = TRUE, fill = TRUE)
es_paths  <- rbindlist(paths_es, use.names = TRUE, fill = TRUE)
diag_tbl  <- rbindlist(diags, use.names = TRUE, fill = TRUE)

# -------------------- write outputs --------------------
fwrite(dat[1:10],            file.path(OUT_DIR, "data_head.csv"))
fwrite(main_twfe,            file.path(OUT_DIR, "did_main_twfe.csv"))
fwrite(es_pre,               file.path(OUT_DIR, "did_eventstudy_pretrend_tests.csv"))
if (exists("es_paths") && length(es_paths)) {
  fwrite(es_paths,           file.path(OUT_DIR, "did_eventstudy_paths.csv"))
} else {
  if (file.exists(file.path(OUT_DIR, "did_eventstudy_paths.csv"))) file.remove(file.path(OUT_DIR, "did_eventstudy_paths.csv"))
}
fwrite(diag_tbl,             file.path(OUT_DIR, "diagnostics.csv"))

cat("\nDone. Outputs in ", OUT_DIR, ":\n",
    "- data_head.csv\n",
    "- did_main_twfe.csv (coef/se/t/p_naive + p_wcb for D)\n",
    "- did_eventstudy_pretrend_tests.csv (joint pretrend)\n",
    "- did_eventstudy_paths.csv (if ES estimated; else removed)\n",
    "- diagnostics.csv\n", sep = "")
