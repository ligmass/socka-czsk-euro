#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
})

set.seed(42)

DATA_PATH <- Sys.getenv("DID_DATA", "data/monthly_panel_clean.csv")
OUT_DIR   <- file.path("result tables baseline", "ukraine", "ukraine_quad"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTCOMES  <- c("hicp_yoy","unemp_rate","hicp","imports_world_meur","exports_world_meur","log_imp","log_exp")
TREAT_DATE <- as.IDate("2022-02-01")

TAU_MIN <- -12L; TAU_MAX <- 18L
BIN_EDGES <- list(pre12_7 = c(-12L,-7L), pre6_1 = c(-6L,-1L), t0 = c(0L,0L), post1_6 = c(1L,6L), post7_12= c(7L,12L), post13_18=c(13L,18L))
REF_BIN <- "pre6_1"

auto_lag <- function(T){ L <- floor(4 * (T/100)^(2/9)); max(0L, min(L, 12L)) }
bin_tau <- function(tau){ out <- rep(NA_character_, length(tau)); for (nm in names(BIN_EDGES)){ lo <- BIN_EDGES[[nm]][1]; hi <- BIN_EDGES[[nm]][2]; out[tau >= lo & tau <= hi] <- nm }; out }
vcov_failsoft <- function(fit, L){ V <- tryCatch(NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE), error = function(e) e); if (inherits(V, "error") || any(!is.finite(V))) V <- tryCatch(vcov(fit), error = function(e) e); if (inherits(V, "error") || any(!is.finite(V))) { X <- model.matrix(fit); V <- matrix(NA_real_, ncol(X), ncol(X)) }; V }

wald_zero <- function(fit, V, which_coefs){ if (length(which_coefs) == 0) return(list(K=0, stat=NA_real_, p=NA_real_)); b <- coef(fit)[which_coefs]; Vsub <- V[which_coefs, which_coefs, drop = FALSE]; stat <- tryCatch(as.numeric(t(b) %*% solve(Vsub, b)), error = function(e) NA_real_); pval <- if (is.finite(stat)) pchisq(stat, df = length(which_coefs), lower.tail = FALSE) else NA_real_; list(K = length(which_coefs), stat = stat, p = pval) }

print_block <- function(lbl, dt, L, Tn, gof_stats = NULL, pre_W = NULL, post_W = NULL){ cat("\n========================\n"); cat(sprintf("%s\n", lbl)); cat(sprintf("Newey–West lag L = %d (T=%d)\n", L, Tn)); if (!is.null(gof_stats)) cat(sprintf("R-squared: %.4f, Adj. R-squared: %.4f\n", gof_stats$r2, gof_stats$ar2)); print(dt); if (!is.null(pre_W)) cat(sprintf("\n-- Joint pretrend --\nK=%d, chi2=%.4f, p=%g\n", pre_W$K, pre_W$stat, pre_W$p)); if (!is.null(post_W)) cat(sprintf("\n-- Joint post period --\nK=%d, chi2=%.4f, p=%g\n", post_W$K, post_W$stat, post_W$p)) }

run_one_series <- function(dat, y, break_date){
  need <- c("date","country", y); if (!all(need %in% names(dat))) return(NULL)
  d0 <- dat[, ..need]; d0 <- dcast(d0, date ~ country, value.var = y); setnames(d0, c("CZ","SK"), c("cz","sk"), skip_absent = TRUE); d0[, diff := cz - sk]
  ## Controls removed for baseline: no policy_diff or fx_vol
  
  d0 <- d0
  d0 <- d0[order(date)]; d0[, time_id := as.integer(factor(date, levels = sort(unique(date))))]; t0 <- unique(d0[date == break_date, time_id]); if (length(t0) != 1L) return(NULL); d0[, tau := time_id - t0]; d0[, trend := time_id]
  d0 <- d0[tau >= TAU_MIN & tau <= TAU_MAX]; if (!nrow(d0)) return(NULL); d0[, bin := bin_tau(tau)]
  keep_cols <- c("diff", "bin", "trend"); d0 <- d0[!is.na(bin)]; d0 <- na.omit(d0, cols = keep_cols)
  if (!(REF_BIN %in% d0$bin)) return(NULL)
  d0[, bin := factor(bin, levels = unique(c(REF_BIN, setdiff(names(BIN_EDGES), REF_BIN))))]; d0[, bin := stats::relevel(bin, ref = REF_BIN)]
  if (nlevels(d0$bin) < 2L) return(NULL)

  # MODIFIED: Added I(trend^2)
  fit <- lm(diff ~ bin + trend + I(trend^2), data = d0)
  
  s <- summary(fit); gof <- data.table(outcome = y, n = nobs(fit), r2 = s$r.squared, ar2 = s$adj.r.squared)
  Tn <- nrow(d0); L  <- auto_lag(Tn); V  <- vcov_failsoft(fit, L)
  cf <- coef(fit); cn <- names(cf); keep <- grep("^bin", cn, value = TRUE)
  if (!length(keep)) return(NULL)
  se <- sqrt(diag(V))[keep]; beta <- cf[keep]; z <- beta / se; p  <- 2 * pnorm(-abs(z))
  rep_tau <- function(nm){ rng <- BIN_EDGES[[sub("^bin","",nm)]]; if (is.null(rng)) rng <- BIN_EDGES[[sub("^bin","", sub("^bin","",nm))]]; as.integer(floor(mean(rng))) }
  taus <- vapply(sub("^bin","", keep), function(s) { nm <- sub("^bin", "", s); if (nzchar(nm)) rep_tau(nm) else rep_tau(s) }, FUN.VALUE = integer(1))
  out <- data.table(outcome = y, bin = keep, tau  = taus, beta = as.numeric(beta), se = as.numeric(se), p = as.numeric(p))[order(tau)]
  pre_bins <- names(BIN_EDGES)[vapply(BIN_EDGES, function(r) r[2] < 0L, logical(1))]; pre_bins <- setdiff(pre_bins, REF_BIN); pre_keep <- paste0("bin", pre_bins); pre_ix <- which(cn %in% pre_keep); pre_W <- if (length(pre_ix)) wald_zero(fit, V, pre_ix) else list(K=0, stat=NA_real_, p=NA_real_)
  post_bins <- names(BIN_EDGES)[vapply(BIN_EDGES, function(r) r[1] >= 0L, logical(1))]; post_keep <- paste0("bin", post_bins); post_ix <- which(cn %in% post_keep); post_W <- if (length(post_ix)) wald_zero(fit, V, post_ix) else list(K=0, stat=NA_real_, p=NA_real_)
  lbl <- sprintf("Outcome: %s (diff = CZ − SK) [UKRAINE QUAD %s]", y, as.character(break_date))
  print_block(lbl, out, L, Tn, gof, pre_W, post_W)
  list(coefs = out, gof = gof)
}

dat <- fread(DATA_PATH); if (!inherits(dat$date, "IDate")) dat[, date := as.IDate(date)]
es_list <- list(); gof_list <- list(); cat("\n--- STARTING UKRAINE QUAD TEST ---\n")
for (y in OUTCOMES[OUTCOMES %in% names(dat)]) { res <- run_one_series(dat, y, TREAT_DATE); if (!is.null(res) && nrow(res$coefs) > 0) { es_list[[y]] <- res$coefs; gof_list[[y]] <- res$gof } }
es_main <- if (length(es_list)) rbindlist(es_list, use.names = TRUE) else data.table()
gof_main <- if (length(gof_list)) rbindlist(gof_list, use.names = TRUE) else data.table()
if (nrow(gof_main) && nrow(es_main)) { es_main <- merge(es_main, gof_main, by = "outcome", all.x = TRUE) }
fwrite(es_main, file.path(OUT_DIR, "es_ukraine_quad.csv")); cat("\nDone. Outputs in ", OUT_DIR, ":\n- es_ukraine_quad.csv (coefficients with GOF columns)\n", sep = "")