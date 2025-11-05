#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)   # Newey–West / HAC
  library(lmtest)     # coeftest / waldtest
})

set.seed(42)

# ---------- config ----------
VERBOSE   <- identical(Sys.getenv("VERBOSE","1"), "1")
DATA_PATH <- Sys.getenv("DID_DATA", "data/monthly_panel_clean.csv")
OUT_DIR   <- "tables5"; dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTCOMES  <- c("hicp_yoy","unemp_rate","hicp","imports_world_meur","exports_world_meur","log_imp","log_exp")
TREAT_DATE <- as.IDate("2009-01-01")   # main break
TAU_MIN <- -12L; TAU_MAX <- 18L        # ES window
# bins: pool months to get residual dof
BIN_EDGES <- list(
  pre12_7 = c(-12L,-7L),
  pre6_1  = c(-6L,-1L),   # reference
  t0      = c(0L,0L),
  post1_6 = c(1L,6L),
  post7_12= c(7L,12L),
  post13_18=c(13L,18L)
)
REF_BIN <- "pre6_1"

# ---------- helpers ----------
quiet <- function(x) suppressWarnings(suppressMessages(x))

auto_lag <- function(T){
  # Andrews–Monahan-ish small-sample rule; keep >=0
  L <- floor(4 * (T/100)^(2/9))
  max(0L, min(L, 12L))
}

bin_tau <- function(tau){
  # returns a character vector of bin names or NA if tau outside all bins
  out <- rep(NA_character_, length(tau))
  for (nm in names(BIN_EDGES)){
    lo <- BIN_EDGES[[nm]][1]; hi <- BIN_EDGES[[nm]][2]
    out[tau >= lo & tau <= hi] <- nm
  }
  out
}

vcov_failsoft <- function(fit, L){
  # try HAC; fall back to OLS; last resort add tiny ridge
  V <- tryCatch(NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE),
                error = function(e) e)
  if (inherits(V, "error") || any(!is.finite(V))){
    V <- tryCatch(vcov(fit),
                  error = function(e) e)
  }
  if (inherits(V, "error") || any(!is.finite(V))){
    # ridge the design slightly
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

print_block <- function(lbl, dt, L, Tn, pre_W = NULL, post_W = NULL){
  cat("\n========================\n")
  cat(sprintf("%s\n", lbl))
  cat(sprintf("Newey–West lag L = %d (T=%d)\n", L, Tn))
  print(dt)
  if (!is.null(pre_W)){
    cat("\n-- Joint pretrend (bins wholly <0) --\n")
    cat(sprintf("K=%d, chi2=%.4f, p=%g\n", pre_W$K, pre_W$stat, pre_W$p))
  }
  if (!is.null(post_W)){
    cat("\n-- Joint post period (bins >=0) --\n")
    cat(sprintf("K=%d, chi2=%.4f, p=%g\n", post_W$K, post_W$stat, post_W$p))
  }
}

wald_zero <- function(fit, V, which_coefs){
  if (length(which_coefs) == 0) return(list(K=0, stat=NA_real_, p=NA_real_))
  b <- coef(fit)[which_coefs]
  M <- diag(length(which_coefs))
  Vsub <- V[which_coefs, which_coefs, drop = FALSE]
  stat <- tryCatch(as.numeric(t(b) %*% solve(Vsub, b)), error = function(e) NA_real_)
  pval <- if (is.finite(stat)) pchisq(stat, df = length(which_coefs), lower.tail = FALSE) else NA_real_
  list(K = length(which_coefs), stat = stat, p = pval)
}

run_one_series <- function(dat, y, break_date){
  # build diff series: CZ - SK
  need <- c("date","country", y)
  if (!all(need %in% names(dat))) return(NULL)
  d0 <- dat[, ..need]
  d0 <- dcast(d0, date ~ country, value.var = y)
  setnames(d0, c("CZ","SK"), c("cz","sk"), skip_absent = TRUE)
  d0[, diff := cz - sk]
  d0 <- d0[!is.na(diff)]

  # time_id & tau
  if (!("time_id" %in% names(dat))) {
    # create from sorted unique dates
    d0 <- d0[order(date)]
    d0[, time_id := as.integer(factor(date, levels = sort(unique(date))))]
  } else {
    # map provided time_id to this frame
    ti_map <- unique(dat[, .(date, time_id)])
    d0 <- ti_map[d0, on = "date"]
  }
  t0 <- unique(d0[date == break_date, time_id])
  if (length(t0) != 1L) return(NULL)
  d0[, tau := time_id - t0]

  # restrict window & drop out-of-bin tau
  d0 <- d0[tau >= TAU_MIN & tau <= TAU_MAX]
  if (!nrow(d0)) return(NULL)

  # binning
  d0[, bin := bin_tau(tau)]
  d0 <- d0[!is.na(bin)]
  # ensure the reference bin exists and has obs
  if (!(REF_BIN %in% d0$bin)) return(NULL)

  # relevel factor (ref = [-6,-1])
  d0[, bin := factor(bin, levels = unique(c(REF_BIN, setdiff(names(BIN_EDGES), REF_BIN))))]
  d0[, bin := stats::relevel(bin, ref = REF_BIN)]

  # if only ref level present (shouldn’t happen after above), skip
  if (nlevels(d0$bin) < 2L) return(NULL)

  # fit lm(diff ~ bin)
  fit <- lm(diff ~ bin, data = d0)

  # HAC
  Tn <- nrow(d0)
  L  <- auto_lag(Tn)
  V  <- vcov_failsoft(fit, L)

  # tidy coefs: we only report non-ref
  cf <- coef(fit)
  cn <- names(cf)
  keep <- grep("^bin", cn, value = TRUE)
  if (!length(keep)) return(NULL)

  se <- sqrt(diag(V))[keep]
  beta <- cf[keep]
  z <- beta / se
  p  <- 2 * pnorm(-abs(z))

  # map back to representative tau per bin (mid-point)
  rep_tau <- function(nm){
    rng <- BIN_EDGES[[sub("^bin","",nm)]]
    if (is.null(rng)) rng <- BIN_EDGES[[sub("^bin","", sub("^bin","",nm))]]
    as.integer(floor(mean(rng)))
  }
  taus <- vapply(sub("^bin","", keep), function(s) {
    nm <- sub("^bin", "", s)
    if (nzchar(nm)) rep_tau(nm) else rep_tau(s)
  }, FUN.VALUE = integer(1))

  out <- data.table(
    outcome = y,
    bin = keep,
    tau  = taus,
    beta = as.numeric(beta),
    se   = as.numeric(se),
    p    = as.numeric(p)
  )[order(tau)]

  # joint tests
  # bins wholly <0: all negative bins except the ref (pre12_7)
  pre_bins  <- names(BIN_EDGES)[vapply(BIN_EDGES, function(r) r[2] < 0L, logical(1))]
  pre_bins  <- setdiff(pre_bins, REF_BIN)
  pre_keep  <- paste0("bin", pre_bins)
  pre_ix    <- which(cn %in% pre_keep)
  pre_W <- if (length(pre_ix)) {
    tmp <- wald_zero(fit, V, pre_ix)
  } else list(K=0, stat=NA_real_, p=NA_real_)

  # bins >=0: t0, post1_6, post7_12, post13_18 (whatever exists)
  post_bins <- names(BIN_EDGES)[vapply(BIN_EDGES, function(r) r[1] >= 0L, logical(1))]
  post_keep <- paste0("bin", post_bins)
  post_ix   <- which(cn %in% post_keep)
  post_W <- if (length(post_ix)) {
    tmp <- wald_zero(fit, V, post_ix)
  } else list(K=0, stat=NA_real_, p=NA_real_)

  if (VERBOSE) {
    lbl <- sprintf("Outcome: %s (diff = CZ − SK) [main %s]", y, as.character(break_date))
    print_block(lbl, out, L, Tn, pre_W, post_W)
  }

  out
}

# ---------- load ----------
dat <- fread(DATA_PATH)
if (!("date" %in% names(dat))) stop("expected 'date' column")
if (!inherits(dat$date, "IDate")) dat[, date := as.IDate(date)]
cat("\n-- columns present --\n"); print(names(dat))

# ---------- run ----------
# ---------- run (DEBUG MODE) ----------
es_list <- list()
cat("\n--- STARTING DEBUG RUN ---\n")
for (y in OUTCOMES[OUTCOMES %in% names(dat)]) {
  cat(sprintf("\nProcessing outcome: %s\n", y))
  
  # Run *without* tryCatch to see the real error
  res <- run_one_series(dat, y, TREAT_DATE) 
  
  if (is.null(res)) {
    cat(sprintf("Function returned NULL for %s.\n", y))
  } else if (!nrow(res)) {
    cat(sprintf("Function returned an empty table for %s.\n", y))
  } else {
    cat(sprintf("Successfully processed %s, found %d coefficients.\n", y, nrow(res)))
    es_list[[y]] <- res
  }
}
cat("\n--- DEBUG RUN COMPLETE ---\n")

es_main <- if (length(es_list)) rbindlist(es_list, use.names = TRUE) else data.table()
if (!nrow(es_main)) {
  warning("No ES coefficients were produced. Consider widening bins or window.")
}

# ---------- write ----------
fwrite(es_main, file.path(OUT_DIR, "es_main.csv"))
cat("\nDone. Outputs in ", OUT_DIR, ":\n- es_main.csv (binned ES with HAC SE)\n", sep = "")
