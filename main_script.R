library(quantmod)
library(fredr)
library(xts)
library(bvarsv)
library(vars)
library(coda)

# Set API key
fredr_set_key("ed46cc40a9dc36b9a59844502490a2d3")

start_date <- as.Date("2010-07-18")
end_date   <- as.Date("2025-05-27")


cat("Step 1: Downloading data from Yahoo Finance and FRED...\n")
btc_ohlc   <- getSymbols("BTC-USD", src="yahoo", from=start_date, to=end_date, auto.assign=FALSE)
gold_ohlc  <- getSymbols("GLD",     src="yahoo", from=start_date, to=end_date, auto.assign=FALSE)
sp500_ohlc <- getSymbols("^GSPC",   src="yahoo", from=start_date, to=end_date, auto.assign=FALSE)
vix_ohlc   <- getSymbols("^VIX",    src="yahoo", from=start_date, to=end_date, auto.assign=FALSE)
dxy_ohlc   <- getSymbols("DX-Y.NYB",src="yahoo", from=start_date, to=end_date, auto.assign=FALSE)

epu_raw <- fredr("USEPUINDXD", observation_start=start_date, observation_end=end_date)
epu     <- xts(epu_raw$value, order.by=as.Date(epu_raw$date))

# Fill missing values
btc_ohlc   <- na.locf(btc_ohlc)
gold_ohlc  <- na.locf(gold_ohlc)
sp500_ohlc <- na.locf(sp500_ohlc)
vix_ohlc   <- na.locf(vix_ohlc)
dxy_ohlc   <- na.locf(dxy_ohlc)
epu        <- na.locf(epu)

# Fix indexes explicitly
index(btc_ohlc)   <- as.Date(index(btc_ohlc))
index(gold_ohlc)  <- as.Date(index(gold_ohlc))
index(sp500_ohlc) <- as.Date(index(sp500_ohlc))
index(vix_ohlc)   <- as.Date(index(vix_ohlc))
index(dxy_ohlc)   <- as.Date(index(dxy_ohlc))
index(epu)        <- as.Date(index(epu))

#STEPS 2+3 COMBINED: daily merge → weekly sample 

# Merge all daily series pairwise
daily_all <- Cl(btc_ohlc)
daily_all <- merge(daily_all, Cl(gold_ohlc),  all = TRUE)
daily_all <- merge(daily_all, Cl(sp500_ohlc), all = TRUE)
daily_all <- merge(daily_all, Cl(vix_ohlc),   all = TRUE)
daily_all <- merge(daily_all, Cl(dxy_ohlc),   all = TRUE)
daily_all <- merge(daily_all, epu,             all = TRUE)
colnames(daily_all) <- c("BTC", "GOLD", "SP500", "VIX", "DXY", "EPU")

# Forward-fill then drop leading NAs
daily_all <- na.locf(daily_all, na.rm = TRUE)
daily_all <- na.omit(daily_all)

cat("Daily rows:", nrow(daily_all), "\n")
cat("Date range:", as.character(start(daily_all)), "to", as.character(end(daily_all)), "\n")

# Sample last observation of each week  all rows now share same weekday
ep           <- endpoints(daily_all, on = "weeks", k = 1)
merged_weekly <- period.apply(daily_all, ep, last)

cat("Weekly rows:", nrow(merged_weekly), "\n")
cat("Sample weekdays:", weekdays(head(index(merged_weekly), 6)), "\n")

# STEP 4: returns 
data_merged <- diff(log(merged_weekly))
data_merged <- na.omit(data_merged)

tvp_data <- as.matrix(data_merged)[, c("EPU", "VIX", "DXY", "SP500", "GOLD", "BTC")]

cat("Observations:", nrow(tvp_data), "\n")
cat("Zero proportions:\n")
print(apply(tvp_data, 2, function(x) mean(x == 0)))

# 5. LAG SELECTION

lag_select <- VARselect(tvp_data, lag.max=24, type="const")
print(lag_select$selection)

# STEP 6: ESTIMATE TVP-VAR
p_lag <- 2

tvp_data_scaled <- scale(tvp_data)

fit_tvp <- bvar.sv.tvp(
  Y      = tvp_data_scaled,
  p      = p_lag,
  tau    = 120,    # ~2.3 years of weekly data as training sample
  nrep   = 30000,
  nburn  = 15000
)

saveRDS(fit_tvp, "fit_tvp_6var.rds")
cat("Model saved.\n")

# Quick sanity check
cat("Beta.postmean dims:", dim(fit_tvp$Beta.postmean), "\n")
cat("H.postmean dims:",    dim(fit_tvp$H.postmean), "\n")
cat("NaN in Beta:", any(is.nan(fit_tvp$Beta.postmean)), "\n")
cat("NaN in H:",    any(is.nan(fit_tvp$H.postmean)), "\n")


# 7. DIMENSION CHECKS

cat("Beta.postmean:", dim(fit_tvp$Beta.postmean), "\n")  # [K, n, T]
cat("H.postmean:",    dim(fit_tvp$H.postmean), "\n")     # [T, n]

cat("NaN in Beta.postmean:", any(is.nan(fit_tvp$Beta.postmean)), "\n")
cat("NaN in H.postmean:",    any(is.nan(fit_tvp$H.postmean)), "\n")


# 8. DIAGNOSTICS — GEWEKE & INEFFICIENCY FACTOR

cat("\nStep 8: Calculating convergence diagnostics...\n")
n       <- ncol(tvp_data)
tau     <- 120
offset  <- tau + p_lag   # = 122

cat("Beta.draws dims:", dim(fit_tvp$Beta.draws), "\n")
K_total <- dim(fit_tvp$Beta.draws)[2]
T_obs   <- dim(fit_tvp$Beta.draws)[3]

ineff_factor <- function(x) length(x) / effectiveSize(as.mcmc(x))

geweke_beta <- numeric(K_total)
ineff_beta  <- numeric(K_total)

for (k in 1:K_total) {
  draws_k        <- fit_tvp$Beta.draws[1, k, ]
  geweke_beta[k] <- coda::geweke.diag(as.mcmc(draws_k))$z
  ineff_beta[k]  <- ineff_factor(draws_k)
}

cat("=== Beta Diagnostics (last t) ===\n")
print(data.frame(
  param        = paste0("Beta_", 1:K_total),
  geweke_z     = round(geweke_beta, 3),
  ineff_factor = round(ineff_beta, 3)
))


# 9. COMPUTE TIME-VARYING FEVD

cat("\nStep 9: Computing structural FEVD values...\n")
H_ahead     <- 24
T_total     <- dim(fit_tvp$Beta.postmean)[3]   # Total weeks estimated (436)
tau_offset  <- 120
p_lag       <- 2                               # Added to match your model specs


model_dates <- index(data_merged)[(tau_offset + p_lag + 1):(tau_offset + p_lag + T_total)]

fevd_at_t <- function(t, H) {
  beta_t <- fit_tvp$Beta.postmean[, , t]        # n x K = 6 x 13
  A      <- beta_t[, 2:(n*p_lag + 1)]           # n x (n*p) = 6 x 12
  
  companion <- matrix(0, n*p_lag, n*p_lag)
  companion[1:n, ] <- A
  if (p_lag > 1)
    companion[(n+1):(n*p_lag), 1:(n*(p_lag-1))] <- diag(n*(p_lag-1))
  
  sigma_structural <- fit_tvp$H.postmean[, , t] # Safe 3D extraction layout [M x M x T]
  sigma_structural <- (sigma_structural + t(sigma_structural)) / 2
  P <- t(chol(sigma_structural))
  
  J <- cbind(diag(n), matrix(0, n, n*(p_lag-1)))
  
  fevd_num    <- matrix(0, n, n)
  companion_h <- diag(n*p_lag)
  
  for (h in 0:(H-1)) {
    Phi_h       <- J %*% companion_h %*% t(J)
    Theta_h     <- Phi_h %*% P
    fevd_num    <- fevd_num + Theta_h^2
    companion_h <- companion_h %*% companion
  }
  
  fevd <- fevd_num / rowSums(fevd_num)
  rownames(fevd) <- colnames(tvp_data)
  colnames(fevd) <- paste0("Shock:", colnames(tvp_data))
  return(fevd)
}

# Target dates — Aligned with the post-estimation matrix range
target_dates <- as.Date(c("2017-06-15", "2018-01-05", "2020-03-15",
                          "2021-11-12", "2022-11-11", "2024-03-15"))

t_indices_model <- sapply(target_dates, function(d) {
  data_idx <- which.min(abs(index(data_merged) - d))
  model_t  <- data_idx - tau_offset
  max(1, min(model_t, T_total))
})

cat("Adjusted t_indices:", t_indices_model, "\n")

cat("Corresponding dates:", as.character(index(data_merged)[t_indices_model + tau_offset + p_lag]), "\n")

for (i in seq_along(target_dates)) {
  cat(sprintf("\n=== FEVD at %s (model t=%d) ===\n", target_dates[i], t_indices_model[i]))
  fevd <- fevd_at_t(t_indices_model[i], H_ahead)
  print(round(fevd[6, ] * 100, 2)) # 6 is BTC variable position
}

# Continuous FEVD timeline
fevd_time <- matrix(NA, nrow=T_total, ncol=n)
colnames(fevd_time) <- colnames(tvp_data)

for (t in 1:T_total) {
  fevd_t        <- fevd_at_t(t, H_ahead)
  fevd_time[t,] <- fevd_t[6, ]
}


fevd_xts <- xts(fevd_time, order.by=model_dates)
colors   <- c("darkgreen", "darkred", "purple", "blue", "darkgoldenrod", "black")

for (i in 1:n) {
  df <- data.frame(date=index(fevd_xts), value=as.numeric(fevd_xts[,i])*100)
  plot(df$date, df$value, type="l",
       main  = paste("Share of BTC FEVD Explained by", colnames(tvp_data)[i]),
       col   = colors[i], lwd=2, ylab="Percent (%)", xlab="")
  abline(h=0,                      col="black",    lty=2)
  abline(v=as.Date("2017-06-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2020-03-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2022-11-11"), col="darkgray", lty=3, lwd=1.5)
}

# DEFINITIONS 

n            <- ncol(tvp_data)
tau          <- 120
offset       <- tau + p_lag
T_total      <- dim(fit_tvp$Beta.postmean)[3]
model_dates  <- index(data_merged)[(offset + 1):(offset + T_total)]
shock_vars   <- c(1, 2, 3, 4, 5)
shock_labels <- c("EPU", "VIX", "DXY", "S&P500", "Gold")
H_irf        <- 12

irf_draws_at_t <- function(t, shock_var, response_var, H, n_draws = 1000) {
  total_draws <- dim(fit_tvp$Beta.draws)[3]   # 3000
  draw_idx    <- sample(1:total_draws, n_draws)
  irf_mat     <- matrix(NA, nrow=n_draws, ncol=H)
  
  for (d in seq_along(draw_idx)) {
    # Beta.draws is [78, 436, 3000] = [n*K, T, nrep]
    # Extract draw d at time t: vector of length 78
    beta_vec <- fit_tvp$Beta.draws[, t, draw_idx[d]]   # length 78
    
    # Reshape to [n x K] = [6 x 13]
    beta_d <- matrix(beta_vec, nrow=n, ncol=n*p_lag+1, byrow=FALSE)
    A_d    <- beta_d[, 2:(n*p_lag + 1)]                # n x (n*p) = 6 x 12
    
    companion <- matrix(0, n*p_lag, n*p_lag)
    companion[1:n, ] <- A_d
    if (p_lag > 1)
      companion[(n+1):(n*p_lag), 1:(n*(p_lag-1))] <- diag(n*(p_lag-1))
    
    sigma_d <- fit_tvp$H.postmean[, , t]
    sigma_d <- (sigma_d + t(sigma_d)) / 2
    sigma_d <- sigma_d + diag(1e-10, n)
    
    tryCatch({
      P <- t(chol(sigma_d))
      J <- cbind(diag(n), matrix(0, n, n*(p_lag-1)))
      companion_h <- diag(n*p_lag)
      for (h in 1:H) {
        Phi_h         <- J %*% companion_h %*% t(J)
        Theta_h       <- Phi_h %*% P
        irf_mat[d, h] <- Theta_h[response_var, shock_var]
        companion_h   <- companion_h %*% companion
      }
    }, error = function(e) NULL)
  }
  
  list(
    median = apply(irf_mat, 2, median,          na.rm=TRUE),
    low95  = apply(irf_mat, 2, quantile, 0.025, na.rm=TRUE),
    high95 = apply(irf_mat, 2, quantile, 0.975, na.rm=TRUE),
    low68  = apply(irf_mat, 2, quantile, 0.16,  na.rm=TRUE),
    high68 = apply(irf_mat, 2, quantile, 0.84,  na.rm=TRUE)
  )
}

plot_irf_bands <- function(irf_result, title, H) {
  ylim <- range(c(irf_result$low95, irf_result$high95), na.rm=TRUE)
  
  plot(1:H, irf_result$median, type="l", col="black", lwd=2,
       ylim=ylim, main=title, xlab="Weeks Ahead", ylab="Impulse Response")
  polygon(c(1:H, rev(1:H)),
          c(irf_result$low95, rev(irf_result$high95)),
          col=adjustcolor("steelblue", alpha.f=0.2), border=NA)
  polygon(c(1:H, rev(1:H)),
          c(irf_result$low68, rev(irf_result$high68)),
          col=adjustcolor("steelblue", alpha.f=0.4), border=NA)
  lines(1:H, irf_result$median, col="black", lwd=2)
  abline(h=0, col="red", lty=2)
  legend("topright",
         legend=c("Median", "68% CI", "95% CI"),
         col=c("black",
               adjustcolor("steelblue", 0.4),
               adjustcolor("steelblue", 0.2)),
         lwd=c(2, 8, 8), bty="n")
}


# 10. POINT-IN-TIME IRF WITH CONFIDENCE BANDS (NATIVE PACKAGE EXTRACTION)

# Explicit index generation using target dates that fall inside estimated timeline
abs_idx_2017 <- which.min(abs(index(data_merged) - as.Date("2017-06-15")))
abs_idx_2020 <- which.min(abs(index(data_merged) - as.Date("2020-03-15")))
abs_idx_2024 <- which.min(abs(index(data_merged) - as.Date("2024-03-15")))

t_idx_2017   <- abs_idx_2017 - tau_offset
t_idx_2020   <- abs_idx_2020 - tau_offset
t_idx_2024   <- abs_idx_2024 - tau_offset

era_list     <- list(
  list(t = t_idx_2017, label = "Retail Era (2017)"),
  list(t = t_idx_2020, label = "COVID Era (2020)"),
  list(t = t_idx_2024, label = "Maturity Era (2024)")
)

for (s in seq_along(shock_vars)) {
  for (era in era_list) {
    cat(sprintf("Computing IRF via package draws: %s shock → BTC, %s\n", shock_labels[s], era$label))
    
    irf_result <- irf_draws_at_t(
      t            = era$t,
      shock_var    = shock_vars[s],
      response_var = 6, 
      H            = H_irf
    )
    
    plot_irf_bands(
      irf_result,
      title = paste0("BTC response to ", shock_labels[s], " shock — ", era$label),
      H     = H_irf
    )
  }
}




# 11. CONTINUOUS IRF TIMELINES (with 68% bands)

cat("\nStep 11: Computing Continuous IRF Timelines with bands...\n")

shock_cols <- list(EPU="darkgreen", VIX="darkred", DXY="purple",
                   SP500="blue", GOLD="darkgoldenrod")

irf_timeline_med  <- matrix(NA, nrow=T_total, ncol=5)
irf_timeline_low  <- matrix(NA, nrow=T_total, ncol=5)
irf_timeline_high <- matrix(NA, nrow=T_total, ncol=5)
colnames(irf_timeline_med)  <- names(shock_cols)
colnames(irf_timeline_low)  <- names(shock_cols)
colnames(irf_timeline_high) <- names(shock_cols)

for (s in seq_along(shock_vars)) {
  cat(sprintf("Computing continuous IRF with bands: %s\n", shock_labels[s]))
  for (t in 1:T_total) {
    irf_t <- irf_draws_at_t(
      t            = t,
      shock_var    = shock_vars[s],
      response_var = 6,
      H            = 1,
      n_draws      = 200
    )
    irf_timeline_med[t, s]  <- irf_t$median[1]
    irf_timeline_low[t, s]  <- irf_t$low68[1]
    irf_timeline_high[t, s] <- irf_t$high68[1]
  }
}

cat("\nStep 11 Complete!\n")


# 12. CUMULATIVE IRFS

cat("\nStep 12: Computing Cumulative IRFs...\n")

cirf_draws_at_t <- function(t, shock_var, response_var, H, n_draws = 1000) {
  total_draws <- dim(fit_tvp$Beta.draws)[3]
  draw_idx    <- sample(1:total_draws, n_draws)
  cirf_mat    <- matrix(NA, nrow=n_draws, ncol=H)
  
  for (d in seq_along(draw_idx)) {
    beta_vec <- fit_tvp$Beta.draws[, t, draw_idx[d]]
    beta_d   <- matrix(beta_vec, nrow=n, ncol=n*p_lag+1, byrow=FALSE)
    A_d      <- beta_d[, 2:(n*p_lag + 1)]
    
    companion <- matrix(0, n*p_lag, n*p_lag)
    companion[1:n, ] <- A_d
    if (p_lag > 1)
      companion[(n+1):(n*p_lag), 1:(n*(p_lag-1))] <- diag(n*(p_lag-1))
    
    sigma_d <- fit_tvp$H.postmean[, , t]
    sigma_d <- (sigma_d + t(sigma_d)) / 2
    sigma_d <- sigma_d + diag(1e-10, n)
    
    tryCatch({
      P <- t(chol(sigma_d))
      J <- cbind(diag(n), matrix(0, n, n*(p_lag-1)))
      companion_h <- diag(n*p_lag)
      cumulative  <- 0
      for (h in 1:H) {
        Phi_h          <- J %*% companion_h %*% t(J)
        Theta_h        <- Phi_h %*% P
        cumulative     <- cumulative + Theta_h[response_var, shock_var]
        cirf_mat[d, h] <- cumulative
        companion_h    <- companion_h %*% companion
      }
    }, error = function(e) NULL)
  }
  
  list(
    median = apply(cirf_mat, 2, median,          na.rm=TRUE),
    low95  = apply(cirf_mat, 2, quantile, 0.025, na.rm=TRUE),
    high95 = apply(cirf_mat, 2, quantile, 0.975, na.rm=TRUE),
    low68  = apply(cirf_mat, 2, quantile, 0.16,  na.rm=TRUE),
    high68 = apply(cirf_mat, 2, quantile, 0.84,  na.rm=TRUE)
  )
}

plot_cirf_bands <- function(cirf_result, title, H) {
  ylim <- range(c(cirf_result$low95, cirf_result$high95), na.rm=TRUE)
  plot(1:H, cirf_result$median, type="l", col="black", lwd=2,
       ylim=ylim, main=title, xlab="Weeks Ahead", ylab="Cumulative Response")
  polygon(c(1:H, rev(1:H)),
          c(cirf_result$low95, rev(cirf_result$high95)),
          col=adjustcolor("firebrick", alpha.f=0.15), border=NA)
  polygon(c(1:H, rev(1:H)),
          c(cirf_result$low68, rev(cirf_result$high68)),
          col=adjustcolor("firebrick", alpha.f=0.35), border=NA)
  lines(1:H, cirf_result$median, col="black", lwd=2)
  abline(h=0, col="blue", lty=2)
  legend("topright",
         legend=c("Median", "68% CI", "95% CI"),
         col=c("black",
               adjustcolor("firebrick", 0.35),
               adjustcolor("firebrick", 0.15)),
         lwd=c(2, 8, 8), bty="n")
}

# Point-in-time cumulative IRFs across eras
for (s in seq_along(shock_vars)) {
  for (era in era_list) {
    cat(sprintf("Computing CIRF: %s shock → BTC, %s\n", shock_labels[s], era$label))
    cirf_result <- cirf_draws_at_t(
      t            = era$t,
      shock_var    = shock_vars[s],
      response_var = 6,
      H            = H_irf
    )
    plot_cirf_bands(
      cirf_result,
      title = paste0("Cumulative IRF: BTC Response to ", shock_labels[s], " — ", era$label),
      H     = H_irf
    )
  }
}

# Continuous cumulative IRF timeline (total effect at H=12)
cirf_total <- matrix(NA, nrow=T_total, ncol=5)
colnames(cirf_total) <- names(shock_cols)

for (s in seq_along(shock_vars)) {
  cat(sprintf("Computing continuous CIRF timeline: %s\n", shock_labels[s]))
  for (t in 1:T_total) {
    cirf_t           <- cirf_draws_at_t(t, shock_vars[s], 6, H=H_irf, n_draws=100)
    cirf_total[t, s] <- cirf_t$median[H_irf]
  }
}

cirf_xts <- xts(cirf_total, order.by=model_dates)

cat("\nStep 12 Complete!\n")


# SAVE ALL PLOTS TO PDF

pdf("BTC_TVP_VAR_Plots.pdf", width=10, height=7)

# FEVD plots
for (i in 1:n) {
  df <- data.frame(date=index(fevd_xts), value=as.numeric(fevd_xts[,i])*100)
  plot(df$date, df$value, type="l",
       main  = paste("Share of BTC FEVD Explained by", colnames(tvp_data)[i]),
       col   = colors[i], lwd=2, ylab="Percent (%)", xlab="")
  abline(h=0,                     col="black",    lty=2)
  abline(v=as.Date("2017-06-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2020-03-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2022-11-11"), col="darkgray", lty=3, lwd=1.5)
}

# Point-in-time IRF plots
for (s in seq_along(shock_vars)) {
  for (era in era_list) {
    irf_result <- irf_draws_at_t(
      t            = era$t,
      shock_var    = shock_vars[s],
      response_var = 6,
      H            = H_irf
    )
    plot_irf_bands(
      irf_result,
      title = paste0("Point-in-Time IRF: BTC Response to ", shock_labels[s], " — ", era$label),
      H     = H_irf
    )
  }
}

# Continuous IRF plots with 68% bands
for (s in seq_along(shock_vars)) {
  var_name <- names(shock_cols)[s]
  dates    <- model_dates
  med      <- irf_timeline_med[, s]
  lo       <- irf_timeline_low[, s]
  hi       <- irf_timeline_high[, s]
  ylim     <- range(c(lo, hi), na.rm=TRUE)
  
  plot(dates, med, type="l", col=shock_cols[[s]], lwd=2,
       ylim=ylim,
       main = paste("Continuous IRF: BTC Response to", var_name, "Shock (h=1)"),
       ylab = "Impulse Response", xlab="")
  polygon(c(dates, rev(dates)),
          c(lo, rev(hi)),
          col=adjustcolor(shock_cols[[s]], alpha.f=0.2), border=NA)
  lines(dates, med, col=shock_cols[[s]], lwd=2)
  abline(h=0,                     col="black",    lty=2)
  abline(v=as.Date("2017-06-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2020-03-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2022-11-11"), col="darkgray", lty=3, lwd=1.5)
  legend("topright", legend=c("Median", "68% CI"),
         col=c(shock_cols[[s]], adjustcolor(shock_cols[[s]], 0.2)),
         lwd=c(2, 8), bty="n")
}

# Point-in-time cumulative IRF plots
for (s in seq_along(shock_vars)) {
  for (era in era_list) {
    cirf_result <- cirf_draws_at_t(
      t            = era$t,
      shock_var    = shock_vars[s],
      response_var = 6,
      H            = H_irf
    )
    plot_cirf_bands(
      cirf_result,
      title = paste0("Cumulative IRF: BTC Response to ", shock_labels[s], " — ", era$label),
      H     = H_irf
    )
  }
}

# Continuous cumulative IRF timeline plots
for (s in seq_along(shock_vars)) {
  var_name <- names(shock_cols)[s]
  df <- data.frame(date=index(cirf_xts), value=as.numeric(cirf_xts[, s]))
  plot(df$date, df$value, type="l",
       main = paste("Continuous Cumulative IRF: BTC Response to", var_name, "Shock (H=12)"),
       col  = shock_cols[[s]], lwd=2, ylab="Cumulative Response", xlab="")
  abline(h=0,                     col="black",    lty=2)
  abline(v=as.Date("2017-06-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2020-03-15"), col="darkgray", lty=3, lwd=1.5)
  abline(v=as.Date("2022-11-11"), col="darkgray", lty=3, lwd=1.5)
}

dev.off()
cat("All plots saved to BTC_TVP_VAR_Plots.pdf\n")