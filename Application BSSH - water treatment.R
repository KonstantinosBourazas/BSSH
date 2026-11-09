# R code for the water treatment plant application of BSSH
# Reproduces the monitoring results, post-alarm diagnosis,
# and the figure reported in the Supplementary Material.
# The file water-treatment.data (UCI Machine Learning Repository)
# must be in the working directory.

# If needed, install the required packages:
# install.packages("mice")
# install.packages("bestNormalize")

library(mice)
library(bestNormalize)

#------------------------------------------------------------------------------
# DATA
#------------------------------------------------------------------------------

var_labels <- c(
  "Q-E", "ZN-E", "PH-E", "DBO-E", "DQO-E", "SS-E", "SSV-E", "SED-E", "COND-E",
  "PH-P", "DBO-P", "SS-P", "SSV-P", "SED-P", "COND-P",
  "PH-D", "DBO-D", "DQO-D", "SS-D", "SSV-D", "SED-D", "COND-D",
  "PH-S", "DBO-S", "DQO-S", "SS-S", "SSV-S", "SED-S", "COND-S",
  "RD-DBO-P", "RD-SS-P", "RD-SED-P", "RD-DBO-S", "RD-DQO-S",
  "RD-DBO-G", "RD-DQO-G", "RD-SS-G", "RD-SED-G"
)

raw <- read.csv("water-treatment.data", header = FALSE, na.strings = "?")
colnames(raw) <- c("date", make.names(var_labels))
raw$date <- as.Date(sub("D-", "", raw$date), format = "%d/%m/%y")
raw <- raw[order(raw$date), ]

# Days with a documented malfunction
fault_dates <- as.Date(c("1990-03-13", "1990-03-14", "1990-03-15", "1990-04-29",
                         "1990-06-05", "1990-08-12", "1990-09-14", "1990-10-22",
                         "1991-05-24", "1991-05-28", "1991-05-31",
                         "1991-07-17", "1991-07-18", "1991-07-19"))
raw$fault <- raw$date %in% fault_dates

# Variables 30-38 are performance indices computed from the direct process measurements and are removed
d <- 29
vars <- make.names(var_labels[1:d])

# Historical data: 1990 without the fault-labelled days. Monitoring data: 1991.
year <- format(raw$date, "%Y")
hist_rows <- year == "1990" & !raw$fault
mon_rows <- year == "1991"

H_raw <- raw[hist_rows, vars]
M_raw <- raw[mon_rows, vars]
hist_dates <- raw$date[hist_rows]
mon_dates <- raw$date[mon_rows]
mon_fault <- raw$fault[mon_rows]

cat("\n===== DATA =====\n")
cat("Observations:", nrow(raw), "from", format(min(raw$date)), "to", format(max(raw$date)), "\n")
cat("Historical data (1990):", sum(year == "1990"), "days,",
    sum(year == "1990" & raw$fault), "fault-labelled days excluded, n0 =", nrow(H_raw), "\n")
cat("Monitoring data (1991):", nrow(M_raw), "days,", sum(mon_fault), "fault-labelled\n")
cat("Number of variables: d =", d, "\n")
cat(sprintf("Missing values: at most %.1f%% per variable and %.1f%% overall in the historical data, %.1f%% overall in the monitoring data\n",
            100 * max(colMeans(is.na(H_raw))), 100 * mean(is.na(H_raw)), 100 * mean(is.na(M_raw))))

#==============================================================================
# MISSING VALUES: PREDICTIVE MEAN MATCHING (mice)
#==============================================================================

# Each missing entry is replaced by the average of 20 imputations
impute_pmm <- function(data, where, ignore, maxit, seed) {
  meth <- rep("", ncol(data))
  names(meth) <- names(data)
  meth[colSums(where) > 0] <- "pmm"
  pred <- matrix(1, ncol(data), ncol(data), dimnames = list(names(data), names(data)))
  diag(pred) <- 0
  imp <- mice(data, m = 20, method = meth, predictorMatrix = pred, where = where,
              ignore = ignore, maxit = maxit, printFlag = FALSE, seed = seed)
  draws <- matrix(NA, sum(where), 20)
  for (m in 1:20) draws[, m] <- as.matrix(complete(imp, m))[where]
  out <- as.matrix(data)
  out[where] <- rowMeans(draws)
  as.data.frame(out)
}

seed <- 20260824

# Historical data
H <- impute_pmm(H_raw, is.na(H_raw), rep(FALSE, nrow(H_raw)), maxit = 20, seed = seed + 5000)

# Each incomplete monitoring day is imputed separately, with the imputation
# models fitted to the completed historical data only
M <- M_raw
for (i in which(!complete.cases(M_raw))) {
  data <- rbind(H, M_raw[i, ])
  n <- nrow(data)
  where <- matrix(FALSE, n, d)
  where[n, ] <- is.na(M_raw[i, ])
  ignore <- c(rep(FALSE, n - 1), TRUE)
  M[i, ] <- impute_pmm(data, where, ignore, maxit = 5, seed = seed + 100000 + i)[n, ]
}

#==============================================================================
# ORQ NORMALIZATION: FITTED TO THE HISTORICAL DATA, APPLIED TO THE MONITORING DATA
#==============================================================================

Y_hist <- matrix(NA, nrow(H), d)
Y_mon <- matrix(NA, nrow(M), d)
for (j in 1:d) {
  orq <- orderNorm(H[, j], warn = FALSE)
  Y_hist[, j] <- orq$x.t
  Y_mon[, j] <- predict(orq, newdata = M[, j], warn = FALSE)
}

# Mardia tests for multivariate skewness and kurtosis
n0 <- nrow(Y_hist)
Z <- scale(Y_hist, scale = FALSE)
G <- Z %*% solve(t(Z) %*% Z / n0) %*% t(Z)
b1 <- mean(G^3)
b2 <- mean(diag(G)^2)
p_skew <- pchisq(n0 * b1 / 6, d * (d + 1) * (d + 2) / 6, lower.tail = FALSE)
p_kurt <- 2 * pnorm(abs(b2 - d * (d + 2)) / sqrt(8 * d * (d + 2) / n0), lower.tail = FALSE)

# Correlation at a lag of one calendar day
lag1 <- which(diff(hist_dates) == 1)
rho1 <- rep(NA, d)
for (j in 1:d) rho1[j] <- cor(Y_hist[lag1, j], Y_hist[lag1 + 1, j])

cat("\n===== TRANSFORMED HISTORICAL DATA =====\n")
cat("Mardia test p-values: skewness", signif(p_skew, 3), ", kurtosis", signif(p_kurt, 3), "\n")
cat(sprintf("Absolute correlation at a lag of one calendar day: mean %.2f, maximum %.2f\n",
            mean(abs(rho1)), max(abs(rho1))))

#######################
#### PRIOR SETTING ####
#######################

# Independence-Jeffreys initial prior
M0n <- rep(0, d)
L0n <- 0
NU0n <- -1
PSIn <- matrix(0, d, d)

# Power prior weights giving a historical effective sample size a0 * n0 = 0, d, 2d, 4d
a0 <- c(0, d, 2 * d, 4 * d) / n0
names(a0) <- c("BSSH_0", "BSSH_d", "BSSH_2d", "BSSH_4d")

# FALSE ALARM RATE: ARL_0 = 370
far <- 1/370

#==============================================================================
# RECURSIVE NIW UPDATE
#==============================================================================

niw_update_one <- function(x, mu, lambda, nu, Psi) {
  x <- as.numeric(x)
  w <- lambda / (lambda + 1)
  diff <- x - mu
  mu_new  <- w * mu + (1 - w) * x
  Psi_new <- Psi + w * tcrossprod(diff)
  list(mu = mu_new, lambda = lambda + 1, nu = nu + 1, Psi = Psi_new)
}

#==============================================================================
# applicationBSSH (every alarm is excluded from the posterior update)
#==============================================================================

applicationBSSH <- function(X, Y, FAR, a0, m0i, l0i, nu0i, psi0i) {
  d <- ncol(X)
  N <- nrow(X)
  
  if (a0 != 0) {
    n0 <- nrow(Y)
    meanY <- colMeans(Y)
    Yc <- scale(Y, scale = FALSE)
    C0 <- t(Yc) %*% Yc
    
    m0 <- (l0i * m0i + a0 * n0 * meanY) / (l0i + a0 * n0)
    l0 <- l0i + a0 * n0
    nu0 <- nu0i + a0 * n0
    psi0 <- psi0i + a0 * C0 + l0i * m0i %*% t(m0i) +
      a0 * n0 * meanY %*% t(meanY) - l0 * m0 %*% t(m0)
  } else {
    m0 <- m0i; l0 <- l0i; nu0 <- nu0i; psi0 <- psi0i
  }
  
  start_n <- max(1, floor(d - nu0))
  
  mu <- m0; lambda <- l0; nu <- nu0; Psi <- psi0
  for (i in 1:start_n) {
    upd <- niw_update_one(X[i, ], mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
  }
  
  stats <- matrix(NA, N, 2)
  colnames(stats) <- c("t2", "qf")
  post <- vector("list", N)
  
  for (n in start_n:(N - 1)) {
    x_next <- X[n + 1, ]
    nu_star <- nu - d + 1
    
    A <- ((lambda + 1) * d) / (lambda * nu_star) * Psi
    diff <- matrix(x_next - mu, ncol = 1)
    
    t2 <- as.numeric(t(diff) %*% solve(A) %*% diff)
    q  <- qf(1 - FAR, d, nu_star)
    
    stats[n + 1, 1] <- t2
    stats[n + 1, 2] <- q
    post[[n + 1]] <- list(mu = mu, lambda = lambda, nu = nu, Psi = Psi)
    
    if (t2 <= q) {
      upd <- niw_update_one(x_next, mu, lambda, nu, Psi)
      mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
    }
  }
  
  list(stats = stats, post = post)
}

#==============================================================================
# POST-ALARM DIAGNOSIS
#==============================================================================

# p-value of the directed test of variable j given the variables in L
directed_test <- function(x, mu, A, nu_star, j, L) {
  d <- length(x)
  if (length(L) == 0) {
    t2 <- (x[j] - mu[j])^2 / A[j, j]
    U <- d * t2
  } else {
    ALL_inv <- solve(A[L, L, drop = FALSE])
    diffL <- x[L] - mu[L]
    mu_jgL <- mu[j] + A[j, L] %*% ALL_inv %*% diffL
    A_jgL <- A[j, j] - A[j, L] %*% ALL_inv %*% A[L, j]
    DeltaL <- d * (t(diffL) %*% ALL_inv %*% diffL)
    t2 <- (x[j] - mu_jgL)^2 / A_jgL
    U <- ((nu_star + length(L)) / (nu_star + DeltaL)) * d * t2
  }
  pf(as.numeric(U), 1, nu_star + length(L), lower.tail = FALSE)
}

# Classical MYT decomposition in the original variable order, alpha = 0.05 for each term
myt <- function(x, post, alpha = 0.05) {
  d <- length(x)
  nu_star <- post$nu - d + 1
  A <- ((post$lambda + 1) * d) / (post$lambda * nu_star) * post$Psi
  p <- rep(NA, d)
  for (j in 1:d) p[j] <- directed_test(x, post$mu, A, nu_star, j, seq_len(j - 1))
  flagged <- which(p < alpha)
  flagged[order(p[flagged])]
}

# Order-free diagnostic: all directed tests with up to kmax conditioning variables,
# BH adjustment at level q. A variable is flagged if at least one of its tests
# is significant.
order_free <- function(x, post, kmax = 2, q = 0.05) {
  d <- length(x)
  nu_star <- post$nu - d + 1
  A <- ((post$lambda + 1) * d) / (post$lambda * nu_star) * post$Psi
  
  M <- d * sum(choose(d - 1, 0:kmax))
  target <- rep(NA, M)
  p <- rep(NA, M)
  r <- 0
  for (j in 1:d) {
    r <- r + 1
    target[r] <- j
    p[r] <- directed_test(x, post$mu, A, nu_star, j, c())
    others <- setdiff(1:d, j)
    for (k in 1:kmax) {
      sets <- combn(others, k)
      for (s in 1:ncol(sets)) {
        r <- r + 1
        target[r] <- j
        p[r] <- directed_test(x, post$mu, A, nu_star, j, sets[, s])
      }
    }
  }
  
  p_adj <- p.adjust(p, method = "BH")
  min_p <- rep(NA, d)
  for (j in 1:d) min_p[j] <- min(p_adj[target == j])
  flagged <- which(min_p <= q)
  flagged[order(min_p[flagged])]
}

#==============================================================================
# RUN BSSH MONITORING
#==============================================================================

results <- list()
for (k in 1:4) results[[k]] <- applicationBSSH(Y_mon, Y_hist, far, a0[k], M0n, L0n, NU0n, PSIn)
names(results) <- names(a0)

alarms <- list()
monitor <- data.frame(Method = names(a0), Tested_days = NA, Total_alarms = NA,
                      Fault_days_detected = NA, Other_alarms = NA, Other_alarm_rate = NA)
for (k in 1:4) {
  st <- results[[k]]$stats
  tested <- !is.na(st[, 1])
  alarm <- tested & st[, 1] > st[, 2]
  alarms[[k]] <- which(alarm)
  monitor$Tested_days[k] <- sum(tested)
  monitor$Total_alarms[k] <- sum(alarm)
  monitor$Fault_days_detected[k] <- paste0(sum(alarm & mon_fault), "/", sum(tested & mon_fault))
  monitor$Other_alarms[k] <- sum(alarm & !mon_fault)
  monitor$Other_alarm_rate[k] <- round(100 * sum(alarm & !mon_fault) / sum(tested & !mon_fault), 1)
}
names(alarms) <- names(a0)

cat("\n===== MONITORING RESULTS =====\n")
print(monitor, row.names = FALSE)

for (k in 1:4) {
  st <- results[[k]]$stats
  al <- alarms[[k]]
  cat("\n===== ALARMS:", names(a0)[k], "=====\n")
  print(data.frame(
    Date = format(mon_dates[al]),
    Fault = ifelse(mon_fault[al], "YES", "NO"),
    Statistic = round(st[al, 1], 4),
    Threshold = round(st[al, 2], 4)
  ), row.names = FALSE)
}

cat("\nBSSH_2d and BSSH_4d give the same alarm dates:",
    identical(alarms$BSSH_2d, alarms$BSSH_4d), "\n")
cat("Alarm dates in BSSH_d but not in BSSH_2d:",
    format(mon_dates[setdiff(alarms$BSSH_d, alarms$BSSH_2d)]), "\n")
cat("Alarm dates in BSSH_2d but not in BSSH_d:",
    format(mon_dates[setdiff(alarms$BSSH_2d, alarms$BSSH_d)]), "\n")

#==============================================================================
# POST-ALARM DIAGNOSIS FOR BSSH_2d (MYT AND ORDER-FREE DIAGNOSTIC)
#==============================================================================

al <- alarms$BSSH_2d
diag_tab <- data.frame(Date = format(mon_dates[al]), Fault = ifelse(mon_fault[al], "YES", "NO"),
                       MYT = NA, Order_free = NA, MYT_variables = "", Order_free_variables = "")
for (k in 1:length(al)) {
  i <- al[k]
  v_myt <- myt(Y_mon[i, ], results$BSSH_2d$post[[i]])
  v_prop <- order_free(Y_mon[i, ], results$BSSH_2d$post[[i]])
  diag_tab$MYT[k] <- length(v_myt)
  diag_tab$Order_free[k] <- length(v_prop)
  diag_tab$MYT_variables[k] <- paste(var_labels[v_myt], collapse = ", ")
  diag_tab$Order_free_variables[k] <- paste(var_labels[v_prop], collapse = ", ")
}

cat("\n===== FLAGGED VARIABLES AT EACH ALARM: BSSH_2d =====\n")
print(diag_tab, row.names = FALSE, right = FALSE)

fault <- diag_tab$Fault == "YES"
cat("\n===== AVERAGE NUMBER OF FLAGGED VARIABLES: BSSH_2d =====\n")
print(data.frame(
  Alarm_type = c("Fault labelled", "Nonfault"),
  Number_of_alarms = c(sum(fault), sum(!fault)),
  MYT = round(c(mean(diag_tab$MYT[fault]), mean(diag_tab$MYT[!fault])), 2),
  Order_free = round(c(mean(diag_tab$Order_free[fault]), mean(diag_tab$Order_free[!fault])), 2)
), row.names = FALSE)

#==============================================================================
# PLOT
#==============================================================================

month_start <- as.Date(c("1991-01-01", "1991-03-01", "1991-05-01",
                         "1991-07-01", "1991-08-01", "1991-10-01"))
month_at <- rep(NA, 6)
for (k in 1:6) month_at[k] <- which.min(abs(as.numeric(mon_dates - month_start[k])))
month_lab <- c("Jan", "Mar", "May", "Jul", "Aug", "Oct")

bssh_panel <- function(st, title) {
  tested <- which(!is.na(st[, 1]))
  alarm <- tested[st[tested, 1] > st[tested, 2]]
  fault <- tested[mon_fault[tested]]
  plot(NA, xlim = c(1, nrow(st)), ylim = c(0, 20), xlab = "", ylab = "",
       axes = FALSE, xaxs = "i", yaxs = "i")
  polygon(c(tested, rev(tested)), c(st[tested, 2], rep(0, length(tested))),
          col = rgb(0.05, 0.85, 0.05, 0.16), border = NA)
  lines(tested, st[tested, 2], col = "firebrick2", lwd = 2.8)
  lines(tested, st[tested, 1], type = "b", pch = 16, cex = 0.68, lwd = 1, col = "dodgerblue3")
  points(fault, st[fault, 1], pch = 17, cex = 1.3, col = "firebrick4")
  points(alarm, st[alarm, 1], pch = 1, cex = 1.75, lwd = 2.6, col = "firebrick1")
  axis(1, at = month_at, labels = month_lab)
  axis(2, at = seq(0, 20, by = 5))
  mtext(title, side = 3, line = 0.05, cex = 1.45)
}

pdf(file = "BSSH_application_water.pdf", width = 13, height = 8.4)

layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.28))
par(bg = "white", oma = c(0.12, 0.55, 1.25, 0.45))

par(mar = c(3.2, 5.2, 1.45, 1))
bssh_panel(results$BSSH_0$stats, expression(paste(BSSH[0], " / FSSH")))
mtext(expression(t[n+1]^2), side = 2, line = 2.05, cex = 1.18)

par(mar = c(3.2, 3.4, 1.45, 1))
bssh_panel(results$BSSH_d$stats, expression(BSSH[d]))

par(mar = c(4.1, 5.2, 1.65, 1))
bssh_panel(results$BSSH_2d$stats, expression(BSSH[2 * d]))
mtext(expression(t[n+1]^2), side = 2, line = 2.05, cex = 1.18)
mtext("1991 Monitoring Period", side = 1, line = 3.15, cex = 1.1)

par(mar = c(4.1, 3.4, 1.65, 1))
bssh_panel(results$BSSH_4d$stats, expression(BSSH[4 * d]))
mtext("1991 Monitoring Period", side = 1, line = 3.15, cex = 1.1)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Non-fault-labelled", "Fault-labelled", "UCL", "Alarm"), ncol = 4,
       col = c("dodgerblue3", "firebrick4", "firebrick2", "firebrick1"),
       pch = c(16, 17, NA, 1), lty = c(1, NA, 1, NA), lwd = c(1.6, NA, 3.2, NA),
       pt.cex = c(1.18, 1.42, NA, 1.7), bty = "n", cex = 1.35, seg.len = 3.25,
       x.intersp = 1.1, y.intersp = 1.1, xpd = NA)

dev.off()

cat("\nPDF saved to: BSSH_application_water.pdf\n")
