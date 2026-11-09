# R code for the transportation cost application of BSSH
# Reproduces the monitoring results, post-alarm diagnosis,
# and the figure reported in the manuscript.

# If needed, install the required packages:
# install.packages("plotrix")
# install.packages("scatterplot3d")

library(plotrix)
library(scatterplot3d)

#------------------------------------------------------------------------------
# DATA
#------------------------------------------------------------------------------

Fuel <- c(16.44, 7.19, 9.92, 4.24, 11.20, 14.25, 13.50, 13.32, 29.11, 12.68,
          7.51, 9.90, 10.25, 11.11, 12.17, 10.24, 10.18, 8.88, 12.34, 8.51,
          26.16, 12.95, 16.93, 14.70, 10.32, 8.98, 9.70, 12.72, 9.49, 8.22,
          13.70, 8.21, 15.86, 9.18, 12.49, 17.32)

Repair <- c(12.43, 2.70, 1.35, 5.78, 5.05, 5.78, 10.98, 14.27, 15.09, 7.61,
            5.80, 3.63, 5.07, 6.15, 14.26, 2.59, 6.05, 2.70, 7.73, 14.02,
            17.44, 8.24, 13.37, 10.78, 5.16, 4.49, 11.59, 8.63, 2.16, 7.95,
            11.22, 9.85, 11.42, 9.18, 4.67, 6.86)

Capital <- c(11.23, 3.92, 9.75, 7.78, 10.67, 9.88, 10.60, 9.45, 3.28, 10.23,
             8.13, 9.13, 10.17, 7.61, 14.39, 6.09, 12.14, 12.23, 11.68, 12.01,
             16.89, 7.18, 17.59, 14.58, 17.00, 4.26, 6.83, 5.59, 6.23, 6.72,
             4.91, 8.17, 13.06, 9.49, 11.49, 4.44)

#######################
#### PRIOR SETTING ####
#######################

M0n <- rep(0, 3)
L0n <- 0
NU0n <- -1
PSIn <- matrix(0, 3, 3)
a0 <- 0

# FALSE ALARM RATE: ARL_0 = 200
far <- 1/200

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

posterior_upto_n_recursive <- function(X, n, m0, l0, nu0, psi0) {
  mu <- m0; lambda <- l0; nu <- nu0; Psi <- psi0
  if (n >= 1) {
    for (i in 1:n) {
      upd <- niw_update_one(X[i,], mu, lambda, nu, Psi)
      mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
    }
  }
  list(mu = mu, lambda = lambda, nu = nu, Psi = Psi)
}

#==============================================================================
# applicationBSSH 
#==============================================================================

applicationBSSH <- function(CX, CY, CZ, PX, PY, PZ, FAR, a0, m0i, l0i, nu0i, psi0i) {
  stats <- matrix(NA, length(CX), 2)
  colnames(stats) <- c("t2", "qf")
  
  if (a0 != 0) {
    Y <- cbind(PX, PY, PZ)
    n0 <- length(PX)
    sumY <- apply(Y, 2, sum)
    meanY <- sumY / n0
    CC0 <- lapply(1:n0, function(x) (Y[x,] - meanY) %*% t(Y[x,] - meanY))
    C0 <- Reduce("+", CC0)
    
    m0 <- (l0i * m0i + a0 * sumY) / (l0i + a0 * n0)
    l0 <- l0i + a0 * n0
    nu0 <- nu0i + a0 * n0
    psi0 <- psi0i + a0 * C0 + l0i * m0i %*% t(m0i) +
      a0 * n0 * meanY %*% t(meanY) - l0 * m0 %*% t(m0)
  } else {
    m0 <- m0i; l0 <- l0i; nu0 <- nu0i; psi0 <- psi0i
  }
  
  X <- cbind(CX, CY, CZ)
  d <- ncol(X); N <- nrow(X)
  start_n <- d + 1
  
  mu <- m0; lambda <- l0; nu <- nu0; Psi <- psi0
  for (i in 1:start_n) {
    upd <- niw_update_one(X[i,], mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
  }
  
  for (n in start_n:(N - 1)) {
    x_next <- X[n + 1,]
    nu_star <- nu - d + 1
    
    A <- ((lambda + 1) * d) / (lambda * nu_star) * Psi
    diff <- matrix(x_next - mu, ncol = 1)
    
    t2 <- as.numeric(t(diff) %*% solve(A) %*% diff)
    q  <- qf(1 - FAR, d, nu_star)
    
    stats[n + 1, 1] <- t2
    stats[n + 1, 2] <- q
    
    upd <- niw_update_one(x_next, mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
  }
  
  return(stats)
}

#==============================================================================
# posteriormeanBSSH
#==============================================================================

prmeanBSSH <- function(CX, CY, CZ, PX, PY, PZ, FAR, a0, m0i, l0i, nu0i, psi0i) {
  MN <- NULL
  
  if (a0 != 0) {
    Y <- cbind(PX, PY, PZ)
    n0 <- length(PX)
    sumY <- apply(Y, 2, sum)
    meanY <- sumY / n0
    CC0 <- lapply(1:n0, function(x) (Y[x,] - meanY) %*% t(Y[x,] - meanY))
    C0 <- Reduce("+", CC0)
    
    m0 <- (l0i * m0i + a0 * sumY) / (l0i + a0 * n0)
    l0 <- l0i + a0 * n0
    nu0 <- nu0i + a0 * n0
    psi0 <- psi0i + a0 * C0 + l0i * m0i %*% t(m0i) +
      a0 * n0 * meanY %*% t(meanY) - l0 * m0 %*% t(m0)
  } else {
    m0 <- m0i; l0 <- l0i; nu0 <- nu0i; psi0 <- psi0i
  }
  
  X <- cbind(CX, CY, CZ)
  d <- ncol(X); N <- nrow(X)
  start_n <- d + 1
  
  mu <- m0; lambda <- l0; nu <- nu0; Psi <- psi0
  for (i in 1:start_n) {
    upd <- niw_update_one(X[i,], mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
  }
  
  for (n in start_n:(N - 1)) {
    MN <- rbind(MN, mu)
    upd <- niw_update_one(X[n + 1,], mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
  }
  
  colnames(MN) <- c("Fuel", "Repair", "Capital")
  return(MN)
}

#==============================================================================
# EXCLUDING OUTLIERS 
#==============================================================================

applicationBSSHexcluding <- function(CX, CY, CZ, PX, PY, PZ, FAR, a0, m0i, l0i, nu0i, psi0i) {
  stats <- matrix(NA, 4, 2)
  colnames(stats) <- c("t2", "qf")
  
  if (a0 != 0) {
    Y <- cbind(PX, PY, PZ)
    n0 <- length(PX)
    sumY <- apply(Y, 2, sum)
    meanY <- sumY / n0
    CC0 <- lapply(1:n0, function(x) (Y[x,] - meanY) %*% t(Y[x,] - meanY))
    C0 <- Reduce("+", CC0)
    
    m0 <- (l0i * m0i + a0 * sumY) / (l0i + a0 * n0)
    l0 <- l0i + a0 * n0
    nu0 <- nu0i + a0 * n0
    psi0 <- psi0i + a0 * C0 + l0i * m0i %*% t(m0i) +
      a0 * n0 * meanY %*% t(meanY) - l0 * m0 %*% t(m0)
  } else {
    m0 <- m0i; l0 <- l0i; nu0 <- nu0i; psi0 <- psi0i
  }
  
  X <- cbind(CX, CY, CZ)
  d <- ncol(X)
  n <- d + 1
  
  while (n <= (nrow(X) - 1)) {
    post <- posterior_upto_n_recursive(X, n, m0, l0, nu0, psi0)
    mu <- post$mu; lambda <- post$lambda; nu <- post$nu; Psi <- post$Psi
    
    nu_star <- nu - d + 1
    A <- ((lambda + 1) * d) / (lambda * nu_star) * Psi
    
    diff <- matrix(X[n + 1, ] - mu, ncol = 1)
    t2 <- as.numeric(t(diff) %*% solve(A) %*% diff)
    q  <- qf(1 - FAR, d, nu_star)
    
    stats <- rbind(stats, c(t2, q))
    
    if (t2 > q) {
      X <- X[-(n + 1), ]
      n <- n - 1
    }
    
    n <- n + 1
  }
  
  return(stats)
}

prmeanBSSHexcluding <- function(CX, CY, CZ, PX, PY, PZ, FAR, a0, m0i, l0i, nu0i, psi0i) {
  MN <- NULL
  
  if (a0 != 0) {
    Y <- cbind(PX, PY, PZ)
    n0 <- length(PX)
    sumY <- apply(Y, 2, sum)
    meanY <- sumY / n0
    CC0 <- lapply(1:n0, function(x) (Y[x,] - meanY) %*% t(Y[x,] - meanY))
    C0 <- Reduce("+", CC0)
    
    m0 <- (l0i * m0i + a0 * sumY) / (l0i + a0 * n0)
    l0 <- l0i + a0 * n0
    nu0 <- nu0i + a0 * n0
    psi0 <- psi0i + a0 * C0 + l0i * m0i %*% t(m0i) +
      a0 * n0 * meanY %*% t(meanY) - l0 * m0 %*% t(m0)
  } else {
    m0 <- m0i; l0 <- l0i; nu0 <- nu0i; psi0 <- psi0i
  }
  
  X <- cbind(CX, CY, CZ)
  d <- ncol(X)
  n <- d + 1
  
  while (n <= (nrow(X) - 1)) {
    post <- posterior_upto_n_recursive(X, n, m0, l0, nu0, psi0)
    mu <- post$mu; lambda <- post$lambda; nu <- post$nu; Psi <- post$Psi
    
    nu_star <- nu - d + 1
    A <- ((lambda + 1) * d) / (lambda * nu_star) * Psi
    
    diff <- matrix(X[n + 1, ] - mu, ncol = 1)
    t2 <- as.numeric(t(diff) %*% solve(A) %*% diff)
    q  <- qf(1 - FAR, d, nu_star)
    
    MN <- rbind(MN, mu)
    
    if (t2 > q) {
      X <- X[-(n + 1), ]
      n <- n - 1
    }
    
    n <- n + 1
  }
  
  colnames(MN) <- c("Fuel", "Repair", "Capital")
  return(MN)
}

#==============================================================================
# POST-ALARM DIAGNOSIS 
#==============================================================================

all_subsets_two <- function(vec) {
  list(integer(0), vec[1], vec[2], vec)
}

diag_all_3d <- function(x_alarm, post, time, q_fdr = 0.05) {
  mu <- post$mu; ln <- post$lambda; nun <- post$nu; Psi <- post$Psi
  d <- length(mu)
  nu_star <- nun - d + 1
  A <- ((ln + 1) * d) / (ln * nu_star) * Psi
  
  rows <- list()
  for (j in 1:d) {
    others <- setdiff(1:d, j)
    subsets <- all_subsets_two(others)
    
    for (L in subsets) {
      L <- as.integer(L)
      Lsize <- length(L)
      
      if (Lsize == 0) {
        t2 <- (x_alarm[j] - mu[j])^2 / A[j, j]
        U  <- d * t2
        df2 <- nu_star
        p <- 1 - pf(U, df1 = 1, df2 = df2)
        test <- sprintf("t^2_{%d,%d}", j, time)
        code <- sprintf("t%d", j)
      } else {
        ALL <- A[L, L, drop = FALSE]
        invALL <- solve(ALL)
        diffL <- matrix(x_alarm[L] - mu[L], ncol = 1)
        
        DeltaL <- as.numeric(d * t(diffL) %*% invALL %*% diffL)
        
        AjL <- matrix(A[j, L], nrow = 1)
        ALj <- matrix(A[L, j], ncol = 1)
        
        mu_jgL <- as.numeric(mu[j] + AjL %*% invALL %*% diffL)
        A_jgL  <- as.numeric(A[j, j] - AjL %*% invALL %*% ALj)
        
        t2 <- (x_alarm[j] - mu_jgL)^2 / A_jgL
        U  <- ((nu_star + Lsize) / (nu_star + DeltaL)) * d * t2
        df2 <- nu_star + Lsize
        p <- 1 - pf(U, df1 = 1, df2 = df2)
        
        Llab <- paste(L, collapse = ",")
        test <- sprintf("t^2_{%d|{%s},%d}", j, Llab, time)
        code <- sprintf("t%dg%s", j, paste(L, collapse = ""))
      }
      
      rows[[length(rows) + 1]] <- data.frame(
        Test = test,
        code = code,
        t2 = t2,
        U = U,
        df2 = df2,
        p = p,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df <- do.call(rbind, rows)
  df$p_adj <- p.adjust(df$p, method = "BH")
  df$Significant <- df$p_adj <= q_fdr
  df <- df[order(df$p), ]
  
  j_from_code <- function(code) as.integer(sub("^t([0-9]).*$", "\\1", code))
  suspects <- sort(unique(j_from_code(df$code[df$Significant])))
  
  list(table = df, suspects = suspects)
}

#==============================================================================
# RUN BSSH MONITORING
#==============================================================================

results   <- applicationBSSH(Fuel, Repair, Capital, NA, NA, NA, far, a0, M0n, L0n, NU0n, PSIn)
resultsex <- applicationBSSHexcluding(Fuel, Repair, Capital, NA, NA, NA, far, a0, M0n, L0n, NU0n, PSIn)

prmean   <- prmeanBSSH(Fuel, Repair, Capital, NA, NA, NA, far, a0, M0n, L0n, NU0n, PSIn)
prmeanex <- prmeanBSSHexcluding(Fuel, Repair, Capital, NA, NA, NA, far, a0, M0n, L0n, NU0n, PSIn)

monitor_all <- data.frame(
  Time = 1:length(Fuel),
  Statistic = round(results[,1], 4),
  Threshold = round(results[,2], 4),
  Alarm = ifelse(!is.na(results[,1]) & results[,1] > results[,2], "YES", "NO")
)

monitor_ex <- data.frame(
  Time = 1:nrow(resultsex),
  Statistic = round(resultsex[,1], 4),
  Threshold = round(resultsex[,2], 4),
  Alarm = ifelse(!is.na(resultsex[,1]) & resultsex[,1] > resultsex[,2], "YES", "NO")
)

alarm_all <- which(!is.na(results[,1]) & results[,1] > results[,2])
alarm_ex  <- which(!is.na(resultsex[,1]) & resultsex[,1] > resultsex[,2])

cat("\n===== STATISTICS AND THRESHOLDS: ALL OBSERVATIONS =====\n")
print(monitor_all, row.names = FALSE)

cat("\n===== ALARM VALUES: ALL OBSERVATIONS =====\n")
if (length(alarm_all) == 0) {
  cat("No alarms.\n")
} else {
  print(data.frame(
    Time = alarm_all,
    Fuel = round(Fuel[alarm_all], 4),
    Repair = round(Repair[alarm_all], 4),
    Capital = round(Capital[alarm_all], 4)
  ), row.names = FALSE)
}

cat("\n===== STATISTICS AND THRESHOLDS: EXCLUDING DETECTED OUTLIERS =====\n")
print(monitor_ex, row.names = FALSE)

cat("\n===== ALARM VALUES: EXCLUDING DETECTED OUTLIERS =====\n")
if (length(alarm_ex) == 0) {
  cat("No alarms.\n")
} else {
  print(data.frame(
    Time = alarm_ex,
    Fuel = round(Fuel[alarm_ex], 4),
    Repair = round(Repair[alarm_ex], 4),
    Capital = round(Capital[alarm_ex], 4)
  ), row.names = FALSE)
}

#==============================================================================
# POST-ALARM DIAGNOSIS (BH-FDR=0.05) FOR t=9 AND t=21
#==============================================================================

Xfull <- cbind(Fuel, Repair, Capital)

t_alarm1 <- 9
post8 <- posterior_upto_n_recursive(Xfull, n = t_alarm1 - 1, m0 = M0n, l0 = L0n, nu0 = NU0n, psi0 = PSIn)
diag9 <- diag_all_3d(Xfull[t_alarm1, ], post8, time = t_alarm1, q_fdr = 0.05)

tab9 <- within(diag9$table, {
  t2 <- round(t2, 6); U <- round(U, 6)
  p <- signif(p, 6); p_adj <- signif(p_adj, 6)
})

Xno9 <- Xfull[-9, , drop = FALSE]
t_alarm2_original <- 21
t_alarm2_stream <- t_alarm2_original - 1

post_before21 <- posterior_upto_n_recursive(Xno9, n = t_alarm2_stream - 1, m0 = M0n, l0 = L0n, nu0 = NU0n, psi0 = PSIn)
diag21 <- diag_all_3d(Xno9[t_alarm2_stream, ], post_before21, time = t_alarm2_original, q_fdr = 0.05)

tab21 <- within(diag21$table, {
  t2 <- round(t2, 6); U <- round(U, 6)
  p <- signif(p, 6); p_adj <- signif(p_adj, 6)
})

cat("\n===== POST-ALARM P-VALUES =====\n")

cat("\n----- Alarm at time 9 -----\n")
cat("Fuel = ", round(Xfull[t_alarm1, 1], 4),
    ", Repair = ", round(Xfull[t_alarm1, 2], 4),
    ", Capital = ", round(Xfull[t_alarm1, 3], 4), "\n", sep = "")
print(tab9, row.names = FALSE)

cat("\n----- Alarm at time 21 -----\n")
cat("Fuel = ", round(Xno9[t_alarm2_stream, 1], 4),
    ", Repair = ", round(Xno9[t_alarm2_stream, 2], 4),
    ", Capital = ", round(Xno9[t_alarm2_stream, 3], 4), "\n", sep = "")
print(tab21, row.names = FALSE)

#==============================================================================
# PLOT
#==============================================================================

test_to_plotmath <- function(test_str) {
  inside <- sub("^t\\^2_\\{", "", test_str)
  inside <- sub("\\}$", "", inside)
  time <- sub("^.*,", "", inside)
  rest <- sub(",[^,]*$", "", inside)
  
  if (!grepl("\\|", rest)) {
    return(parse(text = sprintf("t[%s*','*%s]^2", rest, time))[[1]])
  }
  
  parts <- strsplit(rest, "\\|")[[1]]
  j <- parts[1]
  Lparts <- strsplit(gsub("[{}]", "", parts[2]), ",")[[1]]
  Lstr <- paste(Lparts, collapse = "*','*")
  
  parse(text = sprintf("t[%s*'|'*'{'*%s*'}'*','*%s]^2", j, Lstr, time))[[1]]
}

build_inline_expr <- function(expr_list) {
  if (length(expr_list) == 0) return("none")
  if (length(expr_list) == 1) return(as.expression(expr_list[[1]]))
  
  call_list <- list(as.name("paste"))
  for (k in seq_along(expr_list)) {
    call_list <- c(call_list, list(expr_list[[k]]))
    if (k < length(expr_list)) call_list <- c(call_list, list(", "))
  }
  as.expression(as.call(call_list))
}

sig9_tbl  <- diag9$table[diag9$table$Significant, ]
sig21_tbl <- diag21$table[diag21$table$Significant, ]

sig9_tbl  <- sig9_tbl[order(sig9_tbl$p), ]
sig21_tbl <- sig21_tbl[order(sig21_tbl$p), ]

sig9_tests  <- sig9_tbl$Test
sig21_tests <- sig21_tbl$Test

sig9_inline  <- build_inline_expr(lapply(sig9_tests,  test_to_plotmath))
sig21_inline <- build_inline_expr(lapply(sig21_tests, test_to_plotmath))

header_h <- 0.18
body_top <- 1 - header_h

pdf(file = "BSSH_application_cost.pdf",
    width = 12.75,
    height = 9 / body_top)

par(fig = c(0, 1/3, body_top, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = FALSE)
plot(NA, NA, axes = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))

legend(
  x = 0.05, y = 0.95, xjust = 0, yjust = 1,
  legend = c(
    expression(t[n+1]^2*" (excluding detected outliers)"),
    expression(t[n+1]^2*" (all observations)"),
    "UCL (excluding detected outliers)",
    "UCL (all observations)",
    expression(mu[n]*" (excluding detected outliers)"),
    expression(mu[n]*" (all observations)"),
    "OOC at observation 9",
    "OOC at observation 21"
  ),
  pch = c(16, 1, NA, NA, NA, NA, 8, 4),
  lty = c(1, 2, 1, 2, 1, 2, NA, NA),
  lwd = c(1, 1, 2, 2, 2, 2, NA, NA),
  col = c("dodgerblue4", "dodgerblue2", "firebrick3", "firebrick1",
          "dodgerblue4", "dodgerblue2", "firebrick4", "firebrick1"),
  pt.cex = c(1.2, 1.2, NA, NA, NA, NA, 1.2, 1.2),
  cex = 0.95,
  bty = "n",
  y.intersp = 1.05
)

#################################################################################################

par(fig = c(1/3, 2/3, body_top, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(NA, NA, axes = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.72, "Significant components at observation 9", font = 2, cex = 1.05, adj = c(0.5, 0.5))
text(0.5, 0.40, sig9_inline, cex = 0.95, adj = c(0.5, 0.5))

#################################################################################################

par(fig = c(2/3, 1, body_top, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(NA, NA, axes = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.72, "Significant components at observation 21", font = 2, cex = 1.05, adj = c(0.5, 0.5))
text(0.5, 0.40, sig21_inline, cex = 0.95, adj = c(0.5, 0.5))

#################################################################################################

bord <- c(5:length(Fuel), length(Fuel):5)
bord[1] <- bord[64] <- 4.80
bord[32] <- bord[33] <- 36.20
obs <- 1:length(Fuel)
obs[5] <- 1.80
obs[33] <- 36.20

#################################################################################################

par(fig = c(0, 2/3, 0.55*body_top, 0.95*body_top),
    oma = c(1.5, 1, 1.5, 1), mar = c(2, 2.5, 2.5, 1.5), new = TRUE)

plot(resultsex[,1], type = "b", lty = 1, xlab = "", ylab = "", col = "dodgerblue4",
     lwd = 1, pch = 16, axes = FALSE, yaxs = "i", ylim = c(-0.5, 20), xlim = c(5, 36))

points(results[,1], type = "b", lty = 2, col = "dodgerblue2", lwd = 1, pch = 1)

points(resultsex[,2], type = "l", lty = 1, lwd = 2, col = "firebrick3")
points(results[,2], type = "l", lty = 2, lwd = 2, col = "firebrick1")

polygon(bord, c(resultsex[5:36,2], rep(0, 32)),
        col = rgb(0.05, 0.95, 0.05, 0.2), border = NA)

points(9,  resultsex[9,1],  pch = 8, col = "firebrick4", lwd = 2, cex = 1.5)
points(21, resultsex[21,1], pch = 4, col = "firebrick1", lwd = 2, cex = 1.5)

draw.circle(9,  resultsex[9,1],  2/5, border = "firebrick4", lwd = 2)
draw.circle(21, resultsex[21,1], 2/5, border = "firebrick1", lwd = 2)

axis(1); axis(2)
mtext(expression(t[n+1]^2), side = 2, line = 2.1, cex = 1.1)
mtext("Time",           side = 1, line = 2.1, cex = 1.1)
mtext("BSSH",           side = 3, line = 0.4, cex = 2)

par(fig = c(0, 2/3, 0*body_top, 0.5*body_top),
    oma = c(0, 1, 1.5, 1), mar = c(0, 2.5, 4, 1.5), new = TRUE)

spl <- scatterplot3d(cbind(Fuel, Repair, Capital),
                     pch = c(rep(16, 8), 8, rep(16, 11), 4, rep(16, 15)),
                     color = c(rep(1, 8), "firebrick4", rep(1, 11),
                               "firebrick1", rep(1, 15)),
                     angle = 45, box = FALSE, type = "h")

spl$points3d(Fuel[9],  Repair[9],  Capital[9],  pch = 1, cex = 1.6, col = "firebrick4")
spl$points3d(Fuel[21], Repair[21], Capital[21], pch = 1, cex = 1.6, col = "firebrick1")

mtext("3D scatterplot", side = 3, line = 0.4, cex = 2)

#################################################################################################

par(fig = c(2/3, 1, (2/3)*body_top, 1*body_top),
    oma = c(1.5, 1, 1.5, 1), mar = c(2.5, 2.5, 3, 1.5), new = TRUE)

plot(NA, xlab = "", ylab = "", axes = FALSE, yaxs = "i", xlim = c(1, 36), ylim = c(3, 32))

points(5:36, prmeanex[,1], type = "l", lty = 1, col = "dodgerblue4")
points(5:36, prmean[,1],   type = "l", lty = 2, col = "dodgerblue2")

points(Fuel[1:36], type = "b", lty = 1, lwd = 1,
       pch = c(rep(16, 8), 8, rep(16, 11), 4, rep(16, 15)),
       col = c(rep(1, 8), "firebrick4", rep(1, 11),
               "firebrick1", rep(1, 15)))

draw.circle(9,  Fuel[9],  3/5, border = "firebrick4")
draw.circle(21, Fuel[21], 3/5, border = "firebrick1")

axis(1); axis(2)
mtext("Values", side = 2, line = 2.1, cex = 1.1)
mtext("Time",   side = 1, line = 2.1, cex = 1.1)
mtext("Fuel",   side = 3, line = 0.4, cex = 1.4)

#################################################################################################

par(fig = c(2/3, 1, (1/3)*body_top, (2/3)*body_top),
    oma = c(1.5, 1, 1.5, 1), mar = c(2.5, 2.5, 3, 1.5), new = TRUE)

plot(NA, xlab = "", ylab = "", axes = FALSE, yaxs = "i", xlim = c(1, 36), ylim = c(0.5, 19.5))

points(5:36, prmeanex[,2], type = "l", lty = 1, col = "dodgerblue4")
points(5:36, prmean[,2],   type = "l", lty = 2, col = "dodgerblue2")

points(Repair[1:36], type = "b", lty = 1, lwd = 1,
       pch = c(rep(16, 8), 8, rep(16, 11), 4, rep(16, 15)),
       col = c(rep(1, 8), "firebrick4", rep(1, 11),
               "firebrick1", rep(1, 15)))

draw.circle(9,  Repair[9],  3/5, border = "firebrick4")
draw.circle(21, Repair[21], 3/5, border = "firebrick1")

axis(1); axis(2)
mtext("Values", side = 2, line = 2.1, cex = 1.1)
mtext("Time",   side = 1, line = 2.1, cex = 1.1)
mtext("Repair", side = 3, line = 0.4, cex = 1.4)

#################################################################################################

par(fig = c(2/3, 1, 0*body_top, (1/3)*body_top),
    oma = c(1.5, 1, 1.5, 1), mar = c(2.5, 2.5, 3, 1.5), new = TRUE)

plot(NA, xlab = "", ylab = "", axes = FALSE, yaxs = "i", xlim = c(1, 36), ylim = c(1.5, 19))

points(5:36, prmeanex[,3], type = "l", lty = 1, col = "dodgerblue4")
points(5:36, prmean[,3],   type = "l", lty = 2, col = "dodgerblue2")

points(Capital[1:36], type = "b", lty = 1, lwd = 1,
       pch = c(rep(16, 8), 8, rep(16, 11), 4, rep(16, 15)),
       col = c(rep(1, 8), "firebrick4", rep(1, 11),
               "firebrick1", rep(1, 15)))

draw.circle(9,  Capital[9],  3/5, border = "firebrick4")
draw.circle(21, Capital[21], 3/5, border = "firebrick1")

axis(1); axis(2)
mtext("Values",  side = 2, line = 2.1, cex = 1.1)
mtext("Time",    side = 1, line = 2.1, cex = 1.1)
mtext("Capital", side = 3, line = 0.4, cex = 1.4)

dev.off()

cat("\nPDF saved to: BSSH_application_cost.pdf\n")