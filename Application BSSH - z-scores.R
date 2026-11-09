# R code for the medical laboratory application of BSSH
# Reproduces the monitoring results, posterior summaries,
# and the figure reported in the manuscript.

# If needed, install the required packages:
# install.packages("plotrix")
# install.packages("emdbook")
# install.packages("LaplacesDemon")
# install.packages("mnormt")

library(plotrix)
library(emdbook)
library(LaplacesDemon)
library(mnormt)

#------------------------------------------------------------------------------
# DATA
#------------------------------------------------------------------------------

z1=c(0.51, 0.49, 0.28, -0.01, 0.81, 0.08, 1.27, 2.91)
z2=c(2.22, -0.92, 0.99, 1.07, 1.18, 0.95, 1.15, -0.2)

Z1=c(0.44, 0.33, 0.29, -0.5, 0.57, 0.17, 0.93, 1.41, 0.85, -0.16, 1.81, 0.66,
     -0.68, -0.15, 0.16, 0.51, 0.15, -0.4, 2.26, 0.21, 1.82, 0.98, -0.09,
     1.24, 2.37, 1.45, 0.52, 0.51, 0.49, 0.28, -0.01, 0.81, 0.08, 1.27, 2.91)
Z2=c(0.34, 1.08, 0.77, 0.03, 1.24, 0.56, 1.45, 1.4, 0.32, 0.7, 1, 0.59,
     0.63, 1.2, 0.11, -0.63, 0.64, -0.11, 1.45, 1.76, -0.01, 1.46, 1.87,
     2.7, 0.74, 0.99, 1.45, 2.22, -0.92, 0.99, 1.07, 1.18, 0.95, 1.15, -0.2)

#######################
#### PRIOR SETTING ####
#######################

M0n=c(0,0)
L0n=10
NU0n=10
PSIn=9*matrix(c(1,0.6,0.6,1),2,2)

M0i=c(0,0)
L0i=0
NU0i=-1
PSIi=matrix(0,2,2)

a0=0
far=1/370

#==============================================================================
# RECURSIVE NIW UPDATE 
#==============================================================================

niw_update_one <- function(x, mu, lambda, nu, Psi){
  x <- as.numeric(x)
  w <- lambda/(lambda+1)
  diff <- x - mu
  mu_new  <- w*mu + (1-w)*x
  Psi_new <- Psi + w*tcrossprod(diff)
  list(mu=mu_new, lambda=lambda+1, nu=nu+1, Psi=Psi_new)
}

#==============================================================================
# applicationBSSH 
#==============================================================================

applicationBSSH=function(CX, CY, PX, PY, FAR, a0, m0i, l0i, nu0i, psi0i){
  stats=matrix(NA,length(CX),2)
  colnames(stats)=c("t2","qf")
  
  if (a0!=0) {
    Y=cbind(PX,PY)
    n0=length(PX)
    sumY=apply(Y,2,sum)
    meanY=sumY/n0
    CC0=lapply( 1:n0, function(x) (Y[x,]-meanY)%*%t(Y[x,]-meanY))
    C0=Reduce("+",CC0)
    
    m0 = (l0i*m0i+a0*sumY)/(l0i+a0*n0)
    l0 = l0i + a0*n0
    nu0 = nu0i + a0*n0
    psi0 = psi0i + a0*C0  + l0i*m0i%*%t(m0i) +  a0*n0*meanY%*%t(meanY)  - l0*m0%*%t(m0)
  } else{
    m0=m0i ; l0=l0i ; nu0=nu0i ; psi0=psi0i
  }
  
  X=cbind(CX,CY)
  d=ncol(X)
  N=nrow(X)
  
  start_n <- if (nu0==-1) 3 else 1
  
  mu <- m0; lambda <- l0; nu <- nu0; Psi <- psi0
  
  for(i in 1:start_n){
    upd <- niw_update_one(X[i,], mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
  }
  
  for(n in start_n:(N-1)){
    xnext <- X[n+1,]
    nu_star <- nu - d + 1
    
    A <- (lambda+1)*d/(lambda*nu_star) * Psi
    diff <- matrix(xnext - mu, ncol=1)
    
    t2 <- as.numeric(t(diff) %*% solve(A) %*% diff)
    q  <- qf(1-FAR, d, nu_star)
    
    stats[n+1,1]=t2
    stats[n+1,2]=q
    
    upd <- niw_update_one(xnext, mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
  }
  
  return(stats)
}

#==============================================================================
# posteriormeanBSSH
#==============================================================================

posteriormeanBSSH=function(CX, CY, PX, PY, FAR, a0, m0i, l0i, nu0i, psi0i){
  MU <- matrix(NA, length(CX), 2)
  colnames(MU) <- c("mu1", "mu2")
  
  if (a0!=0) {
    Y=cbind(PX,PY)
    n0=length(PX)
    sumY=apply(Y,2,sum)
    meanY=sumY/n0
    CC0=lapply( 1:n0, function(x) (Y[x,]-meanY)%*%t(Y[x,]-meanY))
    C0=Reduce("+",CC0)
    
    m0 = (l0i*m0i+a0*sumY)/(l0i+a0*n0)
    l0 = l0i + a0*n0
    nu0 = nu0i + a0*n0
    psi0 = psi0i + a0*C0  + l0i*m0i%*%t(m0i) +  a0*n0*meanY%*%t(meanY)  - l0*m0%*%t(m0)
  } else{
    m0=m0i ; l0=l0i ; nu0=nu0i ; psi0=psi0i
  }
  
  MU[1,]=m0
  
  X=cbind(CX,CY)
  N=nrow(X)
  
  start_n <- if (nu0==-1) 3 else 1
  
  mu <- m0; lambda <- l0; nu <- nu0; Psi <- psi0
  
  for(n in 1:(N-1)){
    upd <- niw_update_one(X[n,], mu, lambda, nu, Psi)
    mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
    if(n >= start_n) MU[n+1,]=mu
  }
  
  return(MU)
}

postmeanV=posteriormeanBSSH(Z1, Z2, NA, NA, far, 0, M0n, L0n, NU0n, PSIn)
resultsV=applicationBSSH(Z1, Z2, NA, NA, far, 0, M0n, L0n, NU0n, PSIn)
resultsVI=applicationBSSH(Z1, Z2, NA, NA, far, 0, M0i, L0i, NU0i, PSIi)

monitor_V <- data.frame(
  Time = 1:length(Z1),
  Statistic = round(resultsV[,1], 4),
  Threshold = round(resultsV[,2], 4),
  Alarm = ifelse(!is.na(resultsV[,1]) & resultsV[,1] > resultsV[,2], "YES", "NO")
)

monitor_VI <- data.frame(
  Time = 1:length(Z1),
  Statistic = round(resultsVI[,1], 4),
  Threshold = round(resultsVI[,2], 4),
  Alarm = ifelse(!is.na(resultsVI[,1]) & resultsVI[,1] > resultsVI[,2], "YES", "NO")
)

alarm_V  <- which(!is.na(resultsV[,1])  & resultsV[,1]  > resultsV[,2])
alarm_VI <- which(!is.na(resultsVI[,1]) & resultsVI[,1] > resultsVI[,2])

cat("\n===== STATISTICS AND THRESHOLDS: BSSH =====\n")
print(monitor_V, row.names = FALSE)

cat("\n===== ALARM VALUES: BSSH =====\n")
if (length(alarm_V) == 0) {
  cat("No alarms.\n")
} else {
  print(data.frame(
    Time = alarm_V,
    z1 = round(Z1[alarm_V], 4),
    z2 = round(Z2[alarm_V], 4)
  ), row.names = FALSE)
}

cat("\n===== STATISTICS AND THRESHOLDS: FSSH =====\n")
print(monitor_VI, row.names = FALSE)

cat("\n===== ALARM VALUES: FSSH =====\n")
if (length(alarm_VI) == 0) {
  cat("No alarms.\n")
} else {
  print(data.frame(
    Time = alarm_VI,
    z1 = round(Z1[alarm_VI], 4),
    z2 = round(Z2[alarm_VI], 4)
  ), row.names = FALSE)
}

#==============================================================================
# POST-ALARM DIAGNOSIS 
#==============================================================================

posterior_upto_n <- function(X, n, mu0, lambda0, nu0, Psi0){
  mu <- mu0; lambda <- lambda0; nu <- nu0; Psi <- Psi0
  if(n>=1){
    for(i in 1:n){
      upd <- niw_update_one(X[i,], mu, lambda, nu, Psi)
      mu <- upd$mu; lambda <- upd$lambda; nu <- upd$nu; Psi <- upd$Psi
    }
  }
  list(mu=mu, lambda=lambda, nu=nu, Psi=Psi)
}

build_sig_label <- function(sig_codes){
  
  term <- function(code){
    if(code=="t1")   return(quote(t[1*','*35]^2))
    if(code=="t2")   return(quote(t[2*','*35]^2))
    if(code=="t2g1") return(quote(t[2*'|'*'{'*1*'}'*','*35]^2))
    if(code=="t1g2") return(quote(t[1*'|'*'{'*2*'}'*','*35]^2))
    quote(NA)
  }
  if(length(sig_codes)==0) return("none")
  terms <- lapply(sig_codes, term)
  if(length(terms)==1) return(as.expression(terms[[1]]))
  call_list <- list(as.name("paste"))
  for(k in seq_along(terms)){
    call_list <- c(call_list, list(terms[[k]]))
    if(k < length(terms)) call_list <- c(call_list, list(", "))
  }
  as.expression(as.call(call_list))
}

Xall <- cbind(Z1,Z2)
post34 <- posterior_upto_n(Xall, n=34, mu0=M0n, lambda0=L0n, nu0=NU0n, Psi0=PSIn)

mn  <- post34$mu
ln  <- post34$lambda
nun <- post34$nu
Psin <- post34$Psi
d <- 2
nu_star <- nun - d + 1
An <- (ln+1)*d/(ln*nu_star) * Psin

x_alarm <- c(Z1[35], Z2[35])

# t^2_1
t1 <- (x_alarm[1]-mn[1])^2 / An[1,1]
U1 <- d*t1
p1 <- 1 - pf(U1, df1=1, df2=nu_star)

# t^2_2
t2 <- (x_alarm[2]-mn[2])^2 / An[2,2]
U2 <- d*t2
p2 <- 1 - pf(U2, df1=1, df2=nu_star)

# t^2_{2|1}  
m21 <- mn[2] + (An[2,1]/An[1,1]) * (x_alarm[1]-mn[1])
a21 <- An[2,2] - (An[2,1]^2)/An[1,1]
t21 <- (x_alarm[2]-m21)^2 / a21
Delta1 <- d*t1
U21 <- ((nu_star+1)/(nu_star+Delta1)) * d*t21
p21 <- 1 - pf(U21, df1=1, df2=nu_star+1)

# t^2_{1|2}
m12 <- mn[1] + (An[1,2]/An[2,2]) * (x_alarm[2]-mn[2])
a12 <- An[1,1] - (An[1,2]^2)/An[2,2]
t12 <- (x_alarm[1]-m12)^2 / a12
Delta2 <- d*t2
U12 <- ((nu_star+1)/(nu_star+Delta2)) * d*t12
p12 <- 1 - pf(U12, df1=1, df2=nu_star+1)

diag_tbl <- data.frame(
  Test=c("t^2_{1,35}", "t^2_{2|{1},35}", "t^2_{2,35}", "t^2_{1|{2},35}"),
  code=c("t1","t2g1","t2","t1g2"),
  t2=c(t1,t21,t2,t12),
  U=c(U1,U21,U2,U12),
  df2=c(nu_star,nu_star+1,nu_star,nu_star+1),
  p=c(p1,p21,p2,p12),
  stringsAsFactors=FALSE
)
diag_tbl$p_adj <- p.adjust(diag_tbl$p, method="BH")
diag_tbl$Significant <- diag_tbl$p_adj <= 0.05

sig_tbl <- diag_tbl[diag_tbl$Significant, ]
sig_tbl <- sig_tbl[order(sig_tbl$p), ]   
sig_codes <- sig_tbl$code
sig_expr <- build_sig_label(sig_codes)

cat("\n===== POST-ALARM P-VALUES =====\n")
cat("\n----- Alarm at time 35 -----\n")
cat("z1 = ", round(x_alarm[1], 4),
    ", z2 = ", round(x_alarm[2], 4), "\n", sep = "")
print(diag_tbl[order(diag_tbl$p), ], row.names = FALSE)

#==============================================================================
# For the prior/posterior contour panels
#==============================================================================

mu_post_plot  <- round(mn, 3)
Psi_post_plot <- matrix(round(Psin, 2), 2, 2)
lambda_post <- ln
nu_post <- nun
nu_star_post <- nu_post - d + 1

# prior
mu0_plot <- c(0,0)
Psi0_plot <- PSIn
lambda0 <- L0n
nu0 <- NU0n
nu_star0 <- nu0 - d + 1

#==============================================================================
# PLOT
#==============================================================================

pdf(file = "BSSH_application-z-scores.pdf", width = 12, height = 17)

par(fig=c(0,3/4,13/17+0.01,1), oma=c(1.5,1,1.5,1), mar=c(2,2.5,2.5,1.5))

plot(resultsV[,1],type="b",lty=1,xlab="",ylab="",col="dodgerblue3", lwd=1, pch=16,
     axes=F, yaxs = "i", ylim=c(-0.5, 9), cex=2)
points(c(1:34,35.5), resultsV[,2],type="l",lty=1,lwd=3, col="firebrick2")

polygon( c( 1.5,3:34,35.5, 35.5, 34:3,1.5 ), c(resultsV[2:35,2],rep(0,34)),
         col = rgb(0.05,0.95,0.05,0.2), border = NA)

draw.circle(35,resultsV[35,1],2.2/5,border="firebrick1",lwd=4)

axis(1) ; axis(2)

mtext(expression(t[n+1]^2),side=2,  line=1.8, cex=1.3)
mtext("BSSH",side=3,  line=0.2, cex=2)

#################################################################################################

par(fig=c(0,3/4,9/17+0.02,11/17+0.01), oma=c(1.5,1,1,1), mar=c(2,2.5,1,1.5), new=T)

plot(resultsVI[,1],type="b",lty=1,xlab="",ylab="",col="dodgerblue3", lwd=1,
     pch=16, axes=F, yaxs = "i", ylim=c(-0.5, 9), cex=2)
points(c(4:34,35.5), resultsVI[4:35,2],type="l",lty=1,lwd=2, col="firebrick2")

polygon( c( 3.25,5:34,35.5, 35.5, 34:5,3.25 ), c(resultsVI[4:35,2],rep(0,32)),
         col = rgb(0.05,0.95,0.05,0.2), border = NA)

axis(1) ; axis(2)

mtext("x",  side=2,  at=9.8, cex=1.5)
mtext(expression(t[n+1]^2),side=2,  at=9.8,  line=1.8, cex=1.3)
mtext("Time",side=1,  line=2.3, cex=1.3)

#################################################################################################

par(fig=c(0,3/4,11/17+0.01,13/17+0.01), oma=c(1,1,1.5,1), mar=c(1,2.5,2.5,1.5), new=T)

plot(resultsVI[,1],type="b",lty=1,xlab="",ylab="",col="dodgerblue3", lwd=1,
     pch=16, axes=F, yaxs = "i", ylim=c(70, 74), cex=2)
points(resultsVI[,2],type="l",lty=1,lwd=2, col="firebrick2")
polygon( c( 3.25,5:34,35.5, 35.5, 34:5,3.25 ), c(resultsVI[4:35,2],rep(0,32)),
         col = rgb(0.05,0.95,0.05,0.2), border = NA)

axis(2, at=c(70,72,74), cex.axis=0.9)
mtext("FSSH",side=3,  line=0.2, cex=2)

#################################################################################################

par(fig=c(3/4,1,11/17+0.05,13/17+0.05), oma=c(0,1,1.5,1), mar=c(0,2.5,4,1.5), new=T)

plot(Z1,type="b",lty=1,xlab="",ylab="", lwd=1,
     pch=16, axes=F, yaxs = "i", ylim=c(-1.1, 3.1),
     col=c(rep(1,34),"firebrick1"), cex=c(rep(1,34),1.5) )
axis(1, cex.axis=0.8) ; axis(2, cex.axis=0.9)

points(postmeanV[1:35,1], col = rgb(107, 182, 255,maxColorValue = 255),
       type="l", lty=2, lwd=2)

mtext(expression(z[1]),side=2,  line=2.1, cex=1.3)

#################################################################################################

par(fig=c(3/4,1,9/17+0.05,11/17+0.05), oma=c(0,1,1.5,1), mar=c(0,2.5,4,1.5), new=T)

plot(Z2,type="b",lty=1,xlab="",ylab="", lwd=1,
     pch=16, axes=F, yaxs = "i", ylim=c(-1.1, 3.1),
     col=c(rep(1,34),"firebrick1"), cex=c(rep(1,34),1.5) )
axis(1, cex.axis=0.8) ; axis(2, cex.axis=0.9)

points(postmeanV[1:35,2], col = rgb(107, 182, 255,maxColorValue = 255),
       type="l", lty=2, lwd=2)

mtext(expression(z[2]),side=2,  line=2.1, cex=1.3)
mtext("Time",side=1,  line=2.3, cex=1.3)

#################################################################################################

par(fig=c(0,3/4,0,9/17), oma=c(1.5,1,1.5,1), mar=c(2.5,2.5,3,1.5), new=T)

plot(NA, NA,lty=1,xlab="",ylab="",col=c(rep(1,34),"firebrick1"), lwd=1, pch=16, axes=F,
     yaxs = "i", ylim=c(-3,3), xlim=c(-2.775,3))

segments(c(-3,-2,-1,0,1,2,3),-3,c(-3,-2,-1,0,1,2,3),3, col="gray75")
segments(-3,c(-3,-2,-1,0,1,2,3),3,c(-3,-2,-1,0,1,2,3), col="gray75")

points(Z1, Z2,lty=1,xlab="",ylab="",col=c(rep(1,34),"firebrick1"), cex=c(rep(2,34),3), lwd=1, pch=16,
       yaxs = "i", ylim=c(-1,3), xlim=c(-1,3))

axis(1) ; axis(2)

mtext(expression(z[2]),side=2,  line=2, cex=1.4)
mtext(expression(z[1]),side=1,  line=2, cex=1.4)

#################################################################################################
# PRIOR/POSTERIOR CONTOURS: MEAN VECTOR (mu1,mu2)
#################################################################################################

par(fig=c(3/4,1,6/17,9/17), oma=c(1.5,1,1.5,1), mar=c(2.5,2.5,3,1.5), new=T)

plot(NA, NA,lty=1,xlab="",ylab="",col=c(rep(1,34),2), lwd=1, pch=16, axes=F,
     yaxs = "i", ylim=c(-1,1), xlim=c(-0.93,1))

segments(c(-1,-0.5,0,0.5,1),-1.1,c(-1,-0.5,0,0.5,1),1, col="gray75")
segments(-1.1,c(-1,-0.5,0,0.5,1),1,c(-1,-0.5,0,0.5,1), col="gray75")

axis(1, at=c(-1,-0.5,0,0.5,1), cex.axis=0.9) ; axis(2, c(-1,-0.5,0,0.5,1), cex.axis=0.9)

x = y = seq(-0.8,1,by=0.005)

mu    <- mu0_plot
sigma <- Psi0_plot/(lambda0*nu_star0)
f     <- function(x, y) dmt(cbind(x, y), mu, sigma, df=nu_star0)
z     <- outer(x, y, f)
contour(x, y, z, nlevels = 8, drawlabels = FALSE, add = T,
        col = rgb(255,218,79,maxColorValue = 255), lty=1, lwd=2)

mu    <- mu_post_plot
sigma <- Psi_post_plot/(lambda_post*nu_star_post)
f     <- function(x, y) dmt(cbind(x, y), mu, sigma, df=nu_star_post)
z     <- outer(x, y, f)
contour(x, y, z, nlevels = 8, drawlabels = FALSE, add = T,
        col = rgb(107, 182, 255,maxColorValue = 255), lty=1, lwd=2)

mtext(expression(mu[2]),side=2,  line=2.1, cex=1.1)
mtext(expression(mu[1]),side=1,  line=2.1, cex=1.1)

#################################################################################################

par(fig=c(3/4,1,3/17,6/17), oma=c(1.5,1,1.5,1), mar=c(2.5,2.5,3,1.5), new=T)

plot(NA, NA,lty=1,xlab="",ylab="",col=c(rep(1,34),2), lwd=1, pch=16, axes=F, yaxs = "i",
     ylim=c(0,2.5), xlim=c(0.06,2.5))

segments(c(0.5,1,1.5,2,2.5),0,c(0.5,1,1.5,2,2.5),2.5, col="gray75")
segments(-0.1,c(0.5,1,1.5,2,2.5),2.5,c(0.5,1,1.5,2,2.5), col="gray75")

axis(1, at=c(0,0.5,1,1.5,2,2.5), cex.axis=0.9) ; axis(2, at=c(0,0.5,1,1.5,2,2.5), cex.axis=0.9)

x = y = seq(0.2,2.5,by=0.01)

sigma0 <- Psi0_plot
sigma1 <- Psi_post_plot

ndraw <- 500000
rho0 <- numeric(ndraw)
rho1 <- numeric(ndraw)

v_prior <- matrix(NA_real_, ndraw, 2)
v_post  <- matrix(NA_real_, ndraw, 2)

set.seed(1)
for (i in 1:ndraw) {
  A=rinvwishart(nu0, sigma0)
  B=rinvwishart(nu_post, sigma1)
  rho0[i] <- A[2,1]/sqrt(A[1,1]*A[2,2])
  rho1[i] <- B[2,1]/sqrt(B[1,1]*B[2,2])
  v_prior[i,] <- diag(A)
  v_post[i,]  <- diag(B)
}

# prior IW (df from prior nu0)
sigma <- Psi0_plot
v <- v_prior
z <- MASS::kde2d(v[,1], v[,2],
                 h=2.5*c(MASS::bandwidth.nrd(v[,1]),
                       MASS::bandwidth.nrd(v[,2])),
                 n=length(x), lims=c(range(x),range(y)))$z

contour(x, y, z, nlevels = 8, drawlabels = FALSE, add = T,
        col = rgb(255,218,79,maxColorValue = 255), lty=1, lwd=2)

# posterior IW (df from posterior nu_post)
sigma <- Psi_post_plot
v <- v_post
z <- MASS::kde2d(v[,1], v[,2],
                 h=2.5*c(MASS::bandwidth.nrd(v[,1]),
                       MASS::bandwidth.nrd(v[,2])),
                 n=length(x), lims=c(range(x),range(y)))$z

contour(x, y, z, nlevels = 8, drawlabels = FALSE, add = T,
        col = rgb(107,182,255,maxColorValue = 255), lty=1, lwd=2)

mtext(expression(sigma[2]^2),side=2,  line=2.1, cex=1.1)
mtext(expression(sigma[1]^2),side=1,  line=2.5, cex=1.1)

#################################################################################################

par(fig=c(3/4,1,0/17,3/17), oma=c(1.5,1,1.5,1), mar=c(2.5,2.5,3,1.5), new=T)

plot(NA, NA,lty=1,xlab="",ylab="",col=c(rep(1,34),2), lwd=1, pch=16, axes=F, yaxs = "i",
     ylim=c(0,3.5), xlim=c(-0.45,1))

polygon( density(rho0, adjust=1.3), col=rgb(255,218,79,63.75, maxColorValue = 255), border = NA)
polygon( density(rho1, adjust=1.3), col=rgb(107, 182, 255,63.75, maxColorValue = 255), border = NA)

points(density(rho0, adjust=1.3), type = "l", lwd=2, col=rgb(255,218,79,maxColorValue = 255))
points(density(rho1, adjust=1.3), type = "l", lwd=2, col = rgb(107, 182, 255,maxColorValue = 255))

axis(1, at=c(-0.5,0,0.5,1), cex.axis=0.9) ; axis(2, at=c(0,1,2,3), cex.axis=0.9)

mtext(expression(rho),side=1,  line=2.1, cex=1.1)
mtext("density", side=2,  line=2.1, cex=1.1)

#################################################################################################

par(fig=c(3/4,1,14.5/17,16.5/17),oma=c(0,0,0.2,0), mar=c(0,0,0.2,0), new=T)
plot.new()

Legend=c("OOC observation","prior","posterior",expression(mu[n]))

legend(x="center",y=1, inset=0, legend = Legend , cex=1.25,
       col=c("firebrick1", rgb(255,218,79, maxColorValue = 255),
             rgb(107, 182, 255, maxColorValue = 255),
             rgb(107, 182, 255, maxColorValue = 255)), bg="white",bty ="n", y.intersp = 1.2,
       pch = c(16,-1,-1,-1), pt.cex=2, lty=c(0,1,1,2), lwd=2)

#################################################################################################

par(fig=c(3/4,1,13/17,15/17),oma=c(0,0,0.2,0), mar=c(0,0,0.2,0), new=T)
plot.new()

text(0.45,0.6,substitute(paste(bold("Significant components:"))), cex=1.25 )
text(0.45,0.4, sig_expr, cex=1.25)

dev.off()

cat("\nPDF saved to: BSSH_application-z-scores.pdf\n")