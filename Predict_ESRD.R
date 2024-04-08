
Surv.fit.fun <- function(dat, para) {
  survival <- vector("list", unique(dat$A0))
  for (a in unique(dat$A0)) {
	h0 <- with(subset(dat, Group==1 & A0==a), Y)
	h1 <- with(subset(dat, Group %in% c(2,3)  & A0==a), W)
	h2 <- with(subset(dat, Group==2  & A0==a), Y-W)
	survival[[a]] <- Predict.fun(dat, para, unique(h0, h1+h2), a)
  }
  return(survival)
}

Predict.fun <- function(dat, para, t0, a, h0, h1, h2) {
		 x1 <- subset(dat, A0==a, select=c("x1"))
		 m0 <- length(t0)
     n <- length(x1)
     X <- cbind(rep(1, n), rep(a,n), x1)
     Xd <- cbind(rep(a,n), x1)
     Xe <- cbind(rep(a,n), x1)
     Exbeta0 <- exp(Xd%*%para$beta0vec)
     Exbeta1 <- exp(Xe%*%para$beta1vec)
     Exalp <- exp(X%*%para$alpvec)
     h0 <- with(subset(dat, Group==1), Y)
     h1 <- with(subset(dat, Group %in% c(2,3)), W)
     h2 <- with(subset(dat, Group==2), Y-W)
     n0 <- length(h0)
     n1 <- length(h1)
     n2 <- length(h2)
     indH0 <- (matrix(h0, n0, m0, byrow=FALSE) >= matrix(t0, n0, m0, byrow=TRUE))
     H0 <- para$h0%*%indH0
     
     S1 <- ColSum(exp(-Exbeta0%*%(para$h0%*%indH0))*matrix(1/(1+Exalp), n, m0, byrow=FALSE))
     
     
     PE0 <- 1/(1+Exalp)
     PE1 <- 1-PE0     

  return(dat)
}                                              
  