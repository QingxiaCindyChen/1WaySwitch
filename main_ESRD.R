
source("EM_ESRD3.R")
source("Covest_ESRD.R")

library(Design)
	
nsim <- 1000
n <- 1000
set.seed(1234767)

# set parameters
beta0vec <- c(-1, 1)
beta1vec <- c(-1, 1)
beta2vec <- c(-5, 0.5, 1, 1, 1)
alpvec <- c(-0.5, 1, 1)
lambvec <- c(-0.5, -1, 1)
tau <- 3
istep <- 10
npara <- length(beta0vec)+length(beta1vec)+length(beta2vec)+length(alpvec)
parameters <- matrix(0, istep+1, npara)

A0 <- V <- X1 <- Td <- Te <- Y <- W <- DeltaD <- DeltaE <- TGap <- Cens <- E <- Omega <- NULL

Beta0Vec <- matrix(NA, nsim, length(beta0vec))
Beta1Vec <- matrix(NA, nsim, length(beta1vec))
Beta2Vec <- matrix(NA, nsim, length(beta2vec))
AlpVec <- matrix(NA, nsim, length(alpvec))
Beta0Sd <- matrix(NA, nsim, length(beta0vec))
Beta1Sd <- matrix(NA, nsim, length(beta1vec))
Beta2Sd <- matrix(NA, nsim, length(beta2vec))
AlpSd <- matrix(NA, nsim, length(alpvec))
Beta0CP <- matrix(NA, nsim, length(beta0vec))
Beta1CP <- matrix(NA, nsim, length(beta1vec))
Beta2CP <- matrix(NA, nsim, length(beta2vec))
AlpCP <- matrix(NA, nsim, length(alpvec))
censE <- matrix(NA, nsim, 1)
censD <- matrix(NA, nsim, 1)
G <- matrix(NA, nsim, 4)
seed.mat <- matrix(NA, nsim*0.5, 626)
seed.ind <- 1

for (m in 1:nsim) {

 # simulate data and saved them into simdat
  seed.mat[seed.ind,] <- .Random.seed 
      
  A0 <-  rbinom(n,1,0.5)
  x1 <- rnorm(n, 1, 1)
  lagt <- rexp(n, 1)
  z <- rnorm(n, 1, 2)
  X <- cbind(rep(1,n), A0, x1)
  Xd <- cbind(A0, x1)
  Xe <- cbind(A0, x1)
  E <- rbinom(n, 1, 1/(1+exp(-X%*%alpvec)))
  Te <- -log(runif(n))*exp(-Xe%*%beta1vec) 
  # assuming the baseline hazard function is exponential(1) with density function exp(-t)
  Cens <- 2*runif(n)+2
  Xv <- cbind(rep(1,n), Te, x1)
  V <- ifelse(A0==0 & E==1 & Te <=Cens, rbinom(n, 1, 1/(1+exp(-Xv%*%lambvec))), 0)
  Xg <- cbind(A0, V*(1-A0), Te, x1, z)  
  TGap <- -log(runif(n))*exp(-Xg%*%beta2vec)
  # assuming the baseline hazard function is exponential(1) with density function exp(-t)
  Td <- ifelse(E==0, -log(runif(n))*exp(-Xd%*%beta0vec), Te+TGap)
  Group <- ifelse(E==0 & Td <= Cens, 1, 
             ifelse(E==1 & Te <= Cens & Td <= Cens, 2, 
                ifelse(E==1 & Te <= Cens & Td > Cens, 3, 4)))
  DeltaD <- ifelse(Td <= Cens, 1, 0) 
  DeltaE <- ifelse(Te <= Cens, 1, 0)
  Y <- ifelse(DeltaD==1, Td, Cens)
  W <- ifelse(DeltaE==1, Te, Cens)
  W <- ifelse(Group %in%c(1,4), Cens, W)
  # We don't know if patients in Group 1 or 4 would experience ESRD if they were not censored or dead
  DeltaE[Group==1] <- NA 
  DeltaE[Group==4] <- NA
  censE[m] <- sum(DeltaE, na.rm=TRUE)/n
  censD[m] <- sum(DeltaD, na.rm=TRUE)/n
  G[m,] <- table(Group)/n
  W[Group %in% c(1,4)] <- NA
 
  simdat <- data.frame(A0, V, x1, z, E, DeltaD, DeltaE, Y, W, Group) # all the observed information

  h0 <- with(subset(simdat, Group==1), Y)
  h1 <- with(subset(simdat, Group %in% c(2,3)), W)
  h2 <- with(subset(simdat, Group==2), Y-W)
  n0 <- length(h0)
  n1 <- length(h1)
  n2 <- length(h2)
  indH0 <- (matrix(simdat$Y, n, n0, byrow=FALSE) >= matrix(h0, n, n0, byrow=TRUE))
  indH1 <- (matrix(simdat$W, n, n1, byrow=FALSE) >= matrix(h1, n, n1, byrow=TRUE))
  indH2 <- (matrix(with(simdat, Y-W), n, n2, byrow=FALSE) >= matrix(h2, n, n2, byrow=TRUE))
  indH1[Group %in% c(1,4), ] <- 0
  indH2 <- indH2[Group %in% c(2,3), ]

  oldh0 <- rep(1/n0, n0)
  oldh1 <- rep(1/n1, n1)
  oldh2 <- rep(1/n2, n2)

  oldbeta0vec <- c(0, 0)
  oldbeta1vec <- c(0, 0)
  oldbeta2vec <- c(0, 0, 0, 0, 0)
  oldalpvec <- c(0, 0, 0) 
  oldpara <- list(beta0vec=oldbeta0vec, h0=oldh0, beta1vec=oldbeta1vec, h1=oldh1, beta2vec=oldbeta2vec, h2=oldh2, alpvec=oldalpvec)

  epsilon <- 0.001
  epsilon2 <- 1.0e-10
  maxiter <- 1000
  absdiff <- 1
  diffdist <- 1
  iter <- 0
  diffvec <- rep(1, length(oldpara))
  vec <- c(3,5,7)
  
  #cat(c(beta0vec, beta1vec, beta2vec, ntheta, alpvec, lambvec), "\n") 
  cat("sim=", m, "\n")

  while((absdiff > epsilon | diffdist > epsilon2) & iter < maxiter) {
     simdat <- Estep(simdat, oldpara, indH0, indH1, indH2)
     newpara <- Mstep(simdat, oldpara,indH0, indH1, indH2)
     for (k in 1:istep) parameters[k,] <- parameters[k+1,]
     temp <- oldpara[[1]]
     for (k in vec) temp <- c(temp, oldpara[[k]])
     parameters[istep+1,] <- temp 
     diffdist <- sum((parameters[istep+1,]-parameters[1,])^2)
     if (iter < istep+1) diffdist <- 1
     for (k in 1:length(oldpara)) diffvec[k] <- max(abs(oldpara[[k]]-newpara[[k]]))
	   absdiff <- max(diffvec)
	   oldpara <- newpara
	   # cat(c(iter, absdiff, diffdist), "\n")
	   # cat(c(oldpara$beta0vec, oldpara$beta1vec, oldpara$beta2vec, oldpara$alpvec), "\n") 
	   iter <- iter+1
	}
  # cat(c(oldpara$beta0vec, oldpara$beta1vec, oldpara$beta2vec, oldpara$alpvec), "\n")   
  sdest <- Covest(simdat, oldpara, indH0, indH1, indH2)
  surv.fit <- Predict(simdat, oldpara, h0, h1, h2)  
	
  Beta0Vec[m, ] <- newpara$beta0vec	
  Beta1Vec[m, ] <- newpara$beta1vec	
  Beta2Vec[m, ] <- newpara$beta2vec	
  AlpVec[m, ] <-   newpara$alpvec	
  Beta0Sd[m, ] <- sdest$beta0sd	
  Beta1Sd[m, ] <- sdest$beta1sd	
  Beta2Sd[m, ] <- sdest$beta2sd	
  AlpSd[m, ] <-   sdest$alpsd	

}

Beta0CP <- ifelse(abs((Beta0Vec-matrix(beta0vec, nsim, length(beta0vec), byrow=TRUE))/Beta0Sd) <= 1.96, 1, 0)
Beta1CP <- ifelse(abs((Beta1Vec-matrix(beta1vec, nsim, length(beta1vec), byrow=TRUE))/Beta1Sd) <= 1.96, 1, 0)
Beta2CP <- ifelse(abs((Beta2Vec-matrix(beta2vec, nsim, length(beta2vec), byrow=TRUE))/Beta2Sd) <= 1.96, 1, 0)
AlpCP <- ifelse(abs((AlpVec-matrix(alpvec, nsim, length(alpvec), byrow=TRUE))/AlpSd) <= 1.96, 1, 0)

sim.data <- data.frame(cbind(seq(1,nsim), censD, censE, Beta0Vec, Beta1Vec, Beta2Vec, AlpVec, Beta0Sd, 
 Beta1Sd, Beta2Sd, AlpSd, Beta0CP, Beta1CP, Beta2CP, AlpCP, G))
save(sim.data, file="sim_result_ESRD_Trans.Rdata")
