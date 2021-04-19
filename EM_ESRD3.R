
Estep <- function(dat, oldpara, indH0, indH1, indH2) {
  n <- dim(dat)[1]
  X <- with(dat, cbind(rep(1, n), A0, x1))
  Xd <- with(dat, cbind(A0, x1))
  Xe <- with(dat, cbind(A0, x1))
  Xbeta0 <- Xd%*%oldpara$beta0vec
  Xbeta1 <- Xe%*%oldpara$beta1vec
  Xalp <- X%*%oldpara$alpvec
  Exbeta0 <- exp(Xbeta0)
  Exbeta1 <- exp(Xbeta1)
  Exalp <- exp(Xalp)
  H0 <- indH0%*%oldpara$h0
  H1 <- indH1%*%oldpara$h1
  HExbeta0 <- H0*Exbeta0
  HExbeta1 <- H1*Exbeta1
  PE0 <- 1/(1+Exalp)
  PE1 <- 1-PE0
  PE1.g4 <- PE1*exp(-HExbeta1)/(PE0*exp(-HExbeta0)+PE1*exp(-HExbeta1))   
  dat$Ee <- with(dat, ifelse(Group==1, 0, ifelse(Group %in% c(2,3), 1, PE1.g4)))

  return(dat)
}                                              
  
###################################################################################
  
Mstep <- function(dat, oldpara,indH0, indH1, indH2) {
     n <- dim(dat)[1]
     ind <- with(dat, Group %in% c(2,3))
     X <- with(dat, cbind(rep(1, n), A0, x1))
     Xd <- with(dat, cbind(A0, x1))
     Xe <- with(dat, cbind(A0, x1))
     Xg <- with(dat[ind,], cbind(A0, V*(1-A0), W, x1, z))
     Xg2 <- with(subset(dat, Group==2), cbind(A0, V*(1-A0), W, x1, z))
     pd <- dim(Xd)[2]
     pe <- dim(Xe)[2]
     pg <- dim(Xg)[2]
     n0 <- sum(dat$Group==1)
     n1 <- sum(dat$Group %in% c(2,3))       
     n2 <- sum(dat$Group==2)
     dl2beta0 <- matrix(0, pd, pd)
     dl2beta1 <- matrix(0, pe, pe)
     dl2beta2 <- matrix(0, pg, pg)


     PE1 <- 1/(1+exp(-X%*%oldpara$alpvec))  
     newalpvec <- oldpara$alpvec + solve((t(X)*matrix(PE1*(1-PE1),dim(X)[2],n,byrow=TRUE))%*%X)%*%colSums(matrix(dat$Ee-PE1, n, 3, byrow=FALSE)*X)   
     
     b0t1 <- as.vector(exp(Xd%*%oldpara$beta0vec))*indH0*(1-dat$Ee)
     for (k1 in 1:pd) {
     	for (k2 in 1:k1) {
     		dl2beta0[k1,k2] <- sum(-t(b0t1)%*%(Xd[,k1]*Xd[,k2])/(colSums(b0t1))+(t(b0t1)%*%Xd[,k1])*(t(b0t1)%*%Xd[,k2])/(colSums(b0t1))^2)
     		dl2beta0[k2,k1] <- dl2beta0[k1,k2]
     	}
     }		
     dl1beta0 <- colSums(Xd[dat$Group==1,]-t(b0t1)%*%Xd/(colSums(b0t1)))
  
     b1t1 <- as.vector(exp(Xe%*%oldpara$beta1vec))*indH1*dat$Ee
     for (k1 in 1:pe) {
     	for (k2 in 1:k1) {
     		dl2beta1[k1,k2] <- sum(-t(b1t1)%*%(Xe[,k1]*Xe[,k2])/(colSums(b1t1))+(t(b1t1)%*%Xe[,k1])*(t(b1t1)%*%Xe[,k2])/(colSums(b1t1))^2)
     		dl2beta1[k2,k1] <- dl2beta1[k1,k2]
     	}
     }		
     dl1beta1 <- colSums(Xe[dat$Group %in% c(2,3),]-t(b1t1)%*%Xe/(colSums(b1t1)))
  
     b2t1 <- as.vector(exp(Xg%*%oldpara$beta2vec))*indH2
     for (k1 in 1:pg) {
     	for (k2 in 1:k1) {
     		dl2beta2[k1,k2] <- sum(-t(b2t1)%*%(Xg[,k1]*Xg[,k2])/(colSums(b2t1))+(t(b2t1)%*%Xg[,k1])*(t(b2t1)%*%Xg[,k2])/(colSums(b2t1))^2)
     		dl2beta2[k2,k1] <- dl2beta2[k1,k2]
     	}
     }		
     dl1beta2 <- colSums(Xg2-t(b2t1)%*%Xg/(colSums(b2t1)))
         
     newbeta0vec <- oldpara$beta0vec - solve(dl2beta0)%*%dl1beta0
     newbeta1vec <- oldpara$beta1vec - solve(dl2beta1)%*%dl1beta1
     newbeta2vec <- oldpara$beta2vec - solve(dl2beta2)%*%dl1beta2        
	   newh0 <- 1/colSums(indH0*matrix(exp(Xd%*%oldpara$beta0vec)*(1-dat$Ee), n, n0, byrow=FALSE)) 
	   newh1 <- 1/colSums(indH1*matrix(exp(Xe%*%oldpara$beta1vec)*dat$Ee, n, n1, byrow=FALSE)) 
	   newh2 <- 1/colSums(indH2*matrix(exp(Xg%*%oldpara$beta2vec), sum(ind), n2, byrow=FALSE)) 
          
     newpara <- list(beta0vec=newbeta0vec, h0=newh0, beta1vec=newbeta1vec, h1=newh1, beta2vec=newbeta2vec, h2=newh2, alpvec=newalpvec)
     return(newpara)
}
