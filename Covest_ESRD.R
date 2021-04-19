Covest <- function(dat, para,indH0, indH1, indH2) {
     n <- dim(dat)[1]
     ind0 <- with(dat, Group %in% c(2,3))
     X <- with(dat, cbind(rep(1, n), A0, x1))
     Xd <- with(dat, cbind(A0, x1))
     Xe <- with(dat, cbind(A0, x1))
     Xg <- with(dat[ind0,], cbind(A0, V*(1-A0), W, x1, z))
     Xbeta0 <- Xd%*%para$beta0vec
     Xbeta1 <- Xe%*%para$beta1vec
     Xbeta2 <- Xg%*%para$beta2vec
     Xalp <- X%*%para$alpvec
     pd <- dim(Xd)[2]
     pe <- dim(Xe)[2]
     pg <- dim(Xg)[2]
     Exbeta0 <- exp(Xbeta0)
     Exbeta1 <- exp(Xbeta1)
     Exbeta2 <- exp(Xbeta2)
     Exalp <- exp(Xalp)
     n0 <- length(para$h0)
     n1 <- length(para$h1)       
	   n2 <- length(para$h2)
	   H0 <- indH0%*%para$h0
     H1 <- indH1%*%para$h1
     H2 <- indH2%*%para$h2
     HExbeta0 <- H0*Exbeta0
     HExbeta1 <- H1*Exbeta1
     HExbeta2 <- H2*Exbeta2
     PE0 <- 1/(1+Exalp)
     PE1 <- 1-PE0
     ind <- with(dat, Group %in% c(2,3) & A0==0)
     
     nb0 <- length(para$beta0vec)
     nb1 <- length(para$beta1vec)
     nb2 <- length(para$beta2vec)
     na0 <- length(para$alpvec)
     cn0 <- nb0+n0
     cnb1 <- cn0+nb1
     cn1 <- cnb1+n1
     cnb2 <- cn1+nb2
     cn2 <- cnb2+n2
     ntot <- cn2+na0
     d2Q <- matrix(0, ntot, ntot)
     Ess <- matrix(0, ntot, ntot)
     
     # c(nb0, cn0, cnb1, cn1, cnb2, cn2, ntot)

     # d2(Q)/d(para)d(para)' matrix
     
     d2Q[1:nb0,1:nb0] <- -t(Xd)%*%diag(as.vector(HExbeta0*(1-dat$Ee)))%*%Xd
     d2Q[1:nb0,(nb0+1):cn0] <- -t(Xd)%*%diag(as.vector(Exbeta0*(1-dat$Ee)))%*%indH0
     d2Q[(nb0+1):cn0, 1:nb0] <- t(d2Q[1:nb0,(nb0+1):cn0])
     d2Q[(nb0+1):cn0, (nb0+1):cn0] <- -diag(1/(para$h0)^2)

     d2Q[(cn0+1):cnb1,(cn0+1):(cnb1)] <- -t(Xe)%*%diag(as.vector(HExbeta1*dat$Ee))%*%Xe
     d2Q[(cn0+1):cnb1,(cnb1+1):(cn1)] <- -t(Xe)%*%diag(as.vector(Exbeta1*dat$Ee))%*%indH1
     d2Q[(cnb1+1):cn1,(cn0+1):(cnb1)] <- t(d2Q[(cn0+1):(cnb1),(cnb1+1):cn1])
     d2Q[(cnb1+1):cn1, (cnb1+1):(cn1)] <- -diag(1/(para$h1)^2)

     d2Q[(cn1+1):cnb2,(cn1+1):cnb2] <- -t(Xg)%*%diag(as.vector(HExbeta2))%*%Xg
     d2Q[(cn1+1):cnb2,(cnb2+1):cn2] <- -t(Xg)%*%diag(as.vector(Exbeta2))%*%indH2
     d2Q[(cnb2+1):cn2,(cn1+1):cnb2] <- t(d2Q[(cn1+1):(cnb2),(cnb2+1):cn2])
     d2Q[(cnb2+1):(cn2), (cnb2+1):(cn2)] <- -diag(1/(para$h2)^2)

     d2Q[(cn2+1):(ntot),(cn2+1):(ntot)] <- -t(X)%*%diag(as.vector(PE1*PE0))%*%X

     # E(S_para S'_para | D_obs)
     
     indY0 <- (matrix(dat$Y, n, n0, byrow=FALSE) == matrix(para$h0, n, n0, byrow=TRUE)) 
     indY1 <- (matrix(dat$W, n, n1, byrow=FALSE) == matrix(para$h1, n, n1, byrow=TRUE)) 
     indY2 <- (matrix(with(dat, Y-W), n, n2, byrow=FALSE) == matrix(para$h2, n, n2, byrow=TRUE)) 
     indY1[dat$Group %in% c(1,4), ] <- 0
     indY2 <- indY2[ind0, ]
     h0.mat <- indY0*matrix(1/para$h0,n,n0,byrow=TRUE)
     h1.mat <- indY1*matrix(1/para$h1,n,n1,byrow=TRUE)
     h2.mat <- indY2*matrix(1/para$h2,sum(ind0),n2,byrow=TRUE)
     G <- dat$Group
     G1 <- (G==1)
     G2 <- (G==2)
     G2s <- with(subset(dat, Group %in% c(2,3)), Group==2)
     G23 <- (G %in% c(2,3))
     G14 <- (G %in% c(1,4))
     G234 <- (G %in% c(2,3,4))

     # beta0 by beta0 
     Ess[1:nb0,1:nb0] <- t(Xd[G1,])%*%Xd[G1,]-2*t(Xd[G1,])%*%(as.vector(HExbeta0)*Xd)[G1,]+
        t((as.vector(HExbeta0)*Xd)[G14,])%*%(as.vector((1-dat$Ee)*HExbeta0)*Xd)[G14,]
     # beta0 by h0 
     Ess[1:nb0,(nb0+1):cn0] <- t(Xd[G1,])%*%h0.mat[G1,]-t(Xd[G1,])%*%(as.vector(Exbeta0)*indH0)[G1,]-
        t((as.vector(HExbeta0)*Xd)[G1,])%*%h0.mat[G1,]+
        t((as.vector(HExbeta0)*Xd)[G14,])%*%(as.vector((1-dat$Ee)*Exbeta0)*indH0)[G14,]
     # beta0 by alpha
     Ess[1:nb0,(cn2+1):ntot] <- t(Xd[G1,])%*%(as.vector(dat$Ee-PE1)*X)[G1,]+
        t((as.vector(HExbeta0)*Xd)[G14,])%*%(as.vector((1-dat$Ee)*PE1)*X)[G14,]
       
     # h0 by h0  
     Ess[(nb0+1):cn0,(nb0+1):cn0] <- t(h0.mat[G1,])%*%h0.mat[G1,]-2*t(h0.mat[G1,])%*%(as.vector(Exbeta0)*indH0)[G1,]+
        t((as.vector(Exbeta0)*indH0)[G14,])%*%(as.vector((1-dat$Ee)*Exbeta0)*indH0)[G14,]
     # h0 by alpha
     Ess[(nb0+1):cn0, (cn2+1):ntot] <- t(h0.mat[G1,])%*%(as.vector(dat$Ee-PE1)*X)[G1,]+
        t((as.vector(Exbeta0)*indH0)[G14,])%*%(as.vector((1-dat$Ee)*PE1)*X)[G14,]

     # beta1 by beta1 
     Ess[(cn0+1):(cnb1),(cn0+1):cnb1] <- t(Xe[G23,])%*%Xe[G23,]-
        2*t(Xe[G23,])%*%(as.vector(HExbeta1)*Xe)[G23,]+
        t((as.vector(HExbeta1)*Xe)[G234,])%*%(as.vector(dat$Ee*HExbeta1)*Xe)[G234,]
     # beta1 by h1 
     Ess[(cn0+1):(cnb1),(cnb1+1):(cn1)] <- t(Xe[G23,])%*%h1.mat[G23,]-
        t(Xe[G23,])%*%(as.vector(Exbeta1)*indH1)[G23,]-
        t((as.vector(HExbeta1)*Xe)[G23,])%*%h1.mat[G23,]+
        t((as.vector(HExbeta1)*Xe)[G234,])%*%(as.vector(dat$Ee*Exbeta1)*indH1)[G234,]
     # beta1 by beta2
     Ess[(cn0+1):(cnb1),(cn1+1):(cnb2)] <- t(Xe[G2,])%*%Xg[G2s,]-t((as.vector(HExbeta1)*Xe)[G2,])%*%Xg[G2s,]-
        t(Xe[G23,])%*%(as.vector(dat$Ee[ind0]*HExbeta2)*Xg)+
        t((as.vector(HExbeta1)*Xe)[G23,])%*%(as.vector(dat$Ee[ind0]*HExbeta2)*Xg)
     # beta1 by h2
     Ess[(cn0+1):(cnb1),(cnb2+1):(cn2)] <- t(Xe[G2,])%*%h2.mat[G2s,]-
        t(Xe[G23,])%*%(as.vector(Exbeta2)*indH2)-
        t((as.vector(HExbeta1)*Xe)[G2,])%*%h2.mat[G2s,]+
        t((as.vector(HExbeta1)*Xe)[G23,])%*%(as.vector(dat$Ee[ind0]*Exbeta2)*indH2)
     # beta1 by alpha
     Ess[(cn0+1):(cnb1),(cn2+1):(ntot)] <- t(Xe[G23,])%*%(as.vector(dat$Ee-PE1)*X)[G23,]-
        t((as.vector(HExbeta1)*Xe)[G234,])%*%(as.vector(dat$Ee*(1-PE1))*X)[G234,]

     # h1 by h1
     Ess[(cnb1+1):(cn1),(cnb1+1):(cn1)] <- t(h1.mat[G23,])%*%h1.mat[G23,]-
        2*t(h1.mat[G23,])%*%(as.vector(Exbeta1)*indH1)[G23,]+
        t((as.vector(Exbeta1)*indH1)[G234,])%*%(as.vector(dat$Ee*Exbeta1)*indH1)[G234,]
     # h1 by beta2
     Ess[(cnb1+1):(cn1), (cn1+1):(cnb2)] <- t(h1.mat[G2,])%*%Xg[G2s,]-
        t(h1.mat[G23,])%*%(as.vector(dat$Ee[ind0]*HExbeta2)*Xg)-
        t((as.vector(dat$Ee*Exbeta1)*indH1)[G2,])%*%Xg[G2s,]+
        t((as.vector(Exbeta1)*indH1)[G23,])%*%(as.vector(dat$Ee[ind0]*HExbeta2)*Xg) 
     # h1 by h2
     Ess[(cnb1+1):(cn1),(cnb2+1):(cn2)] <- t(h1.mat[G2,])%*%h2.mat[G2s,]-
        t(h1.mat[G23,])%*%(as.vector(dat$Ee[ind0]*Exbeta2)*indH2)-
        t((as.vector(dat$Ee*Exbeta1)*indH1)[G2,])%*%h2.mat[G2s,]+
        t((as.vector(Exbeta1)*indH1)[G23,])%*%(as.vector(dat$Ee[ind0]*Exbeta2)*indH2) 
     # h1 by alpha
     Ess[(cnb1+1):(cn1), (cn2+1):(ntot)] <- t(h1.mat[G23,])%*%(as.vector(dat$Ee-PE1)*X)[G23,]-
        t((as.vector(Exbeta1)*indH1)[G234,])%*%(as.vector(dat$Ee*(1-PE1))*X)[G234,]

     # beta2 by beta2
     Ess[(cn1+1):(cnb2),(cn1+1):(cnb2)] <- t(Xg[G2s,])%*%Xg[G2s,]-2*t(Xg[G2s,])%*%(as.vector(dat$Ee[ind0]*HExbeta2)*Xg)[G2s,]+
        t((as.vector(HExbeta2)*Xg))%*%(as.vector(dat$Ee[ind0]*HExbeta2)*Xg)
     # beta2 by h2
     Ess[(cn1+1):(cnb2),(cnb2+1):(cn2)] <- t(Xg[G2s,])%*%h2.mat[G2s,]-t(Xg[G2s,])%*%(as.vector(dat$Ee[ind0]*Exbeta2)*indH2)[G2s,]-
        t((as.vector(dat$Ee[ind0]*HExbeta2)*Xg)[G2s,])%*%h2.mat[G2s,]+
        t((as.vector(HExbeta2)*Xg))%*%(as.vector(dat$Ee[ind0]*Exbeta2)*indH2)
     # beta2 by alpha
     Ess[(cn1+1):(cnb2),(cn2+1):(ntot)] <- t(Xg[G2s,])%*%(as.vector(dat$Ee-PE1)*X)[G2,]-
        t((as.vector(HExbeta2)*Xg))%*%(as.vector(dat$Ee*(1-PE1))*X)[G23,]

     # h2 by h2
     Ess[(cnb2+1):(cn2),(cnb2+1):(cn2)] <- t(h2.mat[G2s,])%*%h2.mat[G2s,]-2*t(h2.mat[G2s,])%*%(as.vector(Exbeta2)*indH2)[G2s,]+
        t((as.vector(Exbeta2)*indH2))%*%(as.vector(dat$Ee[ind0]*Exbeta2)*indH2)
     # h2 by alpha
     Ess[(cnb2+1):(cn2), (cn2+1):(ntot)] <- t(h2.mat[G2s,])%*%(as.vector(dat$Ee-PE1)*X)[G2,]-
        t((as.vector(Exbeta2)*indH2))%*%(as.vector(1-PE1)*X)[G23,]

     # alpha by alpha
     Ess[(cn2+1):(ntot), (cn2+1):(ntot)] <- t(X)%*%(as.vector(dat$Ee-2*PE1*dat$Ee+PE1^2)*X) 
     
     Ess[(nb0+1):ntot, 1:nb0] <- t(Ess[1:nb0, (nb0+1):ntot])
     Ess[(cn0+1):ntot, (nb0+1):cn0] <- t(Ess[(nb0+1):cn0, (cn0+1):ntot])
     Ess[(cnb1+1):ntot, (cn0+1):cnb1] <- t(Ess[(cn0+1):cnb1, (cnb1+1):ntot])
     Ess[(cn1+1):ntot, (cnb1+1):(cn1)] <- t(Ess[(cnb1+1):(cn1), (cn1+1):ntot])
     Ess[(cnb2+1):ntot, (cn1+1):(cnb2)] <- t(Ess[(cn1+1):(cnb2), (cnb2+1):ntot])
     Ess[(cn2+1):ntot, (cnb2+1):(cn2)] <- t(Ess[(cnb2+1):(cn2), (cn2+1):ntot])

     tmp1 <- matrix(0, n, nb2)
     tmp2 <- matrix(0, n, n2)
     tmp1[G %in% c(2,3),] <- - as.vector(HExbeta2*dat$Ee[ind0])*Xg
     tmp1[G2, ] <- Xg[G2s,] 
		 tmp2[G %in% c(2,3),] <- - as.vector(Exbeta2*dat$Ee[ind0])*indH2
		 tmp2[G2, ] <- indY2[G2s,]
       
		 Es <- cbind((G1)*Xd-(G%in%c(1,4))*as.vector(HExbeta0*(1-dat$Ee))*Xd, 
		   (G1)*indY0-(G%in%c(1,4))*as.vector(Exbeta0*(1-dat$Ee))*indH0, 
		   (G23)*Xe-(G%in%c(2,3,4))*as.vector(HExbeta1*dat$Ee)*Xe, 
		   (G23)*indY1-(G%in%c(2,3,4))*as.vector(Exbeta1*dat$Ee)*indH1,
		    tmp1, 
		    tmp2,
        as.vector(dat$Ee-PE1)*X)

     tmp <- sqrt(diag(solve(-d2Q-(Ess-t(Es)%*%Es))))

     # tmp <- solve(-d2Q-(Ess-t(Es)%*%Es))

     sdest <- list(beta0sd=tmp[1:nb0], h0sd=tmp[(nb0+1):cn0], beta1sd=tmp[(cn0+1):cnb1], h1sd=tmp[(cnb1+1):cn1], 
        beta2sd=tmp[(cn1+1):cnb2], h2sd=tmp[(cnb2+1):cn2], alpsd=tmp[(cn2+1):ntot])
        
     return(sdest)   
}

