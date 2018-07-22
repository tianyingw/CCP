Allint_logistic<-function(y, W_int, C=NULL, print.summary=TRUE, standardize = TRUE){

  #  W_int is matrix containing replicates, C=vector of cut points

  #####################################
  # check dimensions and replicates
  #####################################

  n <- length(W_int[, 1])
  k <- length(W_int[1, ])

  ####################################
  # standardize the rwo data W
  ####################################

  rowmean_w <- apply(W_int, 1, mean)
  mux_hat <- mean(rowmean_w)
  if(standardize == TRUE){
    a0 <- 1 / sqrt(mean((rowmean_w - mux_hat)^2))
    b0 <- - mux_hat * a0
      }else{
        a0 <- 1
        b0 <- 0}
  ww <- a0 * W_int + b0

  mean_w<-apply(ww,1,mean)
  w<-mean_w
  s2w<-apply(ww,1,var)
  su2_hat<-mean(s2w)
  su2e<- su2_hat/k
  mux_hat<-mean(w)
  s2x_hat<-max(mean((mean_w-mux_hat)^2)-su2e,0.2*(mean((mean_w-mux_hat)^2)))
  lambda_hat=max(s2x_hat/(s2x_hat+su2e),0.2)
  Lambda.par=c(mux_hat,s2x_hat)

  ###########################################
  #Cut points
  if(is.null(C)){
    C=rep(0,4)
    C[1]=qnorm(0.2,mean =mux_hat , sd = sqrt(s2x_hat))
    C[2]=qnorm(0.4,mean =mux_hat , sd = sqrt(s2x_hat))
    C[3]=qnorm(0.6,mean =mux_hat , sd = sqrt(s2x_hat))
    C[4]=qnorm(0.8,mean =mux_hat , sd = sqrt(s2x_hat))
  }else{C = a0*C + b0}
  sizeC<-length(C)
  try(if(sizeC != 4) stop("the size of cutting points is different from 4!"))
  J=length(C)+1   #Number of sets
  ##################
  # Basic functions
  ##################
  fMx<-function(x){
    Mx<-vector()
    Mx[1]=ifelse(x<C[1],1,0)
    Mx[2]=ifelse((C[1]<=x)& (x<C[2]),1,0)
    Mx[3]=ifelse((C[2]<=x)&(x<C[3]),1,0)
    Mx[4]=ifelse((C[3]<=x)&(x<C[4]),1,0)
    Mx[5]=ifelse(x>=C[4],1,0)
    return (Mx)
  }
  fHM<-function(x,theta){
    1/(1+exp(-(fMx(x)%*%theta)))
  }
  fH1<-function(x){1/(1+exp(-x))*(1-1/(1+exp(-x)))}
  fH<-function(x){1/(1+exp(-x))}

  K=c(-Inf,C,Inf)
  TmpIntegrate1 <- function(x){
    n_W = length(w)
    TMP = 0
    for (i in 1:n_W){
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e
      f = dnorm(x,muxw,sqrt(varxw))
      Hmx=1/(1+exp(-(betac[1]+betac[2]*x)))
      TMP = TMP + Hmx* f
    }
    return(TMP)}

  TmpIntegrate2 <- function(x){
    n_W = length(w)
    TMP = 0
    for (i in 1:n_W){
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e
      f = dnorm(x,muxw,sqrt(varxw))
      TMP = TMP + f
    }
    return(TMP)}



  #####################################################
  ## Standard error of the estimate of the parameters
  #####################################################
  S.E.Sigma_thetac_All_Int<-function(beta.par, theta.par, Lambda.par,su2e){#####   Function to get optimized S.E.Sigma_thetac
    mux_hat=Lambda.par[1]
    s2x_hat=Lambda.par[2]

    lambda_hat=s2x_hat/(s2x_hat+su2e)
    npsi<-length(beta.par)+length(theta.par)+length(Lambda.par)+length(su2e)
    J<-length(theta.par)
    ###########################  Find An
    #### we need to find:  E[derivative Phi with respect Lambda],  E[derivative Phi with respect beta]
    ##     E[derivative Q with respect Lambda], E[derivative Q with respect beta],  E[derivative Q with respect theta]
    #####

    An<-matrix(0,nrow = npsi,ncol = npsi)


    #### 1st  E[derivative Phi with respect to beta]

    dphi_i<-function(beta.par){
      lambda_hat=s2x_hat/(s2x_hat+su2e)

      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e
      fxw<-function(x){ dnorm(x,muxw,sqrt(varxw))  } # density of X/W


      Integrant1<-function(x){  # integrant for E(y/w)=P(beta,W)
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ1=(Hmx)*fxw(x)
        return(Integ1)
      }
      Integrant2<-function(x){ # integrant for the derivative of E(y/w) with respect beta_1
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ2=x*Hmx*(1-Hmx)*fxw(x)
        return(Integ2)
      }
      Integrant3<-function(x){ # integrant for the derivative of E(y/w) with respect beta_0
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ2=Hmx*(1-Hmx)*fxw(x)
        return(Integ2)
      }
      Integrant4<-function(x){ # integrant for the derivative of E(y/w) with respect beta_0
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Hmx1=Hmx*(1-Hmx)
        Integ4=Hmx1*(1-2*Hmx)*fxw(x)
        return(Integ4)
      }
      Integrant5<-function(x){ # integrant for the derivative of E(y/w) with respect beta_0
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Hmx1=Hmx*(1-Hmx)
        Integ5=x*Hmx1*(1-2*Hmx)*fxw(x)
        return(Integ5)
      }
      Integrant6<-function(x){ # integrant for the derivative of E(y/w) with respect beta_0
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Hmx1=Hmx*(1-Hmx)
        Integ6=(x^2)*Hmx1*(1-2*Hmx)*fxw(x)
        return(Integ6)
      }
      fHm<-function(x){
        1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
      }

      Int1=integrate(Integrant1, lower = -Inf, upper = Inf)$value
      Int2=integrate(Integrant2, lower = -Inf, upper = Inf)$value
      Int3=integrate(Integrant3, lower = -Inf, upper = Inf)$value
      Int4=integrate(Integrant4, lower = -Inf, upper = Inf)$value
      Int5=integrate(Integrant5, lower = -Inf, upper = Inf)$value
      Int6=integrate(Integrant6, lower = -Inf, upper = Inf)$value

      dPbeta<-c( Int3,  Int2)
      d2Pbeta<-matrix(0,ncol=2, nrow = 2)
      d2Pbeta[1,1]<-Int4
      d2Pbeta[1,2]<- d2Pbeta[2,1]<-Int5
      d2Pbeta[2,2]<-Int6

      Denominator=(Int1*(1-Int1))^2
      #Sum=(d2Pbeta*(y[i]-Int1)+dPbeta%*%t(dPbeta))*(Int1*(1-Int1)-(dPbeta-2*dPbeta*Int1))*dPbeta*(y[i]-Int1)
      Sum=(Int1*(1-Int1))*(d2Pbeta*(y[i]-Int1)-dPbeta%*%t(dPbeta))-dPbeta%*%t(dPbeta)*(y[i]-Int1)*(1-2*Int1)

      return(Sum/Denominator)
    }

    # Output  E[derivative Phi with respect beta]
    dphi_beta<-0
    for (i in 1:n) {
      dphi_beta<- dphi_beta+ dphi_i(beta.par)}
    An[4:5,4:5]<-dphi_beta/(n)

    ##################
    ## 2nd   E[derivative Q with respect beta]

    fQ_i_beta<-function(beta.par ){
      #mux_hat=Lambda.par[1]
      #sx2_hat=Lambda.par[2]
      lambda_hat=s2x_hat/(s2x_hat+su2e)

      fQi<-vector()
      for (j in 1:J){
        muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*su2e
        fxw<-function(x){ dnorm(x,muxw,sqrt(varxw))  } # density of X/W

        fHf<-function(x){
          Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
          return(Hmx*fxw(x))
        }
        l1=integrate(fHf,lower=K[j], upper=K[j+1])$value
        l2=integrate(fxw,lower=K[j], upper=K[j+1])$value

        fQi[j]<-l1-(1/(1+exp(-theta.par[j])))*l2
      }
      return(fQi)
    }

    sumdphi<-0
    for(i in 1:n){
      sumdphi<-sumdphi+numDeriv::jacobian(fQ_i_beta,beta.par)}
    An[6:(J+5),4:5]<-sumdphi/(n)

    #############################
    ##### E[derivative Q with respect theta]
    EfQ_theta<-function(theta.par){

      fx_hat<-function(x){ dnorm(x,mux_hat,sqrt(s2x_hat))  }
      EQtheta<-matrix(0,ncol=J, nrow=J)
      for (j in 1:J){
        IntAn3=integrate(fx_hat, lower = K[j], upper = K[j+1])$value
        EQtheta[j,j]<- -fH1(theta.par[j])*IntAn3}
      return(EQtheta)
    }

    An[6:(J+5),6:(J+5)]<-EfQ_theta(theta.par)

    ################  Aproximate E[derivative of Phi with respect Lambda]
    ##########################################

    phi_ilambda<-function(Lambda.par1){
      lambda_hat=Lambda.par1[2]/(Lambda.par1[2]+Lambda.par1[3])
      fxw<-function(x){
        muxw=Lambda.par1[1]*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*Lambda.par1[3]
        dnorm(x,muxw,sqrt(varxw))  } # density of X/W

      Integrant1<-function(x){  # integrant for E(y/w)
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ1=(Hmx)*fxw(x)
        return(Integ1)
      }
      Integrant2<-function(x){ # integrant for the derivative of E(y/w) with respect beta_1
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ2=x*Hmx*(1-Hmx)*fxw(x)
        return(Integ2)
      }
      Integrant3<-function(x){ # integrant for the derivative of E(y/w) with respect beta_0
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ2=Hmx*(1-Hmx)*fxw(x)
        return(Integ2)
      }
      Int1=integrate(Integrant1, lower = -Inf, upper = Inf)$value
      Int2=integrate(Integrant2, lower = -Inf, upper = Inf)$value
      Int3=integrate(Integrant3, lower = -Inf, upper = Inf)$value

      Sumbeta0=Int3*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      Sumbeta1=Int2*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      return(c( Sumbeta0, Sumbeta1))
    }
    # Aproximate
    sumdphi_lambda<-0
    for(i in 1:n){
      sumdphi_lambda<-sumdphi_lambda+numDeriv::jacobian(phi_ilambda,c(Lambda.par,su2e)) }#
    AEdphilambda<-sumdphi_lambda/(n)

    An[4:5,1:3]<-AEdphilambda

    ###   E[derivative of Q with respect Lambda]

    fQ_i_Lambda<-function(Lambda.par1 ){
      #mux_hat=Lambda.par1[1]
      #sx2_hat=Lambda.par1[2]
      lambda_hat=Lambda.par1[2]/(Lambda.par1[2]+Lambda.par1[3])

      fQi<-vector()
      for (j in 1:J){
        muxw=Lambda.par1[1]*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*Lambda.par1[3]
        fxw<-function(x){ dnorm(x,muxw,sqrt(varxw))  } # density of X/W

        fHf<-function(x){
          Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
          return(Hmx*fxw(x))
        }
        l1=integrate(fHf,lower=K[j], upper=K[j+1])$value
        l2=integrate(fxw,lower=K[j], upper=K[j+1])$value

        fQi[j]<-l1-(1/(1+exp(-theta.par[j])))*l2
      }
      return(fQi)
    }

    sumdphi<-0
    for(i in 1:n){
      sumdphi<-sumdphi+numDeriv::jacobian(fQ_i_Lambda,c(Lambda.par,su2e))}
    EdQ_Lambda<- sumdphi/(n)
    An[6:(J+5),1:3]<-EdQ_Lambda


    #####################   E[derivative of V with respect Lambda]

    An[1:3,1:3]<--diag(3)
    An[2,3]<--1

    ############################   %%%  %%%%%%%%  Complete An %%%%%%%%  ############################

    ############################ ############################ ############################ ############################
    ###########################  Find Bn
    #### we need to find:  Sigma, cov(Psi), E[derivative Phi with respect Lambda],  E[derivative Q with respect Lambda]


    ############################################################################
    ################### cov(Psi)
    ### Aproximate of cov(psi)
    ###################### cov(psi) %%%%%%%%%  Term  : E(Psi Psi^T)
    phi_i_parameter<-function(beta.par,Lambda.par){
      mux_hat=Lambda.par[1]
      s2x_hat=Lambda.par[2]

      lambda_hat=s2x_hat/(s2x_hat+su2e)

      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e
      fxw<-function(x){ dnorm(x,muxw,sqrt(varxw))  } # density of X/W

      Integrant1<-function(x){  # integrant for E(y/w)
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ1=(Hmx)*fxw(x)
        return(Integ1)
      }
      Integrant2<-function(x){ # integrant for the derivative of E(y/w) with respect beta_1
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ2=x*Hmx*(1-Hmx)*fxw(x)
        return(Integ2)
      }
      Integrant3<-function(x){ # integrant for the derivative of E(y/w) with respect beta_0
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        Integ2=Hmx*(1-Hmx)*fxw(x)
        return(Integ2)
      }
      Int1=integrate(Integrant1, lower = -10, upper = 10)$value
      Int2=integrate(Integrant2, lower = -Inf, upper = Inf)$value
      Int3=integrate(Integrant3, lower = -Inf, upper = Inf)$value

      Sumbeta0=Int3*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      Sumbeta1=Int2*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      return(c( Sumbeta0, Sumbeta1))
    }


    fQ_i_parameter<-function(beta.par,theta.par,Lambda.par){
      mux_hat=Lambda.par[1]
      s2x_hat=Lambda.par[2]

      lambda_hat=s2x_hat/(s2x_hat+su2e)
      fQi<-vector()
      for (j in 1:J){
        muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*su2e
        fxw<-function(x){ dnorm(x,muxw,sqrt(varxw))  } # density of X/W

        fHf<-function(x){
          Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
          return(Hmx*fxw(x))
        }
        l1=integrate(fHf,lower=K[j], upper=K[j+1])$value
        l2=integrate(fxw,lower=K[j], upper=K[j+1])$value

        fQi[j]<-l1-(1/(1+exp(-theta.par[j])))*l2
      }
      return(fQi)
    }


    V_int_i<-function(Lambda.par){
      c(w[i]-Lambda.par[1], (w[i]-Lambda.par[1])^2-Lambda.par[2]-su2e,s2w[i]/k-su2e)    # s2w[i] is defined on the sample
    }

    psi_i_parameter<- function(beta.par,theta.par,Lambda.par){
      c(t( V_int_i(Lambda.par)),t(phi_i_parameter(beta.par,Lambda.par)), t(fQ_i_parameter(beta.par,theta.par,Lambda.par)))
    }

    sumdpsi<-0
    for(i in 1:n){

      sumdpsi<-sumdpsi+ psi_i_parameter(beta.par,theta.par,Lambda.par)%*%t( psi_i_parameter(beta.par,theta.par,Lambda.par)) }#
    Epsi_psi<-sumdpsi/(n)  # Aproximate E(Psi Psi^T) ==> aprox cov(psi)


    ######################################## cov(psi) %%%%%%%%%    Finally we calculate cov(Psi)
    covPsi<- Epsi_psi

    ###############  Then Bn is:
    Bn<-matrix()
    Bn<-covPsi

    ############################     %%%%%%%%  Complete Bn %%%%%%%%  ############################
    Sigma_whole = solve(An)%*%Bn%*%t(solve(An))

    mat01<-cbind(matrix(0,nrow = J, ncol = 5), diag(J))

    Sigma_theta<-mat01%*%Sigma_whole%*%t(mat01)
    VarSigma_thetac<-Sigma_theta/n
    s.e_thetac<-sqrt(diag(VarSigma_thetac))

    s.e_thetac_J1<- sqrt((VarSigma_thetac[J,J]+VarSigma_thetac[1,1]-2*VarSigma_thetac[1,J]))

    s.e_thetacfull<-c(s.e_thetac,  s.e_thetac_J1)  # SE of the estimators of theta_1,...,theta_J and (theta_J-theta_1)


    #####################
    mat10<-cbind( diag(5), matrix(0,nrow = J, ncol = J))

    Sigma_par<-mat10%*%Sigma_whole%*%t(mat10)
    VarSigma_par<-Sigma_par/n
    s.e_par<-sqrt(diag(VarSigma_par))
    covab = VarSigma_par[4,5]
    ##################

    return(list(s.e_par,s.e_thetacfull,covab) )
  }
  ##########################################
  ### store estimates related to \theta
  thetac_w_e=vector()
  s.e_thetac_w_e=vector()  # Store s.e of \theta_c by using MLE of beta
  thetac.width.cover_w_e=vector()

  ##################
  ## beta.par
  ##################
  ###  Log likelihood to find estimate of  \beta
  logitLik<-function(betac_0){
    sumlik=0
    for (i in 1:n){
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e
      fxw<-function(x){ dnorm(x,muxw,sqrt(varxw)) } # density of X/W
      Integrant<-function(x){ # integrant for E(y/w)
        Hmx=1/(1+exp(-(betac_0[1]+betac_0[2]*x)))
        Integ=(Hmx)*fxw(x)
        return(Integ)
      }
      Int=integrate(Integrant, lower = -Inf, upper = Inf)$value
      liki=y[i]*log(Int)+(1-y[i])*log(1- Int)
      sumlik=sumlik+liki
    }
    return(-sumlik)
  }## end logitLik
  betac_0<-glm(y ~ w,family=binomial(link="logit"))$coefficients # initial point.
  output.logLk<-optim(betac_0, logitLik, gr=NULL, method="L-BFGS-B",control=list(maxit=500),  hessian=TRUE)
  betac=output.logLk$par
  beta.par=betac

  ##################
  ## theta.par
  ##################
  ###  Log likelihood to find estimate of \theta
  for (j in 1:J){

    l1=integrate(TmpIntegrate1, lower=K[j], upper=K[j+1])$value
    l2=integrate(TmpIntegrate2, lower=K[j], upper=K[j+1])$value

    L=l1/l2
    thetac_w_e[j]=log(L/(1-L))
  }
  theta.par=thetac_w_e
  ###############################################################
  ##                    Assymptotic variance and CI
  ###############################################################
  s.e_thetac_w_e<-S.E.Sigma_thetac_All_Int(beta.par, theta.par, Lambda.par,su2e)
  UL.thetac<-c(theta.par,theta.par[5]-theta.par[1])+1.96*s.e_thetac_w_e[[2]]
  LL.thetac<-c(theta.par,theta.par[5]-theta.par[1])-1.96*s.e_thetac_w_e[[2]]
  I.C<-cbind( LL.thetac,UL.thetac)
  ## Estimate
  thetac_output_w_e<-list(thetac_w_e,  thetac_w_e[J]-thetac_w_e[1] )

  ##############################################################
  ##here we transfer all estimates back
  ##############################################################
  mux_hat = (mux_hat - b0)/a0
  s2x_hat = s2x_hat/(a0^2)
  su2_hat = su2_hat/(a0^2)
  beta.par[2] = a0*betac[2]
  beta.par[1] = b0*betac[2] + betac[1]
  #########################################
  #transfer the variance
  #########################################
  se_Lambda = (s.e_thetac_w_e[[1]])[1:3]
  se_Lambda[1] = se_Lambda[1]/(a0)
  se_Lambda[2] = se_Lambda[2]/(a0^2)
  se_Lambda[3] = se_Lambda[3]/(a0^2)
  se_beta = (s.e_thetac_w_e[[1]])[4:5]
  se_beta[1] = sqrt((se_beta[1])^2+(b0^2)*(se_beta[2])^2+2*(b0)*s.e_thetac_w_e[[3]])
  se_beta[2] = (a0)*se_beta[2]
  ##########################  Output
  #
  #          SUMMARY
  cname <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  if (print.summary == TRUE) {
    all_par<-c(mux_hat, s2x_hat,su2_hat,beta.par, theta.par)
    s.e_all_par<- c(se_Lambda,se_beta, s.e_thetac_w_e[[2]][1:J])
    tval <-all_par/s.e_all_par
    pval <- 2 * (1 - pnorm(abs(tval)))
    outmat <- cbind(all_par,  s.e_all_par, tval, pval)
    name_Lam_beta<-c("mu.x", "sigma^2.x","sigma^2.u", "alpha", "beta")
    name_theta<-vector()
    for(i in 1:J){name_theta[i]<-paste("theta",i)}
    rownames(outmat) <- c(name_Lam_beta, name_theta)

    colnames(outmat) <- cname
    cat("Summary", "\n")
    cat(" ", "\n")
    print(round(outmat, 5))
    cat(" ", "\n")
    dmat<-theta.par[J]-theta.par[1]
    smat<-s.e_thetac_w_e[[2]][J+1]
    tmat<-dmat/smat
    outmat1 <- cbind(dmat, smat, tmat, 2 * (1 - pnorm(abs(tmat))))
    rownames(outmat1) <- paste("theta",J,"- theta 1:")
    colnames(outmat1) <-  cname
    print(round(outmat1, 5))
    cat(" ", "\n")

  }

  lambda_beta<-c(mux_hat, s2x_hat,su2_hat,beta.par)
  names(lambda_beta)<-name_Lam_beta
  out1 <- list(outmat1, theta.par,  lambda_beta, s.e_thetac_w_e[[2]], c(se_Lambda,se_beta))
  names(out1) <- c( "theta5-theta1", "theta","nuisance", "se.theta","se.nuisance")
  return(out1)
}
