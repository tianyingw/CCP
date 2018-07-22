Extint_logistic<-function(y, W_int, W_ext ,C=NULL, print.summary=TRUE, standardize = TRUE){  #  ww is matrix containing replicates, C=vector of cut points
  ####################################
  ## Lambda.par=(mux, s2x) and su2
  #####################################
  n<-dim(W_int)[1]
  r<- 1 # if W_int has replicates, this indicating the number of replicates

  ####################################
  # first we standardize the raw data W
  ####################################
  mux_hat<-mean(W_int)
  if(standardize == TRUE){
    a0 = 1/sd(W_int)
    b0 = -mux_hat*a0}else{
      a0 = 1
      b0 = 0}
  w = a0*W_int + b0
  W_ext = a0*W_ext + b0


  qq = apply(W_ext,1,var)
  su2_ext=qq
  su2e= mean(qq)

  mux_hat=mean(w)
  s2x_hat=max((var(w)-su2e/r),0.2*var(w))
  lambda_hat<-max(s2x_hat/(s2x_hat+su2e/r),0.2)

  #############################################
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


  ###############################################
  #Estimate \Theta

  TmpIntegrate1 <- function(x){
    n_W = length(w)
    TMP = 0
    for (i in 1:n_W){
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e/r
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
      varxw=lambda_hat*su2e/r
      f = dnorm(x,muxw,sqrt(varxw))
      TMP = TMP + f
    }
    return(TMP)}


  ################################################
  # Estimate function for Lambda(mux,sx2)
  V_int<-function(Lambda.par){
    c(mean(w-Lambda.par[1]), mean((w-Lambda.par[1])^2-Lambda.par[2]-su2e/r))
  }

  #####################################################
  ## Standard error of the estimate of the parameters
  #####################################################
  S.E.Sigma_thetac_Int_Ext<-function(beta.par, theta.par, Lambda.par,su2e){
    mux_hat=Lambda.par[1]
    s2x_hat=Lambda.par[2]
    lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)
    npsi<-length(beta.par)+length(theta.par)+length(Lambda.par)
    J<-length(theta.par)
    ###########################  Find An
    #### we need to find:  E[derivative Phi with respect Lambda],  E[derivative Phi with respect beta]
    ##     E[derivative Q with respect Lambda], E[derivative Q with respect beta],  E[derivative Q with respect theta]
    ##### E[derivative V with respect Lambda]

    An<-matrix(0,nrow = npsi,ncol = npsi)


    #### 1st  E[derivative Phi with respect to beta]

    dphi_i<-function(beta.par){
      lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)

      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e/r
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
        Integ3=Hmx*(1-Hmx)*fxw(x)
        return(Integ3)
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
    for (i in 1:n) {dphi_beta<- dphi_beta+ dphi_i(beta.par)}
    An[3:4,3:4]<-dphi_beta/n

    ##################
    ## 2nd   E[derivative Q with respect beta]

    fQ_i_beta<-function(beta.par ){
      #mux_hat=Lambda.par[1]
      #sx2_hat=Lambda.par[2]
      #cl=lambda/3
      lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)
      fQi<-vector()
      for (j in 1:J){
        muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*su2e/r
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
    for(i in 1:n){sumdphi<-sumdphi+numDeriv::jacobian(fQ_i_beta,beta.par)}
    An[5:(J+4),3:4]<-sumdphi/n

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

    An[5:(J+4),5:(J+4)]<-EfQ_theta(theta.par)

    ################  Aproximate E[derivative of Phi with respect Lambda]
    ##########################################

    phi_i_lambda<-function(Lambda.par){
      lambda_hat=max(Lambda.par[2]/(Lambda.par[2]+su2e/r),0.2)
      fxw<-function(x){
        muxw=Lambda.par[1]*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*su2e/r
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
    for(i in 1:n){sumdphi_lambda<-sumdphi_lambda+numDeriv::jacobian(phi_i_lambda,Lambda.par)} #
    AEdphilambda<-sumdphi_lambda/n

    An[3:4,1:2]<-AEdphilambda

    ###   E[derivative of Q with respect Lambda]

    ##########################################
    #####################  derivative of Q with respect to Lambda

    dfQi_lambda<-function(Lambda.par){
      mux_hat=Lambda.par[1]
      sx2_hat=Lambda.par[2]
      lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)

      dfQi_lambda<-matrix(0,ncol = 2,nrow = J)
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i] # mean of X/W
      varxw=lambda_hat*su2e/r   # variance of X/W
      fxw<-function(x){
        dnorm(x,muxw,sqrt(varxw))  } # density of X/W

      Hmx<-function(x){
        Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
        return(Hmx)
      }
      Int1a<-function(x){Hmx(x)*(x-muxw)*(1-lambda_hat)/varxw*fxw(x)}
      Int1b<-function(x){(x-muxw)*(1-lambda_hat)/varxw*fxw(x)}
      Int2a<-function(x){Hmx(x)*((1-lambda_hat)^2/(2*varxw)*fxw(x)*(-1+2*(x-muxw)*(w[i]-mux_hat)/(su2e/r)+(x-muxw)^2/varxw))}
      Int2b<-function(x){(1-lambda_hat)^2/(2*varxw)*fxw(x)*(-1+2*(x-muxw)*(w[i]-mux_hat)/(su2e/r)+(x-muxw)^2/varxw)}
      K0=c(-100,C,100)
      for (j in 1:J){
        l1a=integrate( Int1a,lower=K0[j], upper=K0[j+1])$value
        l1b=integrate( Int1b,lower=K0[j], upper=K0[j+1])$value
        l2a=integrate( Int2a,lower=K0[j], upper=K0[j+1])$value
        l2b=integrate( Int2b,lower=K0[j], upper=K0[j+1])$value
        dfQi_lambda[j,]<-c(l1a-(1/(1+exp(-theta.par[j])))*l1b, l2a-(1/(1+exp(-theta.par[j])))*l2b)
      }
      return(dfQi_lambda)
    }


    #### output
    sumdQi_lambda<-0
    for(i in 1:n){sumdQi_lambda<-sumdQi_lambda+dfQi_lambda(Lambda.par)} #
    EdQ_Lambda<-sumdQi_lambda/n

    An[5:(J+4),1:2]<-EdQ_Lambda


    #####################   E[derivative of V with respect Lambda]

    An[1:2,1:2]<--diag(2)

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
      lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)    #max(s2x_hat/(s2x_hat+su2e),cl)
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e/r
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
      Int1=integrate(Integrant1, lower = -Inf, upper = Inf)$value
      Int2=integrate(Integrant2, lower = -Inf, upper = Inf)$value
      Int3=integrate(Integrant3, lower = -Inf, upper = Inf)$value

      Sumbeta0=Int3*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      Sumbeta1=Int2*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      return(c( Sumbeta0, Sumbeta1))
    }

    fQ_i_parameter<-function(beta.par,theta.par,Lambda.par){
      mux_hat=Lambda.par[1]
      s2x_hat=Lambda.par[2]
      #cl=lambda/3
      lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)    #max(s2x_hat/(s2x_hat+su2e),cl)
      fQi<-vector()
      for (j in 1:J){
        muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*su2e/r
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
      c(w[i]-Lambda.par[1], (w[i]-Lambda.par[1])^2-Lambda.par[2]-su2e/r)
    }

    psi_i_parameter<- function(beta.par,theta.par,Lambda.par){
      c(t( V_int_i(Lambda.par)),t(phi_i_parameter(beta.par,Lambda.par)), t(fQ_i_parameter(beta.par,theta.par,Lambda.par)))
    }

    sumdpsi<-0
    for(i in 1:n){sumdpsi<-sumdpsi+ psi_i_parameter(beta.par,theta.par,Lambda.par)%*%t( psi_i_parameter(beta.par,theta.par,Lambda.par)) }#
    Epsi_psi<-sumdpsi/n  # Aproximate E(Psi Psi^T) ==> aprox cov(psi)


    ######################################## cov(psi) %%%%%%%%%    Finally we calculate cov(Psi)
    covPsi<- Epsi_psi

    ######################################## Sigma %%%%%%%%%

    Sigma<-mean((su2_ext-su2e)^2)  # su2_ext has to be defined out of the function S.E.Sigma_thetac_Int_Ext --- doesn't divide by r!!
    ######################################## End Sigma

    ######################################## E(omega psi respect to su2) %%%%%%%%%
    phi_i_parameter_su2e<-function(su2e){
      mux_hat=Lambda.par[1]
      s2x_hat=Lambda.par[2]
      #cl=lambda/3
      lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)     #max(s2x_hat/(s2x_hat+su2e),cl)

      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e/r
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
      Int1=integrate(Integrant1, lower = -Inf, upper = Inf)$value
      Int2=integrate(Integrant2, lower = -Inf, upper = Inf)$value
      Int3=integrate(Integrant3, lower = -Inf, upper = Inf)$value

      Sumbeta0=Int3*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      Sumbeta1=Int2*(Int1*(1-Int1))^(-1)*(y[i]-Int1)
      return(c( Sumbeta0, Sumbeta1))
    }

    fQ_i_parameter_su2e<-function(su2e){
      mux_hat=Lambda.par[1]
      s2x_hat=Lambda.par[2]
      #cl=lambda/3
      lambda_hat=max(s2x_hat/(s2x_hat+su2e/r),0.2)     #max(s2x_hat/(s2x_hat+su2e),cl)
      fQi<-vector()
      for (j in 1:J){
        muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*su2e/r
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

    V_int_su2e_i<-function(su2e){
      c(w[i]-Lambda.par[1], (w[i]-Lambda.par[1])^2-Lambda.par[2]-su2e/r)
    }

    psi_i_parameter_su2e<- function(su2e){
      c(t( V_int_su2e_i(su2e)),t(phi_i_parameter_su2e(su2e)), t(fQ_i_parameter_su2e(su2e)))
    }

    sumdpsi<-0
    for(i in 1:n){sumdpsi<-sumdpsi+numDeriv::jacobian(psi_i_parameter_su2e,su2e)}#
    EOmegaPsi_su2e<-sumdpsi/n  # Aproximate E(Psi Psi^T) ==> aprox cov(psi)


    ###############  Then Bn is:
    Bn<-matrix()
    m<-length(W_ext[,1])  # size external dataset
    Bn<-covPsi+(n/m)*EOmegaPsi_su2e%*%Sigma%*%t(EOmegaPsi_su2e)

    ############################     %%%%%%%%  Complete Bn %%%%%%%%  ############################
    Sigma_whole = solve(An)%*%Bn%*%t(solve(An))

    Sigma_theta<-Sigma_whole[5:9,5:9]
    VarSigma_thetac<-Sigma_theta/n
    s.e_thetac<-sqrt(diag(VarSigma_thetac))

    s.e_thetac_J1<- sqrt((VarSigma_thetac[J,J]+VarSigma_thetac[1,1]-2*VarSigma_thetac[1,J]))

    s.e_thetacfull<-c(s.e_thetac,  s.e_thetac_J1)  # SE of the estimators of theta_1,...,theta_J and (theta_J-theta_1)
    #####################

    Sigma_par<-Sigma_whole[1:4,1:4]
    VarSigma_par<-Sigma_par/n
    s.e_par<-sqrt(diag(VarSigma_par))
    covab = VarSigma_par[3,4]

    ##################
    return(list(s.e_par,s.e_thetacfull,covab,sqrt(Sigma/n) ))
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
  Lambda.par=c(mux_hat,s2x_hat)
  logitLik<-function(betac_0){
    sumlik=0
    for (i in 1:n){
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e/r
      fxw<-function(x){ dnorm(x,muxw,sqrt(varxw)) } # density of X/W
      Integrant<-function(x){ # integrant for E(y/w)
        Hmx=1/(1+exp(-(betac_0[1]+betac_0[2]*x)))
        Integ=(Hmx)*fxw(x)
        return(Integ)
      }
      Int=integrate(Integrant, lower = -10, upper = 10)$value
      liki=y[i]*log(Int)+(1-y[i])*log(1- Int)
      sumlik=sumlik+liki
    }
    return(-sumlik)
  }

  betac_0<-glm(y ~ w,family=binomial(link="logit"))$coefficients # initial beta from the NaiveLik.
  output.logLk<-optim(betac_0, logitLik, gr=NULL, method="L-BFGS-B",control=list(maxit=2000),  hessian=TRUE)
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
  ###############################################################
  ##                    Assymptotic variance and CI
  ###############################################################
  theta.par=thetac_w_e
  s.e_thetac_w_e<-S.E.Sigma_thetac_Int_Ext(beta.par, theta.par, Lambda.par,su2e)
  UL.thetac<-c(theta.par,theta.par[5]-theta.par[1])+1.96*s.e_thetac_w_e[[2]]
  LL.thetac<-c(theta.par,theta.par[5]-theta.par[1])-1.96*s.e_thetac_w_e[[2]]
  I.C<-cbind( LL.thetac,UL.thetac)
  ## Estimate
  thetac_output_w_e<-cbind(thetac_w_e, thetac_w_e[J]-thetac_w_e[1] )
  ########################## ########################## ##########################

  ##############################################################
  ##here we transfer all estimates back
  ##############################################################
  mux_hat = (mux_hat - b0)/a0
  s2x_hat = s2x_hat/(a0^2)
  su2_hat = su2e/(a0^2)
  beta.par[2] = a0*betac[2]
  beta.par[1] = b0*betac[2] + betac[1]
  #########################################
  #transfer the variance
  #########################################
  se_Lambda = c((s.e_thetac_w_e[[1]])[1:2], s.e_thetac_w_e[[4]])
  se_Lambda[1] = se_Lambda[1]/(a0)
  se_Lambda[2] = se_Lambda[2]/(a0^2)
  se_Lambda[3] = se_Lambda[3]/(a0^2)
  se_beta = (s.e_thetac_w_e[[1]])[3:4]
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
  out1 <- list(outmat1, theta.par,  lambda_beta, s.e_thetac_w_e[[2]], c(se_Lambda, se_beta))

  names(out1) <- c( "theta5-theta1", "theta","nuisance", "se.theta","se.nuisance")
  return(out1)
}
