Allint_linear<-function(y, W_int, C=NULL, print.summary=TRUE, standardize = TRUE){

  # W_int contains replicates, C=vector of cut points


  #####################################
  ## check dimensions and replicates
  #####################################

  n <- nrow(W_int)
  k <- ncol(W_int)

  ####################################
  ## standardize the raw data W
  ####################################
  # calculate row mean of w
  rowmean_w <- apply(W_int, 1, mean)
  mux_hat <- mean(rowmean_w)

  if(standardize == TRUE){
    a0 <- 1 / sqrt(mean((rowmean_w - mux_hat)^2))
    b0 <- - mux_hat * a0}else{
      a0 <- 1
      b0 <- 0}
  ww <- a0 * W_int + b0

  #####################################
  ## calculate nuisance parameters
  #####################################

  # calculate row mean of w
  mean_w <- apply(ww, 1, mean)
  mux_hat <- wmean <- mean(ww)
  w <- mean_w
  # calculate sigma_u

  s2w <- apply(W_int, 1, var)
  su2_hat <- mean(s2w)
  su2e <- su2_hat / k

  # calculate sigma_x
  s2x_hat <- max(mean((mean_w-mux_hat)^2)-su2e,0.2*(mean((mean_w-mux_hat)^2)))
  lambda_hat=max(s2x_hat/(s2x_hat+su2e),0.2)

  # form nuisance parameter vector Lambda
  Lambda.par=c(mux_hat,s2x_hat)

  ###########################################
  #Cut points
  ###########################################

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
    fMx(x)%*%theta
  }


  K=c(-Inf,C,Inf)
  TmpIntegrate1 <- function(x){
    n_W = length(w)
    TMP = 0
    for (i in 1:n_W){
      muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
      varxw=lambda_hat*su2e
      f = dnorm(x,muxw,sqrt(varxw))
      Hmx=betac[1]+betac[2]*x
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

    An[4:5,4:5]<-matrix(c(-1,-mux_hat,-mux_hat,-s2x_hat-mux_hat^2),nrow = 2,ncol = 2)

    ##################
    ## 2nd   E[derivative Q with respect beta]

    fQ_i_beta<-function(beta.par ){
      #mux_hat=Lambda.par[1]
      #sx2_hat=Lambda.par[2]
      lambda_hat=s2x_hat/(s2x_hat+su2e)

      fQi<-matrix(0, ncol = 2, nrow = J)

      for (j in 1:J){
        muxw=mux_hat*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*su2e
        fxw<-function(x){ dnorm(x,muxw,sqrt(varxw))  } # density of X/W

        fxf<-function(x){
          #Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
          return(x*fxw(x))
        }

        l2=integrate(fxf,lower=K[j], upper=K[j+1])$value
        l1=integrate(fxw,lower=K[j], upper=K[j+1])$value

        fQi[j, ]<-c(l1, l2)
      }
      return(fQi)
    }

    sumdphi<-matrix(0, ncol = 2, nrow = J)
    for(i in 1:n){
      sumdphi<-sumdphi+fQ_i_beta(beta.par)}
    An[6:(J+5),4:5]<-sumdphi/(n)

    ##################
    ##### E[derivative Q with respect to theta]
    An[6:(J+5),6:(J+5)]<-diag(-sumdphi[,1]/(n))

    ################ E[derivative of Phi with respect Lambda]
    ##########################################

    An[5,3]<-beta.par[2]
    ###   E[derivative of Q with respect Lambda]

    ##########################################
    #####################  derivative of Q with respect to Lambda
    fQ_i_Lambda<-function(Lambda.par1 ){
      #mux_hat=Lambda.par1[1]
      #sx2_hat=Lambda.par1[2]
      lambda_hat=Lambda.par1[2]/(Lambda.par1[2]+Lambda.par1[3])

      fQi<-vector()
      for (j in 1:J){
        muxw=Lambda.par1[1]*(1-lambda_hat)+lambda_hat*w[i]
        varxw=lambda_hat*Lambda.par1[3]
        fxw<-function(x){ dnorm(x,muxw,sqrt(varxw))  } # density of X/W

        fxf<-function(x){
          #Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
          mx=beta.par[1]+beta.par[2]*x
          return(mx*fxw(x))
        }
        l1=integrate(fxf,lower=K[j], upper=K[j+1])$value
        l2=integrate(fxw,lower=K[j], upper=K[j+1])$value

        fQi[j]<-l1-theta.par[j]*l2
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
    phi_i_parameter<-function(beta.par,su2e){
      phi_i<-c(y[i]-beta.par[1]-beta.par[2]*w[i], w[i]*y[i]-beta.par[1]*w[i]-beta.par[2]*w[i]^2+beta.par[2]*su2e)
      return( phi_i)
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

        fxf<-function(x){
          #Hmx=1/(1+exp(-(beta.par[1]+beta.par[2]*x)))
          Hmx=beta.par[1]+beta.par[2]*x
          return(Hmx*fxw(x))
        }
        l1=integrate(fxf,lower=K[j], upper=K[j+1])$value
        l2=integrate(fxw,lower=K[j], upper=K[j+1])$value

        fQi[j]<-l1-theta.par[j]*l2
      }
      return(fQi)
    }


    V_int_i<-function(Lambda.par){
      c(w[i]-Lambda.par[1], (w[i]-Lambda.par[1])^2-Lambda.par[2]-su2e,s2w[i]/k-su2e)    # s2w[i] is defined on the sample
    }

    psi_i_parameter<- function(beta.par,theta.par,Lambda.par){
      c(t( V_int_i(Lambda.par)),t(phi_i_parameter(beta.par,su2e)), t(fQ_i_parameter(beta.par,theta.par,Lambda.par)))
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


    Sigma_theta<-Sigma_whole[6:10,6:10]
    VarSigma_thetac<-Sigma_theta/n
    s.e_thetac<-sqrt(diag(VarSigma_thetac))

    s.e_thetac_J1<- sqrt((VarSigma_thetac[J,J]+VarSigma_thetac[1,1]-2*VarSigma_thetac[1,J]))

    s.e_thetacfull<-c(s.e_thetac,  s.e_thetac_J1)  # SE of the estimators of theta_1,...,theta_J and (theta_J-theta_1)


    Sigma_par<-Sigma_whole[1:5,1:5]
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
  ###  using estimating equation to find  \beta
  ymean<-mean(y)
  swy<-(y%*%w)/n-ymean*mux_hat
  sww<-w%*%w/n - mux_hat^2
  betac= c(0,0)
  betac[2] = swy/(sww-su2e)
  betac[1] = ymean - betac[2]*mux_hat
  beta.par = betac

  ##################
  ## theta.par
  ##################
  ###  Log likelihood to find estimate of \theta
  for (j in 1:J){

    l1=integrate(TmpIntegrate1, lower=K[j], upper=K[j+1])$value
    l2=integrate(TmpIntegrate2, lower=K[j], upper=K[j+1])$value

    thetac_w_e[j]=l1/l2
    # thetac_w_e[j]=log(L/(1-L))
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
    name_Lam_beta<-c("mu.x", "sigma^2.x","sigma^2.u", "alpha","beta")
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
  return(out1)  }
