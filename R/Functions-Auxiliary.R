

###############################################################################
# Auxiliary functions with introductions
###############################################################################

# --- Chi-squared statistic for linear hypothesis testing --- #
PH.Test <- function(
    object,
    testinfo = NULL
){

  ## preparations
  A   <- as.matrix(testinfo$A)
  b   <- as.matrix(testinfo$b)
  N   <- object$N
  alp <- object$coefficient[,1]
  alp.VarCov <- object$alp.VarCov

  ## do linear hypothesis testing
  value  <- as.numeric(N*t(A%*%alp-b)%*%solve(A%*%alp.VarCov%*%t(A))%*%(A%*%alp-b))
  pvalue <- 1-pchisq(value,df=nrow(A))

  ## output
  out <- list(
    value  = value,
    pvalue = pvalue
  )
  return(out)

}


# --- Cox proportional hazards model (PH) without any auxiliary information --- #
PH <- function(yobs,delta,X,Var=TRUE,Var.Robust=FALSE,eps=1e-6,maxit=5e4){

  ### Preparations
  N <- length(yobs)
  yobs[yobs<=0] <- 1e-6
  pbet <- ncol(X)

  ### calculate initial values for bet
  bet.init <- rep(0,pbet)

  ### calculate MLE for beta
  numit <- 1
  bet.old <- bet.init
  repeat{
    InfoMScore <- PH.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
    Score <- InfoMScore$Score
    InfoM <- InfoMScore$InfoM
    dev <- MASS::ginv(InfoM)%*%Score
    bet <- bet.old + as.vector(dev)
    if( max(abs(bet-bet.old))>eps & numit<maxit ){
      bet.old <- bet
      numit <- numit + 1
    }else{
      break
    }
  }
  convergence <- (numit<maxit)

  ### calculate SEs or not (and tidy them)
  if(Var==TRUE){

    ## calculate SEs for bet using explicit formula !!!
    InfoMScore <- PH.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
    InfoM <- InfoMScore$InfoM
    Score <- InfoMScore$Score
    InfoM.inv <- MASS::ginv(InfoM)
    if(Var.Robust==TRUE){
      Score.Influs <- PH.Influence.EE(yobs=yobs,delta=delta,X=X,bet=bet)
      Score.Influs.Cov <- t(Score.Influs)%*%Score.Influs/N
      VarCov <- InfoM.inv %*% Score.Influs.Cov %*% InfoM.inv
    }else{
      VarCov <- InfoM.inv
    }
    bet.se <- sqrt(diag(VarCov)/N)

    ### tidy the results: inference
    zvalue.bet <- bet/bet.se
    pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
    res <- data.frame(Est=bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                      row.names=colnames(X))

  }else{

    # tidy the results directly
    VarCov <- NULL
    res <- data.frame(Est=bet,row.names=colnames(X))

  }

  ### output
  out <- list(
    info = list(
      convergence = convergence,
      bet.init = bet.init
    ),
    res=res,
    Score=Score,
    InfoM=InfoM,
    VarCov=VarCov
  )
  return(out)

}

# --- Kaplan-Meier (KM) estimator
KM.fit <- function(tm,yobs,delta,Var=FALSE,type="right"){

  ### preparation ###
  time.jump <- sort(yobs[delta==1])
  d         <- sapply(time.jump,function(x){sum(yobs==x)})
  R         <- sapply(time.jump,function(x){sum(yobs>=x)})
  prods     <- 1-d/R

  # calculate the values of KM at pre-specified time points tm
  if(type=="right"){
    Stm <- sapply(tm,function(tmi){prod(prods[time.jump<=tmi])}) # right-continuous
  }else{
    Stm <- sapply(tm,function(tmi){prod(prods[time.jump<tmi])}) # left-continuous
  }

  # prepare variances
  if(Var == FALSE){
    out <- Stm
  }else{
    sums     <- d/(R*(R-d))
    VAR.part <- sapply(tm,function(tmi){sum(sums[time.jump<=tmi])})
    out      <- list(
      sprob = Stm,
      se    = Stm*sqrt(VAR.part)
    )
  }

  # output
  return(out)

}




###############################################################################
# Other auxiliary functions
#   - They are not needed to be displayed explicitly
###############################################################################

PH.SP.Process <- function(yobs,delta,X,auxinfo,bet.init,Score,InfoM,control){

  # preparations
  N    <- length(yobs)
  pbet <- ncol(X)

  # convert auxiliary information into tractable form in coding

  # - prepare the auxframe matrix
  tstar.unique <- sort(unique(unlist(lapply(auxinfo$aux,function(sourcei){lapply(sourcei,function(auxj){auxj$tstar})}))))
  auxframe <- as.numeric()
  for(isource in 1:length(auxinfo$aux)){
    for(itime in 1:length(auxinfo$aux[[isource]])){ # isource <- itime <- 1
      aux.c    <- auxinfo$aux[[isource]][[itime]]
      K.c      <- length(aux.c$sprob)
      auxframe <- rbind(
        auxframe,
        cbind(rep(isource,K.c),
              rep(auxinfo$M[isource],K.c),
              rep(aux.c$tstar,K.c),
              rep(which(tstar.unique==aux.c$tstar),K.c),
              aux.c$sprob,
              aux.c$subidx*1))
    }
  }
  colnames(auxframe) <- c('source','M','tstar','tstaridx','sprob',paste('ind',1:N,sep="")) # rename this matrix
  auxnum <- nrow(auxframe)

  # - solve for the optimal direction needed in de-correlate the original control variate
  Psi.D1 <- PH.SP.EE.D1(bet=bet.init,yobs=yobs,delta=delta,X=X,auxframe=auxframe)
  if(control$penalty$inv == TRUE){

    h <- solve(InfoM)%*%Psi.D1

  }else{

    if(control$penalty$do == FALSE){

      h <- solve(InfoM)%*%Psi.D1

    }else{

      tunparahs <- exp(seq(log(0.005),log(5),length.out=50))*sqrt(log(pbet)/N)
      h         <- matrix(0,nrow=pbet,ncol=auxnum)
      for(k in 1:auxnum){

        # - obtain solutions at all tuning parameters
        h.k.all <- array(0,dim=c(pbet,length(tunparahs)))
        for(itunparah in 1:length(tunparahs)){
          tunparah <- tunparahs[itunparah]
          numit <- 1
          h.old <- rep(0,pbet)
          repeat{
            grad   <- as.vector(InfoM%*%h.old-Psi.D1[,k])
            h.temp <- h.old - 0.1*grad
            h.k    <- Soft.Threshold(theta=h.temp,lam=0.1*tunparah)
            if( max(abs(h.k-h.old))>control$eps & numit<control$maxit ){
              h.old <- h.k
              numit <- numit + 1
            }else{
              break
            }
          }
          h.k.all[,itunparah] <- h.k
        }

        # - get the optimal solution based on the information criterion
        ICs <- sapply(1:length(tunparahs),function(itunparah){
          h.k.itunparah <- h.k.all[,itunparah]
          df.k.itunparah <- sum(abs(h.k.itunparah)>control$eps)
          val <- t(h.k.itunparah)%*%InfoM%*%h.k.itunparah-2*Psi.D1[,k]%*%h.k.itunparah +
            2*df.k.itunparah/N
        })
        h[,k] <- as.matrix(h.k.all[,which.min(ICs)])

      }

    }

  }

  # - prepare for the control variate
  Psi.Individual <- PH.SP.EE(bet=bet.init,yobs=yobs,delta=delta,X=X,auxframe=auxframe,individual=TRUE)
  Psi            <- apply(Psi.Individual,2,mean)
  if(control$penalty$do == FALSE){
    cv <- Psi
  }else{
    cv <- as.vector(Psi+t(h)%*%Score)
  }

  # - influence function of the control variate
  Score.Influs   <- PH.Influence.EE(yobs=yobs,delta=delta,X=X,bet=bet.init)
  Psi.Du1        <- PH.SP.EE.Du1(bet=bet.init,yobs=yobs,delta=delta,X=X,auxframe=auxframe)
  Lamstar.Influs <- PH.Influence.Lam(tm=auxframe[,'tstar'],yobs=yobs,delta=delta,X=X,bet=bet.init,betGiven=TRUE)
  cv.Influs.1    <- Psi.Individual
  cv.Influs.2    <- t(t(Lamstar.Influs)*Psi.Du1)
  cv.Influs.3    <- Score.Influs%*%h
  cv.Influs      <- cv.Influs.1 + cv.Influs.2 + cv.Influs.3
  pai            <- apply(auxframe[,-c(1:5),drop=FALSE],1,mean)

  # variance-covariance matrix and its adjustment concerning uncertainty
  V.aux <- array(0,dim=c(auxnum,auxnum))
  if(any(auxframe[,'M']<Inf)){
    sprob.Influs <- do.call(cbind,lapply(1:auxnum,function(k){
      KM.Influence(tm=auxframe[k,'tstar'],yobs=yobs,delta=delta,subidx=(auxframe[k,-c(1:5)]==1))
    }))
    sprob.VarCov <- t(sprob.Influs)%*%sprob.Influs/N # cov(sprob.Influs) #
    for(ksource in sort(unique(auxframe[,'source']))){
      idxk             <- which(auxframe[,'source']==ksource)
      pai.idxk         <- diag(pai[idxk],nrow=length(idxk))
      M.idxk           <- auxframe[idxk,'M'][1]
      V.aux[idxk,idxk] <- pai.idxk%*%(sprob.VarCov[idxk,idxk])%*%pai.idxk*(N/M.idxk)
    }
  }
  V <- t(cv.Influs)%*%cv.Influs/N + V.aux

  # adjust the control-variate and obtain the estimates
  V.inv      <- solve(V)
  if(auxinfo$hetero==TRUE){
    V.inv.SVD  <- svd(V.inv)
    V.inv.root <- V.inv.SVD$u%*%diag(sqrt(V.inv.SVD$d),nrow=length(cv))%*%t(V.inv.SVD$v)
    tau <- AuxSP.Pen(cv=cv,V.inv.root=V.inv.root,pai=pai,N=N)
    cva <- cv-pai*tau
  }else{
    tau <- rep(0,length(cv))
    cva <- cv
  }

  # output
  out <- list(
    cva       = cva,
    cv.Influs = cv.Influs,
    V.inv     = V.inv,
    V         = V,
    tau       = tau
  )
  return(out)

}

PH.CE.Process <- function(yobs,delta,X,auxinfo,control){

  # preparations
  N    <- length(yobs)

  # convert auxiliary information into tractable form in coding

  # - basic elements
  K            <- length(auxinfo$aux) # number of external studies
  source.index <- do.call(c,lapply(1:K,function(ksource){rep(ksource,length(auxinfo$aux[[ksource]]$ce))})) # a vector to indicate the study number of ce

  # - expand the aux: ce.inn and Influs.inn
  for(ksource in 1:K){
    aux.idx.k    <- auxinfo$aux[[ksource]]$idx
    ce.inn.k     <- PH(yobs,delta,X[,aux.idx.k,drop=F],Var=FALSE,maxit=control$maxit,eps=control$eps)$res[,1]
    Influs.inn.k <- PH.Influence(yobs,delta,X[,aux.idx.k,drop=F],bet=ce.inn.k)
    auxinfo$aux[[ksource]]$ce.inn     <- ce.inn.k
    auxinfo$aux[[ksource]]$Influs.inn <- Influs.inn.k
  }

  # - prepare for the control variate
  ce.ext <- do.call(c,lapply(1:K,function(ksource){auxinfo$aux[[ksource]]$ce}))
  ce.inn <- do.call(c,lapply(1:K,function(ksource){auxinfo$aux[[ksource]]$ce.inn}))
  cv     <- ce.inn - ce.ext

  # - influence function of the control variate
  cv.Influs <- do.call(cbind,lapply(1:K,function(ksource){auxinfo$aux[[ksource]]$Influs.inn}))

  # variance-covariance matrix and its adjustment concerning uncertainty
  V.inn <- t(cv.Influs)%*%cv.Influs/N
  for(ksource in 1:K){
    if(is.null(auxinfo$aux[[ksource]]$V)==TRUE){
      idxk <- (source.index==ksource)
      auxinfo$aux[[ksource]]$V <- V.inn[idxk,idxk]
    }
  }
  V <- V.inn + as.matrix(do.call(Matrix::bdiag,lapply(1:K,function(ksource){auxinfo$aux[[ksource]]$V*N/auxinfo$aux[[ksource]]$M})))

  # adjust the control-variate and obtain the estimates
  V.inv <- solve(V)
  if(auxinfo$hetero==TRUE){
    V.inv.SVD  <- svd(V.inv)
    V.inv.root <- V.inv.SVD$u%*%diag(sqrt(V.inv.SVD$d),nrow=length(cv))%*%t(V.inv.SVD$v)
    tau <- AuxCE.Pen(cv=cv,V.inv.root=V.inv.root,N=N)
    cva <- cv-tau
  }else{
    tau <- rep(0,length(cv))
    cva <- cv
  }

  # output
  out <- list(
    cva       = cva,
    cv.Influs = cv.Influs,
    V.inv     = V.inv,
    V         = V,
    tau       = tau
  )
  return(out)

}

AuxSP.Pen <- function(cv,V.inv.root,pai,N){

  y.tilde <- as.vector( V.inv.root%*%cv )
  X.tilde <- V.inv.root%*%diag(pai,nrow=length(cv))

  # solve the penalized problem: solve adaptive lasso using lars
  w             <-  (1/abs(cv))*(pai)
  X.tilde.star  <- t(t(X.tilde)/w)
  sol.lars      <- lars::lars(X.tilde.star,y.tilde,trace=FALSE,normalize=FALSE,intercept=FALSE)
  tau.path      <- t(as.matrix(sol.lars$beta))/w # each
  tau.path.RSS  <- apply(X.tilde %*% tau.path - y.tilde,2,function(x){sum(x^2)})
  tau.path.Card <- apply(tau.path,2,function(x){sum(x!=0)})
  IC.all        <- as.numeric( tau.path.RSS + log(N)/N * tau.path.Card ) # BIC type criterion
  min_IC.idx    <- which.min( IC.all  )
  tau           <- tau.path[,min_IC.idx]

  # output: return final estimates
  return(tau)

}

AuxCE.Pen <- function(cv,V.inv.root,N){

  # preparations
  y.tilde <- as.vector( V.inv.root%*%cv )
  X.tilde <- V.inv.root

  # solve adaptive lasso using lars
  w <- 1/abs(cv)
  X.tilde.star <- t(t(X.tilde)/w)
  sol.lars <- lars::lars(X.tilde.star,y.tilde,trace=FALSE,normalize=FALSE,intercept=FALSE)
  tau.path <- t(as.matrix(sol.lars$beta))/w # each
  tau.path.RSS <- apply(X.tilde %*% tau.path - y.tilde,2,function(x){sum(x^2)})
  tau.path.Card <- apply(tau.path,2,function(x){sum(x!=0)})
  IC.all <- as.numeric( tau.path.RSS + log(N)/N * tau.path.Card ) # BIC type criterion
  min_IC.idx <- which.min( IC.all  )
  tau <- tau.path[,min_IC.idx]

  # output: return final estimates
  return(tau)

}

PH.SP.EE <- function(bet,yobs,delta,X,auxframe,individual=TRUE){

  # the conditional survival rates
  N <- length(yobs)
  tm      <- auxframe[,'tstar']
  expXbet <- as.vector(exp(X%*%bet))
  SS0 <- sapply(yobs,function(Yi){sum(expXbet*(yobs>=Yi))})
  Lamt <- sapply(tm,function(tmj){sum( (yobs<=tmj)[delta==1]/SS0[delta==1] )})
  St.baseline <- exp(-Lamt)
  Stx <- outer(expXbet,St.baseline,function(x,y){y^x})

  # the estimating equations in individual levels
  Psi <- t((t(Stx)-auxframe[,'sprob'])*auxframe[,-c(1:5),drop=FALSE])
  if(individual==FALSE){Psi <- apply(Psi,2,mean)}

  # output
  return(Psi)
}

PH.SP.EE.Du1 <- function(bet,yobs,delta,X,auxframe){

  # the conditional survival rates
  tm      <- auxframe[,'tstar']
  expXbet <- as.vector(exp(X%*%bet))
  SS0 <- sapply(yobs,function(Yi){sum(expXbet*(yobs>=Yi))})
  Lamt <- sapply(tm,function(tmj){sum( (yobs<=tmj)[delta==1]/SS0[delta==1] ) })
  St.baseline <- exp(-Lamt)
  Stx <- outer(expXbet,St.baseline,function(x,y){y^x})

  # the estimating equations in individual levels
  Psi.Du1 <- - apply(Stx*t(auxframe[,-c(1:5),drop=FALSE])*expXbet,2,mean)

  # output
  return(Psi.Du1)
}

PH.SP.EE.D1 <- function(bet,yobs,delta,X,auxframe){

  # the conditional survival rates
  tm       <- auxframe[,'tstar']
  expXbet  <- as.vector(exp(X%*%bet))
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet*(yobs>=Yi))})
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))
  Lamt <- sapply(tm,function(tmj){sum( (yobs<=tmj)[delta==1]/SS0[delta==1] )})
  Lamt.D1 <- do.call(cbind,lapply(tm,function(tmi){apply(SS1*delta*(yobs<=tmi)/SS0^2,2,sum)}))
  St.baseline <- exp(-Lamt)
  Stx <- outer(expXbet,St.baseline,function(x,y){y^x})
  IStx <- Stx*t(auxframe[,-c(1:5),drop=FALSE])
  Psi.D1 <- do.call(cbind,lapply(1:nrow(auxframe),function(k){
    apply((outer(expXbet,Lamt.D1[,k],FUN="*")-XexpXbet*Lamt[k])*IStx[,k],2,mean)
  }))

  # output
  return(Psi.D1)
}

PH.Beta.LogLik <- function(bet,yobs,delta,X){
  # for maximization

  ## prepare
  N <- length(yobs)
  Xbet <- as.vector(X %*% bet)
  log.SS0 <- sapply(yobs,function(Yi){log(mean(exp(Xbet)*(yobs>=Yi)))})

  ## calculate the partial likelihood function for beta (log form)
  val.log <- sum((Xbet-log.SS0)*delta)

  # output
  return(val.log)

}

PH.ScoreInfoM <-  function(bet,yobs,delta,X,IsScore=FALSE,IsInfoM=TRUE){
  # Score is the first  derivative of [positive] log-likelihood
  # InfoM is the second derivative of [negative] log-likelihood

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])})/N # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))/N # [S(1)(t,bet)] at pre-specified yobs and beta
  out <- list()

  ## prepare information matrix
  if(IsInfoM==TRUE){
    SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
      yobsGYi <- yobs>=Yi
      as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
    }))/N
    I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
    I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
    InfoM <- I1-I2
    out <- c(out,list(InfoM=InfoM))
  }

  ## prepare score vector
  if(IsScore==TRUE){
    U <- (X - SS1/SS0)*delta
    Score <- apply(U,2,mean)
    out <- c(out,list(Score=Score))
  }


  ## output
  return(out)

}

PH.Score.Individual <-  function(bet,yobs,delta,X,Iscov=TRUE){
  # Score is the first  derivative of [positive] log-likelihood

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])})/N # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))/N # [S(1)(t,bet)] at pre-specified yobs and beta
  U <- (X - SS1/SS0)*delta
  if(Iscov==TRUE){
    out <- t(U)%*%U/N
  }else{
    out <- U
  }

  ## output
  return(out)

}

PH.Lam <- function(tm,yobs,delta,X,bet=NULL,type="right"){

  # prepare
  if(is.null(bet)){
    bet <- PH(yobs,delta,X)$res[,1]
  }
  N <- length(yobs)
  expXbet <- as.vector(exp(X %*% bet))
  SS0 <- sapply(yobs,function(Yi){mean(expXbet*(yobs>=Yi))}) # [S(0)(t,bet)] at pre-specified yobs and beta

  # calculate the Lam(t) at specified time points tm
  if(type=="right"){
    Lam <- sapply(tm,function(tmj){
      sum( (yobs<=tmj)[delta==1]/SS0[delta==1] ) / N
    })
  }else if(type=="left"){
    Lam <- sapply(tm,function(tmj){
      sum( (yobs< tmj)[delta==1]/SS0[delta==1] ) / N
    })
  }

  # output
  return(Lam)

}

PH.Influence <- function(yobs,delta,X,bet=NULL){

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## calculate main part (except an inverse of information matrix)
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet

  ## prepare information matrix
  SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
    yobsGYi <- yobs>=Yi
    as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
  }))
  I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
  I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
  InfoM <- I1-I2

  ## final influence functions (individual level) and output
  Influs <- U%*%MASS::ginv(InfoM)
  return(Influs)

}

PH.Influence.EE <- function(yobs,delta,X,bet=NULL){
  # The incluence function for estimating equations

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## calculate main part
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet

  ## final influence functions (individual level) and output
  return(U)

}

PH.Influence.Lam <- function(tm,yobs,delta,X,bet,betGiven=FALSE){

  # tm <- c(0.2,0.4,0.6,0.8,1)

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta

  ## specific calculation
  if(betGiven==FALSE){

    XexpXbet <- X*expXbet
    SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

    ## prepare information matrix
    SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
      yobsGYi <- yobs>=Yi
      as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
    }))
    I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
    I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
    InfoM <- I1-I2

    ## calculate main part 1 (same as bet)
    UU1 <- (X - SS1/SS0)*delta
    UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
    UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
    U   <- UU1-(UU2-UU3)*expXbet
    Influs.bet <- U%*%MASS::ginv(InfoM)

    ## calculate main part 2 (with tm vary)
    LL1.tm <- do.call(cbind,lapply(tm,function(tmi){(yobs<=tmi)*delta/(SS0/N)}))
    LL2.tm <- do.call(cbind,lapply(tm,function(tmi){
      sapply(yobs,function(Yj){sum(delta*(yobs<=min(tmi,Yj))/(SS0^2/N))})
    }))*expXbet
    LL3.pre <- do.call(cbind,lapply(tm,function(tmi){apply(SS1*delta*(yobs<=tmi)/SS0^2,2,sum)}))
    LL3.tm <- Influs.bet %*% LL3.pre

    ## influence function for Lam (individual level) and output
    Influs.Lam <- LL1.tm - LL2.tm - LL3.tm


  }else{

    ## calculate main part 2 (with tm vary)
    LL1.tm <- do.call(cbind,lapply(tm,function(tmi){(yobs<=tmi)*delta/(SS0/N)}))
    LL2.tm <- do.call(cbind,lapply(tm,function(tmi){
      sapply(yobs,function(Yj){sum(delta*(yobs<=min(tmi,Yj))/(SS0^2/N))})
    }))*expXbet

    ## influence function for Lam (individual level) and output
    Influs.Lam <- LL1.tm - LL2.tm

  }



  ## final influence functions
  return(Influs.Lam)

}

PH.Lam.D1 <- function(tm,yobs,delta,X,bet){

  expXbet  <- exp(as.vector(X %*% bet))
  XexpXbet <- X*expXbet
  SS0      <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1      <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta
  Lam.D1   <- -do.call(cbind,lapply(tm,function(tmi){apply(SS1*delta*(yobs<=tmi)/SS0^2,2,sum)}))
  return(Lam.D1)

}


PH.InfoM.Half <-  function(bet,yobs,delta,X){

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  idx.delta <- which(delta==1)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])})/N # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))/N # [S(1)(t,bet)] at pre-specified yobs and beta

  ## prepare information matrix (half)
  InfoM.Half <- array(0,dim=c(N^2,pbet))
  for(j in 1:sum(delta)){
    InfoM.Half[(idx.delta[j]-1)*N+(1:N),] <-
      t(t(X)-SS1[idx.delta[j],]/SS0[idx.delta[j]]) *
      (yobs>=yobs[idx.delta[j]]) * sqrt(expXbet/SS0[idx.delta[j]]) / N
  }

  ## output
  return(InfoM.Half)

}

KM.Cov <- function(tm1,tm2,yobs,delta){

  ### preparation ###
  time.jump <- sort(yobs[delta==1])
  d         <- sapply(time.jump,function(x){sum(yobs==x)})
  R         <- sapply(time.jump,function(x){sum(yobs>=x)})
  prods     <- 1-d/R

  # calculate the values of KM at pre-specified time points tm
  Stm1 <- sapply(tm1,function(tmi){prod(prods[time.jump<=tmi])})
  Stm2 <- sapply(tm2,function(tmi){prod(prods[time.jump<=tmi])})

  # prepare variances
  sums <- d/(R*(R-d))
  Mat.Cov <- array(0,dim=c(length(tm1),length(tm2)))
  for(i in 1:length(tm1)){
    for(j in 1:length(tm2)){
      Mat.Cov[i,j] <- Stm1[i]*Stm2[j]*sum(sums[time.jump<=min(tm1[i],tm2[j])])
    }
  }

  # output
  return(Mat.Cov)

}

KM.Influence <- function(tm,yobs,delta,subidx=NULL){

  N <- length(yobs)

  if(is.null(subidx)){

    Hyobs <- sapply(yobs,function(yobsi){sum(yobs<=yobsi)})/N
    St <- KM.fit(tm=tm,yobs=yobs,delta=delta)
    KM.Influs <- do.call(cbind,lapply(1:length(tm),function(itm){
      tmi <- tm[itm]
      P1 <- delta*(yobs<=tmi)/(1-Hyobs); P1[is.na(P1)] <- 0
      P2 <- sapply(yobs,function(yobsi){sum(delta*(yobs<=min(tmi,yobsi))/((1-Hyobs)^2),na.rm=TRUE)/N})
      (P2-P1)*St[itm]
    }))

  }else{

    subidx <- (1:N) %in% ((1:N)[subidx])
    prop <- mean(subidx)

    Hyobs <- sapply(yobs,function(yobsi){sum((yobs<=yobsi)*subidx)})/N/prop
    St <- KM.fit(tm=tm,yobs=yobs[subidx],delta=delta[subidx])
    KM.Influs <- do.call(cbind,lapply(1:length(tm),function(itm){
      tmi <- tm[itm]
      P1 <- delta*(yobs<=tmi)/(1-Hyobs); P1[is.na(P1)] <- 0
      P2 <- sapply(yobs,function(yobsi){sum(subidx*delta*(yobs<=min(tmi,yobsi))/((1-Hyobs)^2),na.rm=TRUE)/N})/prop
      subidx*(P2-P1)*St[itm]/prop
    }))

  }

  # KM.Influs[is.na(KM.Influs)] <- 0

  # output
  return(KM.Influs)

}

KM.SP.Process <- function(yobs,delta,X=NULL,auxinfo){

  # preparations
  N        <- length(yobs)
  subgroup <- any(unlist(lapply(auxinfo$aux,function(sourcei){lapply(sourcei,function(auxj){!is.null(auxj$gfunc)})})))

  # convert auxiliary information into tractable form in coding
  tstar.unique <- sort(unique(unlist(lapply(auxinfo$aux,function(sourcei){lapply(sourcei,function(auxj){auxj$tstar})}))))
  auxframe <- as.numeric()
  for(isource in 1:length(auxinfo$aux)){
    for(itime in 1:length(auxinfo$aux[[isource]])){ # isource <- itime <- 1
      aux.c    <- auxinfo$aux[[isource]][[itime]]
      K.c      <- length(aux.c$sprob)
      if(is.null(aux.c$gfunc)==TRUE){
        ind.idx.c <- array(1,dim=c(K.c,N))
      }else{
        ind.idx.c <- aux.c$gfunc(X)*1
      }
      auxframe <- rbind(
        auxframe,
        cbind(rep(isource,K.c),
              rep(auxinfo$M[isource],K.c),
              rep(aux.c$tstar,K.c),
              rep(which(tstar.unique==aux.c$tstar),K.c),
              aux.c$sprob,
              ind.idx.c))
    }
  }
  colnames(auxframe) <- c('source','M','tstar','tstaridx','sprob',paste('ind',1:N,sep="")) # rename this matrix
  auxnum <- nrow(auxframe)

  # the control variate and its influence function (homogeneous)
  sprob.inn <- KM.fit(tm=auxframe[,"tstar"],yobs=yobs,delta=delta)
  sprob.ext <- auxframe[,"sprob"]
  cv        <- sprob.inn - sprob.ext
  cv.Influs <- do.call(cbind,lapply(1:auxnum,function(iaux){
    KM.Influence(tm=auxframe[iaux,"tstar"],yobs=yobs,delta=delta,
                 subidx=(auxframe[iaux,-c(1:5)]==1))}))

  # variance-covariance matrix and its adjustment concerning uncertainty
  V.pre <- t(cv.Influs)%*%cv.Influs/N
  V.aux <- array(0,dim=c(auxnum,auxnum))
  if(any(auxframe[,'M']<Inf)){
    for(ksource in sort(unique(auxframe[,'source']))){
      idxk             <- which(auxframe[,'source']==ksource)
      M.idxk           <- auxframe[idxk,'M'][1]
      V.aux[idxk,idxk] <- V.pre[idxk,idxk]*(N/M.idxk)
    }
  }
  V <- V.pre + V.aux

  # adjust the control-variate and obtain the estimates
  V.inv <- solve(V)
  if(auxinfo$hetero==TRUE){
    V.inv.SVD  <- svd(V.inv)
    V.inv.root <- V.inv.SVD$u%*%diag(sqrt(V.inv.SVD$d),nrow=length(cv))%*%t(V.inv.SVD$v)
    tau        <- AuxSP.Pen(cv=cv,V.inv.root=V.inv.root,pai=1,N=N)
    cva        <- cv - tau
  }else{
    tau <- rep(0,length(cv))
    cva <- cv
  }

  # output
  out <- list(
    cva       = cva,
    cv.Influs = cv.Influs,
    V.inv     = V.inv,
    V         = V,
    tau       = tau
  )
  return(out)

}

St.Sub.KM <- function(tstar,yobs,delta,G){

  K <- nrow(G)
  tSP <- rep(0, K)
  for(k in 1:K){
    idx <- (G[k,]==1)
    fit.k <- summary(survival::survfit(survival::Surv(yobs, delta) ~ 1, subset=idx))
    tSP[k] <- min( c(1,fit.k$surv)[ c(0,fit.k$time) <= tstar ] )
  }
  return( tSP )

}

Soft.Threshold <- function(theta,lam){
  # the solution to LASSO penalty under the simplest situation: 2^{-1}(z-theta)^2+penalty
  # this is the equation (2.6) in Fan and Li (2001)
  res <- sign(theta)*pmax(abs(theta)-lam,0)
  return(res)
}

Basis.Bernstein <- function(x,k,d,il,ir){

  # x: evaluated points (can be a vector)
  # k: the kth basis function
  # d: the degree
  # il: the left end point of the interval
  # ir: the right end point of the interval

  x.scale <- (x-il)/(ir-il)
  value <- choose(d,k) * (x.scale^k) * ((1-x.scale)^(d-k))

  return(value)
}

###############################################################################









###############################################################################
# AuxSP-PH using EL: Fit the Cox PH Model with Auxiliary SP using Empirical likelihood method
###############################################################################


#============== The main function that fits the model ==============#
PH.AuxSP.DEL.fit <- function(yobs,delta,X,aux,maxit=30,nboot=100){

  ### Specify the dimension
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data

  ### prepare the auxinfo matrix
  tstar.unique <- unique(unlist(lapply(aux,function(studyi){sapply(studyi,function(timej){timej$tstar})})))
  auxinfo <- as.numeric()
  for(istudy in 1:length(aux)){ # istudy <- 2; itime <- 1
    for(itime in 1:length(aux[[istudy]])){
      aux.c <- aux[[istudy]][[itime]]
      K.c <- length(aux[[istudy]][[itime]]$sprob)
      auxinfo <- rbind(
        auxinfo,
        cbind(rep(istudy,K.c),
              rep(aux.c$tstar,K.c),
              rep(which(tstar.unique==aux.c$tstar),K.c),
              aux.c$sprob,
              aux.c$gfunc(X)))
    }
  }
  colnames(auxinfo) <- c('study','tstar','tstaridx','sprob',paste('ind',1:N,sep="")) # rename this matrix

  ### Initial value for beta and alpha using LY ###
  phfit <- PH.fit(yobs,delta,X)
  bet.init <- phfit$res[,1]
  alp.init <- pmax(PH.Lam(tm=tstar.unique,bet.init,yobs,delta,X),0)
  xi.init <- rep(0, nrow(auxinfo))   # Initial value for xi

  # solve the estimating equations directly
  rot.init <- c(bet.init, alp.init, xi.init)
  sol <- PH.AuxSP.DEL.EE.Solve(yobs=yobs,delta=delta,X=X,auxinfo=auxinfo,
                               rot.init=rot.init,maxit=maxit)
  bet <- sol$bet; alp <- sol$alp
  convergence <- sol$convergence

  ### Variance Estimation###
  # SE <- sqrt(diag(PH.AuxSP.DEL.VCOV(bet,alp,yobs,delta,X,auxinfo))/N)
  bet.boots <- array(0,dim=c(N,p))
  for(iboot in 1:nboot){
    # print(iboot)
    idx <- sample(1:N,N,replace=TRUE)
    sol.iboot <- PH.AuxSP.DEL.EE.Solve(yobs=yobs[idx],delta=delta[idx],X=X[idx,,drop=F],
                                       auxinfo=auxinfo[,c(1:4,idx+4)],
                                       rot.init=rot.init,maxit=maxit)
    bet.boots[iboot,] <- sol.iboot$bet
  }
  SE <- apply(bet.boots,2,sd)

  ### summary the results ###
  zvalue <- bet/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=bet, SE=SE, zvalue=zvalue, pvalue=pvalue,
                     row.names=colnames(X))

  ### Output the Results ###
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X), # collect my data info
    coef=coef,
    convergence=convergence # converge or not
  )
  return(out)

}



#==== score equations for the profiled log-empirical likelihood function ====#
PH.AuxSP.DEL.ProLik.EE <- function(par,yobs,delta,X,auxinfo)
{

  ### Prepare ###
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data
  sumK <- nrow(auxinfo) # number of Groups
  ntime <- length(unique(auxinfo[,'tstar']))
  bet <- par[1:p]
  alp <- par[(p+1):(p+ntime)]
  xi  <- par[-(1:(p+ntime))]
  tstar.unique <- sapply(1:ntime,function(itstaridx){
    auxinfo[auxinfo[,'tstaridx']==itstaridx,'tstar'][1]
  })
  expXbet <- as.vector(exp( X %*% bet ))

  ### calculate [S(k)(t,bet)] at pre-specified yobs and beta ###
  SS0 <- sapply(yobs,function(Yi){mean(expXbet*(yobs>=Yi))})
  SS1 <- t(sapply(yobs,function(Yi){apply(X*expXbet*(yobs>=Yi),2,mean)}))

  ### Prepare values relating to auxiliary information ###
  sur.allG <- t(sapply(1:nrow(auxinfo),function(iG){
    exp(-alp[auxinfo[iG,'tstaridx']]*expXbet)}))
  Psi <- t((sur.allG-auxinfo[,'sprob'])*auxinfo[,-c(1:4)])
  Psi.xi <- as.vector( Psi %*% xi )
  Psidb.xi <- -X*expXbet*as.vector(t(sur.allG*auxinfo[,-c(1:4)]*alp[auxinfo[,'tstaridx']])%*%xi)
  Psida.xi <- - sapply(1:length(alp),function(itstaridx){
    t(sur.allG*auxinfo[,-c(1:4)]*(auxinfo[,'tstaridx']==itstaridx))%*%xi})*expXbet

  ### Estimating Equation U ###
  denom <- 1 + Psi.xi
  v <- apply(Psida.xi/denom,2,mean)
  IYt <- sapply(tstar.unique,function(itstar){yobs <= itstar})
  SS0.vIYt <- SS0 + as.vector(IYt%*%v)
  U1_beta <- apply( (X-SS1/SS0.vIYt)*delta-Psidb.xi/denom, 2, sum )
  U2_xi    <- apply( Psi/denom, 2, sum )
  U3_alpha <- apply( IYt*delta/SS0.vIYt, 2, sum )  - N*alp
  U <- c(U1_beta,U3_alpha,U2_xi)

  ### output
  return(U)

}

#==== numerical second derivative for the profiled log-empirical likelihood function ====#
#' @export
PH.AuxSP.DEL.ProLik.dEE <- function(dpar,yobs,delta,X,auxinfo){
  dh <- 0.00000001
  dl <- length(dpar)
  te<-matrix(0, dl, dl)
  for(i in 1: dl)
  {
    s1 <- s2 <- dpar
    s1[i] <- s1[i] + dh
    s2[i] <- s2[i] - dh
    te[,i]<- ( PH.AuxSP.DEL.ProLik.EE(s1,yobs,delta,X,auxinfo)-
                 PH.AuxSP.DEL.ProLik.EE(s2,yobs,delta,X,auxinfo) ) / (2*dh)
  }
  return(te)
}

PH.AuxSP.DEL.EE.Solve <- function(yobs,delta,X,auxinfo,rot.init,maxit){

  ### prepare
  p <- ncol(X) # number of covariates in internal data

  ### solve the equations
  cytry <- try({
    numit <- 1 # number of iteration
    rot.old <- rot.all <- rot.init
    repeat{ # Newton-Rapson
      EE  <- PH.AuxSP.DEL.ProLik.EE(rot.old,yobs,delta,X,auxinfo)
      dEE <- PH.AuxSP.DEL.ProLik.dEE(rot.old,yobs,delta,X,auxinfo)
      dev <- as.vector(MASS::ginv(dEE)%*%EE) # solve(dEE,EE) #
      rot.new <- rot.old - dev
      # print(sum(dev));print(numit)
      if(sqrt(mean((rot.new-rot.old)^2)) > 1e-6 & numit < maxit){
        rot.old <- rot.new
        rot.all <- rbind(rot.all,rot.new)
        numit <- numit + 1
      }else{
        rot.all <- rbind(rot.all,rot.new)
        break
      }
    }
  },
  silent=T);  rot.new
  convergence <- ifelse(class(cytry) == "try-error" || numit > maxit,F,T)
  bet <- rot.new[1:p]; alp <- rot.new[(p+1):(p+length(unique(auxinfo[,'tstar'])))]

  ### output
  out <- list(
    bet=bet,
    alp=alp,
    convergence=convergence
  )
  return(out)
}


#==== variance-covariance matrix for the method ====#
PH.AuxSP.DEL.VCOV <- function(bet,alp,yobs,delta,X,auxinfo)
{

  ### Prepare ###
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data
  sumK <- nrow(auxinfo) # number of Groups
  ntime <- length(unique(auxinfo[,'tstar']))
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  SS0 <- sapply(yobs,function(Yi){mean(expXbet*(yobs>=Yi))}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- t(sapply(yobs,function(Yi){apply(X*expXbet*(yobs>=Yi),2,mean)})) # [S(1)(t,bet)] at pre-specified yobs and beta
  SS2.vec <- t(sapply(yobs,function(Yi){
    as.vector(t(X)%*%(X*expXbet*(yobs>=Yi)))
  }))/N

  ### Prepare values ###


  KK <- delta/SS0^2



  I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=p,ncol=p)
  I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
  Sigma <- I1-I2



  ### calculate VCOV ###
  VCOV <- solve( Sigma+B%*%solve(Q)%*%t(B) )

  ### output ###
  return(VCOV)

}


###############################################################################




