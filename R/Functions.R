

#=================================#
# The overall package description #
#=================================#

#' aiscv: Auxiliary Information Synthesis via Control Variates.
#'
#' This package provides an information synthesis framework that can transfer auxiliary aggregate
#'    information into the targets via control variates
#' project work
#' @importFrom stats cov model.frame model.matrix model.response optim pchisq pnorm rbinom rexp runif sd uniroot
#' @docType package
#' @name aiscv
NULL




###############################################################################
# Cox proportional hazards model (PH)
#   - with auxiliary subgroup survival probabilities (SP)
#   - Techniques:
#     + Euclidean control variate (ECV) to synthesis information
#     + adaptive penalization (AP) to guarantee transprotability
###############################################################################

#' @title Cox proportional hazards model with auxiliary information in the type
#'   of subgroup survival probabilities at multiple time points
#'
#' @description Fit the Cox proportional hazards (PH) model with auxiliary information
#'   in the type of subgroup survival probabilities (SP) at multiple time points.
#'
#' @param yobs a numeric vector that indicated the observed survival times.
#' @param delta a numeric vector that stores the right-censoring indicators
#'   (1 = event of interest happens, and 0 = censoring).
#' @param X a matrix that includes all interested covariates.
#' @param auxinfo a list used to store information about the inputted auxiliary
#'   information. Specifically, it contains the following elements:
#'   \describe{
#'     \item{\code{aux}}{ a list that includes information about historical
#'       aggregated statistics. It is a list of lists, and each sub-list represents
#'       auxiliary information from the same source study. In each source study,
#'       it contains several time points and each time point contains the following
#'       three elements:
#'       \describe{
#'         \item{\code{tstar}}{ the time point that the auxiliary information
#'           was calculated at.}
#'         \item{\code{sprob}}{ auxiliary subgroup survival rates for each subgroup
#'           at the current time point.}
#'         \item{\code{gfunc}}{ a function used to identify the subgroup.}
#'       } }
#'     \item{\code{M}}{ a vector of integers that indicate the sample sizes of
#'       all external studies. Note that the dimension of this vector should equal
#'       to the number of studies in \code{aux} and the order should be consistent.}
#'     \item{\code{hetero}}{ a logical value. If it is \code{TRUE}, the penalization
#'       will be applied to identify the potential heterogeneous auxiliary subgroup
#'       survival rates and make a refinement to the final estimator.}
#'   }
#' @param iindex indices of our interested coefficients for covariates in \code{X}.
#'   The default is \code{NULL} and all coefficients will be displayed.
#' @param control indicates more detailed control of the underlying model fitting
#'   procedures. It is a list of the following three arguments:
#'   \describe{
#'     \item{\code{maxit}}{ a positive integer that denotes the maximum iteration
#'       number in optimization. The default value is \code{10000}.}
#'     \item{\code{eps}}{ a positive small numeric value that denotes the tolerance
#'       for convergence. The default value is \code{1e-6}.}
#'     \item{\code{penalty}}{ a list that contains information about the decorrelation steps.
#'       The argument \code{do} indicates whether the penalization should be used,
#'       and the argument \code{inv} controls whether the penalization in the decorrelation
#'       setup should be used.}
#'   }
#'
#' @return
#'   A list of fitted results is returned.
#'   Within this outputted list, the following elements can be found:
#'   \describe{
#'     \item{\code{coefficient}}{ estimation and inference of regression coefficients.}
#'     \item{\code{tau}}{ estimated degrees of population discrepency.}
#'   }
#'
#' @export PH.SP

PH.SP <- function(
    yobs,
    delta,
    X,
    auxinfo   = list(aux = NULL, M = NULL, hetero = TRUE),
    iindex    = NULL,
    control   = list(
      maxit   = 1e4,
      eps     = 1e-6,
      penalty = list(do = FALSE, inv = TRUE)
    )
){

  ## Specify several basic elements
  N         <- length(yobs)   # sample size in internal data
  pbet      <- ncol(X)   # number of parameters
  bet.names <- colnames(X)
  if(is.null(iindex)==TRUE){iindex <- 1:pbet}
  gamidx    <- c(1:pbet)[-iindex]
  palp      <- length(iindex)
  pgam      <- pbet-palp

  ## original estimates (without auxiliary information) and other preparations

  # initial estiamtor with different penalties
  if(control$penalty$do == FALSE){

    bet.init <- PH(yobs=yobs,delta=delta,X=X,Var=FALSE,
                   eps=control$eps,maxit=control$maxit)$res[,1]

  }else{

    suppressWarnings({
      bet.init <- as.numeric(coef(glmnet::cv.glmnet(
        X,survival::Surv(yobs,delta),family="cox",
        alpha=1,thresh=control$eps,maxit=control$maxit
      ),s="lambda.min"))
    })

  }

  # basic preparation
  InfoMScore   <- PH.ScoreInfoM(bet=bet.init,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
  Score        <- InfoMScore$Score
  InfoM        <- InfoMScore$InfoM
  InfoM.alp    <- InfoM[iindex,iindex,drop=FALSE]
  InfoM.gamalp <- InfoM[-iindex,iindex,drop=FALSE]

  # solve for the optimal direction
  if(control$penalty$inv == TRUE){

    if(pgam > 0){
      InfoM.gam <- InfoM[-iindex,-iindex,drop=FALSE]
      w         <- solve(InfoM.gam)%*%InfoM.gamalp
    }else{
      w <- matrix(0,nrow=pgam,ncol=palp)
    }

  }else{

    if(control$penalty$do == FALSE & pgam > 0){

      InfoM.gam <- InfoM[-iindex,-iindex,drop=FALSE]
      w         <- solve(InfoM.gam)%*%InfoM.gamalp

    }else{

      w <- matrix(0,nrow=pgam,ncol=palp)
      if(pgam > 0){

        # prepare matrices
        InfoM.gam    <- InfoM[-iindex,-iindex,drop=FALSE]
        InfoM.gamalp <- InfoM[-iindex,iindex,drop=FALSE]

        # get w
        tunparaws <- exp(seq(log(0.005),log(5),length.out=50))*sqrt(log(pbet)/N)
        for(j in 1:palp){

          # - obtain solutions at all tuning parameters
          w.j.all <- array(0,dim=c(pgam,length(tunparaws)))
          for(itunparaw in 1:length(tunparaws)){
            tunparaw <- tunparaws[itunparaw]
            numit <- 1
            w.old <- rep(0,pgam)
            repeat{
              grad <- as.vector(InfoM.gam%*%w.old-InfoM.gamalp[,j])
              w.temp <- w.old - 0.1*grad
              w.j <- Soft.Threshold(theta=w.temp,lam=0.1*tunparaw)
              if( max(abs(w.j-w.old))>control$eps & numit<control$maxit ){
                w.old <- w.j
                numit <- numit + 1
              }else{
                break
              }
            }
            w.j.all[,itunparaw] <- w.j
          }

          # - get the optimal solution based on the information criterion
          ICs <- sapply(1:length(tunparaws),function(itunparaw){
            w.j.itunparaw  <- w.j.all[,itunparaw]
            df.j.itunparaw <- sum(abs(w.j.itunparaw)>control$eps)
            t(w.j.itunparaw)%*%InfoM.gam%*%w.j.itunparaw-2*InfoM.gamalp[,j]%*%w.j.itunparaw +
              2*df.j.itunparaw/N
          })
          w[,j] <- as.matrix(w.j.all[,which.min(ICs)])

        }

      }

    }

  }
  v <- rbind(diag(palp),-w)

  # get original estimator of alpha
  InfoM.inv <- solve(InfoM.alp-t(w)%*%InfoM.gamalp)
  if(control$penalty$do == FALSE){
    alp.ori <- bet.init[iindex]
  }else{
    Score.alp <- Score[iindex]
    Score.gam <- Score[-iindex]
    Score.DC  <- as.vector(Score.alp-t(w)%*%Score.gam)
    alp.ori   <- bet.init[iindex] + as.vector(InfoM.inv%*%Score.DC) # 改成直接求解
  }

  ## fit the model (if auxiliary information exists)
  if(is.null(auxinfo$aux)==TRUE){

    alp        <- alp.ori
    alp.VarCov <- InfoM.inv
    alp.se     <- sqrt(diag(alp.VarCov)/N)
    tau        <- NA
    cvinfo     <- NULL

  }else{

    # process the auxiliary information and extract necessary elements
    aux.process <-  PH.SP.Process(
      yobs=yobs,delta=delta,X=X,auxinfo=auxinfo,bet.init=bet.init,Score=Score,InfoM=InfoM,control=control)
    V         <- aux.process$V
    cva       <- aux.process$cva
    cv.Influs <- aux.process$cv.Influs
    V.inv     <- aux.process$V.inv
    tau       <- aux.process$tau

    # covariance matrix of initial estimator and control variate
    Score.Influs   <- PH.Influence.EE(yobs=yobs,delta=delta,X=X,bet=bet.init)
    alp.ori.Influs <- Score.Influs%*%v%*%t(InfoM.inv)
    Gam            <- t(cv.Influs)%*%alp.ori.Influs/N

    # define the improved estimator formally
    Gam.V.inv <- t(V.inv%*%Gam)
    alp <- as.vector(alp.ori-Gam.V.inv%*%cva)

    # estimator of variance (using oracle properties directly)
    homo.idx <- (tau==0)
    if(any(homo.idx)){
      V.inv.homo <- solve(V[homo.idx,homo.idx,drop=F])
      Gam.homo   <- Gam[homo.idx,,drop=F]
      alp.VarCov <- InfoM.inv - t(Gam.homo)%*%V.inv.homo%*%Gam.homo
      Gam.V.inv  <- t(V.inv.homo%*%Gam.homo)
    }else{
      alp.VarCov <- InfoM.inv
      Gam.V.inv  <- NULL
    }
    alp.se <- sqrt( diag(alp.VarCov)/N )
    alp.se

    # information related to control variates
    cvinfo <- list(
      cv.Influs = cv.Influs,
      Gam.V.inv = Gam.V.inv
    )

  }

  ## do inference and combine results into pre-specified style
  coefficient <- data.frame(Est       = alp,
                            SE        = alp.se,
                            zvalue    = alp/alp.se,
                            pvalue    = 2*(1-pnorm(abs(alp/alp.se))),
                            row.names = bet.names[iindex])

  ## output
  out <- list(
    coefficient = coefficient,
    tau         = tau,
    alp.VarCov  = alp.VarCov,
    control     = control,
    N           = N,
    cvinfo      = cvinfo
  )
  return(out)

}

###############################################################################

###############################################################################
# High-dimensional Cox proportional hazards model (HPH)
#   - with auxiliary subgroup survival probabilities (SP)
#   - Techniques:
#     + Euclidean control variate (ECV) to synthesis information
#     + adaptive penalization (AP) to guarantee transprotability
#     + decorrelation (DC) to stabilize penalized estimators
#     + data splitting (DS) to allow for high dimensional covariates
###############################################################################

#' @title High-Dimensional Cox proportional hazards model with auxiliary information
#'   in the type of subgroup survival probabilities at multiple time points
#'
#' @description Fit the High-dimensional Cox proportional hazards (HPH) model with
#'   auxiliary information in the type of subgroup survival probabilities (SP)
#'   at multiple time points.
#'
#' @param yobs a numeric vector that indicated the observed survival times.
#' @param delta a numeric vector that stores the right-censoring indicators
#'   (1 = event of interest happens, and 0 = censoring).
#' @param X a matrix that includes all interested covariates.
#' @param auxinfo a list used to store information about the inputted auxiliary
#'   information. Specifically, it contains the following elements:
#'   \describe{
#'     \item{\code{aux}}{ a list that includes information about historical
#'       aggregated statistics. It is a list of lists, and each sub-list represents
#'       auxiliary information from the same source study. In each source study,
#'       it contains several time points and each time point contains the following
#'       three elements:
#'       \describe{
#'         \item{\code{tstar}}{ the time point that the auxiliary information
#'           was calculated at.}
#'         \item{\code{sprob}}{ auxiliary subgroup survival rates for each subgroup
#'           at the current time point.}
#'         \item{\code{gfunc}}{ a function used to identify the subgroup.}
#'       } }
#'     \item{\code{M}}{ a vector of integers that indicate the sample sizes of
#'       all external studies. Note that the dimension of this vector should equal
#'       to the number of studies in \code{aux} and the order should be consistent.}
#'     \item{\code{hetero}}{ a logical value. If it is \code{TRUE}, the penalization
#'       will be applied to identify the potential heterogeneous auxiliary subgroup
#'       survival rates and make a refinement to the final estimator.}
#'   }
#' @param dsplit a list used to control the process of data splitting. Specifically,
#'   it contains the following elements:
#'   \describe{
#'     \item{\code{do}}{ a logical that indicate whether the data splitting procedure
#'       should be used. The default value if \code{TRUE}.}
#'     \item{\code{times}}{ an integer that denotes the times of data splitting. The
#'       default value if \code{10}.}
#'     \item{\code{combine}}{ a character that specifies the method of combining
#'       different estimates in each data splitting procedure. It can be either
#'       \code{"median"} (the default) and \code{"mean"}.}
#'   }
#' @param iindex indices of our interested coefficients for covariates in \code{X}.
#'   The default is \code{NULL} and all coefficients will be displayed.
#' @param control indicates more detailed control of the underlying model fitting
#'   procedures. It is a list of the following three arguments:
#'   \describe{
#'     \item{\code{maxit}}{ a positive integer that denotes the maximum iteration
#'       number in optimization. The default value is \code{10000}.}
#'     \item{\code{eps}}{ a positive small numeric value that denotes the tolerance
#'       for convergence. The default value is \code{1e-6}.}
#'   }
#'
#' @return
#'   A list of fitted results is returned.
#'   Within this outputted list, the following elements can be found:
#'   \describe{
#'     \item{\code{coefficient}}{ estimation and inference of regression coefficients.}
#'     \item{\code{tau}}{ estimated degrees of population discrepency.}
#'   }
#'
#' @export HPH.SP

HPH.SP <- function(
    yobs,
    delta,
    X,
    auxinfo = list(aux = NULL, M = NULL, hetero = TRUE),
    iindex  = NULL,
    dsplit  = list(do = TRUE, times = 10, combine = "median"),
    control = list(eps = 1e-6, maxit = 1e4)
){

  # Specify several basic elements
  N         <- length(yobs)   # sample size in internal data
  bet.names <- colnames(X)
  if(is.null(iindex)==TRUE){iindex <- 1:ncol(X)}
  pbet <- ncol(X)
  palp <- length(iindex)

  # the key part of the method
  if(dsplit$do == TRUE){

    # repeat the splitting procedure for multiple times
    alp.list <- alp.VarCov.list <- tau.list <- list()
    for(idsplit in 1:dsplit$times){

      # split the whole data set into two sets
      idx.temp <- sort(sample(1:N,round(N/2)))
      idx.list <- list(idx.temp,(1:N)[-idx.temp])

      # do for all data sets
      for(idx.sample in idx.list){

        # obtain active set from the remaining data

        # - prepare data
        X.remove     <- X[-idx.sample,,drop=FALSE]
        yobs.remove  <- yobs[-idx.sample]
        delta.remove <- delta[-idx.sample]

        # - fit an initial model for weights
        bet.initial <- as.numeric(coef(glmnet::cv.glmnet(
          X.remove,survival::Surv(yobs.remove,delta.remove),
          family="cox",alpha=1,thresh=control$eps,maxit=control$maxit
        ),s="lambda.min"))
        pweights         <- 1/pmax(abs(bet.initial),control$eps)
        pweights[iindex] <- 0

        # - obtain active set
        bet.active <- as.numeric(coef(glmnet::cv.glmnet(
          X.remove,survival::Surv(yobs.remove,delta.remove),
          family="cox",alpha=1,thresh=control$eps,maxit=control$maxit,
          penalty.factor=pweights
        ),s="lambda.min"))
        X.active.idx  <- sort(unique(c(iindex,which(abs(bet.active)>control$eps))))

        # prepare changed inputs
        iindex.refine  <- which(iindex %in% X.active.idx)
        auxinfo.refine <- auxinfo
        for(isource in 1:length(auxinfo$aux)){
          for(itime in 1:length(auxinfo$aux[[isource]])){
            subidx.c <- auxinfo$aux[[isource]][[itime]]$subidx
            auxinfo.refine$aux[[isource]][[itime]]$subidx <- subidx.c[,idx.sample,drop=FALSE]
          }
        }

        # get low-dimensional estimator based on active set and current data
        sol.fit <- PH.SP(
          yobs    = yobs[idx.sample],
          delta   = delta[idx.sample],
          X       = X[idx.sample,X.active.idx,drop=FALSE],
          auxinfo = auxinfo.refine,
          iindex  = iindex.refine,
          control = c(control,list(penalty=list(do=TRUE,inv=TRUE)))
        )

        # store the fitted results
        alp.list        <- c(alp.list,list(sol.fit$coefficient[,1]))
        alp.VarCov.list <- c(alp.VarCov.list,list(sol.fit$alp.VarCov))
        tau.list        <- c(tau.list,list(sol.fit$tau))

      }

    }

    # combined estimator
    alp.all <- matrix(unlist(lapply(1:dsplit$times,function(idsplit){
      (alp.list[[2*idsplit-1]]+alp.list[[2*idsplit]])/2
    })),ncol=dsplit$times)
    tau.all <- matrix(unlist(lapply(1:dsplit$times,function(idsplit){
      (tau.list[[2*idsplit-1]]+tau.list[[2*idsplit]])/2
    })),ncol=dsplit$times)
    alp.VarCov.all <- array(unlist(lapply(1:dsplit$times,function(idsplit){
      (alp.VarCov.list[[2*idsplit-1]]+alp.VarCov.list[[2*idsplit]])/2
    })),dim=c(palp,palp,dsplit$times))
    if(dsplit$combine == "mean"){
      alp        <- apply(alp.all,1,mean,na.rm=TRUE)
      alp.VarCov <- apply(alp.VarCov.all,c(1,2),mean,na.rm=TRUE)
      tau        <- apply(tau.all,1,mean,na.rm=TRUE)
    }else if(dsplit$combine == "median"){
      alp        <- apply(alp.all,1,median,na.rm=TRUE)
      alp.VarCov <- apply(alp.VarCov.all,c(1,2),median,na.rm=TRUE)
      tau        <- apply(tau.all,1,median,na.rm=TRUE)
    }
    alp.se     <- sqrt( diag(alp.VarCov)/N )

  }else{

    # get high-dimensional estimator based on decorrelation
    sol.fit <- PH.SP(
      yobs    = yobs,
      delta   = delta,
      X       = X,
      auxinfo = auxinfo,
      iindex  = iindex,
      control = c(control,list(penalty=list(do=TRUE,inv=FALSE)))
    )

    # store the fitted results
    alp        <- sol.fit$coefficient[,1]
    alp.VarCov <- sol.fit$alp.VarCov
    alp.se     <- sqrt(diag(alp.VarCov)/N )
    tau        <- sol.fit$tau

  }

  ## do inference and combine results into pre-specified style
  coefficient <- data.frame(Est       = alp,
                            SE        = alp.se,
                            zvalue    = alp/alp.se,
                            pvalue    = 2*(1-pnorm(abs(alp/alp.se))),
                            row.names = bet.names[iindex])

  ## output
  out <- list(
    coefficient = coefficient,
    alp.VarCov  = alp.VarCov,
    tau         = tau,
    control     = control,
    N           = N
  )
  return(out)

}

###############################################################################

###############################################################################
# Cox proportional hazards model (PH)
#   - with auxiliary coefficients from external submodels (CE)
#   - Techniques:
#     + Euclidean control variate (ECV) to synthesis information
#     + adaptive penalization (AP) to guarantee transprotability
###############################################################################

#' @title Cox proportional hazards model with auxiliary information in the type
#'   of coefficients from external submodels
#'
#' @description Fit the Cox proportional hazards (PH) model with auxiliary information
#'   in the type of coefficients from external submodels (CE).
#'
#' @param yobs a numeric vector that indicated the observed survival times.
#' @param delta a numeric vector that stores the right-censoring indicators
#'   (1 = event of interest happens, and 0 = censoring).
#' @param X a matrix that includes all interested covariates.
#' @param auxinfo a list used to store information about the inputted auxiliary
#'   information. Specifically, it contains the following elements:
#'   \describe{
#'     \item{\code{aux}}{ a list that includes information about historical
#'       aggregated statistics. It is a list of lists, and each sub-list represents
#'       auxiliary information from the same source study. In each source study,
#'       it contains the following elements:
#'       \describe{
#'         \item{\code{ce}}{ auxiliary regression coefficients.}
#'         \item{\code{ceV}}{ variance-covariance matrix of \code{ce}.}
#'         \item{\code{idx}}{ indices of this submodel's covariates in \code{X}.}
#'         \item{\code{M}}{ an integer that indicate the sample size of the external study.}
#'       } }
#'     \item{\code{hetero}}{ a logical value. If it is \code{TRUE}, the penalization
#'       will be applied to identify the potential heterogeneous auxiliary information
#'       and make a refinement to the final estimator.}
#'   }
#' @param iindex indices of our interested coefficients for covariates in \code{X}.
#'   The default is \code{NULL} and all coefficients will be displayed.
#' @param control indicates more detailed control of the underlying model fitting
#'   procedures. It is a list of the following three arguments:
#'   \describe{
#'     \item{\code{maxit}}{ a positive integer that denotes the maximum iteration
#'       number in optimization. The default value is \code{10000}.}
#'     \item{\code{eps}}{ a positive small numeric value that denotes the tolerance
#'       for convergence. The default value is \code{1e-6}.}
#'     \item{\code{penalty}}{ a list that contains information about the decorrelation steps.
#'       The argument \code{do} indicates whether the penalization should be used,
#'       and the argument \code{inv} controls whether the penalization in the decorrelation
#'       setup should be used.}
#'   }
#'
#' @return
#'   A list of fitted results is returned.
#'   Within this outputted list, the following elements can be found:
#'   \describe{
#'     \item{\code{coefficient}}{ estimation and inference of regression coefficients.}
#'     \item{\code{tau}}{ estimated degrees of population discrepency.}
#'   }
#'
#' @export PH.CE

PH.CE <- function(
    yobs,
    delta,
    X,
    auxinfo   = list(aux = NULL, hetero = TRUE),
    iindex    = NULL,
    control   = list(
      maxit   = 1e4,
      eps     = 1e-6,
      penalty = list(do = FALSE, inv = TRUE)
    )
){

  ## Specify several basic elements
  N         <- length(yobs)   # sample size in internal data
  pbet      <- ncol(X)   # number of parameters
  bet.names <- colnames(X)
  if(is.null(iindex)==TRUE){iindex <- 1:pbet}
  gamidx    <- c(1:pbet)[-iindex]
  palp      <- length(iindex)
  pgam      <- pbet-palp

  ## original estimates (without auxiliary information) and other preparations

  # initial estiamtor with different penalties
  if(control$penalty$do == FALSE){

    bet.init <- PH(yobs=yobs,delta=delta,X=X,Var=FALSE,
                   eps=control$eps,maxit=control$maxit)$res[,1]

  }else{

    bet.init <- as.numeric(coef(glmnet::cv.glmnet(
      X,survival::Surv(yobs,delta),family="cox",
      alpha=1,thresh=control$eps,maxit=control$maxit
    ),s="lambda.min"))

  }

  # basic preparation
  InfoMScore   <- PH.ScoreInfoM(bet=bet.init,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
  Score        <- InfoMScore$Score
  InfoM        <- InfoMScore$InfoM
  InfoM.alp    <- InfoM[iindex,iindex,drop=FALSE]
  InfoM.gamalp <- InfoM[-iindex,iindex,drop=FALSE]

  # solve for the optimal direction
  if(control$penalty$inv == TRUE){

    if(pgam > 0){
      InfoM.gam <- InfoM[-iindex,-iindex,drop=FALSE]
      w         <- solve(InfoM.gam)%*%InfoM.gamalp
    }else{
      w <- matrix(0,nrow=pgam,ncol=palp)
    }

  }else{

    if(control$penalty$do == FALSE & pgam > 0){

      InfoM.gam <- InfoM[-iindex,-iindex,drop=FALSE]
      w         <- solve(InfoM.gam)%*%InfoM.gamalp

    }else{

      w <- matrix(0,nrow=pgam,ncol=palp)
      if(pgam > 0){

        # prepare matrices
        InfoM.gam    <- InfoM[-iindex,-iindex,drop=FALSE]
        InfoM.gamalp <- InfoM[-iindex,iindex,drop=FALSE]

        # get w
        tunparaws <- exp(seq(log(0.005),log(5),length.out=50))*sqrt(log(pbet)/N)
        for(j in 1:palp){

          # - obtain solutions at all tuning parameters
          w.j.all <- array(0,dim=c(pgam,length(tunparaws)))
          for(itunparaw in 1:length(tunparaws)){
            tunparaw <- tunparaws[itunparaw]
            numit <- 1
            w.old <- rep(0,pgam)
            repeat{
              grad <- as.vector(InfoM.gam%*%w.old-InfoM.gamalp[,j])
              w.temp <- w.old - 0.1*grad
              w.j <- Soft.Threshold(theta=w.temp,lam=0.1*tunparaw)
              if( max(abs(w.j-w.old))>control$eps & numit<control$maxit ){
                w.old <- w.j
                numit <- numit + 1
              }else{
                break
              }
            }
            w.j.all[,itunparaw] <- w.j
          }

          # - get the optimal solution based on the information criterion
          ICs <- sapply(1:length(tunparaws),function(itunparaw){
            w.j.itunparaw  <- w.j.all[,itunparaw]
            df.j.itunparaw <- sum(abs(w.j.itunparaw)>control$eps)
            t(w.j.itunparaw)%*%InfoM.gam%*%w.j.itunparaw-2*InfoM.gamalp[,j]%*%w.j.itunparaw +
              2*df.j.itunparaw/N
          })
          w[,j] <- as.matrix(w.j.all[,which.min(ICs)])

        }

      }

    }

  }
  v <- rbind(diag(palp),-w)

  # get original estimator of alpha
  InfoM.inv <- solve(InfoM.alp-t(w)%*%InfoM.gamalp)
  if(control$penalty$do == FALSE){
    alp.ori <- bet.init[iindex]
  }else{
    Score.alp <- Score[iindex]
    Score.gam <- Score[-iindex]
    Score.DC  <- as.vector(Score.alp-t(w)%*%Score.gam)
    alp.ori   <- bet.init[iindex] + as.vector(InfoM.inv%*%Score.DC) # 改成直接求解
  }

  ## fit the model (if auxiliary information exists)
  if(is.null(auxinfo$aux)==TRUE){

    alp        <- alp.ori
    alp.VarCov <- InfoM.inv
    alp.se     <- sqrt(diag(alp.VarCov)/N)
    tau        <- NA
    cvinfo     <- NULL

  }else{

    # process the auxiliary information and extract necessary elements
    aux.process <-  PH.CE.Process(
      yobs=yobs,delta=delta,X=X,auxinfo=auxinfo,control=control)
    V         <- aux.process$V
    cva       <- aux.process$cva
    cv.Influs <- aux.process$cv.Influs
    V.inv     <- aux.process$V.inv
    tau       <- aux.process$tau

    # covariance matrix of initial estimator and control variate
    Score.Influs   <- PH.Influence.EE(yobs=yobs,delta=delta,X=X,bet=bet.init)
    alp.ori.Influs <- Score.Influs%*%v%*%t(InfoM.inv)
    Gam            <- t(cv.Influs)%*%alp.ori.Influs/N

    # define the improved estimator formally
    Gam.V.inv <- t(V.inv%*%Gam)
    alp       <- as.vector(alp.ori-Gam.V.inv%*%cva)

    # estimator of variance (using oracle properties directly)
    homo.idx <- (tau==0)
    if(any(homo.idx)){
      V.inv.homo <- solve(V[homo.idx,homo.idx,drop=F])
      Gam.homo   <- Gam[homo.idx,,drop=F]
      alp.VarCov <- InfoM.inv - t(Gam.homo)%*%V.inv.homo%*%Gam.homo
      Gam.V.inv  <- t(V.inv.homo%*%Gam.homo)
    }else{
      alp.VarCov <- InfoM.inv
      Gam.V.inv  <- NULL
    }
    alp.se <- sqrt( diag(alp.VarCov)/N )
    alp.se

    # information related to control variates
    cvinfo <- list(
      cv.Influs = cv.Influs,
      Gam.V.inv = Gam.V.inv
    )

  }

  ## do inference and combine results into pre-specified style
  coefficient <- data.frame(Est       = alp,
                            SE        = alp.se,
                            zvalue    = alp/alp.se,
                            pvalue    = 2*(1-pnorm(abs(alp/alp.se))),
                            row.names = bet.names[iindex])

  ## output
  out <- list(
    coefficient = coefficient,
    tau         = tau,
    alp.VarCov  = alp.VarCov,
    control     = control,
    N           = N,
    cvinfo      = cvinfo
  )
  return(out)

}

###############################################################################

###############################################################################
# High-dimensional Cox proportional hazards model (HPH)
#   - with auxiliary coefficients from external submodels (CE)
#   - Techniques:
#     + Euclidean control variate (ECV) to synthesis information
#     + adaptive penalization (AP) to guarantee transprotability
#     + decorrelation (DC) to stabilize penalized estimators
#     + data splitting (DS) to allow for high dimensional covariates
###############################################################################

#' @title High-Dimensional Cox proportional hazards model with auxiliary information
#'   in the type of coefficients from external submodels
#'
#' @description Fit the High-dimensional Cox proportional hazards (HPH) model with
#'   auxiliary information in the type of coefficients from external submodels (CE).
#'
#' @param yobs a numeric vector that indicated the observed survival times.
#' @param delta a numeric vector that stores the right-censoring indicators
#'   (1 = event of interest happens, and 0 = censoring).
#' @param X a matrix that includes all interested covariates.
#' @param auxinfo a list used to store information about the inputted auxiliary
#'   information. Specifically, it contains the following elements:
#'   \describe{
#'     \item{\code{aux}}{ a list that includes information about historical
#'       aggregated statistics. It is a list of lists, and each sub-list represents
#'       auxiliary information from the same source study. In each source study,
#'       it contains the following elements:
#'       \describe{
#'         \item{\code{ce}}{ auxiliary regression coefficients.}
#'         \item{\code{ceV}}{ variance-covariance matrix of \code{ce}.}
#'         \item{\code{idx}}{ indices of this submodel's covariates in \code{X}.}
#'         \item{\code{M}}{ an integer that indicate the sample size of the external study.}
#'       } }
#'     \item{\code{hetero}}{ a logical value. If it is \code{TRUE}, the penalization
#'       will be applied to identify the potential heterogeneous auxiliary information
#'       and make a refinement to the final estimator.}
#'   }
#' @param dsplit a list used to control the process of data splitting. Specifically,
#'   it contains the following elements:
#'   \describe{
#'     \item{\code{do}}{ a logical that indicate whether the data splitting procedure
#'       should be used. The default value if \code{TRUE}.}
#'     \item{\code{times}}{ an integer that denotes the times of data splitting. The
#'       default value if \code{10}.}
#'     \item{\code{combine}}{ a character that specifies the method of combining
#'       different estimates in each data splitting procedure. It can be either
#'       \code{"median"} (the default) and \code{"mean"}.}
#'   }
#' @param iindex indices of our interested coefficients for covariates in \code{X}.
#'   The default is \code{NULL} and all coefficients will be displayed.
#' @param control indicates more detailed control of the underlying model fitting
#'   procedures. It is a list of the following three arguments:
#'   \describe{
#'     \item{\code{maxit}}{ a positive integer that denotes the maximum iteration
#'       number in optimization. The default value is \code{10000}.}
#'     \item{\code{eps}}{ a positive small numeric value that denotes the tolerance
#'       for convergence. The default value is \code{1e-6}.}
#'   }
#'
#' @return
#'   A list of fitted results is returned.
#'   Within this outputted list, the following elements can be found:
#'   \describe{
#'     \item{\code{coefficient}}{ estimation and inference of regression coefficients.}
#'     \item{\code{tau}}{ estimated degrees of population discrepency.}
#'   }
#'
#' @export HPH.CE

HPH.CE <- function(
    yobs,
    delta,
    X,
    auxinfo = list(aux = NULL, M = NULL, hetero = TRUE),
    iindex  = NULL,
    dsplit  = list(do = TRUE, times = 10, combine = "median"),
    control = list(eps = 1e-6, maxit = 1e4)
){

  # Specify several basic elements
  N         <- length(yobs)   # sample size in internal data
  bet.names <- colnames(X)
  if(is.null(iindex)==TRUE){iindex <- 1:ncol(X)}
  pbet <- ncol(X)
  palp <- length(iindex)

  # the key part of the method
  if(dsplit$do == TRUE){

    # repeat the splitting procedure for multiple times
    alp.list <- alp.VarCov.list <- tau.list <- list()
    for(idsplit in 1:dsplit$times){

      # split the whole data set into two sets
      idx.temp <- sort(sample(1:N,round(N/2)))
      idx.list <- list(idx.temp,(1:N)[-idx.temp])

      # do for all data sets
      for(idx.sample in idx.list){

        # obtain active set from the remaining data
        bet.active <- as.numeric(coef(glmnet::cv.glmnet(
          X[-idx.sample,,drop=FALSE],survival::Surv(yobs[-idx.sample],delta[-idx.sample]),
          family="cox",alpha=1,thresh=control$eps,maxit=control$maxit
        ),s="lambda.min"))
        X.active.idx  <- sort(unique(c(iindex,which(abs(bet.active)>control$eps))))
        if(is.null(auxinfo$aux)==FALSE){
          X.active.idx <- sort(unique(c(
            X.active.idx,unlist(lapply(auxinfo$aux,function(auxi){auxi$idx})))))
          for(isource in 1:length(auxinfo$aux)){
            auxinfo$aux[[isource]]$idx <- which(auxinfo$aux[[isource]]$idx %in% X.active.idx)
          }
        }
        iindex.refine <- which(iindex %in% X.active.idx)

        # get low-dimensional estimator based on active set and current data
        sol.fit <- PH.CE(
          yobs    = yobs[idx.sample],
          delta   = delta[idx.sample],
          X       = X[idx.sample,X.active.idx,drop=FALSE],
          auxinfo = auxinfo,
          iindex  = iindex.refine,
          control = c(control,list(penalty=list(do=TRUE,inv=TRUE)))
        )

        # store the fitted results
        alp.list        <- c(alp.list,list(sol.fit$coefficient[,1]))
        alp.VarCov.list <- c(alp.VarCov.list,list(sol.fit$alp.VarCov))
        tau.list        <- c(tau.list,list(sol.fit$tau))

      }

    }

    # combined estimator
    alp.all <- matrix(unlist(lapply(1:dsplit$times,function(idsplit){
      (alp.list[[2*idsplit-1]]+alp.list[[2*idsplit]])/2
    })),ncol=dsplit$times)
    tau.all <- matrix(unlist(lapply(1:dsplit$times,function(idsplit){
      (tau.list[[2*idsplit-1]]+tau.list[[2*idsplit]])/2
    })),ncol=dsplit$times)
    alp.VarCov.all <- array(unlist(lapply(1:dsplit$times,function(idsplit){
      (alp.VarCov.list[[2*idsplit-1]]+alp.VarCov.list[[2*idsplit]])/2
    })),dim=c(palp,palp,dsplit$times))
    if(dsplit$combine == "mean"){
      alp        <- apply(alp.all,1,mean)
      alp.VarCov <- apply(alp.VarCov.all,c(1,2),mean)
      tau        <- apply(tau.all,1,mean)
    }else if(dsplit$combine == "median"){
      alp        <- apply(alp.all,1,median)
      alp.VarCov <- apply(alp.VarCov.all,c(1,2),median)
      tau        <- apply(tau.all,1,median)
    }
    alp.se     <- sqrt(diag(alp.VarCov)/N )

  }else{

    # get high-dimensional estimator based on decorrelation
    sol.fit <- PH.CE(
      yobs    = yobs,
      delta   = delta,
      X       = X,
      auxinfo = auxinfo,
      iindex  = iindex,
      control = c(control,list(penalty=list(do=TRUE,inv=FALSE)))
    )

    # store the fitted results
    alp        <- sol.fit$coefficient[,1]
    alp.VarCov <- sol.fit$alp.VarCov
    tau        <- sol.fit$tau

  }

  ## do inference and combine results into pre-specified style
  coefficient <- data.frame(Est       = alp,
                            SE        = alp.se,
                            zvalue    = alp/alp.se,
                            pvalue    = 2*(1-pnorm(abs(alp/alp.se))),
                            row.names = bet.names[iindex])

  ## output
  out <- list(
    coefficient = coefficient,
    alp.VarCov  = alp.VarCov,
    tau         = tau,
    control     = control,
    N           = N
  )
  return(out)

}

###############################################################################








