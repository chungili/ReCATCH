#' ReCatch
#' @description
#' Fit a ReCatch model and predict categorical responses.
#' @param x Input tensor is a list of length n, where is the number of observations. Each element of the list is an array.
#' @param U Input covariate matrix of dimension \eqn{n\times q}, with q>ng, where ng where is the number of observations in class g and q is the number of covariates.
#' @param y Class label. For G class problems, y takes values in 1,...,G.
#' @param testx Input testing tensor of list. Each element of the list is a test case. When it is not provided, the function will predict responses on training tensor x.
#' @param testU Input testing covariate matrix. When it is not provided, the function will predict responses on training covariates U.
#' @param OptMethod A string to specif the method for calculating the tuning parameters. The default method, "pso", uses particle swarm optimization. "kf" refers to K-fold cross-validation. "loo" refers to leave-one-out cross-validation.
#' @param imbalanced Logical; if TRUE the F1 score is used as a measure of predictive performance; otherwise, accuracy is used.
#' @param setseed An integer value that serves as the random seed.
#' @return A list containing the predicted categorical response and the optimized tuning parameters.
#'   \item{pred}{Predicted response of newx for each lambda in model object.}
#'   \item{opt.lambda}{The optimized tuning parameters.}
#' @examples
#' # Generate a training dataset
#' training <- GenCovTensor(q = 50, ng=c(40,40), pk = c(30,36,30), beta2 = 0.8, phi2 = 0.6, alpha2 = 1)
#' testing <- GenCovTensor(q = 50, ng=c(500,500), pk = c(30,36,30), beta2 = 0.8, phi2 = 0.6, alpha2 = 1)
#' x <- training$x; U <- training$U; y <- training$y
#' testx <- testing$x; testU <- testing$U
#' # Fit the model and calculate the predicted categorical response.
#' results <- ReCatch(x=x, U=U, y=y, testx=testx, testU=testU)
#' results$pred
#' @export
#'
ReCatch <- function(x=x, U=U, y=y,
                    testx=NULL, testU=NULL,
                    OptMethod="pso", ...){
  imbalanced <- setseed <- c(...)
  opt.lambda <- switch(OptMethod,
                         kf=kf.ReCatch(x = x, U = U, y = y, ...),
                         loo=loo.ReCatch(x = x, U = U, y = y, ...),
                         pso=pso.ReCatch(x = x, U = U, y = y, ...) )
  fit <- est.ReCatch(x=x, U=U, y=y,
                     lambda1=opt.lambda$lambda1,
                     lambda2=opt.lambda$lambda2,
                     lambda3=opt.lambda$lambda3,
                     lambda4=opt.lambda$lambda4)
  pred <- pred.ReCatch(fit=fit, testx=testx, testU=testU)$pred
  return(list(pred=pred, opt.lambda=opt.lambda))
}

#' est.ReCatch
#' @description
#' Fit a ReCatch model
#'
#' @param x Input tensor is a list of length n, where is the number of observations. Each element of the list is an array.
#' @param U Input covariate is a matrix of dimension \eqn{n\times q}, with q>n, where n where is the number of observations and q is the number of covariates.
#' @param y Class label. For G class problems, y takes values in 1,...,G.
#' @param lambda1 The tuning parameter.
#' @param lambda2 The tuning parameter.
#' @param lambda3 The tuning parameter.
#' @param lambda4 The tuning parameter.
#'
#' @return A list of parameters for the fitted model.
#' @examples
#' # Generate a training dataset
#' training <- GenCovTensor(q = 25, ng=c(500,500), pk = c(30,36,30), beta2 = 0.8, phi2 = 0.6, alpha2 = 1)
#' x <- training$x; U <- training$U; y <- training$y
#' fit <- est.ReCatch(x=x, U=U, y=y, lambda1=0.5, lambda2=0.5, lambda3=1.5, lambda4=2.0)
#' @export
#'
est.ReCatch <- function(x=x, U=U, y=y, lambda1=0, lambda2=0, lambda3=0, lambda4=1){
  library(catch)
  library(tensr)
  library(rTensor)
  #lambda1 <- 0
  #lambda2 <- 0
  #lambda3 <- 0
  #lambda4 <- 1
  #obj <- catch(x=x, z=U, y=y, lambda = lambda4, testx = x, testz = U)
  #objt <- ReAdjten(x=x, z=U, y=y, testx =testx, testz = testU, lambda3 = lambda3)
  pk <- dim(x[[1]])
  objt <- ReAdjten(x=x, z=U, y=y, lambda3 = lambda3)
  invisible(capture.output( objtsda <- tryCatch(catch:::tsda(objt$xres, y, lambda = lambda4),
                                                error=function(e) {NULL} ) ))
  if (is.null(objtsda)) {
    betahat <- matrix(rep(0, prod(pk)), ncol = 1)
  } else{
    betahat <- objtsda$beta$`1`
  }
  # X: tensor, array, c(pk, n)
  # U: Covariates, matrix, n*q
  # y: class labels, vector, n*1
  pk <- dim(x[[1]])
  nvars <- prod(pk)
  M <- length(pk)
  n <- length(y)
  ng <- as.vector(table(y))
  pig<- ng/n
  G <- length(pig)
  q <- dim(U)[2]
  # centering the data
  MatX<-matrix(0, n, nvars) # n*(p1*...pM)

  Ubar <- array(list(), G) # ng*q
  MatX <- t(sapply(x, FUN = as.vector))

  MatXbar <- matrix(NA, G, nvars)
  MatXTilde <- matrix(NA, n, nvars)
  phihat<-array(list(), G)
  Phihat<-array(list(), G)
  Matphi <- matrix(NA, G, q)
  UtildeArray <- array(list(), G)
  Utilde <- c() # n*q
  PhiTemp <- matrix(0, q, q)
  for (g in 1:G){
    # tensor
    MatXbar[g, ]<-colMeans(MatX[y==g,])
    MatXTilde[y==g,] <- sweep( MatX[y==g,], 2, colMeans(MatX[y==g,]) )
    phihat[[g]]<- matrix(apply(U[y==g,], 2, mean), ncol=1)
    Phihat[[g]] <- cov(U[y==g,])
    PhiTemp <- PhiTemp + cov(U[y==g,])*(ng[g]-1)
    Matphi[g,] <- phihat[[g]]
    UtildeArray[[g]] <- matrix(scale(U[y==g,], center = T, scale = F), ncol = q)
    Utilde <- rbind(Utilde, UtildeArray[[g]])
  }
  PhiPool <- PhiTemp/(sum(ng)-1)
  Philambda1 <- array(list(), G)
  Philambda12 <- array(list(), G)
  Philambda12Inv <- array(list(), G)
  for (g in 1:G) {
    Philambda1[[g]] <- (1-lambda1)*Phihat[[g]]+lambda1*PhiPool
    Philambda12[[g]] <- (1-lambda2)*Philambda1[[g]]+(lambda2/G)*sum(diag(Philambda1[[g]]))*diag(q)
    Philambda12Inv[[g]] <- solve(Philambda12[[g]])
  }

  Xbar <- fold(mat = MatXbar, row_idx = (M+1), col_idx = seq(1:M),
               modes = c(pk, G) ) # p1*...*pM*G
  XTilde <- fold(mat = MatXTilde, row_idx = (M+1), col_idx = seq(1:M),
                 modes = c(pk, n) )  # p1*...*pM*n
  X <- fold(mat = MatX, row_idx = (M+1), col_idx = seq(1:M),
            modes = c(pk, n) )  # p1*...*pM*n
  A <- objt$A
  P<-solve(t(Utilde)%*%Utilde+lambda3*diag(q))%*%t(Utilde)
  alphahat <- rTensor::ttm(tnsr =  XTilde, mat = P, m = (M+1) ) # p1*...*pM*q
  Xhat<-rTensor::ttm(tnsr = alphahat, mat = U, m = (M+1) ) # p1*...*pM*n
  adjX <- X-Xhat # p1*...*pM*n
  XhatTilde<-rTensor::ttm(tnsr = alphahat, mat = Utilde, m = (M+1) ) # p1*...*pM*n
  adjXTilde <- XTilde-XhatTilde

  muhat <- Xbar-rTensor::ttm(tnsr = alphahat, mat = Matphi, m = (M+1) ) # p1*...*pM*q
  muhatUnfold <- unfold(tnsr = muhat, row_idx = (M+1), col_idx = c(1:M) )@data # G*(p1***pM)
  return(list(adjX=adjX, Utilde=Utilde, alphahat=alphahat, betahat=betahat, Xhat=Xhat,
              Matphi=Matphi, XTilde=XTilde ,adjXTilde=adjXTilde, Xbar=Xbar, muhat=muhat, muhatUnfold=muhatUnfold,
              RePhi=Philambda12, RePhiInv=Philambda12Inv, phihat=phihat, A=A, pig=pig, x=x, U=U, y=y))
}

#' pred.ReCatch
#' @description
#' Predict categorical responses
#'
#' @param fit A fitted ReCatch model.
#' @param testx Input testing tensor of list.
#' @param testU Input testing covariate matrix.
#'
#' @return pred The predicted categorical responses.
#' @export
#'
#' @examples
#' # Generate a training dataset
#' training <- GenCovTensor(q = 25, ng=c(500,500), pk = c(30,36,30), beta2 = 0.8, phi2 = 0.6, alpha2 = 1)
#' testing <- GenCovTensor(q = 25, ng=c(40,40), pk = c(30,36,30), beta2 = 0.8, phi2 = 0.6, alpha2 = 1)
#' x <- training$x; U <- training$U; y <- training$y
#' testx <- testing$x; testU <- testing$U
#' fit <- est.ReCatch(x=x, U=U, y=y, lambda1=0.5, lambda2=0.5, lambda3=1.5, lambda4=2.0)
#' pred <- pred.ReCatch(fit=fit, testx=testx, testU=testU)$pred
pred.ReCatch <- function(fit=fit, testx=NULL, testU=NULL) {
  pk <- dim(fit$x[[1]])
  nvars <- prod(pk)
  M <- length(pk)
  n <- length(fit$y)
  ng <- as.vector(table(fit$y))
  pig<- ng/n
  G <- length(pig)
  q <- dim(fit$U)[2]
  #testn <- ifelse(!is.matrix(testU), 1, dim(U)[1] )

  # centering the data
  if (!is.null(testx)) {
    if (length(testx)==1) {
      testxn <- 1
      testU <- matrix(testU, nrow = 1, ncol = q)
      testMatX<-matrix(0, testxn, nvars)
      for (i in 1:testxn) testMatX[i,] <- matrix(testx[[i]], nvars, 1)
      testX <- fold(mat = testMatX, row_idx = (M+1), col_idx = seq(1:M),
                    modes = c(pk, testxn) )  # p1*...*pM*n
      testXhat<-Myttv(tnsr = fit$alphahat, v = testU, m = (M+1) ) # p1*...*pM*n
      testX <- fold(mat = matrix(testX@data, nrow = pk[1]), row_idx = 1, col_idx = c(2:M), modes=pk )
      testadjX <- testX-testXhat # p1*...*pM*n
    } else {
      #testxn <- dim(testx)
      testxn <- length(testx)
      testMatX<-matrix(0, testxn, nvars)
      for (i in 1:testxn) testMatX[i,] <- matrix(testx[[i]], nvars, 1)
      testX <- fold(mat = testMatX, row_idx = (M+1), col_idx = seq(1:M),
                    modes = c(pk, testxn) )  # p1*...*pM*n
      testXhat<-rTensor::ttm(tnsr = fit$alphahat, mat = testU, m = (M+1) ) # p1*...*pM*n
      testadjX <- testX-testXhat # p1*...*pM*n
    }
  }
  disScore <- matrix(0, nrow = ifelse(is.null(testU), n, testxn), ncol = G)
  for (g in 2:G ){
    beta <- as.vector(fit$betahat[,g-1])
    gammag <- fit$RePhiInv[[g]]%*%fit$phihat[[g]]- fit$RePhiInv[[1]]%*%fit$phihat[[1]]
    Dg <- fit$RePhiInv[[g]] - fit$RePhiInv[[1]]
    if ( is.null(testU) ) {
      Score1 <- -0.5*diag(fit$U%*%Dg%*%t(fit$U))
      Score2 <- fit$U%*%gammag
      bg <-apply(k_unfold(tnsr = fit$adjX, m = (M+1) )@data, 1, function(x) sum(x*fit$beta))

    } else {
      Score1 <- -0.5*diag(testU%*%Dg%*%t(testU))
      Score2 <- testU%*%gammag
      if (testxn==1) {
        bg <- sum(beta*c(testadjX@data))
      }else{
        bg <-apply(k_unfold(tnsr = testadjX, m = (M+1) )@data, 1, function(x) sum(x*fit$beta))
      }
    }
    mugPlus1 <- fit$muhatUnfold[g,] + fit$muhatUnfold[1,]
    ag.Tensor <- sum( fit$beta*mugPlus1 )
    ag.Cov<-as.numeric( log(det(fit$RePhi[[g]])/det(fit$RePhi[[1]]) )
                        + t(gammag)%*%(fit$phihat[[g]]+fit$phihat[[1]])
                        - t(fit$phihat[[g]])%*%Dg%*%fit$phihat[[1]]
                        -2*log(fit$pig[g]/fit$pig[1]) )
    ag <- ag.Cov + ag.Tensor
    disScore[,g]<- Score1+Score2-0.5*ag+bg
  }
  yhat <- apply(disScore, 1, which.max)
  return(list(pred=yhat))
}

#' loo.ReCatch
#' @description
#' Leave-one-out cross-validation optimization algorithm for the ReCatch model.
#'
#' @param x Input tensor.
#' @param U Input covariate matrix.
#' @param y Class label.
#' @param size.grid A single integer for the maximum number of parameter value combinations returned.
#' @param setseed An integer value that serves as the random seed.
#' @param lambda3Range The parameter interval for \eqn{\lambda_3}.
#' @param imbalanced Logical; if TRUE the F1 score is used as a measure of predictive performance; otherwise, accuracy is used.
#'
#' @return A list of the best set of parameters found.
#' @export
#'
#getLambda <- function(x, U, y, size.grid=100,
loo.ReCatch <- function(x, U, y, size.grid=100,
                      setseed=NULL, lambda3Range=c(0.1, 3.0), imbalanced=FALSE, poslab="2"){
  #####
  library(foreach)
  library(caret)
  #library(doParallel)
  library(parallel)
  library(doSNOW)
  ReAdjten <- function(x, z, y, testx=NULL, testz=NULL, is.centered=FALSE, lambda3=0){
    n <- length(x)
    dimen <- dim(x[[1]])
    nvars <- prod(dimen)
    vec_x <- matrix(0, nvars, n)
    for (i in 1:n){
      vec_x[,i] <- matrix(x[[i]], nvars, 1)
    }
    vec_x = t(vec_x)

    if (!is.null(testx)){
      tesize <- length(testx)
      vec_testx <- matrix(0, nvars, tesize)
      for (i in 1:tesize){
        vec_testx[,i] <- matrix(testx[[i]], nvars,1)
      }
      vec_testx <- t(vec_testx)
    }
    z <- as.matrix(z) #n*q
    cz <- z
    n <- dim(z)[1]
    q <- dim(z)[2]
    cvecx <- vec_x
    nclass <- as.integer(length(unique(y)))
    if (is.centered==FALSE){
      for (i in 1:nclass){
        if (q>1) {
          cz[y == i,] = sweep(z[y == i,], 2, colMeans(z[y == i,]))
        }else {
          cz[y==i]=z[y==i]-mean(z[y==i])
        }
        cvecx[y == i,] = sweep(vec_x[y == i,], 2, colMeans(vec_x[y == i,]))
      }
    }
    c <- solve( (t(cz) %*% cz)+lambda3*diag(q) ) %*% t(cz)
    c = c %*% cvecx
    vec_xres <- vec_x-z%*%c
    A = (t(cz) %*% cz)+lambda3*diag(q)

    xres = array(list(),n)
    for (i in 1:n){
      xres[[i]] = array(vec_xres[i,],dim=dimen)
    }

    if (!is.null(testx)){
      vec_testxres = vec_testx-testz %*% c
      testxres = array(list(),tesize)
      for (i in 1:tesize){
        testxres[[i]] = array(vec_testxres[i,], dim=dimen)
      }
    }else{
      testxres=NULL
    }

    muz = matrix(0, nrow=q, ncol=(nclass-1))    # muz : U_k-U_1  (k>1)
    for (i in 2:nclass){
      if (q>1){muz[,(i-1)] = apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)
      } else {
        muz[i-1] <- mean(z[y==i])-mean(z[y==1])
      }
    }
    #gamma <- solve(cov(z)) %*% muz

    return( list(xres=xres, A=A, testxres=testxres, alpha=c) )
  }

  tsda <- function(x, y, nlambda = 100, lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2, 1E-03), lambda = NULL, dfmax = nobs, pmax = min(dfmax * 2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+05, sml = 1e-06, verbose = FALSE, perturb = NULL) {
    ## data setup
    this.call <- match.call()
    tmp <- prept(x, y)
    ldim <- length(dim(x[[1]])) #m dimensional
    dimen <- dim(x[[1]])
    maxd <- max(dimen)
    sigma <- as.matrix(tmp$sigma)


    sigmalist<-array(list(),ldim)
    for (i in 1:ldim){
      tmpsigma <- matrix(sigma[i,],maxd,maxd)
      if (!is.null(perturb)){
        diag(tmpsigma) <- diag(tmpsigma)+perturb
        sigma[i,] <- matrix(tmpsigma,nrow=1)
      }
      sigmalist[[i]] <- tmpsigma[1:dimen[i],1:dimen[i]]
    }

    delta <- as.matrix(tmp$delta)
    mu <- as.matrix(tmp$mu)
    prior <- tmp$prior

    nvars <- as.integer(prod(dimen))
    nobs <- length(y)
    nclass <- as.integer(length(unique(y)))
    vnames <- colnames(x)
    if (is.null(vnames))
      vnames <- paste("V", seq(nvars), sep = "")
    ng <- as.integer(dim(delta)[1])
    ## parameter setup
    if (length(pf) != nvars)
      stop("The size of penalty factor must be same as the number of input variables")
    maxit <- as.integer(maxit)
    verbose <- as.integer(verbose)
    sml <- as.double(sml)
    pf <- as.double(pf)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax*ng)
    pmax <- as.integer(pmax)
    ## lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
      if (lambda.factor >= 1)
        stop("lambda.factor should be less than 1")
      flmin <- as.double(lambda.factor)
      ulam <- double(1)  #ulam=0 if lambda is missing
    } else {
      # flmin=1 if user define lambda
      flmin <- as.double(1)
      if (any(lambda < 0))
        stop("lambdas should be non-negative")
      ulam <- as.double(rev(sort(lambda)))  #lambda is declining
      nlam <- as.integer(length(lambda))
    }
    ## call Fortran core
    fit <<- .Fortran("catch1", obj = double(nlam), ng, nvars, ldim, dimen,maxd, as.double(sigma), as.double(delta),
                     pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1),
                     theta = double(pmax * ng * nlam), itheta = integer(pmax), ntheta = integer(nlam),
                     alam = double(nlam), npass = integer(1), jerr = integer(1))
    ## output
    outlist <- formatoutputt(fit, maxit, pmax, dfmax, nvars, vnames, ng)
    outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr,
                               sigma <- sigmalist, delta = delta, mu = mu, prior = prior, call = this.call))
    if (is.null(lambda))
      outlist$lambda <- lamfix(outlist$lambda)
    class(outlist) <- c("catchobj")
    outlist

  }

  formatoutputt <- function(fit, maxit, pmax, dfmax, nvars, vnames, ng) {
    nalam <- fit$nalam
    ntheta <- fit$ntheta[seq(nalam)]
    nthetamax <- max(ntheta)
    lam <- fit$alam[seq(nalam)]
    obj <- fit$obj[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    resnames <- paste("delta", seq(ng), sep = "")

    errmsg <- errt(fit$jerr, maxit, pmax,dfmax)  ### error messages from fortran
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))

    dd <- c(nvars, ng)
    df <- rep(0, nalam)
    if (nthetamax > 0) {
      ja <- fit$itheta[seq(nthetamax)]
      oja <- order(ja)
      ja <- rep(ja[oja], ng)
      itheta <- cumsum(c(1, rep(nthetamax, ng)))
      pos <- rep(1:nalam, each = ng * pmax)
      theta <- split(fit$theta[seq(ng * pmax * nalam)], pos)
      for (l in 1:nalam) {
        theta[[l]] <- matrix(theta[[l]], pmax, ng, byrow = TRUE)[seq(nthetamax),
                                                                 , drop = FALSE]
        df[l] <- sum(rowSums(abs(theta[[l]])) != 0)
        theta[[l]] <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames,
                                                                 resnames), x = as.vector(theta[[l]][oja, ]), p = as.integer(itheta -
                                                                                                                               1), i = as.integer(ja - 1))
      }
    } else {
      theta <- list()
      for (l in 1:nalam) {
        theta[[l]] <- zeromat(nvars, ng, vnames, resnames)
      }
      df <- rep(0, nalam)
    }
    list(beta = theta, df = df, dim = dd, lambda = lam, obj = obj)
  }

  # generate sigma, delta and mu from x, y
  prept <- function(x, y) {
    # data setup
    nclass <- as.integer(length(unique(y)))
    prior <- rep(0, nclass)
    for (k in 1:nclass) {
      prior[k] <- mean(y == k)
    }
    ldim <- length(dim(x[[1]])) #m dimensional
    p <- matrix(0,ncol=ldim,nrow=1)
    nvars <- 1
    for (i in 1:ldim){
      p[i] <- as.integer(dim(x[[1]])[i])
      nvars <- nvars*p[i]
    }
    nobs <- length(x)
    nres <- length(y)
    if (nres != nobs)
      stop("x and y have different number of observations")

    dim <- dim(x[[1]])
    num <- matrix(0,nrow=k,ncol=1)
    # estimate mean
    es_M <- array(list(),k)
    for (i in 1:k){
      es_M[[i]] <- array(0,dim=dim)
    }
    ng <- nclass-1
    for (i in 1:nobs){
      es_M[[y[i]]] <- es_M[[y[i]]]+x[[i]]
      num[[y[i]]] <- num[[y[i]]]+1
    }
    for (i in 1:k){
      es_M[[i]] <- es_M[[i]]/num[[i]]
    }

    #estimate sigma
    n <- nobs
    matx <- matrix(list(),n,1)
    es_sigma <- matrix(list(),ldim,1)
    sum_sigma <- matrix(list(),ldim,1)
    matmu <- matrix(list(),k,1)
    tr <- 1
    for (j in 1:ldim){
      for (i in 1:k){
        matmu[[i]] <- mat(es_M[[i]],j)
      }
      sum_sigma[[j]] <- matrix(0,nrow=dim[j],ncol=dim[j])
      for (i in 1:n){
        matx[[i]] <- mat(x[[i]],j)
        sum_sigma[[j]] <- sum_sigma[[j]]+(matx[[i]]-matmu[[y[i]]])%*%t(matx[[i]]-matmu[[y[i]]])
      }
      sum_sigma[[j]] <- sum_sigma[[j]]/(n-k)
      if (j<ldim){
        es_sigma[[j]] <- sum_sigma[[j]]/sum_sigma[[j]][1,1]
        tr <- tr*sum(diag(es_sigma[[j]]))
      }else
      {
        es_sigma[[j]] <- sum_sigma[[j]]/tr
      }
    }
    delta <- matrix(rep(0),nrow=nclass-1,ncol=prod(dim))
    for (i in 1:(nclass-1)){
      delta[i,] <- as.matrix(as.vector(es_M[[i+1]])-as.vector(es_M[[1]]))
    }

    maxd <- max(dim)
    sigma <- matrix(0,nrow=ldim,ncol=(maxd^2))
    for (i in 1:ldim){
      t <- matrix(0,nrow=maxd,ncol=maxd)
      t[1:dim[i],1:dim[i]] <- es_sigma[[i]]
      sigma[i,] <- as.vector(t)
    }
    outlist <- list(sigma=sigma, delta = delta, mu =es_M, prior = prior)
    outlist
  }

  errt <- function(n, maxit, pmax,dfmax) {
    if (n == 0)
      msg <- ""
    if (n > 0) {
      # fatal error
      if (n < 7777)
        msg <- "Memory allocation error; contact package maintainer"
      if (n == 10000)
        msg <- "All penalty factors are <= 0"
      n <- 1
      msg <- paste("in the fortran code -", msg)
    }
    if (n < 0) {
      # non fatal error
      if (n > -10000)
        msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=",
                     maxit, " iterations; solutions for larger lambdas returned.\n",
                     sep = "")
      if (n < -10000)
        msg <- paste("Number of nonzero coefficients along the path exceeds pmax=",
                     pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n",
                     sep = "")
      if (n < -20000)
        msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=",
                     dfmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n",
                     sep = "")
      n <- -1
    }
    list(n = n, msg = msg)
  }

  zeromat <- function(nvars, nalam, vnames, stepnames) {
    ca <- rep(0, nalam)
    ia <- seq(nalam + 1)
    ja <- rep(1, nalam)
    dd <- c(nvars, nalam)
    new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), x = as.vector(ca),
        p = as.integer(ia - 1), i = as.integer(ja - 1))
  }

  lamfix <- function(lam) {
    llam <- log(lam)
    lam[1] <- exp(2 * llam[2] - llam[3])
    lam
  }

  GridPoints<-function(size=100, method="LHC", setseed, lambda3Range){
    #Regular: Regular
    #LH: LatinHypercube
    require(dials)
    lambda1 <- new_quant_param(type ="double",range=c(0.01, 0.9), label = c(lambda1="lambda1"), inclusive = c(FALSE, FALSE))
    lambda2 <- new_quant_param(type ="double",range=c(0.01, 0.9), label = c(lambda2="lambda2"), inclusive = c(FALSE, FALSE))
    lambda3 <- new_quant_param(type ="double",range=lambda3Range, label = c(lambda3="lambda3"), inclusive = c(FALSE, FALSE))
    #lambda4 <- new_quant_param(type ="double",range=c(lambda4Range[1],lambda4Range[2]), label = c(lambda4="lambda4"), inclusive = c(FALSE, FALSE))
    if (!is.null(setseed)) set.seed(setseed)
    if (method=="Regular"){
      Lambda <- grid_regular(x = list(lambda1,lambda2,lambda3), levels=size)
    }else{
      Lambda <- grid_space_filling(list(lambda1,lambda2,lambda3), size=size)
    }
    return(Lambda)
  }

  estRecatch <- function(x=x, U=U, y=y, lambda1=0, lambda2=0, lambda3=0, lambda4=1){
    library(catch)
    library(tensr)
    library(rTensor)
    pk <- dim(x[[1]])
    objt <- ReAdjten(x=x, z=U, y=y, lambda3 = lambda3)
    invisible(capture.output( objtsda <- tryCatch(catch:::tsda(objt$xres, y, lambda = lambda4),
                                                  error=function(e) {NULL} ) ))
    if (is.null(objtsda)) {
      betahat <- matrix(rep(0, prod(pk)), ncol = 1)
    } else{
      #betahat <- objtsda$beta$`1`
      betahat <- objtsda$beta[[1]]
    }
    pk <- dim(x[[1]])
    nvars <- prod(pk)
    M <- length(pk)
    n <- length(y)
    ng <- as.vector(table(y))
    pig<- ng/n
    G <- length(pig)
    q <- dim(U)[2]
    # centering the data
    MatX<-matrix(0, n, nvars) # n*(p1*...pM)

    Ubar <- array(list(), G) # ng*q
    MatX <- t(sapply(x, FUN = as.vector))

    MatXbar <- matrix(NA, G, nvars)
    MatXTilde <- matrix(NA, n, nvars)
    phihat<-array(list(), G)
    Phihat<-array(list(), G)
    Matphi <- matrix(NA, G, q)
    UtildeArray <- array(list(), G)
    Utilde <- c() # n*q
    PhiTemp <- matrix(0, q, q)
    for (g in 1:G){
      # tensor
      MatXbar[g, ]<-colMeans(MatX[y==g,])
      MatXTilde[y==g,] <- sweep( MatX[y==g,], 2, colMeans(MatX[y==g,]) )
      phihat[[g]]<- matrix(apply(U[y==g,], 2, mean), ncol=1)
      Phihat[[g]] <- cov(U[y==g,])
      PhiTemp <- PhiTemp + cov(U[y==g,])*(ng[g]-1)
      Matphi[g,] <- phihat[[g]]
      UtildeArray[[g]] <- matrix(scale(U[y==g,], center = T, scale = F), ncol = q)
      Utilde <- rbind(Utilde, UtildeArray[[g]])
    }
    PhiPool <- PhiTemp/(sum(ng)-1)
    Philambda1 <- array(list(), G)
    Philambda12 <- array(list(), G)
    Philambda12Inv <- array(list(), G)
    for (g in 1:G) {
      Philambda1[[g]] <- (1-lambda1)*Phihat[[g]]+lambda1*PhiPool
      Philambda12[[g]] <- (1-lambda2)*Philambda1[[g]]+(lambda2/G)*sum(diag(Philambda1[[g]]))*diag(q)
      Philambda12Inv[[g]] <- solve(Philambda12[[g]])
    }

    Xbar <- fold(mat = MatXbar, row_idx = (M+1), col_idx = seq(1:M),
                 modes = c(pk, G) ) # p1*...*pM*G
    XTilde <- fold(mat = MatXTilde, row_idx = (M+1), col_idx = seq(1:M),
                   modes = c(pk, n) )  # p1*...*pM*n
    X <- fold(mat = MatX, row_idx = (M+1), col_idx = seq(1:M),
              modes = c(pk, n) )  # p1*...*pM*n
    A <- objt$A
    P<-solve(t(Utilde)%*%Utilde+lambda3*diag(q))%*%t(Utilde)
    alphahat <- rTensor::ttm(tnsr =  XTilde, mat = P, m = (M+1) ) # p1*...*pM*q
    Xhat<-rTensor::ttm(tnsr = alphahat, mat = U, m = (M+1) ) # p1*...*pM*n
    adjX <- X-Xhat # p1*...*pM*n
    XhatTilde<-rTensor::ttm(tnsr = alphahat, mat = Utilde, m = (M+1) ) # p1*...*pM*n
    adjXTilde <- XTilde-XhatTilde

    muhat <- Xbar-rTensor::ttm(tnsr = alphahat, mat = Matphi, m = (M+1) ) # p1*...*pM*q
    muhatUnfold <- unfold(tnsr = muhat, row_idx = (M+1), col_idx = c(1:M) )@data # G*(p1***pM)
    return(list(adjX=adjX, Utilde=Utilde, alphahat=alphahat, betahat=betahat, Xhat=Xhat,
                Matphi=Matphi, XTilde=XTilde ,adjXTilde=adjXTilde, Xbar=Xbar, muhat=muhat, muhatUnfold=muhatUnfold,
                RePhi=Philambda12, RePhiInv=Philambda12Inv, phihat=phihat, A=A, pig=pig, x=x, U=U, y=y))
  }
  estParm <- function(x=x, U=U, y=y, lambda1=0, lambda2=0, lambda3=0){
    library(catch)
    library(tensr)
    library(rTensor)
    # X: tensor, array, c(pk, n)
    # U: Covariates, matrix, n*q
    # y: class labels, vector, n*1
    pk <- dim(x[[1]])
    nvars <- prod(pk)
    M <- length(pk)
    n <- length(y)
    ng <- as.vector(table(y))
    pig<- ng/n
    G <- length(pig)
    q <- dim(U)[2]
    objt <- ReAdjten(x=x, z=U, y=y, lambda3 = lambda3)
    # centering the data
    MatX<-matrix(0, n, nvars) # n*(p1*...pM)

    Ubar <- array(list(), G) # ng*q
    MatX <- t(sapply(x, FUN = as.vector))

    MatXbar <- matrix(NA, G, nvars)
    MatXTilde <- matrix(NA, n, nvars)
    phihat<-array(list(), G)
    Phihat<-array(list(), G)
    Matphi <- matrix(NA, G, q)
    UtildeArray <- array(list(), G)
    Utilde <- c() # n*q
    PhiTemp <- matrix(0, q, q)
    for (g in 1:G){
      # tensor
      MatXbar[g, ]<-colMeans(MatX[y==g,])
      MatXTilde[y==g,] <- sweep( MatX[y==g,], 2, colMeans(MatX[y==g,]) )
      phihat[[g]]<- matrix(apply(U[y==g,], 2, mean), ncol=1)
      Phihat[[g]] <- cov(U[y==g,])
      PhiTemp <- PhiTemp + cov(U[y==g,])*(ng[g]-1)
      Matphi[g,] <- phihat[[g]]
      UtildeArray[[g]] <- matrix(scale(U[y==g,], center = T, scale = F), ncol = q)
      Utilde <- rbind(Utilde, UtildeArray[[g]])
    }
    PhiPool <- PhiTemp/(sum(ng)-1)
    Philambda1 <- array(list(), G)
    Philambda12 <- array(list(), G)
    Philambda12Inv <- array(list(), G)
    for (g in 1:G) {
      Philambda1[[g]] <- (1-lambda1)*Phihat[[g]]+lambda1*PhiPool
      Philambda12[[g]] <- (1-lambda2)*Philambda1[[g]]+(lambda2/G)*sum(diag(Philambda1[[g]]))*diag(q)
      Philambda12Inv[[g]] <- solve(Philambda12[[g]])
    }

    Xbar <- fold(mat = MatXbar, row_idx = (M+1), col_idx = seq(1:M),
                 modes = c(pk, G) ) # p1*...*pM*G
    XTilde <- fold(mat = MatXTilde, row_idx = (M+1), col_idx = seq(1:M),
                   modes = c(pk, n) )  # p1*...*pM*n
    X <- fold(mat = MatX, row_idx = (M+1), col_idx = seq(1:M),
              modes = c(pk, n) )  # p1*...*pM*n
    A <- objt$A
    P<-solve(t(Utilde)%*%Utilde+lambda3*diag(q))%*%t(Utilde)
    alphahat <- rTensor::ttm(tnsr =  XTilde, mat = P, m = (M+1) ) # p1*...*pM*q
    Xhat<-rTensor::ttm(tnsr = alphahat, mat = U, m = (M+1) ) # p1*...*pM*n
    adjX <- X-Xhat # p1*...*pM*n
    XhatTilde<-rTensor::ttm(tnsr = alphahat, mat = Utilde, m = (M+1) ) # p1*...*pM*n
    adjXTilde <- XTilde-XhatTilde

    muhat <- Xbar-rTensor::ttm(tnsr = alphahat, mat = Matphi, m = (M+1) ) # p1*...*pM*q
    muhatUnfold <- unfold(tnsr = muhat, row_idx = (M+1), col_idx = c(1:M) )@data # G*(p1***pM)
    return(list(adjX=adjX, Utilde=Utilde, alphahat=alphahat, Xhat=Xhat,
                Matphi=Matphi, XTilde=XTilde ,adjXTilde=adjXTilde, Xbar=Xbar, muhat=muhat, muhatUnfold=muhatUnfold,
                RePhi=Philambda12, RePhiInv=Philambda12Inv, phihat=phihat, A=A, pig=pig, x=x, U=U, y=y))
  }
  pred.ReCatch <- function(fit=fit, x=x, U=U, y=y, testx=NULL, testU=NULL) {
    pk <- dim(x[[1]])
    nvars <- prod(pk)
    M <- length(pk)
    n <- length(y)
    ng <- as.vector(table(y))
    pig<- ng/n
    G <- length(pig)
    q <- dim(U)[2]
    # centering the data
    if (!is.null(testx)) {
      if (length(testx)==1) {
        testxn <- 1
        testU <- matrix(testU, nrow = 1, ncol = q)
        testMatX<-matrix(0, testxn, nvars)
        for (i in 1:testxn) testMatX[i,] <- matrix(testx[[i]], nvars, 1)
        testX <- fold(mat = testMatX, row_idx = (M+1), col_idx = seq(1:M),
                      modes = c(pk, testxn) )  # p1*...*pM*n
        testXhat<-Myttv(tnsr = fit$alphahat, v = testU, m = (M+1) ) # p1*...*pM*n
        testX <- fold(mat = matrix(testX@data, nrow = pk[1]), row_idx = 1, col_idx = c(2:M), modes=pk )
        testadjX <- testX-testXhat # p1*...*pM*n
      } else {
        testxn <- length(testx)
        testMatX<-matrix(0, testxn, nvars)
        for (i in 1:testxn) testMatX[i,] <- matrix(testx[[i]], nvars, 1)
        testX <- fold(mat = testMatX, row_idx = (M+1), col_idx = seq(1:M),
                      modes = c(pk, testxn) )  # p1*...*pM*n
        testXhat<-rTensor::ttm(tnsr = fit$alphahat, mat = testU, m = (M+1) ) # p1*...*pM*n
        testadjX <- testX-testXhat # p1*...*pM*n
      }
    }
    disScore <- matrix(0, nrow = ifelse(is.null(testU), n, testxn), ncol = G)
    for (g in 2:G ){
      beta <- as.vector(fit$betahat[,g-1])
      gammag <- fit$RePhiInv[[g]]%*%fit$phihat[[g]]- fit$RePhiInv[[1]]%*%fit$phihat[[1]]
      Dg <- fit$RePhiInv[[g]] - fit$RePhiInv[[1]]
      if ( is.null(testU) ) {
        Score1 <- -0.5*diag(U%*%Dg%*%t(U))
        Score2 <- U%*%gammag
        bg <-apply(k_unfold(tnsr = fit$adjX, m = (M+1) )@data, 1, function(x) sum(x*fit$beta))

      } else {
        Score1 <- -0.5*diag(testU%*%Dg%*%t(testU))
        Score2 <- testU%*%gammag
        if (testxn==1) {
          bg <- sum(beta*c(testadjX@data))
        }else{
          bg <-apply(k_unfold(tnsr = testadjX, m = (M+1) )@data, 1, function(x) sum(x*fit$beta))
        }
      }
      mugPlus1 <- fit$muhatUnfold[g,] + fit$muhatUnfold[1,]
      ag.Tensor <- sum( fit$beta*mugPlus1 )
      ag.Cov<-as.numeric( log(det(fit$RePhi[[g]])/det(fit$RePhi[[1]]) )
                          + t(gammag)%*%(fit$phihat[[g]]+fit$phihat[[1]])
                          - t(fit$phihat[[g]])%*%Dg%*%fit$phihat[[1]]
                          -2*log(fit$pig[g]/fit$pig[1]) )
      ag <- ag.Cov + ag.Tensor
      disScore[,g]<- Score1+Score2-0.5*ag+bg
    }
    yhat <- apply(disScore, 1, which.max)
    return(list(pred=yhat))
  }
  #####

  #####
  myLambda <- GridPoints(size = size.grid, lambda3Range=lambda3Range, setseed=setseed)
  library(foreach)
  library(iterators)
  cat("The progress of Leave-One-Out cross-validation\n")
  cat("Imbalanced Data: ", imbalanced,"\n")
  cat("Random seed:", ifelse(is.null(setseed), FALSE, setseed), "\n")
  n_cores <- detectCores() - 2
  workers <- makeCluster(n_cores)
  registerDoSNOW(workers)
  pb <- txtProgressBar(0, size.grid, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  #####
  acc<-foreach(i=myLambda$lambda1, j=myLambda$lambda2,
               k=myLambda$lambda3,
               .options.snow=opts,
               .combine = 'rbind',
               .packages = c("tensr","rTensor","catch","dials","caret") ) %dopar% {
                 objt <- ReAdjten(x=x, z=U, y=y, lambda3 = k)
                 objtsda <- tsda(objt$xres, y)
                 l4 <- objtsda$lambda[objtsda$df>0]
                 n.lam4 <- length(l4)
                 acc.lambda4 <- rep(NA, n.lam4)
                 acc.lambda4<-sapply(1:n.lam4, function(w){
                   l <- l4[w]
                   fit<-estRecatch(x = x, U = U, y=y, lambda1 = i, lambda2 = j, lambda3 = k, lambda4 = l)
                   pred <- pred.ReCatch(fit=fit, x=x, U=U, y=y,
                                        testx=x, testU=U)$pred
                   CMatrix <- suppressWarnings(confusionMatrix(data=as.factor(pred),
                                                               reference = as.factor(y), positive="2"))
                   if (imbalanced) {
                     ifelse( is.na(CMatrix$byClass["F1"]), 0, CMatrix$byClass["F1"])
                   }else{
                     CMatrix$overall["Accuracy"]
                   }
                 })
                 loc.lam4 <- which.max(acc.lambda4)
                 c(l4[loc.lam4], acc.lambda4[loc.lam4])
               }
  #####
  lambda4 <- acc[which.max(acc[,2]),1];names(lambda4) <- "lambda4"
  opt.lambda <- append(as.list(myLambda[which.max(acc[,2]), ]),
                       lambda4, after = 3)
  close(pb)
  stopCluster(workers)
  return(as.vector(opt.lambda))
}

adjten_alpha <- function(x, z, y, testx=NULL, testz=NULL, is.centered=FALSE, lambda3=0){
  n <- length(x)
  dimen <- dim(x[[1]])
  nvars <- prod(dimen)
  vec_x <- matrix(0, nvars, n)
  for (i in 1:n){
    vec_x[,i] <- matrix(x[[i]], nvars, 1)
  }
  vec_x = t(vec_x)

  if (!is.null(testx)){
    tesize <- length(testx)
    vec_testx <- matrix(0, nvars, tesize)
    for (i in 1:tesize){
      vec_testx[,i] <- matrix(testx[[i]], nvars,1)
    }
    vec_testx <- t(vec_testx)
  }
  z <- as.matrix(z) #n*q
  cz <- z
  n <- dim(z)[1]
  q <- dim(z)[2]
  cvecx <- vec_x
  nclass <- as.integer(length(unique(y)))
  if (is.centered==FALSE){
    for (i in 1:nclass){
      if (q>1) {
        cz[y == i,] = sweep(z[y == i,], 2, colMeans(z[y == i,]))
      }else {
        cz[y==i]=z[y==i]-mean(z[y==i])
      }
      cvecx[y == i,] = sweep(vec_x[y == i,], 2, colMeans(vec_x[y == i,]))
    }
  }
  c <- solve(t(cz) %*% cz+lambda3*diag(q)) %*% t(cz)
  c = c %*% cvecx
  A <- t(cz) %*% cz+lambda3*diag(q)

  vec_xres <- vec_x-z%*%c

  xres = array(list(),n)
  for (i in 1:n){
    xres[[i]] = array(vec_xres[i,],dim=dimen)
  }

  if (!is.null(testx)){
    vec_testxres = vec_testx-testz %*% c
    testxres = array(list(),tesize)
    for (i in 1:tesize){
      testxres[[i]] = array(vec_testxres[i,], dim=dimen)
    }
  }else{
    testxres=NULL
  }

  muz = matrix(0, nrow=q, ncol=(nclass-1))    # muz : U_k-U_1  (k>1)
  for (i in 2:nclass){
    if (q>1){muz[,(i-1)] = apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)
    } else {
      muz[i-1] <- mean(z[y==i])-mean(z[y==1])
    }
  }
  gamma <- solve(cov(z)) %*% muz

  return( list(alpha=c, gamma=gamma, xres=xres, A=A, testxres=testxres) )
}

ReAdjten <- function(x, z, y, testx=NULL, testz=NULL, is.centered=FALSE, lambda3=0){
  n <- length(x)
  dimen <- dim(x[[1]])
  nvars <- prod(dimen)
  vec_x <- matrix(0, nvars, n)
  for (i in 1:n){
    vec_x[,i] <- matrix(x[[i]], nvars, 1)
  }
  vec_x = t(vec_x)

  if (!is.null(testx)){
    tesize <- length(testx)
    vec_testx <- matrix(0, nvars, tesize)
    for (i in 1:tesize){
      vec_testx[,i] <- matrix(testx[[i]], nvars,1)
    }
    vec_testx <- t(vec_testx)
  }
  z <- as.matrix(z) #n*q
  cz <- z
  n <- dim(z)[1]
  q <- dim(z)[2]
  cvecx <- vec_x
  nclass <- as.integer(length(unique(y)))
  if (is.centered==FALSE){
    for (i in 1:nclass){
      if (q>1) {
        cz[y == i,] = sweep(z[y == i,], 2, colMeans(z[y == i,]))
      }else {
        cz[y==i]=z[y==i]-mean(z[y==i])
      }
      cvecx[y == i,] = sweep(vec_x[y == i,], 2, colMeans(vec_x[y == i,]))
    }
  }
  c <- solve( (t(cz) %*% cz)+lambda3*diag(q) ) %*% t(cz)
  c = c %*% cvecx
  vec_xres <- vec_x-z%*%c
  A = (t(cz) %*% cz)+lambda3*diag(q)

  xres = array(list(),n)
  for (i in 1:n){
    xres[[i]] = array(vec_xres[i,],dim=dimen)
  }

  if (!is.null(testx)){
    vec_testxres = vec_testx-testz %*% c
    testxres = array(list(),tesize)
    for (i in 1:tesize){
      testxres[[i]] = array(vec_testxres[i,], dim=dimen)
    }
  }else{
    testxres=NULL
  }

  muz = matrix(0, nrow=q, ncol=(nclass-1))    # muz : U_k-U_1  (k>1)
  for (i in 2:nclass){
    if (q>1){muz[,(i-1)] = apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)
    } else {
      muz[i-1] <- mean(z[y==i])-mean(z[y==1])
    }
  }
  #gamma <- solve(cov(z)) %*% muz

  return( list(xres=xres, A=A, testxres=testxres, alpha=c) )
}

estParm <- function(x=x, U=U, y=y, lambda1=0, lambda2=0, lambda3=0){
  library(catch)
  library(tensr)
  library(rTensor)
  pk <- dim(x[[1]])
  nvars <- prod(pk)
  M <- length(pk)
  n <- length(y)
  ng <- as.vector(table(y))
  pig<- ng/n
  G <- length(pig)
  q <- dim(U)[2]
  objt <- ReAdjten(x=x, z=U, y=y, lambda3 = lambda3)
  # centering the data
  MatX<-matrix(0, n, nvars) # n*(p1*...pM)

  Ubar <- array(list(), G) # ng*q
  #for (i in 1:n) MatX[i,] <- matrix(x[[i]], nvars, 1)
  MatX <- t(sapply(x, FUN = as.vector))

  MatXbar <- matrix(NA, G, nvars)
  MatXTilde <- matrix(NA, n, nvars)
  phihat<-array(list(), G)
  Phihat<-array(list(), G)
  #PhihatInv<-array(list(), G)
  Matphi <- matrix(NA, G, q)
  UtildeArray <- array(list(), G)
  Utilde <- c() # n*q
  PhiTemp <- matrix(0, q, q)
  for (g in 1:G){
    # tensor
    MatXbar[g, ]<-colMeans(MatX[y==g,])
    MatXTilde[y==g,] <- sweep( MatX[y==g,], 2, colMeans(MatX[y==g,]) )
    phihat[[g]]<- matrix(apply(U[y==g,], 2, mean), ncol=1)
    Phihat[[g]] <- cov(U[y==g,])
    PhiTemp <- PhiTemp + cov(U[y==g,])*(ng[g]-1)
    Matphi[g,] <- phihat[[g]]
    # covariate
    UtildeArray[[g]] <- matrix(scale(U[y==g,], center = T, scale = F), ncol = q)
    Utilde <- rbind(Utilde, UtildeArray[[g]])
  }
  PhiPool <- PhiTemp/(sum(ng)-1)
  ### Regularization
  Philambda1 <- array(list(), G)
  Philambda12 <- array(list(), G)
  Philambda12Inv <- array(list(), G)
  for (g in 1:G) {
    Philambda1[[g]] <- (1-lambda1)*Phihat[[g]]+lambda1*PhiPool
    Philambda12[[g]] <- (1-lambda2)*Philambda1[[g]]+(lambda2/G)*sum(diag(Philambda1[[g]]))*diag(q)
    Philambda12Inv[[g]] <- solve(Philambda12[[g]])
  }

  Xbar <- fold(mat = MatXbar, row_idx = (M+1), col_idx = seq(1:M),
               modes = c(pk, G) ) # p1*...*pM*G
  XTilde <- fold(mat = MatXTilde, row_idx = (M+1), col_idx = seq(1:M),
                 modes = c(pk, n) )  # p1*...*pM*n
  X <- fold(mat = MatX, row_idx = (M+1), col_idx = seq(1:M),
            modes = c(pk, n) )  # p1*...*pM*n
  A <- objt$A
  P<-solve(t(Utilde)%*%Utilde+lambda3*diag(q))%*%t(Utilde)
  alphahat <- rTensor::ttm(tnsr =  XTilde, mat = P, m = (M+1) ) # p1*...*pM*q
  Xhat<-rTensor::ttm(tnsr = alphahat, mat = U, m = (M+1) ) # p1*...*pM*n
  adjX <- X-Xhat # p1*...*pM*n
  XhatTilde<-rTensor::ttm(tnsr = alphahat, mat = Utilde, m = (M+1) ) # p1*...*pM*n
  adjXTilde <- XTilde-XhatTilde

  muhat <- Xbar-rTensor::ttm(tnsr = alphahat, mat = Matphi, m = (M+1) ) # p1*...*pM*q
  muhatUnfold <- unfold(tnsr = muhat, row_idx = (M+1), col_idx = c(1:M) )@data # G*(p1***pM)
  return(list(adjX=adjX, Utilde=Utilde, alphahat=alphahat, Xhat=Xhat,
              Matphi=Matphi, XTilde=XTilde ,adjXTilde=adjXTilde, Xbar=Xbar, muhat=muhat, muhatUnfold=muhatUnfold,
              RePhi=Philambda12, RePhiInv=Philambda12Inv, phihat=phihat, A=A, pig=pig, x=x, U=U, y=y))
}
#####

#####
Myttv <- function(tnsr, v, m){
  tnsr.dim <- dim(tnsr)
  M <- length(tnsr.dim)
  Mattnsr <- k_unfold(tnsr = tnsr, m = M )@data
  new("Tensor", as.integer(M-1), as.integer(tnsr.dim[-m]),
      data=apply(Mattnsr, 2, FUN = function(x) sum(x*v) ) )
}

tsda <- function(x, y, nlambda = 100, lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2, 1E-03), lambda = NULL, dfmax = nobs, pmax = min(dfmax * 2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+05, sml = 1e-06, verbose = FALSE, perturb = NULL) {
  ## data setup
  this.call <- match.call()
  tmp <- prept(x, y)
  ldim <- length(dim(x[[1]])) #m dimensional
  dimen <- dim(x[[1]])
  maxd <- max(dimen)
  sigma <- as.matrix(tmp$sigma)


  sigmalist<-array(list(),ldim)
  for (i in 1:ldim){
    tmpsigma <- matrix(sigma[i,],maxd,maxd)
    if (!is.null(perturb)){
      diag(tmpsigma) <- diag(tmpsigma)+perturb
      sigma[i,] <- matrix(tmpsigma,nrow=1)
    }
    sigmalist[[i]] <- tmpsigma[1:dimen[i],1:dimen[i]]
  }

  delta <- as.matrix(tmp$delta)
  mu <- as.matrix(tmp$mu)
  prior <- tmp$prior

  nvars <- as.integer(prod(dimen))
  nobs <- length(y)
  nclass <- as.integer(length(unique(y)))
  vnames <- colnames(x)
  if (is.null(vnames))
    vnames <- paste("V", seq(nvars), sep = "")
  ng <- as.integer(dim(delta)[1])
  ## parameter setup
  if (length(pf) != nvars)
    stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit)
  verbose <- as.integer(verbose)
  sml <- as.double(sml)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax*ng)
  pmax <- as.integer(pmax)
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)  #ulam=0 if lambda is missing
  } else {
    # flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))  #lambda is declining
    nlam <- as.integer(length(lambda))
  }
  ## call Fortran core
  fit <<- .Fortran("catch1", obj = double(nlam), ng, nvars, ldim, dimen,maxd, as.double(sigma), as.double(delta),
                   pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1),
                   theta = double(pmax * ng * nlam), itheta = integer(pmax), ntheta = integer(nlam),
                   alam = double(nlam), npass = integer(1), jerr = integer(1))
  ## output
  outlist <- formatoutputt(fit, maxit, pmax, dfmax, nvars, vnames, ng)
  outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr,
                             sigma <- sigmalist, delta = delta, mu = mu, prior = prior, call = this.call))
  if (is.null(lambda))
    outlist$lambda <- lamfix(outlist$lambda)
  class(outlist) <- c("catchobj")
  outlist

}

formatoutputt <- function(fit, maxit, pmax, dfmax, nvars, vnames, ng) {
  nalam <- fit$nalam
  ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta)
  lam <- fit$alam[seq(nalam)]
  obj <- fit$obj[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  resnames <- paste("delta", seq(ng), sep = "")

  errmsg <- errt(fit$jerr, maxit, pmax,dfmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))

  dd <- c(nvars, ng)
  df <- rep(0, nalam)
  if (nthetamax > 0) {
    ja <- fit$itheta[seq(nthetamax)]
    oja <- order(ja)
    ja <- rep(ja[oja], ng)
    itheta <- cumsum(c(1, rep(nthetamax, ng)))
    pos <- rep(1:nalam, each = ng * pmax)
    theta <- split(fit$theta[seq(ng * pmax * nalam)], pos)
    for (l in 1:nalam) {
      theta[[l]] <- matrix(theta[[l]], pmax, ng, byrow = TRUE)[seq(nthetamax),
                                                               , drop = FALSE]
      df[l] <- sum(rowSums(abs(theta[[l]])) != 0)
      theta[[l]] <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames,
                                                               resnames), x = as.vector(theta[[l]][oja, ]), p = as.integer(itheta -
                                                                                                                             1), i = as.integer(ja - 1))
    }
  } else {
    theta <- list()
    for (l in 1:nalam) {
      theta[[l]] <- zeromat(nvars, ng, vnames, resnames)
    }
    df <- rep(0, nalam)
  }
  list(beta = theta, df = df, dim = dd, lambda = lam, obj = obj)
}

prept <- function(x, y) {
  # data setup
  nclass <- as.integer(length(unique(y)))
  prior <- rep(0, nclass)
  for (k in 1:nclass) {
    prior[k] <- mean(y == k)
  }
  ldim <- length(dim(x[[1]])) #m dimensional
  p <- matrix(0,ncol=ldim,nrow=1)
  nvars <- 1
  for (i in 1:ldim){
    p[i] <- as.integer(dim(x[[1]])[i])
    nvars <- nvars*p[i]
  }
  nobs <- length(x)
  nres <- length(y)
  if (nres != nobs)
    stop("x and y have different number of observations")

  dim <- dim(x[[1]])
  num <- matrix(0,nrow=k,ncol=1)
  # estimate mean
  es_M <- array(list(),k)
  for (i in 1:k){
    es_M[[i]] <- array(0,dim=dim)
  }
  ng <- nclass-1
  for (i in 1:nobs){
    es_M[[y[i]]] <- es_M[[y[i]]]+x[[i]]
    num[[y[i]]] <- num[[y[i]]]+1
  }
  for (i in 1:k){
    es_M[[i]] <- es_M[[i]]/num[[i]]
  }

  #estimate sigma
  n <- nobs
  matx <- matrix(list(),n,1)
  es_sigma <- matrix(list(),ldim,1)
  sum_sigma <- matrix(list(),ldim,1)
  matmu <- matrix(list(),k,1)
  tr <- 1
  for (j in 1:ldim){
    for (i in 1:k){
      matmu[[i]] <- mat(es_M[[i]],j)
    }
    sum_sigma[[j]] <- matrix(0,nrow=dim[j],ncol=dim[j])
    for (i in 1:n){
      matx[[i]] <- mat(x[[i]],j)
      sum_sigma[[j]] <- sum_sigma[[j]]+(matx[[i]]-matmu[[y[i]]])%*%t(matx[[i]]-matmu[[y[i]]])
    }
    sum_sigma[[j]] <- sum_sigma[[j]]/(n-k)
    if (j<ldim){
      es_sigma[[j]] <- sum_sigma[[j]]/sum_sigma[[j]][1,1]
      tr <- tr*sum(diag(es_sigma[[j]]))
    }else
    {
      es_sigma[[j]] <- sum_sigma[[j]]/tr
    }
  }
  delta <- matrix(rep(0),nrow=nclass-1,ncol=prod(dim))
  for (i in 1:(nclass-1)){
    delta[i,] <- as.matrix(as.vector(es_M[[i+1]])-as.vector(es_M[[1]]))
  }

  maxd <- max(dim)
  sigma <- matrix(0,nrow=ldim,ncol=(maxd^2))
  for (i in 1:ldim){
    t <- matrix(0,nrow=maxd,ncol=maxd)
    t[1:dim[i],1:dim[i]] <- es_sigma[[i]]
    sigma[i,] <- as.vector(t)
  }
  outlist <- list(sigma=sigma, delta = delta, mu =es_M, prior = prior)
  outlist
}

errt <- function(n, maxit, pmax,dfmax) {
  if (n == 0)
    msg <- ""
  if (n > 0) {
    # fatal error
    if (n < 7777)
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000)
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    # non fatal error
    if (n > -10000)
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=",
                   maxit, " iterations; solutions for larger lambdas returned.\n",
                   sep = "")
    if (n < -10000)
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=",
                   pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n",
                   sep = "")
    if (n < -20000)
      msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=",
                   dfmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n",
                   sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}

zeromat <- function(nvars, nalam, vnames, stepnames) {
  ca <- rep(0, nalam)
  ia <- seq(nalam + 1)
  ja <- rep(1, nalam)
  dd <- c(nvars, nalam)
  new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), x = as.vector(ca),
      p = as.integer(ia - 1), i = as.integer(ja - 1))
}

lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}

#' Generate dataset
#' @description
#' This function is used to generate a dataset with two classes.
#' @param q The number of covariates.
#' @param ng A vector where the g-th component represents the sample size for the g-th class.
#' @param pk Dimension of tensor data.
#' @param phi2 The mean vector of U for the class 2.
#' @param beta2 The tensor coefficients for class 2.
#' @param alpha2 The tensor coefficients.
#' @param CorrType The structure of the correlation matrix. The The default structure, ID, refers to the identity matrix. AR refers to auto-correlation, and CS refers to compound symmetry.
#' @param a It is used to indicate the difference in the covariance of matrices between two classes. When a = 1, it means the variances of covariates are equal. When a > 1, it means one variance is a times larger than the other.
#' @param alpha The coefficient of correlation.
#' @return x tensor data.
#' @return U covariate matrix.
#' @return y Class label.
#' @examples
#' mydt <- GenCovTensor(q=50, ng=c(20,20), pk=c(30,36,30))
#'
#' @export
#'
GenCovTensor <- function(q, ng, pk, phi2=0.3, beta2=0.4, alpha2=1, CorrType="ID", a=1, rho=0.8){
  K <- 2
  library(MASS)
  library(tensr)
  library(parallel)
  library(doSNOW)
  #####
  AR <- function(q, rho) {
    mat <- diag(q)
    return(rho^abs(row(mat)-col(mat)))
  }
  CS <- function(q, rho) {
    mat <- diag(q)
    mat[lower.tri(mat)] <- rho
    mat[upper.tri(mat)] <- rho
    return( mat )
  }
  #####
  sqrtMatrix <- function(M) eigen(M)$vectors %*% sqrt(diag(eigen(M)$values)) %*% t(eigen(M)$vectors)
  if (length(ng)==1) ng <- rep(ng, K)
  nvars <- prod(pk)
  M <- length(pk)
  n <- sum(ng)
  y <- rep(c(1:K), times=ng)
  ###########
  sigma <- array(list(), length(pk) )
  for (k in 1:length(pk)){
    sigma[[k]] <- switch(
      CorrType,
      ID = diag(pk[k]),
      AR = AR(q = pk[k], rho = rho),
      CS = CS(q = pk[k], rho = rho)
    )
  }
  library(MASS)
  # generate covariates
  Phi <- array(list(), K)
  Phi[[1]] <- diag(q)
  Phi[[2]] <- a*diag(q)
  phi <- c(0, phi2)
  U <- matrix(0, nrow=n, ncol=q) # n*q
  for (i in 1:K){
    U[y==i,] <- mvrnorm(n = ng[i], mu = rep(phi[i], q), Sigma = Phi[[i]])# ng*2
  }
  # generate Tensor data
  # Define the coefficients for the mean of Tensor
  CoefForMu <- array(0, dim = pk)
  indi = c(1, 2, 11, 12)
  indj = c(1, 11)
  indk = c(1, 11)
  CoefForMu[indi, indj, indk] <- beta2
  # mu
  muk = array(list(), K)
  muk[[1]] = array(0, dim = pk)  # p1*...*pK
  muk[[2]] = atrans(A = CoefForMu, B = sigma) # Tucker product between an array and a list of matrices
  # alpha star
  alphaStar <- array(0, dim = c(pk, q)) #p1*...*pK*q
  indi<-indj<-indk<-1:15
  alphaStar[indi, indj, indk, 1] <- alpha2
  # tensor
  library(tensr)
  vec_z = matrix(rnorm(nvars*n), ncol=n) # (p1*...*pK)*n
  #z <- array(list(), n)
  sigmaSqrt <- lapply(sigma, sqrtMatrix)
  sigmaSqrt1 <- sigmaSqrt
  sigmaSqrt1[[M+1]] <- diag(q)
  alpha <- atrans(alphaStar, sigmaSqrt1)
  x = array(list(), n)
  library(rTensor)
  library(tensr)
  n_cores <- detectCores() - 2
  workers <- makeCluster(n_cores)
  registerDoSNOW(workers)
  x <- foreach(i=1:n, .packages = c("tensr") ) %dopar% {
    z = array(vec_z[,i], pk)
    muk[[y[i]]] +
      amprod(A=alpha, M=t(U[i,]), k=(M+1) )[,,,1] + # k-mode product dim( amprod(coef, t(z[i,]),4))=30*36*30*1
      atrans(z, sigmaSqrt)
  }
  stopCluster(workers)
  return(list(U=U, x=x, y=y))
}
#####

#' kf.ReCatch
#' @description
#' K-fold cross-validation optimization algorithm for the ReCatch model.
#'
#' @param x Input tensor.
#' @param U Input covariate matrix.
#' @param y Class label.
#' @param size.grid A single integer for the maximum number of parameter value combinations returned.
#' @param n.folds The number of folds. Defaults to 5.
#' @param lambda3Range The parameter interval for \eqn{\lambda_3}.
#' @param setseed An integer value that serves as the random seed.
#' @param imbalanced Logical; if TRUE the F1 score is used as a measure of predictive performance; otherwise, accuracy is used.
#'
#' @return A list of the best set of parameters found.
#' @export
#'
kf.ReCatch <- function(x=x, U=U, y=y, size.grid=100, n.folds=5,
                       lambda3Range = c(0.1, 3.0), setseed=NULL, imbalanced=FALSE) {
  library(foreach)
  library(doSNOW)
  library(caret)
  library(dials)
  library(scales)
  library(splitTools)
  ReAdjten <- function(x, z, y, testx=NULL, testz=NULL, is.centered=FALSE, lambda3=0){
    n <- length(x)
    dimen <- dim(x[[1]])
    nvars <- prod(dimen)
    vec_x <- matrix(0, nvars, n)
    for (i in 1:n){
      vec_x[,i] <- matrix(x[[i]], nvars, 1)
    }
    vec_x = t(vec_x)

    if (!is.null(testx)){
      tesize <- length(testx)
      vec_testx <- matrix(0, nvars, tesize)
      for (i in 1:tesize){
        vec_testx[,i] <- matrix(testx[[i]], nvars,1)
      }
      vec_testx <- t(vec_testx)
    }
    z <- as.matrix(z) #n*q
    cz <- z
    n <- dim(z)[1]
    q <- dim(z)[2]
    cvecx <- vec_x
    nclass <- as.integer(length(unique(y)))
    if (is.centered==FALSE){
      for (i in 1:nclass){
        if (q>1) {
          cz[y == i,] = sweep(z[y == i,], 2, colMeans(z[y == i,]))
        }else {
          cz[y==i]=z[y==i]-mean(z[y==i])
        }
        cvecx[y == i,] = sweep(vec_x[y == i,], 2, colMeans(vec_x[y == i,]))
      }
    }
    c <- solve( (t(cz) %*% cz)+lambda3*diag(q) ) %*% t(cz)
    c = c %*% cvecx
    vec_xres <- vec_x-z%*%c
    A = (t(cz) %*% cz)+lambda3*diag(q)

    xres = array(list(),n)
    for (i in 1:n){
      xres[[i]] = array(vec_xres[i,],dim=dimen)
    }

    if (!is.null(testx)){
      vec_testxres = vec_testx-testz %*% c
      testxres = array(list(),tesize)
      for (i in 1:tesize){
        testxres[[i]] = array(vec_testxres[i,], dim=dimen)
      }
    }else{
      testxres=NULL
    }

    muz = matrix(0, nrow=q, ncol=(nclass-1))    # muz : U_k-U_1  (k>1)
    for (i in 2:nclass){
      if (q>1){muz[,(i-1)] = apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)
      } else {
        muz[i-1] <- mean(z[y==i])-mean(z[y==1])
      }
    }
    #gamma <- solve(cov(z)) %*% muz

    return( list(xres=xres, A=A, testxres=testxres, alpha=c) )
  }

  tsda <- function(x, y, nlambda = 100, lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2, 1E-03), lambda = NULL, dfmax = nobs, pmax = min(dfmax * 2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+05, sml = 1e-06, verbose = FALSE, perturb = NULL) {
    ## data setup
    this.call <- match.call()
    tmp <- prept(x, y)
    ldim <- length(dim(x[[1]])) #m dimensional
    dimen <- dim(x[[1]])
    maxd <- max(dimen)
    sigma <- as.matrix(tmp$sigma)


    sigmalist<-array(list(),ldim)
    for (i in 1:ldim){
      tmpsigma <- matrix(sigma[i,],maxd,maxd)
      if (!is.null(perturb)){
        diag(tmpsigma) <- diag(tmpsigma)+perturb
        sigma[i,] <- matrix(tmpsigma,nrow=1)
      }
      sigmalist[[i]] <- tmpsigma[1:dimen[i],1:dimen[i]]
    }

    delta <- as.matrix(tmp$delta)
    mu <- as.matrix(tmp$mu)
    prior <- tmp$prior

    nvars <- as.integer(prod(dimen))
    nobs <- length(y)
    nclass <- as.integer(length(unique(y)))
    vnames <- colnames(x)
    if (is.null(vnames))
      vnames <- paste("V", seq(nvars), sep = "")
    ng <- as.integer(dim(delta)[1])
    ## parameter setup
    if (length(pf) != nvars)
      stop("The size of penalty factor must be same as the number of input variables")
    maxit <- as.integer(maxit)
    verbose <- as.integer(verbose)
    sml <- as.double(sml)
    pf <- as.double(pf)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax*ng)
    pmax <- as.integer(pmax)
    ## lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
      if (lambda.factor >= 1)
        stop("lambda.factor should be less than 1")
      flmin <- as.double(lambda.factor)
      ulam <- double(1)  #ulam=0 if lambda is missing
    } else {
      # flmin=1 if user define lambda
      flmin <- as.double(1)
      if (any(lambda < 0))
        stop("lambdas should be non-negative")
      ulam <- as.double(rev(sort(lambda)))  #lambda is declining
      nlam <- as.integer(length(lambda))
    }
    ## call Fortran core
    fit <<- .Fortran("catch1", obj = double(nlam), ng, nvars, ldim, dimen,maxd, as.double(sigma), as.double(delta),
                     pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1),
                     theta = double(pmax * ng * nlam), itheta = integer(pmax), ntheta = integer(nlam),
                     alam = double(nlam), npass = integer(1), jerr = integer(1))
    ## output
    outlist <- formatoutputt(fit, maxit, pmax, dfmax, nvars, vnames, ng)
    outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr,
                               sigma <- sigmalist, delta = delta, mu = mu, prior = prior, call = this.call))
    if (is.null(lambda))
      outlist$lambda <- lamfix(outlist$lambda)
    class(outlist) <- c("catchobj")
    outlist

  }

  formatoutputt <- function(fit, maxit, pmax, dfmax, nvars, vnames, ng) {
    nalam <- fit$nalam
    ntheta <- fit$ntheta[seq(nalam)]
    nthetamax <- max(ntheta)
    lam <- fit$alam[seq(nalam)]
    obj <- fit$obj[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    resnames <- paste("delta", seq(ng), sep = "")

    errmsg <- errt(fit$jerr, maxit, pmax,dfmax)  ### error messages from fortran
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))

    dd <- c(nvars, ng)
    df <- rep(0, nalam)
    if (nthetamax > 0) {
      ja <- fit$itheta[seq(nthetamax)]
      oja <- order(ja)
      ja <- rep(ja[oja], ng)
      itheta <- cumsum(c(1, rep(nthetamax, ng)))
      pos <- rep(1:nalam, each = ng * pmax)
      theta <- split(fit$theta[seq(ng * pmax * nalam)], pos)
      for (l in 1:nalam) {
        theta[[l]] <- matrix(theta[[l]], pmax, ng, byrow = TRUE)[seq(nthetamax),
                                                                 , drop = FALSE]
        df[l] <- sum(rowSums(abs(theta[[l]])) != 0)
        theta[[l]] <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames,
                                                                 resnames), x = as.vector(theta[[l]][oja, ]), p = as.integer(itheta -
                                                                                                                               1), i = as.integer(ja - 1))
      }
    } else {
      theta <- list()
      for (l in 1:nalam) {
        theta[[l]] <- zeromat(nvars, ng, vnames, resnames)
      }
      df <- rep(0, nalam)
    }
    list(beta = theta, df = df, dim = dd, lambda = lam, obj = obj)
  }

  # generate sigma, delta and mu from x, y
  prept <- function(x, y) {
    # data setup
    nclass <- as.integer(length(unique(y)))
    prior <- rep(0, nclass)
    for (k in 1:nclass) {
      prior[k] <- mean(y == k)
    }
    ldim <- length(dim(x[[1]])) #m dimensional
    p <- matrix(0,ncol=ldim,nrow=1)
    nvars <- 1
    for (i in 1:ldim){
      p[i] <- as.integer(dim(x[[1]])[i])
      nvars <- nvars*p[i]
    }
    nobs <- length(x)
    nres <- length(y)
    if (nres != nobs)
      stop("x and y have different number of observations")

    dim <- dim(x[[1]])
    num <- matrix(0,nrow=k,ncol=1)
    # estimate mean
    es_M <- array(list(),k)
    for (i in 1:k){
      es_M[[i]] <- array(0,dim=dim)
    }
    ng <- nclass-1
    for (i in 1:nobs){
      es_M[[y[i]]] <- es_M[[y[i]]]+x[[i]]
      num[[y[i]]] <- num[[y[i]]]+1
    }
    for (i in 1:k){
      es_M[[i]] <- es_M[[i]]/num[[i]]
    }

    #estimate sigma
    n <- nobs
    matx <- matrix(list(),n,1)
    es_sigma <- matrix(list(),ldim,1)
    sum_sigma <- matrix(list(),ldim,1)
    matmu <- matrix(list(),k,1)
    tr <- 1
    for (j in 1:ldim){
      for (i in 1:k){
        matmu[[i]] <- mat(es_M[[i]],j)
      }
      sum_sigma[[j]] <- matrix(0,nrow=dim[j],ncol=dim[j])
      for (i in 1:n){
        matx[[i]] <- mat(x[[i]],j)
        sum_sigma[[j]] <- sum_sigma[[j]]+(matx[[i]]-matmu[[y[i]]])%*%t(matx[[i]]-matmu[[y[i]]])
      }
      sum_sigma[[j]] <- sum_sigma[[j]]/(n-k)
      if (j<ldim){
        es_sigma[[j]] <- sum_sigma[[j]]/sum_sigma[[j]][1,1]
        tr <- tr*sum(diag(es_sigma[[j]]))
      }else
      {
        es_sigma[[j]] <- sum_sigma[[j]]/tr
      }
    }
    delta <- matrix(rep(0),nrow=nclass-1,ncol=prod(dim))
    for (i in 1:(nclass-1)){
      delta[i,] <- as.matrix(as.vector(es_M[[i+1]])-as.vector(es_M[[1]]))
    }

    maxd <- max(dim)
    sigma <- matrix(0,nrow=ldim,ncol=(maxd^2))
    for (i in 1:ldim){
      t <- matrix(0,nrow=maxd,ncol=maxd)
      t[1:dim[i],1:dim[i]] <- es_sigma[[i]]
      sigma[i,] <- as.vector(t)
    }
    outlist <- list(sigma=sigma, delta = delta, mu =es_M, prior = prior)
    outlist
  }

  errt <- function(n, maxit, pmax,dfmax) {
    if (n == 0)
      msg <- ""
    if (n > 0) {
      # fatal error
      if (n < 7777)
        msg <- "Memory allocation error; contact package maintainer"
      if (n == 10000)
        msg <- "All penalty factors are <= 0"
      n <- 1
      msg <- paste("in the fortran code -", msg)
    }
    if (n < 0) {
      # non fatal error
      if (n > -10000)
        msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=",
                     maxit, " iterations; solutions for larger lambdas returned.\n",
                     sep = "")
      if (n < -10000)
        msg <- paste("Number of nonzero coefficients along the path exceeds pmax=",
                     pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n",
                     sep = "")
      if (n < -20000)
        msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=",
                     dfmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n",
                     sep = "")
      n <- -1
    }
    list(n = n, msg = msg)
  }

  zeromat <- function(nvars, nalam, vnames, stepnames) {
    ca <- rep(0, nalam)
    ia <- seq(nalam + 1)
    ja <- rep(1, nalam)
    dd <- c(nvars, nalam)
    new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), x = as.vector(ca),
        p = as.integer(ia - 1), i = as.integer(ja - 1))
  }

  lamfix <- function(lam) {
    llam <- log(lam)
    lam[1] <- exp(2 * llam[2] - llam[3])
    lam
  }
  GridPoints<-function(size=100, method="LHC", setseed=NA, lambda3Range=lambda3Range ){
    #Regular: Regular
    #LH: LatinHypercube
    require(dials)
    lambda1 <- new_quant_param(type ="double",range=c(0.1, 0.9), label = c(lambda1="lambda1"), inclusive = c(FALSE, FALSE))
    lambda2 <- new_quant_param(type ="double",range=c(0.1, 0.9), label = c(lambda2="lambda2"), inclusive = c(FALSE, FALSE))
    lambda3 <- new_quant_param(type ="double",range=lambda3Range, label = c(lambda3="lambda3"), inclusive = c(FALSE, FALSE))
    if (!is.na(setseed)) set.seed(setseed)
    if (method=="Regular"){
      Lambda <- grid_regular(x = list(lambda1,lambda2,lambda3), levels=50)
    }else{
      Lambda <- grid_space_filling(list(lambda1,lambda2,lambda3), size=size)
    }
    return(Lambda)
  }

  estRecatch <- function(x=x, U=U, y=y, lambda1=0, lambda2=0, lambda3=0, lambda4=1){
    library(catch)
    library(tensr)
    library(rTensor)
    pk <- dim(x[[1]])
    objt <- ReAdjten(x=x, z=U, y=y, lambda3 = lambda3)
    invisible(capture.output( objtsda <- tryCatch(catch:::tsda(objt$xres, y, lambda = lambda4),
                                                  error=function(e) {NULL} ) ))
    if (is.null(objtsda)) {
      betahat <- matrix(rep(0, prod(pk)), ncol = 1)
    } else{
      betahat <- objtsda$beta$`1`
    }
    nvars <- prod(pk)
    M <- length(pk)
    n <- length(y)
    ng <- as.vector(table(y))
    pig<- ng/n
    G <- length(pig)
    q <- dim(U)[2]
    # centering the data
    MatX<-matrix(0, n, nvars) # n*(p1*...pM)

    Ubar <- array(list(), G) # ng*q
    MatX <- t(sapply(x, FUN = as.vector))

    MatXbar <- matrix(NA, G, nvars)
    MatXTilde <- matrix(NA, n, nvars)
    phihat<-array(list(), G)
    Phihat<-array(list(), G)
    Matphi <- matrix(NA, G, q)
    UtildeArray <- array(list(), G)
    Utilde <- c() # n*q
    PhiTemp <- matrix(0, q, q)
    for (g in 1:G){
      # tensor
      MatXbar[g, ]<-colMeans(MatX[y==g,])
      MatXTilde[y==g,] <- sweep( MatX[y==g,], 2, colMeans(MatX[y==g,]) )
      phihat[[g]]<- matrix(apply(U[y==g,], 2, mean), ncol=1)
      Phihat[[g]] <- cov(U[y==g,])
      PhiTemp <- PhiTemp + cov(U[y==g,])*(ng[g]-1)
      Matphi[g,] <- phihat[[g]]
      UtildeArray[[g]] <- matrix(scale(U[y==g,], center = T, scale = F), ncol = q)
      Utilde <- rbind(Utilde, UtildeArray[[g]])
    }
    PhiPool <- PhiTemp/(sum(ng)-1)
    Philambda1 <- array(list(), G)
    Philambda12 <- array(list(), G)
    Philambda12Inv <- array(list(), G)
    for (g in 1:G) {
      Philambda1[[g]] <- (1-lambda1)*Phihat[[g]]+lambda1*PhiPool
      Philambda12[[g]] <- (1-lambda2)*Philambda1[[g]]+(lambda2/G)*sum(diag(Philambda1[[g]]))*diag(q)
      Philambda12Inv[[g]] <- solve(Philambda12[[g]])
    }

    Xbar <- fold(mat = MatXbar, row_idx = (M+1), col_idx = seq(1:M),
                 modes = c(pk, G) ) # p1*...*pM*G
    XTilde <- fold(mat = MatXTilde, row_idx = (M+1), col_idx = seq(1:M),
                   modes = c(pk, n) )  # p1*...*pM*n
    X <- fold(mat = MatX, row_idx = (M+1), col_idx = seq(1:M),
              modes = c(pk, n) )  # p1*...*pM*n
    A <- objt$A
    P<-solve(t(Utilde)%*%Utilde+lambda3*diag(q))%*%t(Utilde)
    alphahat <- rTensor::ttm(tnsr =  XTilde, mat = P, m = (M+1) ) # p1*...*pM*q
    Xhat<-rTensor::ttm(tnsr = alphahat, mat = U, m = (M+1) ) # p1*...*pM*n
    adjX <- X-Xhat # p1*...*pM*n
    XhatTilde<-rTensor::ttm(tnsr = alphahat, mat = Utilde, m = (M+1) ) # p1*...*pM*n
    adjXTilde <- XTilde-XhatTilde

    muhat <- Xbar-rTensor::ttm(tnsr = alphahat, mat = Matphi, m = (M+1) ) # p1*...*pM*q
    muhatUnfold <- unfold(tnsr = muhat, row_idx = (M+1), col_idx = c(1:M) )@data # G*(p1***pM)
    return(list(adjX=adjX, Utilde=Utilde, alphahat=alphahat, betahat=betahat, Xhat=Xhat,
                Matphi=Matphi, XTilde=XTilde ,adjXTilde=adjXTilde, Xbar=Xbar, muhat=muhat, muhatUnfold=muhatUnfold,
                RePhi=Philambda12, RePhiInv=Philambda12Inv, phihat=phihat, A=A, pig=pig, x=x, U=U, y=y))
  }

  ReCatch.pred <- function(fit=fit, x=x, U=U, y=y, testx=NULL, testU=NULL) {
    pk <- dim(x[[1]])
    nvars <- prod(pk)
    M <- length(pk)
    n <- length(y)
    ng <- as.vector(table(y))
    pig<- ng/n
    G <- length(pig)
    q <- dim(U)[2]
    # centering the data
    if (!is.null(testx)) {
      if (length(testx)==1) {
        testxn <- 1
        testU <- matrix(testU, nrow = 1, ncol = q)
        testMatX<-matrix(0, testxn, nvars)
        for (i in 1:testxn) testMatX[i,] <- matrix(testx[[i]], nvars, 1)
        testX <- fold(mat = testMatX, row_idx = (M+1), col_idx = seq(1:M),
                      modes = c(pk, testxn) )  # p1*...*pM*n
        testXhat<-Myttv(tnsr = fit$alphahat, v = testU, m = (M+1) ) # p1*...*pM*n
        testX <- fold(mat = matrix(testX@data, nrow = pk[1]), row_idx = 1, col_idx = c(2:M), modes=pk )
        testadjX <- testX-testXhat # p1*...*pM*n
      } else {
        testxn <- length(testx)
        testMatX<-matrix(0, testxn, nvars)
        for (i in 1:testxn) testMatX[i,] <- matrix(testx[[i]], nvars, 1)
        testX <- fold(mat = testMatX, row_idx = (M+1), col_idx = seq(1:M),
                      modes = c(pk, testxn) )  # p1*...*pM*n
        testXhat<-rTensor::ttm(tnsr = fit$alphahat, mat = testU, m = (M+1) ) # p1*...*pM*n
        testadjX <- testX-testXhat # p1*...*pM*n
      }
    }
    disScore <- matrix(0, nrow = ifelse(is.null(testU), n, testxn), ncol = G)
    for (g in 2:G ){
      beta <- as.vector(fit$betahat[,g-1])
      gammag <- fit$RePhiInv[[g]]%*%fit$phihat[[g]]- fit$RePhiInv[[1]]%*%fit$phihat[[1]]
      Dg <- fit$RePhiInv[[g]] - fit$RePhiInv[[1]]
      if ( is.null(testU) ) {
        Score1 <- -0.5*diag(U%*%Dg%*%t(U))
        Score2 <- U%*%gammag
        bg <-apply(k_unfold(tnsr = fit$adjX, m = (M+1) )@data, 1, function(x) sum(x*fit$beta))

      } else {
        Score1 <- -0.5*diag(testU%*%Dg%*%t(testU))
        Score2 <- testU%*%gammag
        if (testxn==1) {
          bg <- sum(beta*c(testadjX@data))
        }else{
          bg <-apply(k_unfold(tnsr = testadjX, m = (M+1) )@data, 1, function(x) sum(x*fit$beta))
        }
      }
      mugPlus1 <- fit$muhatUnfold[g,] + fit$muhatUnfold[1,]
      ag.Tensor <- sum( fit$beta*mugPlus1 )
      ag.Cov<-as.numeric( log(det(fit$RePhi[[g]])/det(fit$RePhi[[1]]) )
                          + t(gammag)%*%(fit$phihat[[g]]+fit$phihat[[1]])
                          - t(fit$phihat[[g]])%*%Dg%*%fit$phihat[[1]]
                          -2*log(fit$pig[g]/fit$pig[1]) )
      ag <- ag.Cov + ag.Tensor
      disScore[,g]<- Score1+Score2-0.5*ag+bg
    }
    yhat <- apply(disScore, 1, which.max)
    return(list(pred=yhat))
  }
  if (!(is.null(setseed))) set.seed(setseed)
  myLambda <- GridPoints(size = size.grid, lambda3Range = lambda3Range, method="LHC")
  #apply(as.matrix(myLambda), MARGIN = 1, Leave1out, x = x, U = U, y = y)
  ## generate a sequence of lambda 4 based on all data
  library(iterators)
  library(foreach)
  n_cores <- detectCores() - 2
  workers <- makeCluster(n_cores)
  registerDoSNOW(workers)
  n.lam3 <- length(myLambda$lambda3)
  myLambda4 <- foreach(w = 1:n.lam3, .combine = c,.packages = c("tensr","rTensor","catch","dials","caret")) %dopar% {
    objt <- ReAdjten(x= x, z = U, y= y, lambda3 = myLambda$lambda3[w])
    invisible(capture.output( objtsda <- catch:::tsda(objt$xres, y=y, nlambda = 100)))
    objtsda$lambda
  }
  stopCluster(workers)

  folds <- splitTools::create_folds(y=y, k=n.folds, type="stratified")
  cat("The progress of K-fold cross-validation\n")
  cat("Imbalanced Data: ", imbalanced, "\n")
  cat("Random seed:", ifelse(is.null(setseed), FALSE, setseed), "\n")
  n_cores <- detectCores() - 2
  workers <- makeCluster(n_cores)
  registerDoSNOW(workers)
  pb <- txtProgressBar(0, size.grid, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  acc<-foreach(i=myLambda$lambda1, j=myLambda$lambda2, k=myLambda$lambda3,
               .combine = cbind, .options.snow=opts,
               .packages = c("tensr","rTensor","catch","dials","caret")) %dopar% {
                 loc.lam3 <- which(myLambda$lambda3==k)
                 nlambda4 <- length(myLambda4[[loc.lam3]])
                 lambda4 <- myLambda4[[loc.lam3]]
                 acc1 <- matrix(NA, nrow = n.folds, ncol = nlambda4)
                 for (m in 1:n.folds){
                   trainID <- folds[[m]]
                   trainY <- y[trainID];trainX <- x[trainID];trainU <- U[trainID,]
                   testY <- y[-trainID];testX <- x[-trainID];testU <- U[-trainID,]
                   acc1[m,] <- sapply(1:nlambda4, function(u) {
                     fit <- estRecatch(x = trainX, U = trainU, y=trainY,
                                       lambda1 = i, lambda2 = j, lambda3 = k, lambda4 = lambda4[u] )
                     pred <- ReCatch.pred(fit=fit, x=trainX, U=trainU, y=trainY,
                                          testx=testX, testU=testU)$pred
                     CMatrix <- suppressWarnings(confusionMatrix(data=as.factor(pred),
                                                                 reference = as.factor(testY), positive="2"))
                     ifelse(imbalanced, CMatrix$byClass["F1"], CMatrix$overall["Accuracy"])
                   })
                 }
                 accl4 <- colMeans(acc1)
                 loc.opt.lam4 <- which(accl4==max(accl4, na.rm = T))
                 c(accl4[min(loc.opt.lam4)],i,j,k,lambda4[min(loc.opt.lam4)])
               }
  opt.lambda<-acc[2:5, which.max(acc[1,])]
  close(pb)
  stopCluster(workers)
  list(lambda1=opt.lambda[1],lambda2=opt.lambda[2],
       lambda3=opt.lambda[3],lambda4=opt.lambda[4])
}

#' loo.ReCatch
#' @description
#' Leave-one-out cross-validation optimization algorithm for the ReCatch model.
#'
#' @param x Input tensor.
#' @param U Input covariate matrix.
#' @param y Class label.
#' @param size.grid A single integer for the maximum number of parameter value combinations returned.
#' @param setseed An integer value that serves as the random seed.
#' @param lambda3Range The parameter interval for \eqn{\lambda_3}.
#' @param imbalanced Logical; if TRUE the F1 score is used as a measure of predictive performance; otherwise, accuracy is used.
#'
#' @return A list of the best set of parameters found.
#' @export
#'
loo.ReCatch0 <-function(x=x, U=U, y=y, size.grid=100,
                       setseed=NULL, lambda3Range=c(0.1, 3.0), imbalanced=FALSE){
  ReAdjten <- function(x, z, y, testx=NULL, testz=NULL, is.centered=FALSE, lambda3=0){
    n <- length(x)
    dimen <- dim(x[[1]])
    nvars <- prod(dimen)
    vec_x <- matrix(0, nvars, n)
    for (i in 1:n){
      vec_x[,i] <- matrix(x[[i]], nvars, 1)
    }
    vec_x = t(vec_x)

    if (!is.null(testx)){
      tesize <- length(testx)
      vec_testx <- matrix(0, nvars, tesize)
      for (i in 1:tesize){
        vec_testx[,i] <- matrix(testx[[i]], nvars,1)
      }
      vec_testx <- t(vec_testx)
    }
    z <- as.matrix(z) #n*q
    cz <- z
    n <- dim(z)[1]
    q <- dim(z)[2]
    cvecx <- vec_x
    nclass <- as.integer(length(unique(y)))
    if (is.centered==FALSE){
      for (i in 1:nclass){
        if (q>1) {
          cz[y == i,] = sweep(z[y == i,], 2, colMeans(z[y == i,]))
        }else {
          cz[y==i]=z[y==i]-mean(z[y==i])
        }
        cvecx[y == i,] = sweep(vec_x[y == i,], 2, colMeans(vec_x[y == i,]))
      }
    }
    c <- solve( (t(cz) %*% cz)+lambda3*diag(q) ) %*% t(cz)
    c = c %*% cvecx
    vec_xres <- vec_x-z%*%c
    A = (t(cz) %*% cz)+lambda3*diag(q)

    xres = array(list(),n)
    for (i in 1:n){
      xres[[i]] = array(vec_xres[i,],dim=dimen)
    }

    if (!is.null(testx)){
      vec_testxres = vec_testx-testz %*% c
      testxres = array(list(),tesize)
      for (i in 1:tesize){
        testxres[[i]] = array(vec_testxres[i,], dim=dimen)
      }
    }else{
      testxres=NULL
    }

    muz = matrix(0, nrow=q, ncol=(nclass-1))    # muz : U_k-U_1  (k>1)
    for (i in 2:nclass){
      if (q>1){muz[,(i-1)] = apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)
      } else {
        muz[i-1] <- mean(z[y==i])-mean(z[y==1])
      }
    }
    #gamma <- solve(cov(z)) %*% muz

    return( list(xres=xres, A=A, testxres=testxres, alpha=c) )
  }

  tsda <- function(x, y, nlambda = 100, lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2, 1E-03), lambda = NULL, dfmax = nobs, pmax = min(dfmax * 2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+05, sml = 1e-06, verbose = FALSE, perturb = NULL) {
    ## data setup
    this.call <- match.call()
    tmp <- prept(x, y)
    ldim <- length(dim(x[[1]])) #m dimensional
    dimen <- dim(x[[1]])
    maxd <- max(dimen)
    sigma <- as.matrix(tmp$sigma)


    sigmalist<-array(list(),ldim)
    for (i in 1:ldim){
      tmpsigma <- matrix(sigma[i,],maxd,maxd)
      if (!is.null(perturb)){
        diag(tmpsigma) <- diag(tmpsigma)+perturb
        sigma[i,] <- matrix(tmpsigma,nrow=1)
      }
      sigmalist[[i]] <- tmpsigma[1:dimen[i],1:dimen[i]]
    }

    delta <- as.matrix(tmp$delta)
    mu <- as.matrix(tmp$mu)
    prior <- tmp$prior

    nvars <- as.integer(prod(dimen))
    nobs <- length(y)
    nclass <- as.integer(length(unique(y)))
    vnames <- colnames(x)
    if (is.null(vnames))
      vnames <- paste("V", seq(nvars), sep = "")
    ng <- as.integer(dim(delta)[1])
    ## parameter setup
    if (length(pf) != nvars)
      stop("The size of penalty factor must be same as the number of input variables")
    maxit <- as.integer(maxit)
    verbose <- as.integer(verbose)
    sml <- as.double(sml)
    pf <- as.double(pf)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax*ng)
    pmax <- as.integer(pmax)
    ## lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
      if (lambda.factor >= 1)
        stop("lambda.factor should be less than 1")
      flmin <- as.double(lambda.factor)
      ulam <- double(1)  #ulam=0 if lambda is missing
    } else {
      # flmin=1 if user define lambda
      flmin <- as.double(1)
      if (any(lambda < 0))
        stop("lambdas should be non-negative")
      ulam <- as.double(rev(sort(lambda)))  #lambda is declining
      nlam <- as.integer(length(lambda))
    }
    ## call Fortran core
    fit <<- .Fortran("catch1", obj = double(nlam), ng, nvars, ldim, dimen,maxd, as.double(sigma), as.double(delta),
                     pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1),
                     theta = double(pmax * ng * nlam), itheta = integer(pmax), ntheta = integer(nlam),
                     alam = double(nlam), npass = integer(1), jerr = integer(1))
    ## output
    outlist <- formatoutputt(fit, maxit, pmax, dfmax, nvars, vnames, ng)
    outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr,
                               sigma <- sigmalist, delta = delta, mu = mu, prior = prior, call = this.call))
    if (is.null(lambda))
      outlist$lambda <- lamfix(outlist$lambda)
    class(outlist) <- c("catchobj")
    outlist

  }

  formatoutputt <- function(fit, maxit, pmax, dfmax, nvars, vnames, ng) {
    nalam <- fit$nalam
    ntheta <- fit$ntheta[seq(nalam)]
    nthetamax <- max(ntheta)
    lam <- fit$alam[seq(nalam)]
    obj <- fit$obj[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    resnames <- paste("delta", seq(ng), sep = "")

    errmsg <- errt(fit$jerr, maxit, pmax,dfmax)  ### error messages from fortran
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))

    dd <- c(nvars, ng)
    df <- rep(0, nalam)
    if (nthetamax > 0) {
      ja <- fit$itheta[seq(nthetamax)]
      oja <- order(ja)
      ja <- rep(ja[oja], ng)
      itheta <- cumsum(c(1, rep(nthetamax, ng)))
      pos <- rep(1:nalam, each = ng * pmax)
      theta <- split(fit$theta[seq(ng * pmax * nalam)], pos)
      for (l in 1:nalam) {
        theta[[l]] <- matrix(theta[[l]], pmax, ng, byrow = TRUE)[seq(nthetamax),
                                                                 , drop = FALSE]
        df[l] <- sum(rowSums(abs(theta[[l]])) != 0)
        theta[[l]] <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames,
                                                                 resnames), x = as.vector(theta[[l]][oja, ]), p = as.integer(itheta -
                                                                                                                               1), i = as.integer(ja - 1))
      }
    } else {
      theta <- list()
      for (l in 1:nalam) {
        theta[[l]] <- zeromat(nvars, ng, vnames, resnames)
      }
      df <- rep(0, nalam)
    }
    list(beta = theta, df = df, dim = dd, lambda = lam, obj = obj)
  }

  # generate sigma, delta and mu from x, y
  prept <- function(x, y) {
    # data setup
    nclass <- as.integer(length(unique(y)))
    prior <- rep(0, nclass)
    for (k in 1:nclass) {
      prior[k] <- mean(y == k)
    }
    ldim <- length(dim(x[[1]])) #m dimensional
    p <- matrix(0,ncol=ldim,nrow=1)
    nvars <- 1
    for (i in 1:ldim){
      p[i] <- as.integer(dim(x[[1]])[i])
      nvars <- nvars*p[i]
    }
    nobs <- length(x)
    nres <- length(y)
    if (nres != nobs)
      stop("x and y have different number of observations")

    dim <- dim(x[[1]])
    num <- matrix(0,nrow=k,ncol=1)
    # estimate mean
    es_M <- array(list(),k)
    for (i in 1:k){
      es_M[[i]] <- array(0,dim=dim)
    }
    ng <- nclass-1
    for (i in 1:nobs){
      es_M[[y[i]]] <- es_M[[y[i]]]+x[[i]]
      num[[y[i]]] <- num[[y[i]]]+1
    }
    for (i in 1:k){
      es_M[[i]] <- es_M[[i]]/num[[i]]
    }

    #estimate sigma
    n <- nobs
    matx <- matrix(list(),n,1)
    es_sigma <- matrix(list(),ldim,1)
    sum_sigma <- matrix(list(),ldim,1)
    matmu <- matrix(list(),k,1)
    tr <- 1
    for (j in 1:ldim){
      for (i in 1:k){
        matmu[[i]] <- mat(es_M[[i]],j)
      }
      sum_sigma[[j]] <- matrix(0,nrow=dim[j],ncol=dim[j])
      for (i in 1:n){
        matx[[i]] <- mat(x[[i]],j)
        sum_sigma[[j]] <- sum_sigma[[j]]+(matx[[i]]-matmu[[y[i]]])%*%t(matx[[i]]-matmu[[y[i]]])
      }
      sum_sigma[[j]] <- sum_sigma[[j]]/(n-k)
      if (j<ldim){
        es_sigma[[j]] <- sum_sigma[[j]]/sum_sigma[[j]][1,1]
        tr <- tr*sum(diag(es_sigma[[j]]))
      }else
      {
        es_sigma[[j]] <- sum_sigma[[j]]/tr
      }
    }
    delta <- matrix(rep(0),nrow=nclass-1,ncol=prod(dim))
    for (i in 1:(nclass-1)){
      delta[i,] <- as.matrix(as.vector(es_M[[i+1]])-as.vector(es_M[[1]]))
    }

    maxd <- max(dim)
    sigma <- matrix(0,nrow=ldim,ncol=(maxd^2))
    for (i in 1:ldim){
      t <- matrix(0,nrow=maxd,ncol=maxd)
      t[1:dim[i],1:dim[i]] <- es_sigma[[i]]
      sigma[i,] <- as.vector(t)
    }
    outlist <- list(sigma=sigma, delta = delta, mu =es_M, prior = prior)
    outlist
  }

  errt <- function(n, maxit, pmax,dfmax) {
    if (n == 0)
      msg <- ""
    if (n > 0) {
      # fatal error
      if (n < 7777)
        msg <- "Memory allocation error; contact package maintainer"
      if (n == 10000)
        msg <- "All penalty factors are <= 0"
      n <- 1
      msg <- paste("in the fortran code -", msg)
    }
    if (n < 0) {
      # non fatal error
      if (n > -10000)
        msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=",
                     maxit, " iterations; solutions for larger lambdas returned.\n",
                     sep = "")
      if (n < -10000)
        msg <- paste("Number of nonzero coefficients along the path exceeds pmax=",
                     pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n",
                     sep = "")
      if (n < -20000)
        msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=",
                     dfmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n",
                     sep = "")
      n <- -1
    }
    list(n = n, msg = msg)
  }

  zeromat <- function(nvars, nalam, vnames, stepnames) {
    ca <- rep(0, nalam)
    ia <- seq(nalam + 1)
    ja <- rep(1, nalam)
    dd <- c(nvars, nalam)
    new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), x = as.vector(ca),
        p = as.integer(ia - 1), i = as.integer(ja - 1))
  }

  lamfix <- function(lam) {
    llam <- log(lam)
    lam[1] <- exp(2 * llam[2] - llam[3])
    lam
  }

  Myttv <- function(tnsr, v, m){
    tnsr.dim <- dim(tnsr)
    M <- length(tnsr.dim)
    Mattnsr <- k_unfold(tnsr = tnsr, m = M )@data
    new("Tensor", as.integer(M-1), as.integer(tnsr.dim[-m]),
        data=apply(Mattnsr, 2, FUN = function(x) sum(x*v) ) )
  }

  estRecatch <- function(x=x, U=U, y=y, lambda1=0, lambda2=0, lambda3=0, lambda4=1){
    library(catch)
    library(tensr)
    library(rTensor)
    #lambda1 <- 0
    #lambda2 <- 0
    #lambda3 <- 0
    #lambda4 <- 1
    #obj <- catch(x=x, z=U, y=y, lambda = lambda4, testx = x, testz = U)
    #objt <- ReAdjten(x=x, z=U, y=y, testx =testx, testz = testU, lambda3 = lambda3)
    pk <- dim(x[[1]])
    objt <- ReAdjten(x=x, z=U, y=y, lambda3 = lambda3)
    invisible(capture.output( objtsda <- tryCatch(catch:::tsda(objt$xres, y, lambda = lambda4),
                                                  error=function(e) {NULL} ) ))
    if (is.null(objtsda)) {
      betahat <- matrix(rep(0, prod(pk)), ncol = 1)
    } else{
      betahat <- objtsda$beta$`1`
    }
    # X: tensor, array, c(pk, n)
    # U: Covariates, matrix, n*q
    # y: class labels, vector, n*1
    pk <- dim(x[[1]])
    nvars <- prod(pk)
    M <- length(pk)
    n <- length(y)
    ng <- as.vector(table(y))
    pig<- ng/n
    G <- length(pig)
    q <- dim(U)[2]
    # centering the data
    MatX<-matrix(0, n, nvars) # n*(p1*...pM)

    Ubar <- array(list(), G) # ng*q
    MatX <- t(sapply(x, FUN = as.vector))

    MatXbar <- matrix(NA, G, nvars)
    MatXTilde <- matrix(NA, n, nvars)
    phihat<-array(list(), G)
    Phihat<-array(list(), G)
    Matphi <- matrix(NA, G, q)
    UtildeArray <- array(list(), G)
    Utilde <- c() # n*q
    PhiTemp <- matrix(0, q, q)
    for (g in 1:G){
      # tensor
      MatXbar[g, ]<-colMeans(MatX[y==g,])
      MatXTilde[y==g,] <- sweep( MatX[y==g,], 2, colMeans(MatX[y==g,]) )
      phihat[[g]]<- matrix(apply(U[y==g,], 2, mean), ncol=1)
      Phihat[[g]] <- cov(U[y==g,])
      PhiTemp <- PhiTemp + cov(U[y==g,])*(ng[g]-1)
      Matphi[g,] <- phihat[[g]]
      UtildeArray[[g]] <- matrix(scale(U[y==g,], center = T, scale = F), ncol = q)
      Utilde <- rbind(Utilde, UtildeArray[[g]])
    }
    PhiPool <- PhiTemp/(sum(ng)-1)
    Philambda1 <- array(list(), G)
    Philambda12 <- array(list(), G)
    Philambda12Inv <- array(list(), G)
    for (g in 1:G) {
      Philambda1[[g]] <- (1-lambda1)*Phihat[[g]]+lambda1*PhiPool
      Philambda12[[g]] <- (1-lambda2)*Philambda1[[g]]+(lambda2/G)*sum(diag(Philambda1[[g]]))*diag(q)
      Philambda12Inv[[g]] <- solve(Philambda12[[g]])
    }

    Xbar <- fold(mat = MatXbar, row_idx = (M+1), col_idx = seq(1:M),
                 modes = c(pk, G) ) # p1*...*pM*G
    XTilde <- fold(mat = MatXTilde, row_idx = (M+1), col_idx = seq(1:M),
                   modes = c(pk, n) )  # p1*...*pM*n
    X <- fold(mat = MatX, row_idx = (M+1), col_idx = seq(1:M),
              modes = c(pk, n) )  # p1*...*pM*n
    A <- objt$A
    P<-solve(t(Utilde)%*%Utilde+lambda3*diag(q))%*%t(Utilde)
    alphahat <- rTensor::ttm(tnsr =  XTilde, mat = P, m = (M+1) ) # p1*...*pM*q
    Xhat<-rTensor::ttm(tnsr = alphahat, mat = U, m = (M+1) ) # p1*...*pM*n
    adjX <- X-Xhat # p1*...*pM*n
    XhatTilde<-rTensor::ttm(tnsr = alphahat, mat = Utilde, m = (M+1) ) # p1*...*pM*n
    adjXTilde <- XTilde-XhatTilde

    muhat <- Xbar-rTensor::ttm(tnsr = alphahat, mat = Matphi, m = (M+1) ) # p1*...*pM*q
    muhatUnfold <- unfold(tnsr = muhat, row_idx = (M+1), col_idx = c(1:M) )@data # G*(p1***pM)
    return(list(adjX=adjX, Utilde=Utilde, alphahat=alphahat, betahat=betahat, Xhat=Xhat,
                Matphi=Matphi, XTilde=XTilde ,adjXTilde=adjXTilde, Xbar=Xbar, muhat=muhat, muhatUnfold=muhatUnfold,
                RePhi=Philambda12, RePhiInv=Philambda12Inv, phihat=phihat, A=A, pig=pig, x=x, U=U, y=y))
  }

  GridPoints<-function(size=size.grid, method="LHC", setseed=NA, lambda3Range=lambda3Range ){
    #Regular: Regular
    #LH: LatinHypercube
    require(dials)
    lambda1 <- new_quant_param(type ="double", range=c(0.1, 0.9), label = c(lambda1="lambda1"), inclusive = c(FALSE, FALSE))
    lambda2 <- new_quant_param(type ="double", range=c(0.1, 0.9), label = c(lambda2="lambda2"), inclusive = c(FALSE, FALSE))
    lambda3 <- new_quant_param(type ="double", range=lambda3Range, label = c(lambda3="lambda3"), inclusive = c(FALSE, FALSE))
    lambda4 <- new_quant_param(type ="double", range=c(0.1, 4.0), label = c(lambda4="lambda4"), inclusive = c(FALSE, FALSE))
    if (!is.na(setseed)) set.seed(setseed)
    if (method=="Regular"){
      Lambda <- grid_regular(x = list(lambda1,lambda2,lambda3), levels=50)
    }else{
      Lambda <- grid_space_filling(list(lambda1,lambda2,lambda3, lambda4), size=size)
    }
    return(Lambda)
  }
  if (!(is.null(setseed))) set.seed(setseed)
  myLambda <- GridPoints(size = size.grid, lambda3Range = lambda3Range, method="LHC")
  n <- length(y)
  library(caret)
  library(doSNOW)
  library(foreach)
  library(parallel)
  library(iterators)
  cat("The progress of Leave-One-Out cross-validation\n")
  cat("Imbalanced Data: ", imbalanced,"\n")
  cat("Random seed:", ifelse(is.null(setseed), FALSE, setseed), "\n")
  n_cores <- detectCores() - 2
  workers <- makeCluster(n_cores)
  registerDoSNOW(workers)
  pb <- txtProgressBar(0, size.grid*n, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  opts1 <- list(chunkSize=5)
  acc<-foreach(i=myLambda$lambda1, j=myLambda$lambda2, k=myLambda$lambda3, l=myLambda$lambda4,
               .combine = cbind, .options.snow=opts, .options.nws=opts1,
               .packages = c("tensr","rTensor","catch","dials","caret")) %:%
    foreach(ni=1:n, .combine = c, .packages = c("tensr","rTensor","catch","dials","caret")) %dopar% {
      #disScore <- matrix(0, nrow = n, ncol = G)
      ###
      fit <- estRecatch(x=x, U=U, y=y, lambda1=i, lambda2=j, lambda3=k, lambda4=l)
      Xbar <- fit$Xbar
      Utilde <- fit$Utilde
      phihat <- fit$phihat
      Phihat <- fit$RePhi
      alphahat <- fit$alphahat
      PhihatInv <- fit$RePhiInv
      XTilde <- fit$XTilde
      Matphi <- fit$Matphi
      adjXTilde <- fit$adjXTilde
      A <- fit$A
      ng <- as.numeric(table(y))
      n <- sum(ng)
      G <- length(table(y))
      q <- dim(U)[2]
      pk <- dim(x[[1]])
      M <- length(dim(x[[1]]))
      #adjX_Mp1 <-  unfold(tnsr = fit$adjX, row_idx = (M+1), col_idx = c(1:M))@data
      adjX_Mp1 <-  unfold(tnsr = adjXTilde, row_idx = (M+1), col_idx = c(1:M))@data
      Xbar_Mp1 <- unfold(tnsr = Xbar, row_idx = (M+1), col_idx = c(1:M))@data
      AInv <- solve(A)
      phi_gi <- array(list(), G) # ng*q
      detPhi <- unlist( lapply(Phihat, function(x) det(x)))
      #for (g in 1:G) detPhi_gi[g] <- det(fit$RePhi[[g]])
      ##
      Utilde1 <- Utilde/(rep(ng, times=ng)-1)
      UUT <- t(apply(Utilde, 1, function(x) x%o%x ))
      Hi <- diag(Utilde%*%AInv%*%t(Utilde))
      ng.all <- rep(ng, times=ng)
      Utildestar2 <- Utilde%*%AInv # A^(-1)%*%Utilde_i

      hi <- c()
      Utildestar <- c()
      for (g in 1:G){
        hi <- c(hi, diag(Utilde[y==g,] %*% PhihatInv[[g]] %*% t(Utilde[y==g,])) )
        Utildestar <- rbind(Utildestar, Utilde[y==g,]%*%PhihatInv[[g]])
      }
      coef.detPhi <- ( (ng.all-2)/(ng.all-1) )^q*(1-(hi*ng.all)/(ng.all-1)^2)
      phi_gi <- phihat
      Phi_gi <- Phihat
      PhiInv_gi <- PhihatInv
      detPhi_gi <- detPhi
      ###
      pig_gi <- as.numeric(table(y[-ni])/n)
      pig <- as.vector(table(y[-ni]))/sum(table(y[-ni]))
      objt <- ReAdjten(x=x[-ni], z=U[-ni,], y=y[-ni], lambda3 = k)
      rangelam4 <- range(tsda(objt$xres, y[-ni], nlambda = 100)$lambda)
      if (l<rangelam4[1] | l>rangelam4[2]) {
        disScore <- NA
      }else{
        #####
        betahat <- tsda(objt$xres, y[-ni], lambda = l)$beta$'1'
        ind_g <- y[ni];ng_gi<-ng.all[ni]
        phi_gi[[ind_g]] <- phihat[[ ind_g ]] - (1/(ng_gi-1))*( Utilde[ni,] )
        #Phi_gi[[ind_g]] <- ( ng_gi/(ng_gi-1) )*( Phihat[[ ind_g ]] - ( 1/(ng_gi- 1))*matrix(UUT[i,],q,q) )
        Phi_gi[[ind_g]] <- ( (ng_gi-1)/(ng_gi-2) )*( Phihat[[ ind_g ]] - ( ng_gi/(ng_gi- 1)^2)*matrix(UUT[ni,],q,q) )
        PhiInv_gi[[ind_g]] <- ( (ng_gi-2)/(ng_gi-1) )*( PhihatInv[[ind_g]]  + Utildestar[ni,]%o% Utildestar[ni,]/((ng_gi-1)^2/ng_gi-hi[ni]) )
        detPhi_gi[ind_g] <- coef.detPhi[ni]*detPhi[ind_g]
        #cat(detPhi_gi,"\n")

        Utilde_i <- Utilde[ni, ]
        x_i <- x[[ni]]
        U_i <- U[ni,]

        coef.alpha <- 1/(1-1/ng_gi-Hi[ni])
        Matalpha_adji <- coef.alpha*kronecker( as.matrix(Utildestar2[ni,]) , t(adjX_Mp1[ni,]) )
        alpha_adji <- fold(mat = Matalpha_adji,
                           row_idx = (M+1), col_idx = seq(1:M), modes = c(pk, q) )
        alpha_i <- alphahat - alpha_adji
        #
        mu_gi <- Xbar-rTensor::ttm(tnsr = alpha_i, mat = Matphi, m = (M+1) ) # p1*...*pM*q
        Matmu_gi <- unfold(tnsr = mu_gi, row_idx = M+1, col_idx = 1:M)@data

        #UnfoldMu_gi <- Xbar_Mp1
        UnfoldXTilde <- unfold(tnsr = XTilde, row_idx = M+1, col_idx = 1:M)@data  # n*(p1*...*pM)
        Xbar_gi <- Xbar_Mp1[ind_g, ] - (1/(ng_gi-1))*UnfoldXTilde[ni, ]
        Matmu_gi[ind_g,] <- Xbar_gi - c(Myttv(tnsr = alpha_i, v = phi_gi[[ind_g]], m = M+1)@data)
        Mu_gi <- fold(mat = Matmu_gi, row_idx = (M+1), col_idx = seq(1:M), modes = c(pk,G) )
        UnfoldMu_gi <- unfold(tnsr = Mu_gi, row_idx = (M+1), col_idx = c(1:M) )@data # G*(p1***pM)

        xres_gi <- x_i - Myttv(tnsr = alpha_i, v= U_i, m = (M+1) )@data
        ##
        for (g in 2:G ){
          beta <- as.vector(betahat)
          bg <- sum(beta*c(xres_gi))
          gammag <- PhiInv_gi[[g]]%*%phi_gi[[g]]- PhiInv_gi[[1]]%*%phi_gi[[1]]
          Dg <- PhiInv_gi[[g]] - PhiInv_gi[[1]]
          Score1 <- -0.5*diag(t(U_i)%*%Dg%*%U_i)
          Score2 <- t(gammag)%*%U_i
          mugPlus1 <- UnfoldMu_gi[g,] + UnfoldMu_gi[1,]
          ag.Tensor <- sum( beta*mugPlus1 )
          ag.Cov<-as.numeric( suppressWarnings(log(detPhi_gi[g])/detPhi_gi[1])
                              + t(gammag)%*%(phi_gi[[g]]+phi_gi[[1]])
                              - t(phi_gi[[g]])%*%Dg%*%phi_gi[[1]]
                              -2*log(pig_gi[g]/pig_gi[1]) )
          ag <- ag.Cov + ag.Tensor
          disScore<- Score1+Score2-0.5*ag+bg
        }
        #####

        disScore
      }
    }
  close(pb)
  stopCluster(workers)
  loc <- as.numeric(which(colSums(is.na(acc))==0))
  pred0 <- acc[,loc]
  loc
  pred<-apply(pred0, 2, FUN = function(x) ifelse(x>0,2,1) )
  #sum(pred[1,])
  acc1 <- apply(pred, 2, FUN = function(x) {
    CMatrix <- suppressWarnings(confusionMatrix(data=as.factor(x),
                                                reference = as.factor(y), positive="2"))
    ifelse(imbalanced, CMatrix$byClass["F1"], CMatrix$overall["Accuracy"])
  } )
  as.list(myLambda[loc[which.max(acc1)],])
}
#' pso.ReCatch
#' @description
#' Particle swarm optimization algorithm for the ReCatch model.
#'
#' @param x Input tensor.
#' @param U Input covariate matrix.
#' @param y Class label.
#' @param lam3Range The parameter interval for \eqn{\lambda_3}.
#' @param neighbor_number The exponent for calculating number of informants. Defaults to 3.
#' @param maxit_num The maximum number of iterations. Defaults to 100.
#' @param inertia_weight The exploitation constant.
#' @param setseed An integer value that serves as the random seed.
#' @return A list of the best set of parameters found.
#' @export
#'
pso.ReCatch <- function(x, U, y,
                        lam3Range = c(0,15),
                        neighbor_number=3, maxit_num=100,
                        inertia_weight = 1/(2*log(2)), setseed = NULL, imbalanced=FALSE){
  library(pso)
  library(catch)
  library(caret)
  ReAdjten <- function(x, z, y, testx=NULL, testz=NULL, is.centered=FALSE, lambda3=0){
    n <- length(x)
    dimen <- dim(x[[1]])
    nvars <- prod(dimen)
    vec_x <- matrix(0, nvars, n)
    for (i in 1:n){
      vec_x[,i] <- matrix(x[[i]], nvars, 1)
    }
    vec_x = t(vec_x)

    if (!is.null(testx)){
      tesize <- length(testx)
      vec_testx <- matrix(0, nvars, tesize)
      for (i in 1:tesize){
        vec_testx[,i] <- matrix(testx[[i]], nvars,1)
      }
      vec_testx <- t(vec_testx)
    }
    z <- as.matrix(z) #n*q
    cz <- z
    n <- dim(z)[1]
    q <- dim(z)[2]
    cvecx <- vec_x
    nclass <- as.integer(length(unique(y)))
    if (is.centered==FALSE){
      for (i in 1:nclass){
        if (q>1) {
          cz[y == i,] = sweep(z[y == i,], 2, colMeans(z[y == i,]))
        }else {
          cz[y==i]=z[y==i]-mean(z[y==i])
        }
        cvecx[y == i,] = sweep(vec_x[y == i,], 2, colMeans(vec_x[y == i,]))
      }
    }
    c <- solve( (t(cz) %*% cz)+lambda3*diag(q) ) %*% t(cz)
    c = c %*% cvecx
    vec_xres <- vec_x-z%*%c
    A = (t(cz) %*% cz)+lambda3*diag(q)

    xres = array(list(),n)
    for (i in 1:n){
      xres[[i]] = array(vec_xres[i,],dim=dimen)
    }

    if (!is.null(testx)){
      vec_testxres = vec_testx-testz %*% c
      testxres = array(list(),tesize)
      for (i in 1:tesize){
        testxres[[i]] = array(vec_testxres[i,], dim=dimen)
      }
    }else{
      testxres=NULL
    }

    muz = matrix(0, nrow=q, ncol=(nclass-1))    # muz : U_k-U_1  (k>1)
    for (i in 2:nclass){
      if (q>1){muz[,(i-1)] = apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)
      } else {
        muz[i-1] <- mean(z[y==i])-mean(z[y==1])
      }
    }
    #gamma <- solve(cov(z)) %*% muz

    return( list(xres=xres, A=A, testxres=testxres, alpha=c) )
  }
  getlam4<-function(x, U, y, lambda3){
    library(catch)
    objt <- ReAdjten(x=x, z=U, y=y, lambda3 = lambda3)
    #objt <- adjten(x = x,  z = U, y = y)
    invisible(capture.output( objtsda <- catch:::tsda(objt$xres, y=y, nlambda = 100)))
    return(objtsda$lambda)
  }
  fn <- function(lambda, x, U, y, imbalanced=FALSE){
    fit <- est.ReCatch(x=x, U=U, y=y,
                       lambda1=lambda[1], lambda2=lambda[2],
                       lambda3=lambda[3], lambda4=lambda[4])
    #yhat <- ReCatch.pred(fit=fit, x=x, U=U, y=y, testx = x, testU = U)
    yhat <- pred.ReCatch(fit=fit, testx = x, testU = U)
    y <- as.factor(y)
    G <- length(unique(y))
    yhat$pred <- factor(yhat$pred, levels = levels(y))
    #ConfMatrix <- confusionMatrix(yhat$pred, y)
    CMatrix <- suppressWarnings(confusionMatrix(data=as.factor(yhat$pred),
                                                reference = as.factor(y), positive="2"))
    conf.df <- as.data.frame(CMatrix$byClass)
    F1 <- conf.df$F1
    F1 <- ifelse(G==2, conf.df["F1",], conf.df$F1)
    F1 <- ifelse(is.na(F1), 0, F1)
    weight <- (1/table(y))
    normal_weight <- weight/sum(weight)
    wF1 <- ifelse(G==2, F1, normal_weight%*%F1)
    #wF1 <- normal_weight%*%F1
    accuracy <- ifelse(imbalanced, wF1, CMatrix$overall["Accuracy"])
    #accuracy <- ConfMatrix$overall["Accuracy"]
    error <- 1 - accuracy
    return(error)
  }
  update <- TRUE
  cat("Particle swarm optimization is performed.\nPlease wait for the results!\n")
  n.iter <- 1
  control <- list(maxit = maxit_num, trace = FALSE,
                  abstol = 1e-8,
                  w = inertia_weight, k = neighbor_number)
  if (!(is.null(setseed))) set.seed(setseed)
  while (update){
    cat("The", paste(n.iter, "-th iteration is running.",sep = ""),"\n")
    initial_lambda1 <- runif(1, 0, 1)
    initial_lambda2 <- runif(1, 0, 1)
    initial_lambda3 <- runif(1, min(lam3Range), max(lam3Range))
    #trainlam4 <- getlam4(x = x, U = U, y = y, lambda3 = initial_lambda3)
    trainlam4 <- NULL
    while (is.null(trainlam4)){
      trainlam4 <- tryCatch( getlam4(x = x, U = U, y = y, lambda3 = initial_lambda3), error=function(e) {NULL} )
    }
    #testlam4 <- getlam4(x = testx, U = testU, y = testy, initial_lambda3)
    #if (is.null(validx)) {
    min_lambda4 <- min(trainlam4)
    max_lambda4 <- max(trainlam4)
    lambda4_range <- c( min_lambda4, max_lambda4 )
    initial_lambda4 <- runif(1, min(lambda4_range), max(lambda4_range))

    lowerParm <- c(0, 0, min(lam3Range), min(lambda4_range))
    upperParm <- c(1, 1, max(lam3Range), max(lambda4_range))
    initial_Parm <- c(initial_lambda1, initial_lambda2, initial_lambda3, initial_lambda4)
    #imbalanced <- c(...)
    #cat("Imbalanced Data:", ifelse(is.null(imbalanced), FALSE, imbalanced), "\n")
    cat("Imbalanced Data:", imbalanced, "\n")
    cat("Random seed:", ifelse(is.null(setseed), FALSE, setseed), "\n")
    Optlambda <- psoptim(
      fn = fn, x = x, U = U, y = y, imbalanced=imbalanced,
      control = control,
      lower = lowerParm,
      upper = upperParm,
      par = initial_Parm
    )
    lambda1 <- Optlambda$par[1]
    lambda2 <- Optlambda$par[2]
    lambda3 <- Optlambda$par[3]
    lambda4 <- Optlambda$par[4]
    fit=tryCatch(est.ReCatch(x = x, U = U, y = y,
                             lambda1 = lambda1,
                             lambda2 = lambda2,
                             lambda3 = lambda3,
                             lambda4 = lambda4),
                 error=function(e) {NULL} )
    #print(fit)
    betahatlag <- sum(fit$betahat!=0) == 0
    #if (is.null(validx)) validx <- x
    #if (is.null(validU)) validU <- U
    yhat = tryCatch(pred.ReCatch(fit = fit, testx = x, testU = U), error=function(e) {NULL} )
    update <- is.null(fit) | is.null(yhat)  | betahatlag
    n.iter <- n.iter + 1
  }
  return(list(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4))
}


#' CatchPuls
#' @description
#' Fit a CatchPlus model and predict categorical responses.
#' @param x Input tensor.
#' @param U Input covariate matrix.
#' @param y Class label.
#' @param testx Input testing tensor of list.
#' @param testU Input testing covariate matrix.
#' @param testy Input testing class label.
#' @param setseed An integer value that serves as the random seed.
#' @param n.folds The number of folds. Defaults to 5.
#' @param imbalanced Logical; if TRUE the F1 score is used as a measure of predictive performance; otherwise, accuracy is used.
#'
#' @return
#'  \item{pred}{The predicted categorical responses.}
#'  \item{opt.lambda}{The optimized lambda parameter used to calculate the predicted categorical responses.}
#' @export
CatchPlus <- function(x, U, y, testx=x, testU=U, testy=y,
                      setseed=NULL, n.folds=5, imbalanced=TRUE){
  opt.lambda <- kf.CatchPlus(x = x, U = U, y = y, n.folds = n.folds,
                             setseed=setseed, imbalanced = imbalanced)$lambda.min
  pred <- pred.CatchPlus(x = x, U = U, y = y, testx = testx, testU = testU, testy = testy, lambda = opt.lambda)$pred
  return(list(pred=as.numeric(pred), opt.lambda=opt.lambda))
}

#' pred.CatchPlus
#' @description
#' Predict categorical responses
#' @param x Input tensor.
#' @param U Input covariate matrix.
#' @param y Class label.
#' @param testx Input testing tensor of list.
#' @param testU Input testing covariate matrix.
#' @param testy Input testing class label.
#' @param setseed An integer value that serves as the random seed.
#' @param n.folds The number of folds. Defaults to 5.
#' @param imbalanced Logical; if TRUE the F1 score is used as a measure of predictive performance; otherwise, accuracy is used.
#'
#' @return
#'  \item{pred}{Predict categorical responses.}
#'  \item{lambda}{The actual lambda sequence used.}
#' @export
pred.CatchPlus <- function(x, U, y,
                           testx, testU, testy, lambda=NULL){
  obj <- adjten(x = x, z = U, y = y, testx=testx, testz=testU)
  invisible(capture.output( objtsda <- catch:::tsda(obj$xres, y=y, nlambda = 100, lambda = lambda) ))
  nlam <- length(objtsda$beta)

  mybeta <- array(list(), nlam)
  n <- length(objtsda$y)
  pik<-objtsda$prior
  G <- length(pik)
  mubar <- array(list(), G-1 )
  gammak <- array(list(), G-1 )
  Phi <- array(list(), G)
  phi <- array(list(), G)
  Phi_inv <- array(list(), G)
  PhiInvdiff <- array(list(), G-1 )
  for (i in 1:G){
    phi[[i]] <- matrix(apply(U[y==i,], 2, mean), ncol=1)
    Phi[[i]] <- cov(U[y==i,])
    Phi_inv[[i]] <- solve(Phi[[i]])
  }

  for (i in 2:G){
    mubar[[i-1]] <- (objtsda$mu[[1]]+objtsda$mu[[i]])/2 # p1*...*pk, array
    gammak[[i-1]] <- Phi_inv[[i]]%*%phi[[i]]-Phi_inv[[1]]%*%phi[[1]]
    PhiInvdiff[[i-1]] <- Phi_inv[[i]] - Phi_inv[[1]]
  }

  testn <- length(testy)
  alpha <- k_fold(adjten_alpha(x = x, z = U, y = y)$alpha ,
                  m=length(c(dim(x[[1]]), dim(U)[2])),
                  modes=c(dim(x[[1]]), dim(U)[2]) )@data #p1*...pk*q
  #testn<- length(testy)
  adjx <- array(list(), testn)
  adjx<-adjten(x = x, z = U, y = y, testx = testx, testz = testU)$testxres #array
  pred<- matrix(NA, testn, nlam)
  tensor1 <- matrix(NA, nrow = testn, ncol = 2 ) #
  covTerm <- matrix(NA, nrow = testn, ncol = 2 )
  for (w in 1:nlam){
    Beta <- objtsda$beta[[w]]
    #gammak <- Phi2_inv%*%phi[[2]]-Phi1_inv%*%phi[[1]]
    #library(rTensor);library(tensr)
    discScore<-matrix(NA, testn, G )
    discScore[,1] <- rep(0, testn)
    ak <- rep(NA, 2)
    for (j in 2:G ){
      beta <- as.vector(Beta[,j-1])
      for (i in 1:testn) {
        tensor1[i,j-1] <- sum( beta*vec(as.tensor(adjx[[i]])) )
        covTerm[i,j-1] <- -0.5*t(testU[i,])%*%PhiInvdiff[[j-1]]%*%testU[i,]+t(gammak[[j-1]])%*%testU[i,]
      }
      tensor2 <- sum( beta*vec(as.tensor(mubar[[j-1]])) )
      ak[j-1]<-as.numeric( 0.5*t(phi[[j]])%*%PhiInvdiff[[j-1]]%*%phi[[1]]
                           -0.5*log(det(Phi[[j]])/det(Phi[[1]]) )
                           -0.5*t(gammak[[j-1]])%*%(phi[[j]]+phi[[1]])+log(pik[j]/pik[1]) )-tensor2
      discScore[,j] <- covTerm[,j-1]+tensor1[,j-1]+ak[j-1]
    }
    pred[,w] <- apply(discScore, 1, which.max)
  }
  return(list(pred=pred, lambda=objtsda$lambda))
}


#' kf.CatchPlus
#' @description
#' K-fold cross-validation optimization algorithm for the CatchPlus model.
#'
#' @param x Input tensor.
#' @param U Input covariate matrix.
#' @param y Class label.
#' @param n.folds The number of folds. Defaults to 5.
#' @param imbalanced Logical; if TRUE the F1 score is used as a measure of predictive performance; otherwise, accuracy is used.
#' @param setseed An integer value that serves as the random seed.
#' @param lambda.opt The optimal criteria when multiple elements in lambda return the same minimum classification error. "min" will return the smallest lambda with minimum cross validation error. "max" will return the largest lambda with the minimum cross validation error.
#'
#' @return
#'   \item{lambda}{The actual lambda sequence used.}
#'   \item{cvm}{The mean of cross validation errors for each lambda.}
#'   \item{cvsd}{The standard error of cross validaiton errors for each lambda.}
#'   \item{lambda.min}{The lambda with minimum cross validation error. If lambda.opt is min, then returns the smallest lambda with minimum cross validation error. If lambda.opt is max, then returns the largest lambda with minimum cross validation error.}
#'   \item{lambda.1se}{The largest lambda with cross validation error within one standard error of the minimum.}
#' @export
#'
kf.CatchPlus <- function(x, U, y, n.folds, imbalanced, setseed, lambda.opt = "min"){
  library(splitTools)
  library(rTensor)
  y <- drop(y)
  n <- length(x)
  dimen <- dim(x[[1]])
  nvars <- prod(dimen)
  G <- length(unique(y))
  prior <- rep(0, G)
  for (i in 1:G) {
    prior[i] <- mean(y == i)
  }
  U <- as.matrix(U)
  q <- dim(U)[2]
  obj <- adjten(x, U, y)
  obj.tsda <- catch:::tsda(obj$xres, y)

  lambda <- obj.tsda$lambda
  if (n.folds < 3)
    stop("n.folds must be bigger than 3; n.folds=5 recommended")
  if (n.folds > n)
    stop("The number of folds should be smaller than the sample size.")
  residmat <- matrix(NA, n.folds, length(lambda))
  if (!(is.null(setseed))) set.seed(setseed)
  folds <- splitTools::create_folds(y=y, k=n.folds, type="stratified")
  cat(paste0(n.folds,"-fold validation is performed.\nPlease wait for the results!\n"))
  cat("Imbalanced Data: ", imbalanced,"\n")
  cat("Random seed:", ifelse(is.null(setseed), FALSE, setseed), "\n")
  for (m in 1:n.folds) {
    trainID <- folds[[m]]
    trainY <- y[trainID]; trainX <- x[trainID]; trainU <- U[trainID,]
    testY <- y[-trainID]; testX <- x[-trainID]; testU <- U[-trainID,]
    fit <- pred.CatchPlus(x = trainX, U = trainU, y = trainY,
                          testx = testX, testU = testU, testy = testY, lambda = lambda)
    nlami <- length(fit$lambda)
    loc.lambda <- which(lambda %in% fit$lambda)
    accuracy<- rep(NA, nlami)
    yhatMat <- fit$pred
    accuracy <- apply(yhatMat, MARGIN = 2, function(y){
      yhat <- factor(y, levels=levels(as.factor(y) ) )
      CMatrix <- suppressWarnings(caret::confusionMatrix(data=as.factor(yhat),
                                                         reference = as.factor(testY), positive="2"))
      conf.df <- as.data.frame(CMatrix$byClass)
      #F1 <- conf.df$F1
      F1 <- ifelse(G==2, conf.df["F1",], conf.df$F1)
      F1 <- ifelse(is.na(F1), 0, F1)
      weight <- (1/table(testY))
      normal_weight <- weight/sum(weight)
      #wF1 <- normal_weight%*%F1
      wF1 <- ifelse(G==2, F1, normal_weight%*%F1)
      #accuracy[j] <- ifelse(imbalanced, wF1, CMatrix$overall["Accuracy"])
      ifelse(imbalanced, wF1, CMatrix$overall["Accuracy"])
    })
    residmat[m, loc.lambda] <- 1-accuracy
  }
  residmat[is.na(residmat)] <- 1
  cvm <- colMeans(residmat, na.rm = TRUE)
  cvsd <- sqrt(colMeans(scale(residmat, cvm, FALSE)^2, na.rm = TRUE)/(n.folds -
                                                                        1))
  if (lambda.opt == "min") {
    lambda.min <- min(lambda[which(cvm == min(cvm, na.rm = TRUE))])
  } else {
    lambda.min <- max(lambda[which(cvm == min(cvm, na.rm = TRUE))])
  }
  id.min <- which.min(cvm)
  if (cvsd[id.min] == 0) {
    lambda.1se = lambda.min
  }  else {
    lambda.1se <- max(lambda[cvm < min(cvm, na.rm = TRUE) +
                               cvsd[id.min]], na.rm = TRUE)
  }
  return(list(lambda = lambda, cvm = cvm, cvsd = cvsd, lambda.min = lambda.min,
              lambda.1se = lambda.1se))
}
