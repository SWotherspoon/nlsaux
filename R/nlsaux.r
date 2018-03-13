##' Simulate one or more responses from the distribution corresponding to a fitted \code{\link{nls}} object.
##'
##' This is a simple wrapper function for \code{simulate.lm}.
##'
##' @title Simulate Responses.
##' @param object an object representing a fitted model.
##' @param nsim number of response vectors to simulate. Defaults to 1.
##' @param seed an object specifying if and how the random number
##'   generator should be initialized ('seeded').
##' @param ... additional optional arguments.
##' @return A dataframe with attribute "seed", where each column
##'   represents a set of simulated responses.
##' @importFrom stats simulate df.residual deviance fitted rnorm runif
##' @export
simulate.nls <- function (object, nsim = 1, seed = NULL, ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ftd <- fitted(object)
  n <- length(ftd)
  vars <- deviance(object)/df.residual(object)
  if(!is.null(object$weights)) vars <- vars/object$weights
  val <- rnorm(n*nsim,ftd,sqrt(vars))
  dim(val) <- c(n, nsim)
  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  row.names(val) <- names(ftd)
  attr(val, "seed") <- RNGstate
  val
}


##' Use missing value information to extend a vector
##'
##' The data, weights and fitted values stored within an
##' \code{nlsModel} object have had cases with any missing values
##' filtered out.  This function reintroduces missing values at the
##' right locations so that the data, weights and fitted values match
##' the original data source.
##'
##' @title Adjust for missing values
##' @param omit an object creates by an \code{\link{na.action}} function
##' @param x a vector or matrix
##' @return a vector or matrix
na_extend <- function(omit,x) {
  if(!is.null(omit)) {
    if(is.data.frame(x)) {
      keep <- rep.int(NA, nrow(x)+length(omit))
      keep[-omit] <- seq_len(nrow(x))
      x <- x[keep,,drop=FALSE]
    } else {
      keep <- rep.int(NA, length(x)+length(omit))
      keep[-omit] <- seq_len(length(x))
      x <- x[keep]
    }
  }
  x
}

##' Calculate DFBETAS and DFFITS for a fitted nls model
##'
##' Returns DFBETAS as columns labelled by the corresponding parameter
##' names and DFFITS labelled by the name of the response.
##' @title NLS Influence Measures
##' @param object a fitted nls object.
##' @return A matrix where each row corresponds to an observation and
##'   each column a DFBETAS or DFFITS
##' @importFrom stats formula coef setNames update
##' @export
nlsInfluence <- function(object) {
  ## Extract the data and (best) starting values
  params <- names(environment(object$m$getEnv)$ind)
  values <- as.list(object$m$getEnv())
  data <-  as.data.frame(values[!(names(values) %in% params)])
  start <-  values[params]
  ## Model summary and name of response
  smry <- summary(object)
  response <- as.character(formula(object)[[2]])
  ## Extract model weights, fitted values and diagonals of the hat matrix
  w <- object$weights
  if(is.null(w)) w <- rep(1,length(object$m$fitted()))
  f <- object$m$fitted()
  g <- object$m$gradient()
  if(inherits(object$m,"nlsModel.plinear"))
    g <- cbind(apply(g,c(1,3),'%*%',object$m$getAllPars()[-seq_len(dim(g)[3])]),
               eval(formula(object)[[3L]],object$m$getEnv()))
  h <- diag(g%*%smry$cov.unscaled%*%t(g))
  ## Reinsert missing values
  data <- na_extend(object$na.action,data)
  w <- na_extend(object$na.action,w)
  f <- na_extend(object$na.action,f)
  h <- na_extend(object$na.action,h)
  ## Only reweight non-missing cases
  k <- which(!is.na(w))
  w[is.na(w)] <- 0
  dfs <- t(sapply(k,
                  function(k) {
                    ## Weight out current case
                    w[k] <- 0
                    ## Update fit
                    object.k <- update(object,data=data,start=start,weights=w)
                    smry.k <- summary(object.k)
                    f.k <- na_extend(object.k$na.action,object.k$m$fitted())
                    c(## DFBETAS
                      (coef(object)-coef(object.k))/(smry.k$sigma*sqrt(diag(smry$cov.unscaled))),
                      ## DFFITS
                      setNames((f[k]-f.k[k])/(smry.k$sigma*sqrt(h[k])),response))
                  }))
  rownames(dfs) <- k
  dfs
}




##' Perform the Jackknife for a fitted nls model.
##'
##' The user must provide a function to compute the test statistic
##' from a fitted nls object.  The test statistic is computed for each
##' jacknife sample.
##'
##' @title NLS Jackknife
##' @param object a fitted nls object.
##' @param statistic a function that takes a fitted nls object as its first
##'   argument and returns a vector of test statistics.
##' @return Returns an object of class \code{nlsJKnife} containing
##' \item{\code{call}}{the matched call}
##' \item{\code{statistic}}{the test statistic}
##' \item{\code{nls}}{the fitted nls object}
##' \item{\code{jknife}}{an array of jackknife samples}
##' @importFrom stats formula coef update
##' @export
nlsJKnife <- function(object,statistic=coef) {
  cl <- match.call()
  ## Extract the data and (best) starting values
  params <- names(environment(object$m$getEnv)$ind)
  values <- as.list(object$m$getEnv())
  data <-  as.data.frame(values[!(names(values) %in% params)])
  start <-  values[params]
  ## Extract model weights, fitted values and diagonals of the hat matrix
  w <- object$weights
  if(is.null(w)) w <- rep(1,length(object$m$fitted()))
  ## Reinsert missing values
  data <- na_extend(object$na.action,data)
  w <- na_extend(object$na.action,w)
  ## Only reweight non-missing cases
  k <- which(!is.na(w))
  w[is.na(w)] <- 0
  jknife <- sapply(k,
                   function(k) {
                     ## Weight out current case
                     w[k] <- 0
                     ## Update fit
                     statistic(update(object,data=data,start=start,weights=w))
                   })
  jknife <- if(is.null(dim(jknife))) t(jknife) else jknife

  structure(list(call=cl,statistic=statistic,nls=object,jknife=jknife),
            class="nlsJKnife")
}



##' Constructs a table of estimates and confidence intervals for
##' the test statistic.
##'
##' @title Summary for nlsJKnife
##' @param object an object of class \code{nlsJKnife}
##' @param level the confidence level required.
##' @param x an object of class \code{summary.nlsJKnife}
##' @param digits the number of significant digits to use when printing
##' @param ... currently ignored
##' @return Returns an object of class \code{nlsJKnife} containing
##' \item{\code{object}}{the \code{nlsJKnife} object}
##' \item{\code{theta}}{estimate of the test statistic}
##' \item{\code{corrected}}{jackknife bias corrected estimate}
##' \item{\code{bias}}{jackknife estimate of bias}
##' \item{\code{cov}}{jackknife estimate of covariance}
##' \item{\code{level}}{requested confidence level}
##' \item{\code{lwr}}{jackknife lower confidence limit}
##' \item{\code{upr}}{jackknife upper confidence limit}
##' @importFrom stats qt
##' @export
summary.nlsJKnife <- function(object,level=0.95,...) {

  a <- (1-level)/2
  n <- dim(object$jknife)[2]
  npar <- length(coef(object$nls))
  theta <- object$statistic(object$nls)
  jkmean <- apply(object$jknife,1,mean)
  cov <- (n-1)/n*tcrossprod(object$jknife-jkmean)
  bias <- (n-1)*(jkmean-theta)
  corrected <- theta-bias
  se <- sqrt(diag(cov))
  ci <- corrected+se%o%qt(c(a,1-a),n-npar)
  colnames(ci) <- paste(format(100*c(a,1-a),trim=TRUE,scientific=FALSE,digits=3),"%",sep="")
  structure(list(object=object,
                 theta=theta,
                 corrected=corrected,
                 bias=bias,
                 cov=cov,
                 level=level,
                 ci=ci),
            class="summary.nlsJKnife")
}


##' @rdname summary.nlsJKnife
##' @export
print.summary.nlsJKnife <- function(x,digits = max(3L, getOption("digits") - 3L),...) {
  cat("\nCall:\n",
      paste(deparse(x$object$call),sep = "\n",collapse = "\n"),
      "\n\n", sep = "")

  a <- (1-x$level)/2
  tab <- cbind(x$theta,x$bias,x$corrected,sqrt(diag(x$cov)),x$ci)
  colnames(tab) <- c("Estimate","Bias","Corrected","SE",
                     paste(format(100*c(a,1-a),trim=TRUE,scientific=FALSE,digits=3),"%",sep=""))
  rownames(tab) <- names(x$theta)

  printCoefmat(tab,digits=digits,has.Pvalues=FALSE)
}




##' Perform a parametric bootstrap for a fitted nls model.
##'
##' The user must provide a function to compute the test statistic
##' from a fitted nls object.  The test statistic is computed for each
##' bootstrap sample and the results are returned as an array.
##'
##' @title NLS Parametric Bootstrap
##' @param object a fitted nls object.
##' @param nboot number of bootstrap samples.
##' @param statistic a function that takes a fitted nls object as its first
##'   arguement and returns a vector of test statistics.
##' @return Returns an object of class \code{nlsBoot} containing
##' \item{\code{call}}{the matched call}
##' \item{\code{statistic}}{the test statistic}
##' \item{\code{nls}}{the fitted nls object}
##' \item{\code{boot}}{an array of bootrap samples}
##' @importFrom stats formula coef update
##' @export
nlsParBoot <- function(object,nboot=99,statistic=coef) {
  cl <- match.call()
  ## Store seed
  if(!exists(".Random.seed")) set.seed(NULL)
  seed <- .Random.seed
  ## Extract the data and (best) starting values
  params <- names(environment(object$m$getEnv)$ind)
  values <- as.list(object$m$getEnv())
  data <-  as.data.frame(values[!(names(values) %in% params)])
  start <-  values[params]
  ## Replace missing values
  data <- na_extend(object$na.action,data)
  ## Determine response
  resp <- as.character(formula(object)[[2]])
  ## Ensure simulate predicts for the NA values
  object.sim <- object
  if(!is.null(object.sim$na.action) && class(object.sim$na.action)!="exclude")
    class(object.sim$na.action) <- "exclude"
  boot <- lapply(simulate(object.sim,nsim=nboot),
                 function(y) {
                   data[[resp]] <- y
                   tryCatch(statistic(update(object,data=data,start=start)),
                            error=function(e) NULL)
                 })
  structure(list(call=cl,statistic=statistic,nls=object,boot=boot,seed=seed),
            class="nlsBoot")
}





##' Perform a Bayesian bootstrap for a fitted nls model.
##'
##' The user must provide a function to compute the test statistic
##' from a fitted nls object.  The test statistic is computed for each
##' bootstrap sample and the results are returned as an array.
##'
##' @title NLS Bayesian Bootstrap
##' @param object a fitted nls object.
##' @param nboot number of bootstrap samples.
##' @param statistic a function that takes a fitted nls object as its first
##'   argument and returns a vector of test statistics.
##' @return Returns an object of class \code{nlsBoot} containing
##' \item{\code{call}}{the matched call}
##' \item{\code{statistic}}{the test statistic}
##' \item{\code{nls}}{the fitted nls object}
##' \item{\code{boot}}{an array of bootrap samples}
##' @importFrom stats formula coef rexp update
##' @export
nlsBayesBoot <- function(object,nboot=99,statistic=coef) {
  cl <- match.call()
  ## Store seed
  if(!exists(".Random.seed")) set.seed(NULL)
  seed <- .Random.seed
  ## Extract the data and (best) starting values
  params <- names(environment(object$m$getEnv)$ind)
  values <- as.list(object$m$getEnv())
  data <-  as.data.frame(values[!(names(values) %in% params)])
  start <-  values[params]
  ## Extract model weights
  w <- object$weights
  if(is.null(w)) w <- rep(1,length(object$m$fitted()))
  ## Reinsert missing values
  data <- na_extend(object$na.resid,data)
  w <- na_extend(object$na.action,w)
  ## Only reweight non-missing cases
  k <- which(!is.na(w))
  w[is.na(w)] <- 0
  boot <- sapply(seq_len(nboot),
                 function(.) {
                   ## Reweight by dirichlet weights
                   b <- rexp(length(k),1)
                   w[k] <- w[k]*(b/sum(b))
                   tryCatch(statistic(update(object,data=data,start=start,weights=w)),
                            error=function(e) NULL)
                 })
  structure(list(call=cl,statistic=statistic,nls=object,boot=boot,seed=seed),
            class="nlsBoot")
}



##' Constructs a table of estimates and boostrap confidence intervals for
##' the test statistic.
##'
##' @title Summary for nlsBoot
##' @param object an object of class \code{nlsBoot}
##' @param level the confidence level required.
##' @param method confidence interval type
##' @param x an object of class \code{summary.nlsBoot}
##' @param digits the number of significant digits to use when printing
##' @param ... ignored
##' @importFrom stats quantile
##' @return Returns the estimates as an array.
##' @export
summary.nlsBoot <- function(object,level=0.95,method=c("basic","percentile","Normal"),...) {

  method <- match.arg(method)
  a <- (1-level)/2
  theta <- object$statistic(object$nls)
  failed <- sapply(object$boot,is.null)
  bs <- simplify2array(object$boot[!failed],FALSE)
  if(is.null(dim(bs))) bs <- t(bs)
  ci <- switch(method,
               basic=t(apply(bs,1,quantile,c(a,1-a))),
               percentile=2*theta-t(apply(bs,1,quantile,c(a,1-a)))[,2:1],
               Normal=theta+apply(bs,1,sd)%o%qnorm(c(a,1-a)))

  structure(list(object=object,
                 method=method,
                 level=level,
                 theta=theta,
                 ci=ci,
                 nfailed=sum(failed)),
            class="summary.nlsBoot")
}

##' @rdname summary.nlsBoot
print.summary.nlsBoot <- function(x,digits = max(3L, getOption("digits") - 3L),...) {
  cat("\nCall:\n",
      paste(deparse(x$object$call),sep = "\n",collapse = "\n"),
      "\n\n", sep = "")
  if(nfailed>0) cat("Convergence Failures:",nfailed,"\n\n")

  a <- (1-x$level)/2
  tab <- cbind(x$theta,ci)
  colnames(tab) <- c("Estimate",
                     paste(format(100*c(a,1-a),trim=TRUE,scientific=FALSE,digits=3),"%",sep=""))
  rownames(tab) <- names(x$theta)
  printCoefmat(x$est,digits=digits,has.Pvalues=FALSE)
}


##' Normal QQ plots of the bootstrap samples of an \code{nlsBoot} object.
##'
##' @title Normal QQ Plots for nlsBoot
##' @param y an object of class \code{nlsBoot}
##' @param which a numeric vector indicating with statistics to plot
##' @param ... extra arguments to be passed to \code{qqnorm}
##' @export
qqnorm.nlsBoot <- function(y,which,...) {
  failed <- sapply(y$boot,is.null)
  bt <- simplify2array(object$boot[!failed],FALSE)
  if(missing(which)) which <- seq_len(nrow(bt))
  for(k in which) qqnorm(bt[k,],...)
}

