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





##' Perfrom a parametric bootstrap for a fitted nls model.
##'
##' The user must provide a function to compute the test statistic
##' from a fitted nls object.  The test statistic is computed for each
##' bootstrap sample and the results are returned as an array.
##'
##' @title NLS Parametric Bootstrap
##' @param object a fitted nls object.
##' @param nboot number of bootstrap samples.
##' @param stat a function that takes a fitted nls object and returns
##'   the test statistic.
##' @return an array of bootstrap samples.
##' @importFrom stats formula coef
##' @export
nlsParBoot <- function(object,nboot=99,stat=coef) {
  ## Regenerate the original data and name of response
  d <- eval(object$data,environment(formula(object)))
  resp <- as.character(formula(object)[[2]])
  ## Ensure simulate predicts for the NA values
  object.pr <- object
  if(!is.null(object.pr$na.action) && class(object.pr$na.action)!="exclude")
    class(object.pr$na.action) <- "exclude"
  sapply(simulate(object.pr,nsim=nboot),
         function(y) {
           d[[resp]] <- y
           cl <- eval(substitute(update(object,data=data,start=start,evaluate=FALSE),
                                 list(data=d,start=coef(object))))
           stat(eval(cl,environment(formula(object))))
         })
}



##' Use missing value information to extend a vector
##'
##' Cases with missing values are removed from teh fitted values and
##' weights stored in \code{nls}.  This function reintroduces missing
##' values at the right locations so that the weights and fitted
##' values match the original data source.
##'
##' @title Adjust for missing values
##' @param omit an object creates by an \code{\link{na.action}} function
##' @param x a vector
##' @return  a vector
na_extend <- function(omit,x) {
  if(!is.null(omit)) {
    keep <- rep.int(NA, length(x)+length(omit))
    keep[-omit] <- seq_along(x)
    x <- x[keep]
    if(!is.null(names(x))) names(x)[omit] <- names(omit)
  }
  x
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
##' @param stat a function that takes a fitted nls object and returns
##'   the test statistic.
##' @return an array of bootstrap samples.
##' @importFrom stats formula coef rexp
##' @export
nlsBayesBoot <- function(object,nboot=99,stat=coef) {
  ## Extract model weights
  w <- object$weights
  if(is.null(w)) w <- rep(1,length(object$m$fitted()))
  w <- na_extend(object$na.action,w)
  ## Only non-missing cases get reweighted
  k <- which(!is.na(w))
  w[is.na(w)] <- 0
  sapply(seq_len(nboot),
         function(.) {
           ## Reweight by dirichlet weights
           b <- rexp(length(k),1)
           w[k] <- w[k]*(b/sum(b))
           ## Update call and evaluate in the environment of the model formula
           cl <- eval(substitute(update(object,weights=weights,start=start,evaluate=FALSE),
                                 list(weights=w,start=coef(object))))
           stat(eval(cl,environment(formula(object))))
         })
}


##' Calculate DFBETAS and DFFITS for a fitted nls model
##'
##' Returns DFBETAS as columns labelled by the corresponding parameter
##' names and DFFITS labelled by the name of the response.
##' @title NLS Influence Measures
##' @param object a fitted nls object.
##' @return A matrix where each row corresponds to an observation and
##'   each column a DFBETAS or DFFITS
##' @importFrom stats formula coef setNames
##' @export
nlsInfluence <- function(object) {
  ## Model summary and name of response
  smry <- summary(object)
  resp <- as.character(formula(object)[[2]])
  ## Extract model weights, fitted values and diagonals of the hat matrix
  w <- object$weights
  if(is.null(w)) w <- rep(1,length(object$m$fitted()))
  w <- na_extend(object$na.action,w)
  f <- na_extend(object$na.action,object$m$fitted())
  g <- object$m$gradient()
  h <- na_extend(object$na.action,diag(g%*%smry$cov.unscaled%*%t(g)))
  ## Only non-missing cases get reweighted
  k <- which(!is.na(w))
  w[is.na(w)] <- 0
  dfs <- t(sapply(k,
                  function(k) {
                    ## Weight out current case
                    w[k] <- 0
                    ## Update call and evaluate in the environment of the model formula
                    cl <- eval(substitute(update(object,weights=weights,start=start,evaluate=FALSE),
                                          list(weights=w,start=coef(object))))
                    object.k <- eval(cl,environment(formula(object)))
                    smry.k <- summary(object.k)
                    f.k <- na_extend(object.k$na.action,object.k$m$fitted())
                    c(
                      ## DFBETAS
                      (coef(object)-coef(object.k))/(smry.k$sigma*sqrt(diag(smry$cov.unscaled))),
                      ## DFFITS
                      setNames((f[k]-f.k[k])/(smry.k$sigma*sqrt(h[k])),resp))

                  }))
  rownames(dfs) <- k
  dfs
}

