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
##' @param x a vector
##' @return  a vector
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





##' Perform a parametric bootstrap for a fitted nls model.
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
##' @importFrom stats formula coef update
##' @export
nlsParBoot <- function(object,nboot=99,stat=coef) {
  ## Extract the data and (best) starting values
  params <- names(environment(object$m$getEnv)$ind)
  values <- as.list(object$m$getEnv())
  data <-  as.data.frame(values[!(names(values) %in% params)])
  start <-  as.data.frame(values[(names(values) %in% params)])
  ## Replace missing values
  data <- na_extend(object$na.action,data)
  ## Determine response
  resp <- as.character(formula(object)[[2]])
  ## Ensure simulate predicts for the NA values
  object.sim <- object
  if(!is.null(object.sim$na.action) && class(object.sim$na.action)!="exclude")
    class(object.sim$na.action) <- "exclude"
  sapply(simulate(object.sim,nsim=nboot),
         function(y) {
           data[[resp]] <- y
           stat(update(object,data=data,start=start))
         })
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
##' @importFrom stats formula coef rexp update
##' @export
nlsBayesBoot <- function(object,nboot=99,stat=coef) {
  ## Extract the data and (best) starting values
  params <- names(environment(object$m$getEnv)$ind)
  values <- as.list(object$m$getEnv())
  data <-  as.data.frame(values[!(names(values) %in% params)])
  start <-  as.data.frame(values[(names(values) %in% params)])
  ## Extract model weights
  w <- object$weights
  if(is.null(w)) w <- rep(1,length(object$m$fitted()))
  ## Reinsert missing values
  data <- na_extend(object$na.resid,data)
  w <- na_extend(object$na.action,w)
  ## Only reweight non-missing cases
  k <- which(!is.na(w))
  w[is.na(w)] <- 0
  sapply(seq_len(nboot),
         function(.) {
           ## Reweight by dirichlet weights
           b <- rexp(length(k),1)
           w[k] <- w[k]*(b/sum(b))
           stat(update(object,data=data,start=start,weights=w))
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
##' @importFrom stats formula coef setNames update
##' @export
nlsInfluence <- function(object) {
  ## Model summary and name of response
  smry <- summary(object)
  response <- as.character(formula(object)[[2]])
  ## Extract the data and (best) starting values
  params <- names(environment(object$m$getEnv)$ind)
  values <- as.list(object$m$getEnv())
  data <-  as.data.frame(values[!(names(values) %in% params)])
  start <-  as.data.frame(values[(names(values) %in% params)])
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

