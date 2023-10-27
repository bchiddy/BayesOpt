# Implementation of the Bayesian Optimization Algorithm 
# Author: Ben Chiddy
# Date: 24/10/2023
# Reference: (Gramacy, 2020)

initialize_values <- function(num_par, pmin, pmax, ninit)
{
  ## Sample initial values using a rescaled LHS scheme
  lhs          <- randomLHS(ninit, num_par)     # Get LHS 
  X            <- pmin + (pmax - pmin) * lhs    # Re-scale LHS to param. space
  return(X)
}

## Main wrapper for bayes opt algorithm

bayes_opt <- function(f, X, num_par, ninit, end, pmin=0, pmax=1, noise=TRUE, 
                      names, type="ef", plt=FALSE, resolution=200, lambda=5,
                      log_surfaces=FALSE)
{
  #' Performs Surrogate Assisted (Bayesian) Optimization to find an 'optimal' 
  #' set of parameter values for a noisy (or deterministic) function within a 
  #' specified budget. The Gaussian Process is used as a surrogate function 
  #' which maps parameter values to observed outputs. Uses an active learning 
  #' (sequential design) scheme where the next parameter set is selected by
  #' an acquisition function. Function will return both a data frame containing 
  #' parameters attempted and their outputs as well as a list containing the 
  #' predictive surfaces at each time point for animation the algorithm.
  #' -------------------------------------------------------------------------
  #' @param f           object of type Function. Function to evaluate pars.
  #' @param X           object of type Numeric. Matrix of starting values.
  #' @param ninit       object of type Numeric. Amount of initial points to 
  #'                      sample (using LHS) before starting the BO algorithm. 
  #' @param num_par     object of type Numeric. Number of par values needed.
  #' @param pmin        object of type Numeric. Lower bound for par values.
  #' @param pmax        object of type Numeric. Upper bound for par values.
  #' @param end         object of type Numeric. Maximum number a iterations.        
  #' @param names       object of type Character. Vector of hp names.
  #' @param type        object of type Character. Acquisition function.
  #'                      "ef"  = Mean Critereon
  #'                      "ei"  = Expected Improvement 
  #'                      (Jones et al., 1998)
  #'                      "ucb" = Upper Confidence Bound 
  #'                      (Srinivas et al., 2009)
  #' @param plt         object of type Character. Specify whether or not 
  #'                      plotting should be shown during runtime. Used for
  #'                      debugging purposes.
  #' @param resolution  object of type Numeric. Resolution of surface plots.
  #' @param lambda      object of type Numeric. Lambda param for UCB AF.
  #' -------------------------------------------------------------------------
  #' @return            object of type List. Data frame of parameters and 
  #'                    associated outputs on f as well as a list of surface 
  #'                    states collected during runtime.
  #' -----------------------------------------------------------------------
  
  ## Initialize storage objects
  cols         <- c("darkblue", "blue", "lightblue", "white", "orange", "red")
  colour       <- colorRampPalette(cols)(100)
  prd_mean     <- list()                        # For posterior mean
  prd_var      <- list()                        # For posterior variance 
  step         <- 1                             # Initialize step counter
  kpars        <- list()                        # Stores kernel parameters
  y            <- c() 
  
  ## Run f at ninit initial locations
  for (i in 1:ninit) {
    params     <- X[i,]                    
    y[i]       <- do.call(f, list(params))  # Run f at parameter setting
  }
  
  ## Fits GP using initially sampled locations using laGP (Gramacy, 2020).
  ## Note that if function being evaluated is noise corrupted, a nugget is 
  ## included in the Gaussian Process to account for extra variation.
  if (!noise) {
    
    # Create vague prior for kernel parameters for Bayesian inference
    da         <- darg(list(mle=TRUE, max=pmax), randomLHS(1000, num_par))
    kpars[[1]] <- da # Save prior kernel parameters
    
    # Fit Gaussian Process using MLE
    gpi        <- newGPsep(X, y, d=da$start, g=1e-6, dK=TRUE) 
    mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)  
    
  } else {
    
    # Create vague prior for kernel parameters for Bayesian inference
    d       <- darg(NULL, X)            # Lengthscale
    g       <- garg(list(mle=TRUE), y)  # Nugget
    kpars[[1]] <- d                     # Save prior lengthscale
    kpars[[2]] <- g                     # Save prior nugget
    
    # Fit Gaussian Process using MLE
    gpi     <- newGPsep(X, y, d=d$start, g=1e-6, dK=TRUE) 
    mleGPsep(gpi, param="both", tmin=c(d$min, g$min), 
             tmax=c(d$max, g$max),ab=c(d$ab, g$ab))
    
  }
  
  if (log_surfaces) {
    
    ## Create a grid over parameter space to create predictive surfaces
    xx     <- seq(pmin, pmax, length = resolution)
    XX     <- expand.grid(xx, xx)
    
    ## Create initial predictive surfaces for param space and save in list
    prd               <- predGPsep(gpi, XX, lite = TRUE) 
    prd_mean[[step]]  <- prd$mean     # Posterior mean surface 
    prd_var[[step]]   <- prd$s2       # Posterior variance mean surface
    step              <- step + 1     # Increment timestep
    
  }
  
  ## Plots initial predictive surfaces for param space. Used for Debugging.
  if (plt) {
    par(mfrow = c(1, 2)) 
    
    # Plot posterior mean surface
    Z <- matrix(prd_mean[[step-1]], ncol=length(xx))
    image(xx, xx, Z, col=colour, ylab="", xlab="")
    contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
    title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
    points(X, pch = 1)
    
    # Plot posterior variance surface
    Z <- matrix(prd_var[[step-1]], ncol=length(xx))
    image(xx, xx, Z, col=colour, ylab="", xlab="")
    contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
    title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
    points(X, pch = 1)
  }
  
  ## Start BO algorithm using selected acquisition strategy.
  if (type == "ef") {
    
    # Predicted improvement acquisition strategy 
    res <- mean_critereon(f, gpi, X, y, log_surfaces=log_surfaces, 
                          end=end, prd_mean=prd_mean, prd_var=prd_var,
                          kpars=kpars, plt=plt, ninit=ninit, 
                          noise=noise, pmin=pmin, pmax=pmax, 
                          colour=colour)
    
  } else if (type == "ei") {
    
    # Expected improvement acquisition strategy
    res <- expected_improvement(f, gpi, X, y, log_surfaces=log_surfaces, 
                                prd_mean=prd_mean, prd_var=prd_var,
                                end=end, kpars=kpars, plt=plt, 
                                ninit=ninit, noise=noise, colour=colour,
                                pmin=pmin, pmax=pmax)
    
  } else if (type == "ucb") {
    
    # Upper confidence bound acquisition strategy
    res <- ucb_improvement(f, gpi, X, y, log_surfaces=log_surfaces, 
                           prd_mean=prd_mean, prd_var=prd_var,
                           end=end, kpars=kpars, plt=plt, 
                           ninit=ninit, noise=noise, colour=colour,
                           pmin=pmin, pmax=pmax, lambda=lambda)
  }
  
  ## Combine parameter vectors and outputs on f as well as list of surfaces
  D           <- data.frame(res$data)
  colnames(D) <- c(names, "y")       
  bo_results  <- list(data = D, prd_mean = res$prd_mean, prd_var = res$prd_var)
  
  if (!log_surfaces) { bo_results  <- list(data = D) }
  
  return(bo_results)
}

## Mean criterion acquisition policy (Gramacy, 2020)

predicted_value <- function(x, gpi) 
{
  #' Helper function for Predicted Improvement acquisition strategy. Returns
  #' predicted value at a given coordinate within the surrogate param space.
  #' -------------------------------------------------------------------------
  #' @param x         object of type Numeric. Vector of param coordinates.
  #' @param gpi       object of type laGP. Fitted Gaussian Process.
  #' -------------------------------------------------------------------------
  #' @return          object of type Numeric. Predicted f value in par space.
  #' -------------------------------------------------------------------------
  
  return(predGPsep(gpi, matrix(x, nrow=1), lite=TRUE)$mean)
}

mean_critereon <- function(f, gpi, X, y, log_surfaces, prd_mean, prd_var, 
                           step=2, end, kpars, plt=FALSE, ninit, 
                           noise=TRUE, colour, pmin, pmax, resolution=200) 
{
  #' Runs Bayesian Optimization using the Mean Critereon acquisition 
  #' strategy. Not intended to be run by user, access provided through 
  #' bayes_opt function. Returns list of BO algorithm results as well as the 
  #' predictive surfaces at each timestep.
  #' -----------------------------------------------------------------------
  #' @param f         object of type Function. Function to evaluate hps over.
  #' @param gpi       object of type laGP. Fitted Gaussian Process. 
  #' @param X         object of type Matrix Array. Initial param values tested.
  #' @param y         object of type Numeric. Values for f obtained at X.
  #' @param xx        object of type Numeric. Values to create hp grid.
  #' @param XX        object of type Data Frame. Prediction grid for plotting.
  #' @param prd_mean  object of type List. Posterior mean at each timestep.
  #' @param prd_var   object of type List. Posterior var. at each timestep.
  #' @param ninit     object of type Numeric. Number of initial points.
  #' @param step      object of type Numeric. Current timestep.
  #' @param end       object of type Numeric. Maximum number of iterations.
  #' @param kpars     object of type List. GP Kernel parameters.
  #' @param plt       object of type Character. Specify whether or not 
  #'                    plotting should be done during runtime.
  #' @param colour    object of type Character. Colour palette for plots.
  #' @param pmin      object of type Numeric. Lower bound for par values.
  #' @param pmax      object of type Numeric. Upper bound for par values.
  #' -----------------------------------------------------------------------
  #' @return          object of type List. Data frame of hps and their outputs 
  #'                  on f as well as a list of surface states collected 
  #'                  during runtime.
  #' -----------------------------------------------------------------------
  
  ## Sequential Design (Optimization) Loop
  for(i in (nrow(X)+1):end) {
    
    ## Using L-BFGS-B to obtain the predicted global minimum param values 
    ## based on the GP Surrogate. Initializes at current minimum.
    params    <- optim(X[which.min(y),], predicted_value, lower=pmin, 
                       upper=pmax, method="L-BFGS-B", gpi=gpi)$par
    ynew      <- do.call(f, list(params))  # Run f at proposed hp values.
    
    ## Update GP Surrogate based on new information
    updateGPsep(gpi, matrix(params, nrow=1), ynew) 
    
    ## Re-estimate parameters based on new information
    if (!noise) {
      
      # Unpack kernel parameters
      da <- kpars[[1]] 
      
      # Fit Gaussian Process using MLE
      mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
      
    } else {
      
      # Unpack kernel parameters
      d <- kpars[[1]]    # Unpack prev. estimate of lengthscale 
      g <- kpars[[2]]    # Unpack prev. estimate of nugget
      
      # Fit Gaussian Process using MLE
      mleGPsep(gpi, param="both", tmin=c(d$min, g$min), 
               tmax=c(d$max, g$max), ab=c(d$ab, g$ab))
      
    }
    
    ## Update D (data) based on new evaluation of f
    X            <- rbind(X, params)
    y            <- c(y, ynew)
    rownames(X)  <- NULL # Reset row names. Assists with plotting.
    
    if (log_surfaces) {
      
      ## Create a grid over parameter space to create predictive surfaces
      xx     <- seq(pmin, pmax, length = resolution)
      XX     <- expand.grid(xx, xx)
      
      ## Log predictive surfaces at current timestep.
      prd                <- predGPsep(gpi, XX, lite = TRUE)
      prd_mean[[step]]   <- prd$mean    # Posterior mean surface 
      prd_var[[step]]    <- prd$s2      # Posterior variance mean surface
      step               <- step + 1    # Increment timestep
    }
    
    
    ## Plots current predictive surface for param space. Used for debugging
    if (plt) {
      par(mfrow = c(1, 2))
      
      # Plot posterior mean surface
      Z <- matrix(prd_mean[[step-1]], ncol=length(xx))
      image(xx, xx, Z, col=colour, ylab="", xlab="")
      contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
      title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
      
      # Plot initial points
      points(X[1:ninit,], pch=1)
      
      #' Plots the points selected through predicted improvement. 
      #' Ugly if-else is the only way I could fix a plotting bug.
      #' Will fix in the future.
      if ((ninit+1)==nrow(X)) {
        points(t(X[(ninit+1):nrow(X),]), pch=19, cex=0.7)
      } else {
        points(X[(ninit+1):nrow(X),], pch=19, cex=0.7)
      }
      
      # Plot lines showing where point was selected.
      abline(v=X[nrow(X),1], lty = 2)
      abline(h=X[nrow(X),2], lty = 2)
      
      # Plot posterior variance surface
      Z <- matrix(prd_var[[step-1]], ncol=length(xx))
      image(xx, xx, Z, col=colour, ylab="", xlab="")
      contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
      title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
      
      # Plot the initial points
      points(X[1:ninit,], pch=1)
      
      #' Plots the points selected through predicted improvement. 
      #' Ugly if-else is the only way I could fix a plotting bug.
      #' Will fix in the future.
      if ((ninit+1)==nrow(X)) {
        points(t(X[(ninit+1):nrow(X),]), pch=19, cex=0.7)
      } else {
        points(X[(ninit+1):nrow(X),], pch=19, cex=0.7)
      }
      
      # Plot lines showing where point was selected.
      abline(v=X[nrow(X),1], lty = 2)
      abline(h=X[nrow(X),2], lty = 2)
    }
  }
  
  ## Clean and return results of BO algorithm
  deleteGPsep(gpi) 
  data  <- as.matrix(cbind(X,y))
  ef    <- list(data = data, prd_mean = prd_mean, prd_var = prd_var)
  
  return(ef)
}

## Upper Confidence Bound acquisition policy (Srinivas et al., 2009)

ucb <- function(x, gpi, lambda=1) 
{
  #' Helper function for UCB acquisition strategy. Returns a weighted sum
  #' between posterior mean and variance at a specified x. Lambda parameter
  #' specifies if the algorithm should favor exploration (high lambda values)
  #' or exploitation (low lambda values). Heuristic prevents algorithm
  #' getting stuck in local minima. Lambda should be externally set.
  #' -------------------------------------------------------------------------
  #' @param x         object of type Numeric. Vector of param coordinates.
  #' @param gpi       object of type laGP. Fitted Gaussian Process.
  #' @param lambda    object of type Numeric. Level of exploration needed.
  #' -------------------------------------------------------------------------
  #' @return          object of type Numeric. Weighted sum between posterior 
  #'                  mean and variance at a specified x
  #' -------------------------------------------------------------------------
  
  p <- predGPsep(gpi, matrix(x, nrow=1), lite=TRUE)
  return(  p$mean - (lambda * p$s2) )
}

ucb_improvement <- function(f, gpi, X, y, log_surfaces, prd_mean, prd_var, 
                            step=2, end, kpars, plt=FALSE, ninit, noise=TRUE,
                            colour, pmin, pmax, lambda, resolution=200) 
{
  #' Runs Bayesian Optimization using the UCB acquisition strategy. Not 
  #' intended to be run by user, access provided through bayes_opt function. 
  #' Returns list of BO algorithm results as well as the predictive surfaces 
  #' at each timestep. Based onthe heuristic proposed by 
  #' (Srinivas et al., 2009).
  #' -----------------------------------------------------------------------
  #' @param f         object of type Function. Function to evaluate hps over.
  #' @param gpi       object of type laGP. Fitted Gaussian Process. 
  #' @param X         object of type Matrix Array. Initial param values tested.
  #' @param y         object of type Numeric. Values for f obtained at X.
  #' @param xx        object of type Numeric. Values to create hp grid.
  #' @param XX        object of type Data Frame. Prediction grid for plotting.
  #' @param prd_mean  object of type List. Posterior mean at each timestep.
  #' @param prd_var   object of type List. Posterior var. at each timestep.
  #' @param ninit     object of type Numeric. Number of initial points.
  #' @param step      object of type Numeric. Current timestep.
  #' @param end       object of type Numeric. Maximum number of iterations.
  #' @param kpars     object of type List. GP Kernel parameters.
  #' @param plt       object of type Character. Specify whether or not 
  #'                    plotting should be done during runtime.
  #' @param colour    object of type Character. Colour palette for plots.
  #' @param pmin      object of type Numeric. Lower bound for par values.
  #' @param pmax      object of type Numeric. Upper bound for par values.
  #' @param lambda    object of type Numeric. Lambda param for UCB AF.
  #' -----------------------------------------------------------------------
  #' @return          object of type List. Data frame of hps and their outputs 
  #'                  on f as well as a list of surface states collected 
  #'                  during runtime.
  #' -----------------------------------------------------------------------
  
  ## Sequential Design (Optimization) Loop
  for(i in (nrow(X)+1):end) {
    
    ## Using L-BFGS-B to obtain the next parameter set using the UCB 
    ## acquisition function. Initializes at the current minimum. 
    params   <- optim(X[which.min(y),], ucb, lower=pmin, upper=pmax, 
                      method="L-BFGS-B", gpi=gpi, lambda=lambda)$par
    ynew     <- do.call(f, list(params))  # Run f at proposed hp values.
    
    ## Update GP Surrogate based on new information
    updateGPsep(gpi, matrix(params, nrow=1), ynew) 
    
    ## Re-estimate parameters based on new information
    if (!noise) {
      
      # Unpack kernel parameters
      da <- kpars[[1]] 
      
      # Fit Gaussian Process using MLE
      mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
      
    } else {
      
      # Unpack kernel parameters
      d <- kpars[[1]]    # Unpack prev. estimate of lengthscale 
      g <- kpars[[2]]    # Unpack prev. estimate of nugget
      
      # Fit Gaussian Process using MLE
      mleGPsep(gpi, param="both", tmin=c(d$min, g$min), 
               tmax=c(d$max, g$max), ab=c(d$ab, g$ab))
      
    }
    
    ## Update D (data) based on new evaluation of f
    X            <- rbind(X, params)
    y            <- c(y, ynew)
    rownames(X)  <- NULL # Reset row names. Assists with plotting.
    
    if (log_surfaces) {
      
      ## Create a grid over parameter space to create predictive surfaces
      xx     <- seq(pmin, pmax, length = resolution)
      XX     <- expand.grid(xx, xx)
      
      ## Log predictive surfaces at current timestep.
      prd                <- predGPsep(gpi, XX, lite = TRUE)
      prd_mean[[step]]   <- prd$mean    # Posterior mean surface 
      prd_var[[step]]    <- prd$s2      # Posterior variance mean surface
      step               <- step + 1    # Increment timestep
    }
    
    ## Plots current predictive surface for param space. Used for debugging
    if (plt) {
      par(mfrow = c(1, 2))
      
      # Plot posterior mean surface
      Z <- matrix(prd_mean[[step-1]], ncol=length(xx))
      image(xx, xx, Z, col=colour, ylab="", xlab="")
      contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
      title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
      
      # Plot initial points
      points(X[1:ninit,], pch=1)
      
      #' Plots the points selected through predicted improvement. 
      #' Ugly if-else is the only way I could fix a plotting bug.
      #' Will fix in the future.
      if ((ninit+1)==nrow(X)) {
        points(t(X[(ninit+1):nrow(X),]), pch=19, cex=0.7)
      } else {
        points(X[(ninit+1):nrow(X),], pch=19, cex=0.7)
      }
      
      # Plot lines showing where point was selected.
      abline(v=X[nrow(X),1], lty = 2)
      abline(h=X[nrow(X),2], lty = 2)
      
      # Plot posterior variance surface
      Z <- matrix(prd_var[[step-1]], ncol=length(xx))
      image(xx, xx, Z, col=colour, ylab="", xlab="")
      contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
      title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
      
      # Plot the initial points
      points(X[1:ninit,], pch=1)
      
      #' Plots the points selected through predicted improvement. 
      #' Ugly if-else is the only way I could fix a plotting bug.
      #' Will fix in the future.
      if ((ninit+1)==nrow(X)) {
        points(t(X[(ninit+1):nrow(X),]), pch=19, cex=0.7)
      } else {
        points(X[(ninit+1):nrow(X),], pch=19, cex=0.7)
      }
      
      # Plot lines showing where point was selected.
      abline(v=X[nrow(X),1], lty = 2)
      abline(h=X[nrow(X),2], lty = 2)
    }
  }
  
  ## Clean and return results of BO algorithm
  deleteGPsep(gpi) 
  data <- as.matrix(cbind(X,y))
  ucb <- list(data = data , prd_mean = prd_mean, prd_var = prd_var)
  
  return(ucb)
}

## Expected improvement acquisition policy (Jones et al., 1998)

EI <- function(gpi, x, fmin, pred=predGPsep) 
{
  #' Using a fitted Gaussian Process calculates the expected improvement 
  #' of a given point.
  #' -------------------------------------------------------------------------
  #' @param x         object of type Numeric. Vector of param coordinates.
  #' @param gpi       object of type laGP. Fitted Gaussian Process.
  #' @param pred      object of type laGP. Posterior GP.
  #' -------------------------------------------------------------------------
  #' @return          object of type Numeric. Expected improvement value
  #' -------------------------------------------------------------------------
  
  ## Checks input is valid
  if(is.null(nrow(x))) x <- matrix(x, nrow=1) 
  
  ## Expected Improvement calculations
  p       <- pred(gpi, x, lite=TRUE)  # Posterior GP
  d       <- fmin - p$mean           
  sigma   <- sqrt(p$s2)
  dn      <- d/sigma
  ei      <- d*pnorm(dn) + sigma*dnorm(dn) 
  
  return(ei)
}

obj.EI <- function(x, fmin, gpi, pred=predGPsep) 
{
  #' Helper function for EI acquisition strategy. Wrapper function which 
  #' provides access to EI calculations for single points in the param space. 
  #' -------------------------------------------------------------------------
  #' @param x         object of type Numeric. Vector of param coordinates.
  #' @param fmin      object of type Numeric. Current minimum after n evals.
  #' @param gpi       object of type laGP. Fitted Gaussian Process.
  #' @param pred      object of type laGP. Posterior GP.
  #' -------------------------------------------------------------------------
  #' @return          object of type Numeric. Weighted sum between posterior 
  #'                  mean and variance at a specified x
  #' -------------------------------------------------------------------------
  
  return(- EI(gpi, x, fmin, pred))
}

EI.search <- function(X, y, gpi, pred=predGPsep, multi.start=5, pmin, pmax) 
{
  #' Explores parameter space of the surrogate GP using the expected
  #' improvement acquisition strategy. Starts at multiple locations and
  #' returns information on each search but the maximum should be selected. 
  #' --------------------------------------------------------------------------
  #' @param X             object of type Matrix Array. Param values already 
  #'                        tested.
  #' @param y             object of type Numeric. Values for f obtained at X.
  #' @param pred          object of type laGP. Posterior GP.
  #' @param multi.start   object of type Numeric. Number of locations to start 
  #'                        search at.
  #' @param pmin          object of type Numeric. Lower bound for par values.
  #' @param pmax          object of type Numeric. Upper bound for par values.
  #' --------------------------------------------------------------------------
  #' @return          object of type Numeric. Matrix of EI locations and values
  #' --------------------------------------------------------------------------
  
  ## Set tolerance level
  tol <- sqrt(.Machine$double.eps) # Used for numerical stability
  
  ## Initial calculations
  fmin   <- y[which.min(y)]                  # Current minimum 
  start  <- matrix(X[which.min(y),], nrow=1) # Start location
  
  ## Sets up multiple start locations for EI search
  if(multi.start > 1)
    start <- rbind(start, randomLHS(multi.start - 1, ncol(X))) 
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
  
  ## Runs EI Search starting at specified locations
  for(i in 1:nrow(start)) {
    
    # Checks if improvement is below tolerance level set
    if(EI(gpi, start[i,], fmin) <= tol) { out <- list(value=-Inf); next } 
    
    # Using L-BFGS-B search parameter spaces using the EI wrapper function.
    out <- optim(start[i,], obj.EI, method="L-BFGS-B", lower=pmin, 
                 upper=pmax, gpi=gpi, pred=pred, fmin=fmin) 
    xnew[i,] <- c(out$par, -out$value) # Save output from search
  }
  
  ## Cleanup and return best solutions found during EI search
  solns <- data.frame(cbind(start, xnew)) 
  
  ## Rename columns based on parameter space
  nms <- c()
  for (i in 1:(2*ncol(X)+1)) {
    if (i <= ncol(X)) { nms <- c(nms, paste0("s", i)) }
    if (i > ncol(X) && i!=2*ncol(X)+1) { nms <- c(nms, paste0("v",i-ncol(X))) }
    if(i == (2*ncol(X)+1)) { nms <- c(nms, "val") }
  }
  
  names(solns) <- nms 
  
  ## Get rid of values below tolerance
  solns <- solns[solns$val > tol,]
  
  return(solns)
}

expected_improvement <- function(f, gpi, X, y, log_surfaces, prd_mean, prd_var,
                                 step=2, end, kpars, plt=FALSE, ninit, 
                                 noise=TRUE, colour, pmin, pmax, resolution=200) 
{
  #' Runs Bayesian Optimization using the Expected Improvement strategy. Not 
  #' intended to be run by user, access provided through bayes_opt function. 
  #' Returns list of BO algorithm results as well as the predictive surfaces 
  #' at each timestep. Based on the EI critereon (Jones et al., 1998).
  #' -----------------------------------------------------------------------
  #' @param f         object of type Function. Function to evaluate hps over.
  #' @param gpi       object of type laGP. Fitted Gaussian Process. 
  #' @param X         object of type Matrix Array. Initial param values tested.
  #' @param y         object of type Numeric. Values for f obtained at X.
  #' @param xx        object of type Numeric. Values to create hp grid.
  #' @param XX        object of type Data Frame. Prediction grid for plotting.
  #' @param prd_mean  object of type List. Posterior mean at each timestep.
  #' @param prd_var   object of type List. Posterior var. at each timestep.
  #' @param ninit     object of type Numeric. Number of initial points.
  #' @param step      object of type Numeric. Current timestep.
  #' @param end       object of type Numeric. Maximum number of iterations.
  #' @param kpars     object of type List. GP Kernel parameters.
  #' @param plt       object of type Character. Specify whether or not 
  #'                    plotting should be done during runtime.
  #' @param colour    object of type Character. Colour palette for plots.
  #' @param pmin      object of type Numeric. Lower bound for par values.
  #' @param pmax      object of type Numeric. Upper bound for par values.
  #' -----------------------------------------------------------------------
  #' @return          object of type List. Data frame of hps and their outputs 
  #'                  on f as well as a list of surface states collected 
  #'                  during runtime.
  #' -----------------------------------------------------------------------
  
  ## Sequential Design (Optimization) Loop
  for(i in (nrow(X)+1):end) {
    
    ## Find next set of parameter values which maximizes the 
    ## expected improvement criterion
    solns     <- EI.search(X, y, gpi, pmin=pmin, pmax=pmax)  
    m         <- which.max(solns$val)         # Selects maximum
    params    <- as.matrix(solns[m,(ncol(X)+1):(2*ncol(X))])     
    ynew      <- do.call(f, list(params))     # Run f at proposed hp values.
    
    ## Update GP Surrogate based on new information
    updateGPsep(gpi, matrix(params, nrow=1), ynew) 
    
    ## Re-estimate parameters based on new information
    if (!noise) {
      
      # Unpack kernel parameters
      da <- kpars[[1]] 
      
      # Fit Gaussian Process using MLE
      mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
      
    } else {
      
      # Unpack kernel parameters
      d <- kpars[[1]]    # Unpack prev. estimate of lengthscale 
      g <- kpars[[2]]    # Unpack prev. estimate of nugget
      
      # Fit Gaussian Process using MLE
      mleGPsep(gpi, param="both", tmin=c(d$min, g$min), 
               tmax=c(d$max, g$max), ab=c(d$ab, g$ab))
      
    }
    
    ## Update D (data) based on new evaluation of f
    X            <- rbind(X, params)
    y            <- c(y, ynew)
    rownames(X)  <- NULL # Reset row names. Assists with plotting.
    
    if (log_surfaces) {
      
      ## Create a grid over parameter space to create predictive surfaces
      xx     <- seq(pmin, pmax, length = resolution)
      XX     <- expand.grid(xx, xx)
      
      ## Log predictive surfaces at current timestep.
      prd                <- predGPsep(gpi, XX, lite = TRUE)
      prd_mean[[step]]   <- prd$mean    # Posterior mean surface 
      prd_var[[step]]    <- prd$s2      # Posterior variance mean surface
      step               <- step + 1    # Increment timestep
    }
    
    ## Plots current predictive surface for param space. Used for debugging
    if (plt) {
      par(mfrow = c(1, 2))
      
      # Plot posterior mean surface
      Z <- matrix(prd_mean[[step-1]], ncol=length(xx))
      image(xx, xx, Z, col=colour, ylab="", xlab="")
      contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
      title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
      
      # Plot initial points
      points(X[1:ninit,], pch=1)
      
      #' Plots the points selected through predicted improvement. 
      #' Ugly if-else is the only way I could fix a plotting bug.
      #' Will fix in the future.
      if ((ninit+1)==nrow(X)) {
        points(t(X[(ninit+1):nrow(X),]), pch=19, cex=0.7)
      } else {
        points(X[(ninit+1):nrow(X),], pch=19, cex=0.7)
      }
      
      # Plot lines showing where point was selected.
      abline(v=X[nrow(X),1], lty = 2)
      abline(h=X[nrow(X),2], lty = 2)
      
      # Plot posterior variance surface
      Z <- matrix(prd_var[[step-1]], ncol=length(xx))
      image(xx, xx, Z, col=colour, ylab="", xlab="")
      contour(xx, xx, Z, col=colour, ylab="", xlab="", add=TRUE)
      title(xlab = expression(lambda), ylab = expression(alpha), adj = 0.5)
      
      # Plot the initial points
      points(X[1:ninit,], pch=1)
      
      #' Plots the points selected through predicted improvement. 
      #' Ugly if-else is the only way I could fix a plotting bug.
      #' Will fix in the future.
      if ((ninit+1)==nrow(X)) {
        points(t(X[(ninit+1):nrow(X),]), pch=19, cex=0.7)
      } else {
        points(X[(ninit+1):nrow(X),], pch=19, cex=0.7)
      }
      
      # Plot lines showing where point was selected.
      abline(v=X[nrow(X),1], lty = 2)
      abline(h=X[nrow(X),2], lty = 2)
    }
  }
  
  ## Clean and return results of BO algorithm
  deleteGPsep(gpi) 
  data <- as.matrix(cbind(X,y))
  ei   <- list(data = data, prd_mean = prd_mean, prd_var = prd_var)
  
  return(ei)
}

## Monte Carlo simulation

## Track best objective function value
bov <- function(y, end=length(y)) 
{
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))] 
  for(i in 2:end)
    if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1] 
  return(prog)
}

## Log posterior surfaces
log_surface_information <- function(prd_mean, prd_var, dir, i)
{
  ## Write state information to csv to reduce memory load.
  for ( j in 1:length(prd_mean))  {
    
    ## Predicted mean surface
    
    # Create the directory if it doesn't exist
    dir_name <- paste0(dir, "/mc_run", i)
    if (!dir.exists(dir_name)) { dir.create(dir_name, recursive = TRUE) }
    
    # Write the data frame to CSV
    file_name <- paste0(dir_name, "/", "mean_surface", j, ".csv")
    preds <- data.frame(prd_mean[[j]])
    colnames(preds) <- "pred"
    write.csv(preds, file = file_name)
    
    ## Predictive variance surface
    
    # Create the directory if it doesn't exist
    dir_name <- paste0(dir, "/mc_run", i)
    if (!dir.exists(dir_name)) { dir.create(dir_name, recursive = TRUE) }
    
    # Write the data frame to CSV
    file_name <- paste0(dir_name, "/", "var_surface", j, ".csv")
    preds <- data.frame(prd_var[[j]])
    colnames(preds) <- "pred"
    write.csv(preds, file = file_name)
  }
}

## Get performance data from Bayesian optimization algorithms
bake_off <- function(f, names=c("x1", "x2"), dir, num_par=2, noise=FALSE,
                     ninit=12, end=25, runs=1, mc=TRUE, pmin=0, pmax=1,
                     log_surfaces=FALSE, start)
{
  ## Keep data frame for storage of results
  data <- data.frame()
  
  ## Start Bake-off
  for (i in start:runs) {
    
    ## Initialize MC run
    cat("--------------------------------------------------------\n")
    tic_main <- Sys.time()
    cat("Starting MC simulation", i, "of", runs, ":\n")
    
    ## Bayesian Optimization
    
    data <- data.frame() # Reset data
    
    ## Starting values (ensures all are equal)
    start_X <- initialize_values(num_par, pmin, pmax, ninit)
    
    ## Start EF
    cat("EF starting ...\n")
    tic <- Sys.time()
    
    ## Predicted Improvement
    res <- bayes_opt(f=f, X=start_X, ninit = ninit, end = end, 
                     num_par = num_par, pmin = pmin, pmax = pmax, 
                     names = names, type = "ef", noise=noise,
                     log_surfaces=log_surfaces)
    
    ## Log results to main data frame
    data <- rbind(data, data.frame(
      x1    = res$data$x1,
      x2    = res$data$x2,
      y     = res$data$y,
      bov   = bov(res$data$y),
      type  = "pi",
      run   = i))
    
    if (!mc) {
      log_surface_information(prd_mean = res$prd_mean, 
                              prd_var = res$prd_var,
                              dir = paste0(dir, "/mean_critereon"), 
                              i = i)
    }
    
    toc <- round(difftime(Sys.time(), tic, units = "mins"), 2)
    cat("EF completed in", toc, "minutes.\n")
    
    
    cat("UCB starting ...\n")
    tic <- Sys.time()
    
    ## UCB. Lambda = 10
    res <- bayes_opt(f=f, X=start_X, ninit = ninit, end = end, 
                     num_par = num_par, pmin = pmin, pmax = pmax, 
                     names = names, type = "ucb", noise=noise, 
                     lambda=(pmax-pmin)*10, log_surfaces=log_surfaces)
    
    ## Log results to main data frame
    data <- rbind(data, data.frame(
      x1    = res$data$x1,
      x2    = res$data$x2,
      y     = res$data$y,
      bov   = bov(res$data$y),
      type  = "ucb",
      run   = i))
    
    if (!mc) {
      log_surface_information(prd_mean = res$prd_mean, 
                              prd_var = res$prd_var,
                              dir = paste0(dir, "/upper_cb"), 
                              i = i)
    }
    
    toc <- round(difftime(Sys.time(), tic, units = "mins"), 2)
    cat("UCB completed in", toc, "minutes.\n")
    
    
    cat("EI starting ...\n")
    tic <- Sys.time()
    
    ## Expected Improvement
    res <- bayes_opt(f=f, X=start_X, ninit = ninit, end = end, 
                     num_par = num_par, pmin = pmin, pmax = pmax, 
                     names = names, type = "ei", noise=noise,
                     log_surfaces=log_surfaces)
    
    ## Log results to main data frame
    data <- rbind(data, data.frame(
      x1    = res$data$x1,
      x2    = res$data$x2,
      y     = res$data$y,
      bov   = bov(res$data$y),
      type  = "ei",
      run   = i))
    
    if (!mc) {
      log_surface_information(prd_mean = res$prd_mean, 
                              prd_var = res$prd_var,
                              dir = paste0(dir, "/expected_improvement"), 
                              i = i)
    }
    
    toc <- round(difftime(Sys.time(), tic, units = "mins"), 2)
    cat("EI completed in", toc, "minutes.\n")
    
    ## End MC run
    toc_main <- round(difftime(Sys.time(), tic_main, units = "mins"), 2)
    time_left <- toc_main*(runs-i)
    cat("--------------------------------------------------------\n")
    cat("MC simulation", i, "successfully completed in", toc_main, "minutes.\n")
    cat("Predicted time left:", time_left, "minutes or", 
        round(time_left/60, 2), "hours.\n")
    
    if (i==1) {
      write.table(data, paste0(dir, "/results_1000.csv"), sep=",", 
                  row.names = FALSE)
    } else {
      write.table(data, paste0(dir, "/results_1000.csv"),
                  sep=",", append = TRUE, col.names = FALSE,
                  row.names = FALSE)
    }
    
  }
}
