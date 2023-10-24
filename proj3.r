# Em Belanger, s2409525
# Functions for smoothing and estimating functions using B-spline basis functions, summarizing data from the 
# estimated model, making predictions using a pspline model, and creating several plots from the pspline data



#################################################################################################################
# Data smoothing function using psplines
# Inputs:
# x, y - the data vectors to be smoothed
# k - the number of basis functions to use. If no value specified default is 20
# logsp - the ends of the interval over which to search for the smoothing parameter
# If no value specified default for logsp is (-5, 5)
# bord - the pspline order to use, if no number specified default is 3 (cubic)
# pord - the order of difference to use in the penalty, if no value specified the
# default is 2
# ngrid - the number of smoothing parameter values to try. If no value specified the
# default is 100

# Outputs:
# m - a list of class pspline defining the best fit spline smoother
# List contains:
# original x, y data, estimated coefficients, fitted values, sig2, covariance matrix, 
# knots vector used for spline creation, lambda value used for the model, 
# dk, bord, pord, rsqr, ks, n, k, gcv
#################################################################################################################


pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
  
  # Create the x matrix
  # knot spacing
  dk <- diff(range(x))/(k-bord) 
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  # Spline
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  
  # create the D matrix
  D <- diff(diag(k),differences=pord)
  
  # Do the qr decomposition of X
  qrd <- qr(X)
  
  # Create the R matrix
  R <- qr.R(qrd)
  
  # Create the Q matrix
  Q <- qr.Q(qrd)
  
  # Create the U matrix
  # eigen decompose the covariance matrix
  ec <- eigen(t(X)%*%X/(nrow(X)-1)) 
  # extract eigenvectors and values
  U <- ec$vectors;lam <- ec$values 
  
  # Create the capital lambda matrix
  L <- diag(lam)
  
  # Create the identity matrix
  I <- diag(1, nrow = dim(L)[1], ncol = dim(L)[2])
  
  # Store value of n
  n <- length(y)
  
  # Create a vector of log-lambda values to try out
  lambda_t <- seq(from = logsp[1], to = logsp[2], length.out = ngrid)
  
  # Create a vector to store GCV values in for different values of lambda 
  gcv <-  rep(0, times = ngrid)
  
  # Check all values of lambda iteratively and store values in the GCV vector
  for (i in length(lambda_t)){
    # Estimate betas using current lambda value
   betaH <- solve(R)%*%U%*%solve(I + lambda_t[i]*L)%*%t(U)%*%t(Q)%*%y
   # compute effective degrees of freedom for given lambda
   ks <- sum(diag(solve(I + lambda_t[i]*L)))
   # Computer sigma^s for the given lambda
   sigma2 <-  norm(y - X%*%betaH)^2/(n - ks)
   # Check the GCV for the given lambda
   gcv[i] <- sigma2/(n-ks)
  }
  
  # Find the minimum GCV value and the corresponding lambda value
  lambda <- exp(lambda_t[which.min(gcv)])
  
  # Compute squiggly k (effective degrees of freedom)
  ks <- sum(diag(solve((I + lambda*L))))
  
  # estimate the coefficients, fitted values, and sigma squared
  coef <- solve(R)%*%U%*%solve(I + lambda*L)%*%t(U)%*%t(Q)%*%y
  
  # Estimate the fitted values
  fitted <- t(X%*%coef)
  
  # Estimate sig^2
  sig2 <- (norm(y - fitted))^2/(n-ks)
  
  # Compute the covariance matrix
  covmat <- (solve(t(X)%*%X + lambda*t(D)%*%D))*sig2
  
  # Computer r-squared
  rsqr <- 1 - (n-1)*sig2/sum((y-mean(y))^2)
  
  # make note of the GCV used for the model
  GCV <- sigma2/(n-ks)
  
  # store important information in a list
  m <- list(coef, fitted, sig2, covmat, lambda, dk, bord, pord, rsqr, ks, n, y, x, k, GCV, knots)
  # gives names to the list
  names(m) <- c("coef", "fitted", "sig2", "covmat", "lambda", "dk", "bord", "pord", "rsqr", 
                "ks", "n", "y", "x", "k", "GCV", "knots")

  # The class of the list is "pspline"
  class(m) <- 'pspline'
  
  # Return the list
  return(m)
  
}

#################################################################################################################
## Function which summarizes some of the information stored in a pspline object
## Input:
## Takes m - a pspline object
## Output:
## silently returns the gcv, effective degrees of freedom (edf), and r-squared (r2)
## Prints the order of the function, the order of the pspline, the effective degrees of freedom,
## the residual standard deviation, r-squared, and GCV
#################################################################################################################

print.pspline <- function(m){
  cat("Order", m$bord, "p-spline with order", m$pord, "penalty\n")
  cat("Effective degrees of freedom", m$ks, "Coefficients:", m$coef, '\n')
  cat("residual std dev:", m$sig2^0.5, "r-squared:", m$rsqr, "GCV:", m$GCV)
  
  gcv <- m$gcv
  edf <- m$ks
  r2 <- m$rsqr
  
  list <- list(gcv, edf, r2)
  
  invisible(list)
}

#################################################################################################################
## Function for predicting new values from already estimated pspline model
## Inputs:
## m (an object of class pspline), x (new x values to estimate within the original data's 
## range), and se (TRUE or FALSE depending on whether you want the function to return the standard error or not 
## - default is TRUE)
## Outputs:
## the fitted values of x, and the standard errors of the estimates if se = TRUE
#################################################################################################################


predict.pspline <- function(m, x, se=TRUE){
  
  # Create new x matrix
  dk <- diff(range(x))/(m$k-m$bord) 
  Xp <- splines::splineDesign(m$knots,x,ord=m$bord+1,outer.ok=TRUE)
  
  # compute fitted values of x using m
  fit <- Xp%*%m$coef
  
  # If se == FALSE only return the fitted values
  if(se == FALSE){
    invisible(fit)
  }
  
  # If se==TRUE return the fitted values and the standard errors of the fitted values
  if(se == TRUE){
    # take the covariance matrix from m
    V <- m$covmat
    # Calculate the standard errors of the fitted values
    se <- rowSums(Xp*(Xp%*%V))^0.5
    # create list containing fitted values and standard errors
    li <- list(fit, se)
    # Gives names to elements of the list
    names(li) <- c("fit", "se")
    # Silently return the list
    invisible(li)
  }
  
}


#################################################################################################################
# A function for plotting objects of class pspline
# Inputs:
# m - an object of class pspline


# The function creates 3 plots:
# A plot of the original data, the fitted values, and 95% credible intervals for the fitted values
# A plot of the fitted values and the residuals
# A QQ-Plot of the residuals of the fitted function
# The function returns ll and ul, vectors of the lower credible interval and upper credible interval (respectively),
# as well as x, a vector of x values used to construct the plots
# The credible intervals are symmetric, and they are calculated using the t-value
#################################################################################################################

plot.pspline <- function(m){
  
  # Calculate 95% symmetric credible intervals for each value in X
  ll <- m$fitted-qt(0.025, m$ks)*((m$sig2/m$n)^0.5)
  ul <- m$fitted+qt(0.025, m$ks)*((m$sig2/m$n)^0.5)
  
  # put x-values in a vector for constructing the plots and to return later
  x <- m$x
  
  # The first plot is original data, the fitted values, and the credible interval
  # Start by plotting the original data
  plot(x, m$y, pch = 20, col ="#FFB481", ylab = "y",
       main = "Original Data Vs Fitted Values with Credible Intervals")
  # Add a line for the fitted line
  lines(x, m$fitted, type ='l', col = '#9A71CA')
  # Add a line for the lower limit of the CI
  lines(x, ll, col = '#F55EB5')
  # Add a line for the upper limit of the CI
  lines(x, ul, col = '#4B86E6')
  # Add a legend
  legend("bottomright", legend = c('original data', 'fitted function',
        'lower 95% credible limit', 'upper 95% credible limit'), lty= c(NA, 1, 1, 1), pch = c(20, NA, NA, NA),
        col = c('#FFB481', '#9A71CA', '#F55EB5', '#4B86E6'))
  
  
  ## The second plot is the fitted values and the residuals
  # Plot the fitted values with a line
  plot(x, m$fitted, type = 'l', col ='#9A71CA', ylab = "y", lwd=2,
       main = "Fitted Values Vs. Residuals")
  # Add residual points
  points(x, m$y-m$fitted, col='#00B2FF', bg='#58BD4E', pch=23)
  # Add a legend
  legend("bottomright", legend = c('fitted values', 'residuals'),pch = c(NA, 23),
         col = c("#9A71CA", "#00B2FF"), lty=c(1, NA), pt.bg=c(NA, '#58BD4E'), lwd=c(2, NA))

  # Make a qq plot of the fitted values
  qqnorm(m$fitted)
  # Add a straight line for interpretation
  qqline(m$fitted, col ='#23CBC6')


  # Put ll, ul, and x in a list
  li <- list(ll, ul, x)
  # Add names to the list
  names(li) <- c("ll", "ul", "x")
  # Silently return the list
  invisible(li)
  
}




 
 
 
 
 
 