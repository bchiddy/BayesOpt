# Analytical test functions for benchmarking optimization algorithms
# Author: Ben Chiddy
# Date: 24/10/2023
# Reference: (Surjanovic and Bingham, 2013)

## Scaled Goldstein-Price Function (Picheny et al., 2013) evaluated on [0,1]^2

goldprsc <- function(xx) 
{
  ##########################################################################
  #
  # GOLDSTEIN-PRICE FUNCTION, SCALED
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1bar <- 4*xx[1] - 2
  x2bar <- 4*xx[2] - 2
  
  fact1a <- (x1bar + x2bar + 1)^2
  fact1b <- 19 - 14*x1bar + 3*x1bar^2 - 14*x2bar + 6*x1bar*x2bar + 3*x2bar^2
  fact1 <- 1 + fact1a*fact1b
  
  fact2a <- (2*x1bar - 3*x2bar)^2
  fact2b <- 18 - 32*x1bar + 12*x1bar^2 + 48*x2bar - 36*x1bar*x2bar + 27*x2bar^2
  fact2 <- 30 + fact2a*fact2b
  
  prod <- fact1*fact2
  
  y <- (log(prod) - 8.693) / 2.427  ## Noise added
  
  return(y)
}

## Function for tracking the values tried by optim function.
goldprsc_opt <- function(xx) 
{
  ynew <- goldprsc(xx)
  x1new <- xx[1]
  x2new <- xx[2]
  xx1 <<- c(xx1, x1new)
  xx2 <<- c(xx2, x2new)
  y    <<- c(y, ynew) 
  return(ynew)
} 

## Plot function surface using plotly 
cols       <- c("darkblue", "blue", "lightblue", "white", "orange", "red")
colour     <- colorRampPalette(cols)(1000)
xx         <- seq(0, 1, length=200)
XX         <- expand.grid(xx, xx)
plot_ly(x = XX[1], y = XX[2], showscale=FALSE,
        z = matrix(as.matrix(goldprsc(XX)), ncol=length(xx)), 
        colors = colour) %>% 
  add_surface(contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  ),
  ) 

## Ackley Function (Ackley, 1987) evaluated on [-5,5]^2

ackley <- function(xx) 
{
  x <- xx[1]; y <- xx[2]
  return(-20.0 * exp(-0.2 * sqrt(0.5 * (x**2 + y**2))) - 
           exp(0.5 * (cos(2 * 22/7 * x) + cos(2 * 22/7 * y))) + exp(1) + 20)
}

## Function for tracking the values tried by optim function.
ackley_opt <- function(xx) 
{
  ynew <- ackley(xx)
  x1new <- xx[1]
  x2new <- xx[2]
  xx1 <<- c(xx1, x1new)
  xx2 <<- c(xx2, x2new)
  y <<- c(y, ynew) 
  return(ynew)
} 

## Plot function surface using plotly 
cols       <- c("darkblue", "blue", "lightblue", "white", "orange", "red")
colour     <- colorRampPalette(cols)(1000)
xx         <- seq(-5, 5, length=200)
XX         <- expand.grid(xx, xx)
plot_ly(x = XX[1], y = XX[2], showscale=FALSE,
        z = matrix(as.matrix(ackley(XX)), ncol=length(xx)), 
        colors = colour) %>% 
  add_surface(contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  ),
  ) 

# Bukin Number 6 Function (Silagadze, 2007) evaluated on [-5,5]^2

bukin6 <- function(xx) 
{
  ##########################################################################
  #
  # BUKIN FUNCTION N. 6
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- 100 * sqrt(abs(x2 - 0.01*x1^2))
  term2 <- 0.01 * abs(x1+10)
  
  y <- term1 + term2
  return(y)
}

## Function for tracking the values tried by optim function.
bukin6_opt <- function(xx) 
{
  ynew <- bukin6(xx)
  x1new <- xx[1]
  x2new <- xx[2]
  xx1 <<- c(xx1, x1new)
  xx2 <<- c(xx2, x2new)
  y <<- c(y, ynew) 
  return(ynew)
} 

## Plot function surface using plotly 
cols       <- c("darkblue", "blue", "lightblue", "white", "orange", "red")
colour     <- colorRampPalette(cols)(1000)
xx         <- seq(-5, 5, length=200)
XX         <- expand.grid(xx, xx)
plot_ly(x = XX[1], y = XX[2], showscale=FALSE,
        z = matrix(as.matrix(bukin6(XX)), ncol=length(xx)), 
        colors = colour) %>% 
  add_surface(contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  ),
  ) 
