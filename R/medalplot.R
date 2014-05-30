#' @title Medal plot
#'
#' @description Computes the radii of the inner, middle and outer disks in a `medal plot'.
#' 
#' @author Andrew Zammit-Mangion and Jonathan C. Rougier
#' 
#' @param Q sparse or dense precision matrix of the Gaussian random field \code{X}
#' @param Qo sparse or dense precision of the observation error \code{e}. Typically this is diagonal
#' @param A sparse or dense matrix mapping the observations \code{Z} to the underlying field \code{X} under the model \code{Z = AX + e}
#' @param subset a vector identifying which of the rows in \code{A} to use in constructing the medals
#' @param Sigma_part (optional) the partial or full posterior covariance matrix \code{Var(X|Z)} (if available)
#' @param xsamp (optional)  an \code{n x N} matrix containing the samples of the posterior distribution of \code{X}. \code{N} is the number of samples. 
#' @details The medal plot is a visualisation tool for large-scale latent Gaussian models of the form \code{Z = AX + e} 
#' where \code{Z} are the observations, \code{A} is a mapping matrix, \code{X} is the collection of hidden states and 
#' \code{e} is the additive noise with variance \code{T}. The medal plot reveals the relation of the uncertainty of the update on linear combinations 
#' of \code{Y = LX} (where \code{Y} is a subet of the rows of \code{A}) with respect to the observation uncertainty \code{T} and the prior uncertainty \code{Var(X)}. 
#' The medal also highlights the effect of the field's dependence structure on uncertainty reduction. The plot consists of a set of concentric disks with the following properties:
#' \itemize{
#'    \item The radius of the outer disk is proportional to either the square root of the observation variance \code{T} or the prior variance  
#'    \code{Var(X)}, whichever is the smaller (note that the lower of these two is a lower bound for the posterior variance in a pure Gaussian system).
#'    \item The colour of the outer disk is blue if the disk represents the observation uncertainty, red the prior uncertainty.
#'    \item The radius of the middle disk (in white), corresponding to \code{Var(Yi | Zi)}. i.e. the uncertainty on \code{Yi} if only \code{Zi} were present (or, alternatively, there was no dependence structure in \code{X}).
#'    \item The inner disk is proportional to the updated variance \code{Var(Y|Z)} and is gold in colour.
#'}
#' @return A data frame with five fields, \code{r1}: the radius of the outer disk, \code{r2}: the radius of the middle disk, \code{r3}: the radius of the inner disk, \code{col_outer}: the colour of the outer disk and \code{col_inner}: the colour of the inner disk.
#' @keywords medal plot, uncertainty visualisation
#' @examples 
#' ### Define the Matern function
#' Matern <- function(r=0:100,nu=3/2,var=1,kappa=0.1) {
#' K <- var/((2^(nu-1))*gamma(nu))*(kappa*abs(r))^nu*besselK(kappa*abs(r),nu=nu)
#'     diag(K) = var
#' return(K) }
#'
#'### Construct a valid covariance matrix on a grid of 20 cells
#' n <- 20
#' s <- 1:n
#' S <- Matern(as.matrix(dist(s)), nu = 3/2,var=4, kappa = 0.1)
#' Q <- chol2inv(chol(S))
#'
#' ### Observations and mapping matrix
#' sy <- 1:n
#' sigmav <- rep(1,n)
#' Qo <- diag(1/sigmav^2)
#' ny <- length(sy)
#' A <- diag(rep(1,n))
#' 
#' ### Medal plot function
#' M <- medalplot(Q=Q,Qo=Qo,A=A)
#' @references Jonathan C. Rougier, Andrew Zammit-Mangion and Nana Schoen (2014). Visualisation for large-scale Gaussian updates. \url{http://www.maths.bris.ac.uk/~MAZJCR/rougierVLSGU.pdf}
medalplot <- function(Q,Qo,A,subset=1L:nrow(A),Sigma_part=NULL,xsamp=NULL) {
  ## V1: Prior variance
  ## V2: Updated variance conditioned on observations individually
  ## V3: Updated variance conditioned on all observations
  if (!is.null(Sigma_part) & !(is.null(xsamp))) stop("Cannot specify both Sigma_part and xsamp")
  
  ## Cast matrices into correct types
  Q <-  as(Q,"dgCMatrix")
  Qo <- as(Qo,"dgCMatrix")
  A <- as(A,"dgCMatrix")
  
  if(!(is.null(Sigma_part)))
    Sigma_part <- as(Sigma_part,"dgCMatrix")
  
  
  Asub <- A[subset,]
  Qosub <- Qo[subset,]
  
  # Find the prior variance
  cholQ <- .Matrix_chol(Q)
  PP <- t(.cholsolve(t(Asub),cholQ$Lp,cholQ$P))  # A backsolve is a transposed forward solve ...
  V1 <- rowSums(PP * Asub)
  
  # Find the update conditioned only on the ith observation
  Prior_prec_diag <- 1/diag(PP %*% t(Asub)) # Just take the diagonal elements
  V2 <- 1/(Prior_prec_diag + diag(Qosub))
                          
  # Find the update based on all observations
  if (is.null(xsamp)) { # In this case we evaluate the partial matrix inverse
    if (is.null(Sigma_part)) {
      ## Do Takahashi
      Qstar <- t(Asub)%*%Qo%*%Asub + Q
      cholQ <- .Matrix_chol(Qstar)
      Sigma_part <- .Takahashi_Davis(cholQ$Lp,cholQ$P)
    }
    
    PP <- Asub %*% Sigma_part
    V3 <- rowSums(PP * Asub)
  
  }   else {
    V3 <- rep(0,nrow(Asub))
    for (i in 1:nrow(Asub)) {
      ind <- which(!(Asub[i,] == 0))
      Cov <- (xsamp[ind,] - apply(xsamp[ind,],1,mean)) %*% t(xsamp[ind,] - apply(xsamp[ind,],1,mean))/(ncol(xsamp)-1)
      V3[i] <- C_full[i,ind] %*% Cov %*% (C_full[i,ind])
    }
  }
  
  ## Find colour of outer disk
  col1 <- apply(cbind(1/diag(Qosub),V1),1,function(x) {
    minbnd <- which.min(x)
    if(minbnd == 1) {
      return("#3A3A98FF")
    } else {
      return("#832424FF")
    }})
  
  return(data.frame(r1 = sqrt(pmin(1/diag(Qosub),V1)),
                    r2 = sqrt(V2),
                    r3 = sqrt(V3),
                    col_outer = col1,
                    col_inner = "#FFD700"))
  
}


### deprecated
.spam_chol <- function(Q) {
  Qspam <- as.spam.dgCMatrix(Q)
  X  <- spam::chol(Qspam)
  P <- sparseMatrix(i=X@pivot,j=1:nrow(X),x=1)
  Lp <- as(as.dgCMatrix.spam(t(X)),"dtCMatrix")
  return(list(Lp = Lp,P=P))
}

.Matrix_chol <- function(Q) {
  symchol <- Cholesky(Q)
  j <- 1:n
  i <- symchol@perm + 1
  P <- sparseMatrix(i,j,x=rep(1,n))
  Lp <- t(chol(t(P)%*%Q%*%P)) 
  return(list(Lp = Lp, P=P))
}




.cholsolve <- function(y,Lp,P)  {
    ## Solve Qx = y
    v <- solve(Lp,t(P)%*%y)
    w <- solve(t(Lp),v)
    x <- P%*%w
  return(x)
}


.Takahashi_Davis <- function(Lp,P) {
  
  n <- nrow(Lp)
  
  
  d <- diag (Lp)
  L <- tril(Lp%*%sparseMatrix(i=1:n,j=1:n,x=1/d),-1)
  d <- d^2
  D <- sparseMatrix(i=1:n,j=1:n,x=d)
  
  #ii <- L@i + 1 # in {1,...,n}
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp) # in {1,...,n}, non-decreasing
  
  gc()
  Zpattern <- sparseMatrix(c(L@i + 1,jj,1:n),c(jj,L@i + 1,1:n))
  rm(dp,jj)
  
  gc()
  Z <- .sparseinv_wrapper(L,d,L,Zpattern)
  return(P%*%Z%*%t(P))
}

.sparseinv_wrapper <- function(L,d,U,Zpattern) {
  
  n <- nrow(L)
  Lp <- L@p
  Li <- L@i
  Lx <- L@x
  
  Up <- U@p
  Uj <- U@i
  Ux <- U@x
  
  Zpatp <- Zpattern@p
  Zpati <- Zpattern@i
  znz = Zpatp [n+1]
  
  
  X <- .C("sparseinv",as.integer(n),as.integer(Lp),as.integer(Li),as.double(Lx),as.double(d),as.integer(Up),as.integer(Uj),as.double(Ux),as.integer(Zpatp),as.integer(Zpati),result = double(znz))
  X <- X$result
  
  rm(U,L,Zpattern,Ux,Uj,Up,Lp,Li,Lx)
  Z <- sparseMatrix(p = Zpatp, i =Zpati, x = X,index1=F)
  
  return(Z)
}

.circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

