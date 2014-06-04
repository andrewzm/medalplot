## ----Matern--------------------------------------------------------------
Matern <- function(r=0:100,nu=3/2,var=1,kappa=0.1) {
  r <- kappa * abs(r)
  robj <- r^nu * besselK(r, nu = nu) / (2^(nu - 1) * gamma(nu))
  robj[is.nan(robj)] <- 1
  var * robj
}

## ----init, fig.height=4--------------------------------------------------
n <- 99
s <- 1:n
V <- Matern(as.matrix(dist(s)), nu = 3/2,var=4, kappa = 0.1)

## ----Q-------------------------------------------------------------------
Q <- chol2inv(chol(V))

## ----obs-----------------------------------------------------------------
sy <- c(10,14,18,22,26,30,60,80,85,90)
sigmav <- c(1.6,1.6,1.6,0.8,1.6,1.6,1.5,2.2,2.2,2.2)
Qo <- diag(1/sigmav^2) # Precision matrix

## ----A-------------------------------------------------------------------
ny <- length(sy)
A <- matrix(0,ny,n)
A[cbind(seq(along = sy), sy)] <- 1

## ----medal---------------------------------------------------------------
library(medalplot)
M <- medalplot(Q=Q,Qo=Qo,A=A)

## ----M-------------------------------------------------------------------
print(M)

## ----load_lib2-----------------------------------------------------------
library(ggplot2)
library(plyr)
library(ellipse)

## ----Madd----------------------------------------------------------------
M$x <- sy
M$y <- 2.7
M$sigmav <- sigmav

## ----gmedals,echo=FALSE--------------------------------------------------
g_medals <- function(data,print_middle=T,alpha=1,scale=0.004,
                     clamp_below=F,show_rim=F,g=ggplot()) {
  
  
  .ellipseFun <- function(center=c(0,0),scale=c(1,4),npoints=100) {
    df <- data.frame(ellipse(0,scale=scale,npoints=npoints))
    df$x <- df$x + center[1]
    df$y <- df$y + center[2]
    df
  }
  
  size_outer <-  (data$r1)
  size_inner <- (data$r3)
  size_middle <- (data$r2)
  
  if (show_rim) {
    ## Enlarge outer rim
    foo <- (size_outer - size_middle) / size_outer
    size_outer <- ifelse(foo < 0.1, (1.1 - foo) * size_outer, size_outer)
  }
  
  ## If medals are too small scale them up
  if(clamp_below) {
    min_size <- min(diff(range(data$x)),diff(range(data$y)))/900
    ind <- which(size_outer < min_size)   
    scales <- min_size/size_outer
    size_outer[ind] <- min_size 
    size_inner[ind] <- (size_inner*scales)[ind]
    size_middle[ind] <- (size_middle*scales)[ind]
  }
  
  mobs <- nrow(M)
  
  ## Plot outer medals
  for (i in 1:mobs) {
    g <- g + geom_polygon(data=.ellipseFun(c(data$x[i],data$y[i]),
                                           size_outer[i]*scale,
                                           npoints=100),
                          aes(x,y),fill=data$col_outer[i],alpha=alpha)
  }
  ## Plot middle medals
  if(print_middle) for (i in 1:mobs) {
      g <- g + geom_polygon(data=.ellipseFun(c(data$x[i],data$y[i]),
                                             size_middle[i]*scale,
                                             npoints=100),
                            aes(x,y),fill="white",alpha=alpha)
  }
   
  ## Plot inner medals
  for (i in 1:mobs) {
    g <- g + geom_polygon(data=.ellipseFun(c(data$x[i],data$y[i]),
                                           size_inner[i]*scale,
                                           npoints=100),
                          aes(x,y),fill=data$col_inner[i])
  }
  return(g)
}

## ----Qstar---------------------------------------------------------------
Qstar <- Qtot <- crossprod(chol(Qo) %*% A) + Q
Sigma <- chol2inv(chol(Qstar))
x_std <- sqrt(diag(Sigma))

## ----plot2,fig.height=4,fig.cap='Medals showing the relation between the prior, posterior and observation uncertainty and the effect of vicinity of the observations on each other. The blue bars denote the observation uncertainty, the red shading the prior uncertainty and the yellow shading the posterior uncertainty. For interpretation of the medals see main text.'----
### Prior and posterior uncertainty
X <- data.frame(s=s,std = x_std,prior_std = sqrt(diag(V)))

### Bars for observation uncertainty
Obs_bars <- ddply(M,"s",function(df) { 
  X <- data.frame(s1 = c(df$x-0.5, df$x-0.5, df$x + 0.5, df$x + 0.5),
                  y1 = c(0,df$sigmav,df$sigmav,0))
  })

### Final plot without medals
g <- ggplot() +  
   geom_ribbon(data=X,aes(x=s,ymin = 0,ymax = prior_std),
               colour="black",fill="#FF525A",alpha=1) +
   geom_ribbon(data=X,aes(x=s,ymin = 0,ymax = std),
               fill="#FFEB80",colour="black",alpha=1)   +
   geom_polygon(data=Obs_bars,aes(x=s1,y=y1,group=s),
                fill="blue") +
   theme(panel.background = element_rect(fill='white',colour="black"),
         legend.position="none",
         panel.grid.major =  element_line(colour = "light gray", size = 0.05),
         text = element_text(size=20)) +
   xlab("s") + ylab("") + xlim(0,100) + ylim(0,3.2)


### Add medals to plot 
g <- g_medals(g=g,data=M,print_middle=T,alpha=1,scale=c(0.8,0.05),clamp_below=F)
print(g)

## ----gmedals2------------------------------------------------------------
g_medals <- function(data,print_middle=T,alpha=1,scale=0.004,
                     clamp_below=F,show_rim=F,g=ggplot()) {
  
  
  .ellipseFun <- function(center=c(0,0),scale=c(1,4),npoints=100) {
    df <- data.frame(ellipse(0,scale=scale,npoints=npoints))
    df$x <- df$x + center[1]
    df$y <- df$y + center[2]
    df
  }
  
  size_outer <-  (data$r1)
  size_inner <- (data$r3)
  size_middle <- (data$r2)
  
  if (show_rim) {
    ## Enlarge outer rim
    foo <- (size_outer - size_middle) / size_outer
    size_outer <- ifelse(foo < 0.1, (1.1 - foo) * size_outer, size_outer)
  }
  
  ## If medals are too small scale them up
  if(clamp_below) {
    min_size <- min(diff(range(data$x)),diff(range(data$y)))/900
    ind <- which(size_outer < min_size)   
    scales <- min_size/size_outer
    size_outer[ind] <- min_size 
    size_inner[ind] <- (size_inner*scales)[ind]
    size_middle[ind] <- (size_middle*scales)[ind]
  }
  
  mobs <- nrow(M)
  
  ## Plot outer medals
  for (i in 1:mobs) {
    g <- g + geom_polygon(data=.ellipseFun(c(data$x[i],data$y[i]),
                                           size_outer[i]*scale,
                                           npoints=100),
                          aes(x,y),fill=data$col_outer[i],alpha=alpha)
  }
  ## Plot middle medals
  if(print_middle) for (i in 1:mobs) {
      g <- g + geom_polygon(data=.ellipseFun(c(data$x[i],data$y[i]),
                                             size_middle[i]*scale,
                                             npoints=100),
                            aes(x,y),fill="white",alpha=alpha)
  }
   
  ## Plot inner medals
  for (i in 1:mobs) {
    g <- g + geom_polygon(data=.ellipseFun(c(data$x[i],data$y[i]),
                                           size_inner[i]*scale,
                                           npoints=100),
                          aes(x,y),fill=data$col_inner[i])
  }
  return(g)
}

