data_simulated <- function(T = 50, L = 1, ni = c(9,9,1), theta = c(0.5,0,0,0.4)){
  # This function simulate a dynamic factor model.
  # (ni1,ni2,ni3) is the voxel index.
  # ni1 is always required.
  # the dimension is equal to the number of filled voxel indexes.
  # theta is the vector of the coefficients of the VAR(1).
  # the length must be L^2.
  # Mu is the expectation of the Z. 
  # Required packages:
  # MASS : library(MASS)
  

  # 1- Simulation of the latents and the voxels :
  #============================================
  
  J<-prod(ni)
  d2<-as.numeric(ni[2]>1)
  d3<-as.numeric(ni[3]>1)
  dimension <- 1+d2+d3
  
  data <- data_fMRI(T=T,L=L,ni1=ni[1],ni2=ni[2],ni3=ni[3],theta = theta)
  X <- data$X
  Z <- data$Z
  xi <- cbind(1,Z)

  
  # 2- Building M matrix:
  #===================
  
  M<-matrix(0,nrow=L+1,ncol=J)
  
  # 1 d imension:
  if(dimension==1){
    M[1,]<-cos(2*pi*X)
    if(L>1){M[2,]<-sin(2*pi*X)}
    if(L>2){M[3,]<-cos(4*pi*X)}
    if(L>=3){M[3,]<-sin(4*pi*X)}
  }
  # 2 dimensions
  if(dimension==2){
    M[1,]<-rep(1,J)
    if(L>=1){M[2,]<-1.34*sin(2*pi*X[,1])}
    if(L>=2){M[3,]<-8.32*((X[,1]-0.5)^2 + (X[,2]-0.5)^2-0.136)}
    if(L>=3){M[4,]<-8.17*(5*(X[,2]-0.5)^3-(X[,2]-0.5))}
  }
  # 3 dimensions
  if(dimension==3){
    M[1,]<-sin(2*pi*X[,1])
    if(L>=1){M[2,]<-9.45*((X[,2]-0.5)^2+(X[,3]-0.5)^2)-1.6}
    if(L>=2){M[3,]<-3.46*(X[,2]-0.5)}
  }
  
  # 3- Simulation of Y:
  #===================
  E<-matrix(rnorm(n=T*J,mean= 0,sd=sqrt(J/2)), nrow = T,ncol = J)  # White noise process
  Y=xi%*%M+E 
  
  return(list(Y=Y,M=M,Z=Z,X=X))
}



data_fMRI <- function(T, L, ni1, ni2 = 1, ni3 = 1, theta = rep(0, L^2), MU = rep(0,L)){
  # Simulate the latent variables and the voxel matrix.
  
  theta<-matrix(theta, nrow = L, ncol=L)
  J<-ni1*ni2*ni3
  X <- X.comp(ni)  # Simulation of X:
  
  
  # 2- Simulation of Z:
  #===================
  
  # Z is simulated as a VAR(1)
  # from a white noise process U.
  
  # White noise process:
  SIGMA <- diag(1,L)    # variance of the white noise process
  U <- mvrnorm(n=T,mu=rep(0,L),Sigma=SIGMA) # We simulate U from a multivariate normal
  
  # Z process:
  Z<-U
  for(i in 2:T){
    Z[i,]<-(diag(L)-theta)%*%MU+theta%*%Z[i-1,]+U[i,]
  }
  
  return(list(X=X,Z=Z,theta=theta))
}



X.comp  <-  function(ni){
  J <- prod(ni)   
  d2  <-  as.numeric(ni[2] > 1)
  d3  <-  as.numeric(ni[3] > 1)

  i1 <- (1:ni[1] %x% rep(1, ni[2] * ni[3])) / (ni[1] + 1)
  i2 <- (1:ni[2] %x% rep(1, ni[3])) / (ni[2] + 1)
  i3 <- (1:ni[3]) / (ni[3] + 1)
  X  <-  cbind(i1, i2, i3)    
  X  <-  X[, c(1, d2*2, d3*3)]  
}

