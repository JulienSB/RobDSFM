library(orthogonalsplinebasis)

B_splines  <-  function(X, ni, K_vec) {
  X <- X.comp(ni)
  K1 <- K_vec[1]; K2 <- K_vec[2]; K3 <- K_vec[3]
  K  <-  prod(K_vec)
  splines  <-  B_splines.comp(X, K1 = K_vec[1], K2 = K_vec[2], K3 = K_vec[3], ni[1],ni[2],ni[3])
  psimatrix <- signif(splines$psimatrix,digits = 6)
  return(list(splines.mat = psimatrix, X = X))
}



B_splines.comp <- function(X, K1, K2, K3, ni1, ni2, ni3, degree = 2){
  #Description:
  #The function return the 2 dimensional B_splines matrix
  # evaluated at X.
  # X must be  a matrix with 2 columns.
  #K1 is the number of knots for the 1st voxel index
  #K2 is the number of knots for the 2nd voxel index
  # degree is the degree of the polynomials: 2 is the quadratic B-spline
  
  X  <-  as.matrix(X)
  d2  <-  as.numeric(K2>1)
  d3  <-  as.numeric(K3>1)
  dimension  <-  1+d2+d3
  
  
  # first dimension:
  x1  <-  (1:ni1)/(ni1+1)
  nknots1  <-  K1-degree-1  # nknots is the number of interior points for 1st dimension.
  knots1  <-  seq(min(x1)*0.95,max(x1)*1.05, length.out =nknots1+2)[2:(nknots1+1)]
  fullknots1  <-  c(rep(min(x1)*0.95,degree+1),knots1,rep(max(x1)*1.05,degree+1))
  basis1  <-  SplineBasis(knots = fullknots1,order=degree+1)
  base1  <-  evaluate(basis1,x1)
  
  
  # second dimension:
  if(d2==1){
    x2  <-  (1:ni2)/(ni2+1)
    nknots2  <-  K2-degree-1  # nknots is the number of interior points for 1st dimension.
    knots2  <-  seq(min(x2)*0.95,max(x2)*1.05, length.out =nknots2+2)[2:(nknots2+1)]
    fullknots2 <- c(rep(min(x2)*0.95,degree+1),knots2,rep(max(x2)*1.05,degree+1))
    basis2 <- SplineBasis(knots = fullknots2,order=degree+1)
    base2 <- evaluate(basis2,x2)
  }else{basis2 <- NULL}
  
  # third dimension:
  if(d3==1){
    x3 <- (1:ni3)/(ni3+1)
    nknots3 <- K3-degree-1  # nknots is the number of interior points for 1st dimension.
    knots3 <- seq(min(x3)*0.95,max(x3)*1.05, length.out =nknots3+2)[2:(nknots3+1)]
    fullknots3 <- c(rep(min(x3)*0.95,degree+1),knots3,rep(max(x3)*1.05,degree+1))
    basis3 <- SplineBasis(knots = fullknots3,order=degree+1)
    base3 <- evaluate(basis3,x3)
  }else{basis3 <- NULL}
  
  
  # Tensor product
  if(dim(X)[2]==1){
    base <- base1
  }
  
  if(d2==1){
    base <- base2%x%base1
  }
  
  if(d3==1){
    base <- base3%x%base2%x%base1
  }
  return(list(psimatrix=t(base),basis1=basis1,basis2=basis2,basis3=basis3))
}
