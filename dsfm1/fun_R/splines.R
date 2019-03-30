if(!require(R.matlab)){
  install.packages("R.matlab")
  library(R.matlab)
}
# library(orthogonalsplinebasis)
source("dsfm1/fun_R/lib_simulate_data.R")
source("dsfm1/fun_R/lib_splines.R")

input <- readMat("dsfm1/fun_R/input.mat") # Load input:
X <- X.comp(as.vector(input$ni)) # compute design matrix
out <- B_splines(X, as.vector(input$ni), as.vector(input$K.vec))  # compute splines matrix

writeMat("dsfm1/fun_R/output.mat", splines = out$splines, X = out$X)







