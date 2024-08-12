####### This file contains the row-wise Kronecker product and the corresponding aggregation into matrices, if convenient
######### Author: Zheyuan Li and Arthur Nicolaus Fendrich

## `sm` is the result of `sm <- SmoothCon(...)[[1]]`
## `s.margin` is the ID for spatial margins
## `newdata` is a list (without NA) giving essential covariate values
FastPredictMat2 <- function (sm, newdata, qrQ) {
  ## deal with spatial margins
  n.s <- length(sm$margin)
  Xs <- vector("list", n.s)
  i <- 1
  while (i <= n.s) {
    Xs[[i]] <- Predict.matrix2(sm$margin[[i]], newdata)
    i <- i + 1
  }
  Xs <- tensor.prod.model.matrix(Xs)

  ## then, deal with temporal part
  nbr <- dim(newdata[[sm$by]])[1]
  nbc <- dim(newdata[[sm$by]])[2]
  for(j in 1:nbc){
    st <- (j - 1) * nbr + 1
    en <- j * nbr

    if(j == 1) X <- newdata[[sm$by]][,j] * Xs[st:en,]
    else X <- X + newdata[[sm$by]][,j] * Xs[st:en,]
  }

  ## handles the centering constraints
  X <- X %*% qrQ[,2:ncol(X)]

  ## return X
  X
}


# This versions may even take more time but uses way less memory
FastPredictMat3 <- function (sm, newdata, qrQ) {
  ## stores the temporal vectors to reuse later in a memory-efficient way
  n.s <- length(sm$margin)
  v.s <- lapply(1:n.s, function(i){
  	newdata[[sm$margin[[i]]$term]]
  })

  nbr <- dim(newdata[[sm$by]])[1]
  nbc <- dim(newdata[[sm$by]])[2]
  for(j in 1:nbc){
    Xs <- vector("list", n.s)
    for(i in 1:n.s){
      newdata[[sm$margin[[i]]$term]] <- as.matrix(v.s[[i]][,j])
      Xs[[i]] <- Predict.matrix2(sm$margin[[i]], newdata)
    }

    if(j == 1) X <- tensor.prod.model.matrix(Xs) # in our case, the by variable will always be = 1
    else X <- X + tensor.prod.model.matrix(Xs)
    rm(Xs); gc()
  }
  
  ## handles the centering constraints
  X <- X %*% qrQ[,2:ncol(X)]

  ## return X
  X
}

