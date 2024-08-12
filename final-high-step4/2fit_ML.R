################### Script: 2fit_ML.R ##################################
## This script defines the aggregation function from the fine to the coarse scale ('ag.m'), which will be called repeatedly in every iteration.
## This is very similar to 'final-high-step1', but concerns the model for the variance instead of the mean.

## The following lines load the necessary libraries
rm(list=ls())

library(terra)
library(INLA)
library(mgcv)
library(snow)
library(Rcpp)
library(numDeriv)

#####################################################
######## Step 1: setting up the environment #########
#####################################################
## In this step, we just define some necessary functions
np_var <- 2
ll.t1 <- readRDS('../final-high-step1/output/1gen_mtc-mu.rds')

## Loads the target coarse variables Y to a single dataframe 'df.Y'
for(j in 10:20){
	k = readRDS(paste0('../large_domain_allpixels/input/1df_20', j, '.rds'))
	## Below is a trick to treat block data separately in time by assigning a unique ID. ID bl = 2010, ID bl+2e4 = 2011, ID bl + 2*2e4 = 2012 etc
	k[,1] <- k[,1] + (j-10) * 2e4
	assign(paste0('df.Y', j), k)
}
df.Y <- rbind(df.Y10, df.Y11, df.Y12, df.Y13, df.Y14, df.Y15, df.Y16, df.Y17, df.Y18, df.Y19, df.Y20)

## We define the number of rows and blocks that we will split our processing
nr <- nrow(rast('../large_domain_allpixels/input/1blocks_high.tif'))
blks <- 2000

## We load the likelihood package, and set the link function and its derivative. We also load a faster version of the row-wise Kronecker product.
source('utils/PredictMat_new.R')
library(disag.v8)
g <- function(x) exp(x)
dg <- function(x) exp(x)

## These two functions, respectively: i) help to subdivide the raster in blocks, ii) allow us to calculate the cell numbers for a given row.
bssize <- function(nr, nb) {
	size <- ceiling(nr/nb)
	nb <- ceiling(nr / size)

	row <- (0:(nb-1))*size + 1
	nrows <- rep(size, length(row))
	dif <- nb * size - nr 
	nrows[length(nrows)] = nrows[length(nrows)] - dif
	return(list(row=row, nrows=nrows, n=nb))
}
cellFromRow <- function(object, rownr) {
	cols <- rep(1:ncol(object), times=length(rownr))
	rows <- rep(rownr, each=ncol(object))	
	cellFromRowCol(object, rows, cols)
}

## This is an alternative function to accelerate INLA's creation of model matrix
## Such a function may not be needed after fmesher v0.1.6 - https://github.com/inlabru-org/fmesher/issues/14
inla.spde.make.A_f <- function(obj, tmp, k){
	n <- length(tmp$x1)
	n = round(seq(1, n, length.out = k), 0)
	n[1] <- 0
	n[k] <- length(tmp$x1)
	l = lapply(1:(k-1), function(j){
		v <- inla.spde.make.A(obj, cbind(tmp$x1[(n[j] + 1):n[j+1]], tmp$x2[(n[j] + 1):n[j+1]]))
		v
	})
	l = do.call(rbind, l)
	l
}

## We load the reference smoothers for the mean.
sm = readRDS('../final-high-step1/output/1arable_sm1.rds')
sm1 = sm[[1]]; sm2 = sm[[2]]; sm3 = sm[[3]]; sm4 = sm[[4]]; sm5 = sm[[5]]
bd <- sapply(1:4, function(i) get(paste0('sm', i))[[1]]$df)
bd <- c(bd, rep(sm5[[1]]$df, 3))
bd <- 1 + sum(bd) # 1 = intercept
S1 <- S2 <- S3 <- S4 <- S5 <- S6 <- S7 <- S8 <- S9 <- matrix(0, bd, bd); idx <- 1
S1[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[1]]
S2[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[2]]; idx <- idx + sm1[[1]]$df
S3[idx + (1:sm2[[1]]$df), idx + (1:sm2[[1]]$df)] <- sm2[[1]]$S[[1]]; idx <- idx + sm2[[1]]$df
S4[idx + (1:sm3[[1]]$df), idx + (1:sm3[[1]]$df)] <- sm3[[1]]$S[[1]]
S5[idx + (1:sm3[[1]]$df), idx + (1:sm3[[1]]$df)] <- sm3[[1]]$S[[2]]; idx <- idx + sm3[[1]]$df
S6[idx + (1:sm4[[1]]$df), idx + (1:sm4[[1]]$df)] <- sm4[[1]]$S[[1]]; idx <- idx + sm4[[1]]$df
S7[idx + (1:sm5[[1]]$df), idx + (1:sm5[[1]]$df)] <- sm5[[1]]$S[[1]]; idx <- idx + sm5[[1]]$df
S8[idx + (1:sm5[[1]]$df), idx + (1:sm5[[1]]$df)] <- sm5[[1]]$S[[1]]; idx <- idx + sm5[[1]]$df
S9[idx + (1:sm5[[1]]$df), idx + (1:sm5[[1]]$df)] <- sm5[[1]]$S[[1]]

# We load the reference smoothers for the variance, and construct the corresponding penalty matrices.
sm_var1 = readRDS('mesh.rds')
sm_var2 <- sm4
qrQv2 <- qr.Q(attr(sm_var2[[1]], 'qrc'), complete = TRUE)

bd_var <- sm_var1$n + sm_var2[[1]]$df
Sv1 <- Sv2 <- Sv3 <- Sv4 <- matrix(0, bd_var, bd_var); idx <- 0
Sv1[idx + (1:sm_var1$n), idx + (1:sm_var1$n)] <- as.matrix(inla.mesh.fem(sm_var1)$c1)
Sv2[idx + (1:sm_var1$n), idx + (1:sm_var1$n)] <- as.matrix(inla.mesh.fem(sm_var1)$g1)
Sv3[idx + (1:sm_var1$n), idx + (1:sm_var1$n)] <- as.matrix(inla.mesh.fem(sm_var1)$g2); idx <- idx + sm_var1$n
Sv4[idx + (1:sm_var2[[1]]$df), idx + (1:sm_var2[[1]]$df)] <- sm_var2[[1]]$S[[1]]

## Splits the raster in chunks
bss <- bssize(nr, blks)
print(paste0('bss: ', bss$n))
d <- bss$n

## 'ag.m' is the function that calculates a chunk of the aggregated matrices based on two arguments (i = raster line; b = vector of coefficients)
### The initial vector of coefficients is a random guess, which is iteratively updated
ag.m <- function(i, b, v){
	print(paste0('received, i = ', i, ' : ', Sys.time()))

	cci <- rast('../large_domain_allpixels/input/1CCI_biomass_2010.tif')
	block <- rast('../large_domain_allpixels/input/1blocks_high.tif')
	map <- rast('../large_domain_allpixels/input/1Mapbiomas_2010.tif')
	for(j in 10:20){
		y = as.numeric(paste0('20', j))
		assign(paste0('tmf', j), rast(paste0('../large_domain_allpixels/input/1TMF_extent_20', j, '.tif')))
	}
	readStart(cci)
	readStart(block)
	readStart(map)
	for(j in 10:20){
		readStart(get(paste0('tmf', j)))
	}

	## Then, a dataframe containing their values is created
	bl <- readValues(block, row = bss$row[i], nrows = bss$nrow[i])
	xy <- xyFromCell(block, cellFromRow(block, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1)))
	v.map <- readValues(map, row = bss$row[i], nrows = bss$nrow[i])
	na = v.map == 0 | is.na(bl)

	tmp <- lapply(10:20, function(j){
		if(sum(na) == length(bl)) return(NULL)

		v.tmf <- readValues(get(paste0('tmf', j)), row = bss$row[i], nrows = bss$nrow[i], mat = TRUE)
		## Reminder - TMF:
		# Original TMF: 1 = Undisturbed, 2 = Disturbed, 3 = Deforested, 4 = Regrowth, 6 = OtherCover
		# Modified TMF: 1 = Undisturbed, 2 = Regrowth, 3 = Disturbed, 4 = Deforested, 5 = OtherCover, 8 = NA
		v.tmf[v.tmf == 2] <- 7
		v.tmf[v.tmf == 5] <- 8
		v.tmf[v.tmf == 3] <- 9
		v.tmf[v.tmf == 4] <- 2
		v.tmf[v.tmf == 7] <- 3
		v.tmf[v.tmf == 9] <- 4
		v.tmf[v.tmf == 6] <- 5

		yr <- rep(2000 + j, length(v.tmf))
		mat <- cbind(xy, (j-10)*2e4 + bl, v.map, v.tmf, yr) # (j-10)*2e4 is a trick to have an unique id for each combination (coarse area, year)
		colnames(mat) <- c('x1', 'x2', 'block', 'map', 'tmf', 'year')
		mat
	})
	tmp <- do.call(rbind, tmp)
	if(is.null(tmp)) return(list(NULL, NULL, NULL))

	readStop(cci)
	readStop(block)
	readStop(map)
	for(j in 10:20){
		readStop(get(paste0('tmf', j)))
	}

	## The pixels that are not in MapBiomas are removed, as well as those NA pixels in TMF
	tmp[tmp[,'map'] == 0 | tmp[,'tmf'] %in% c(0, 8), 'block'] <- NA
	na <- getRowIndex_NA(as.matrix(tmp))
	if(length(na) == nrow(tmp)) return(list(NULL, NULL, NULL))

	## I convert the dataframe to a list (for practical reasons)
	if(length(na) > 0){
		tmp <- list(x1 = tmp[,'x1'][-na], x2 = tmp[,'x2'][-na], tmf = tmp[,'tmf'][-na], block = tmp[,'block'][-na])
	} else {
		tmp <- list(x1 = tmp[,'x1'], x2 = tmp[,'x2'], tmf = tmp[,'tmf'][-na], block = tmp[,'block'])
	}
	
	## We save the ordering of groups in the aggregation function, which allows us identify which row belongs to each coarse pixel
	groups <- sort(unique(tmp$block))

	## Then, we evaluate the smoothers at the pixel values extracted, and create the model matrix at the fine scale
	### Based on this model matrix, we derive the matrices at the coarse scale, that allow us to recover the coefficients while approximating the disaggregation constraints
	tmpgam_var1 <- inla.spde.make.A_f(sm_var1, tmp, k = 10)
	tmpgam_var2 <- mgcv:::Predict.matrix3(sm_var2[[1]], tmp)$X
	tmpgam_var2 <- tmpgam_var2 %*% qrQv2[,2:ncol(tmpgam_var2)]
	tmpgam_var <- cbind(tmpgam_var1, tmpgam_var2)
	v. <- exp(tmpgam_var %*% v[[1]])

	#AVA'
	v.out <- rowsum(as.numeric(v.), group = tmp$block)
	tmp_AVA <- v.out

	# A(Xv)A
	blk = fac2sparse(tmp$block)
	tmpgam_var <- as.numeric(v.) * tmpgam_var
	v.out <- blk %*% tmpgam_var
	tmp_AXvA <- v.out

	return(list(tmp_AVA, tmp_AXvA, groups))
}

## Here, the cluster is created and the variables, exported
clusterEvalQ(cl, rm(list=ls()))
clusterExport(cl, list('bss', 'sm_var1', 'sm_var2', 'qrQv2', 'ag.m', 'cellFromRow', 'inla.spde.make.A_f'), envir = nenv)
clusterApply(cl, seq_along(cl), function(i) sink(paste0('output/worker_', i)))
clusterEvalQ(cl, {
	library(terra)
	library(INLA)
	library(mgcv)
	library(Rcpp)
	library(numDeriv)
	source('utils/PredictMat_new.R')
	library(disag.v8)

	return(TRUE)
})

## This is the function to coordinate the generation of model matrices in parallel
out.fn <- function(b, v){
	## All target variables are started with zero values, and then updated after the result of each raster block is received
	b <- list(matrix(b, ncol = 1))
	v <- list(matrix(v, ncol = 1))
	Xv_ <- matrix(0, nrow = dim(df.Y)[1], ncol = bd_var)
	diagAVA <- matrix(0, nrow = dim(df.Y)[1], ncol = 1)

	## Sends the function to the workers, and updates the target variables
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i, b = b, v = v), tag = i)
	}
	for (i in 1:d){
		k <- recvOneData(cl)

		ni <- nodes + i
		if (ni <= d) {
			sendCall(cl[[k$node]], ag.m, list(i = ni, b = b, v = v), tag = ni)
		}
		k <- k$value$value

		if(i == 1 | i %% round(d/5, 0) == 0) print(paste0('received: ', i, ' - ag.m - ', Sys.time()))

		if(!is.null(k[[3]])){
			groups <- k[[3]]

			v.out <- k[[1]]
			v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
			diagAVA <- diagAVA + v.out

			v.out <- k[[2]]
			v.out <- rbind(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1),]
			Xv_ <- Xv_ + as.matrix(v.out)

			rm(v.out, k)
		}
	}
	print(paste0('received: all - ag.m - ', Sys.time()))
	b <- b[[1]]
	v <- v[[1]]

	## After all coarse-scale model matrices are received, we have to derive yet another auxiliary matrices
	### Again, it is necessary to refer to the disaggregation paper's Section 2.2.3 

	## Any linear operation between the fine-scale and the coarse-scale data is supported (sum, mean, weighted sum etc.)
	### In the AGC case, the following block of lines is necessary because the coarse-scale values are the average of the fine-scale ones
	n <- readRDS('../final-high-step1/output/1Nn.rds')[[2]]
	diagAVA <- diagAVA * (1/n)^2
	Xv_ <- Xv_ * (1/n)^2

	return(list(diagAVA, Xv_, v))
}

## Here, we set some constants before effectivelly running the model:
### 'maxit': the maximum number of iterations
#### In practice I'm not using 'maxit' anymore since the model is converging with relatively few iterations, but it can be useful in more complex cases
### 'b0': the initial guess of model coefficients. I modify the first element so that feasiblity is more likely
### 'st': the initial guess of optimization parameters
maxit <- 50
ml <- rep(NA, maxit)
n.dif <- 1
res = FALSE
mtcs <- rep(list(NA), maxit)
it <- 0
b0 = readRDS('../final-high-step1/output/1mtcs_sol.rds')[[5]]
st = readRDS('../final-high-step1/output/1mtcs_sol.rds')[[7]]

if(file.exists('output/1mtcs_sol-var.rds')){
	v0 <- as.numeric(readRDS('output/1mtcs_sol-var.rds')[[3]])
} else {
	v0 <- rep(0, bd_var)
}

st.var <- rep(0, np_var)

source('2fit_ML2.R', local = nenv)

