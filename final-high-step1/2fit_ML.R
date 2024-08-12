################### Script: 2fit_ML.R ##################################
## This script defines the aggregation function from the fine to the coarse scale ('ag.m'), which will be called repeatedly in every iteration.

## The following lines load the necessary libraries
rm(list=ls())
library(terra)
library(mgcv)
library(snow)
library(Rcpp)
library(numDeriv)

#####################################################
######## Step 1: setting up the environment #########
#####################################################
## In step 1, we define some necessary functions and loads the coarse variable

## If we had a previous iteration for the variance parameter, then we load it to incorporate in our likelihood function
if(file.exists('../final-high-step4/output/1mtcs_sol-var.rds')){
	ll.t2 <- as.numeric(readRDS('../final-high-step4/output/1mtcs_sol-var.rds')[[1]])
}

## We define the number of rows and blocks that we will split our processing
nr <- nrow(rast('../large_domain_allpixels/input/1blocks_high.tif'))
blks <- 4000

## We load the likelihood package, and set the link function and its derivative. We also load a faster version of the row-wise Kronecker product.
library(disag.v8)
g <- function(x) exp(x)
dg <- function(x) exp(x)
source('utils/PredictMat_new.R')

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

## We load the target coarse variables Y to a single dataframe 'df.Y'
for(j in 10:20){
	k = readRDS(paste0('../large_domain_allpixels/input/1df_20', j, '.rds'))
	## Below is a trick to treat block data separately in time by assigning a unique ID. ID bl = 2010, ID bl+2e4 = 2011, ID bl + 2*2e4 = 2012 etc
	k[,1] <- k[,1] + (j-10) * 2e4
	assign(paste0('df.Y', j), k)
}
df.Y <- rbind(df.Y10, df.Y11, df.Y12, df.Y13, df.Y14, df.Y15, df.Y16, df.Y17, df.Y18, df.Y19, df.Y20)

## We load the reference smoothers, previously calculated in '1smoothers.R', and fill our penalty matrices.
sm = readRDS('output/1arable_sm1.rds')
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

## Splits the raster in chunks
bss <- bssize(nr, blks)
print(paste0('bss: ', bss$n))
d <- bss$n

## 'ag.m' is the function that calculates a chunk of the aggregated matrices based on two arguments (i = raster line; b = vector of coefficients)
### The initial vector of coefficients is a random guess, which is iteratively updated
ag.m <- function(i, b){
	print(paste0('received, i = ', i, ' : ', Sys.time()))

	age <- function(A){
		A <- as.matrix(A)
		age <- rep(0, nrow(A))
		if(ncol(A) == 1) return(age)
		for(i in 2:ncol(A)){
			tmp <- A[,i] - A[,i-1]

			na <- tmp == 0
			if(sum(na) > 0) age[na] <- age[na] + 1
			
			na <- tmp != 0
			if(sum(na) > 0) age[na] <- 0
		}
		age
	}

	## This block loads all spatial files and set them to read mode
	cci <- rast('../large_domain_allpixels/input/1CCI_biomass_2010.tif')
	block <- rast('../large_domain_allpixels/input/1blocks_high.tif')
	map <- rast('../large_domain_allpixels/input/1Mapbiomas_2010.tif')
	for(j in 10:20){
		y = as.numeric(paste0('20', j))
		
		assign(paste0('tmf', j), rast(paste0('../large_domain_allpixels/input/1TMF_extent_20', j, '.tif')))
		assign(paste0('tem0', j), rast(paste0('../large_domain_allpixels/input/1Tmax_', y-0, '.tif'))[[1]])
		assign(paste0('tem1', j), rast(paste0('../large_domain_allpixels/input/1Tmax_', y-1, '.tif'))[[1]])
		assign(paste0('tem2', j), rast(paste0('../large_domain_allpixels/input/1Tmax_', y-2, '.tif'))[[1]])
		assign(paste0('tem3', j), rast(paste0('../large_domain_allpixels/input/1Tmax_', y-3, '.tif'))[[1]])
		assign(paste0('mcw0', j), rast(paste0('../large_domain_allpixels/input/1MCWD_', y-0, '.tif')))
		assign(paste0('mcw1', j), rast(paste0('../large_domain_allpixels/input/1MCWD_', y-1, '.tif')))
		assign(paste0('mcw2', j), rast(paste0('../large_domain_allpixels/input/1MCWD_', y-2, '.tif')))
		assign(paste0('mcw3', j), rast(paste0('../large_domain_allpixels/input/1MCWD_', y-3, '.tif')))
	}
	readStart(cci)
	readStart(block)
	readStart(map)
	for(j in 10:20){
		readStart(get(paste0('tmf', j)))
		readStart(get(paste0('tem0', j)))
		readStart(get(paste0('tem1', j)))
		readStart(get(paste0('tem2', j)))
		readStart(get(paste0('tem3', j)))
		readStart(get(paste0('mcw0', j)))
		readStart(get(paste0('mcw1', j)))
		readStart(get(paste0('mcw2', j)))
		readStart(get(paste0('mcw3', j)))
	}

	## We also force to close these connections after processing
	on.exit({
		readStop(cci)
		readStop(block)
		readStop(map)
		for(j in 10:20){
			readStop(get(paste0('tmf', j)))
			readStop(get(paste0('tem0', j)))
			readStop(get(paste0('tem1', j)))
			readStop(get(paste0('tem2', j)))
			readStop(get(paste0('tem3', j)))
			readStop(get(paste0('mcw0', j)))
			readStop(get(paste0('mcw1', j)))
			readStop(get(paste0('mcw2', j)))
			readStop(get(paste0('mcw3', j)))
		}
	})

	## Then, a dataframe containing their values is created
	bl <- readValues(block, row = bss$row[i], nrows = bss$nrow[i])
	xy <- xyFromCell(block, cellFromRow(block, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1)))
	v.cci <- readValues(cci, row = bss$row[i], nrows = bss$nrow[i])
	v.map <- readValues(map, row = bss$row[i], nrows = bss$nrow[i])
	na = v.map == 0 | is.na(bl)

	tmp <- lapply(10:20, function(j){
		if(sum(na) == length(bl)) return(NULL)

		v.tmf <- readValues(get(paste0('tmf', j)), row = bss$row[i], nrows = bss$nrow[i], mat = TRUE)
		## Here, the TMF dataset is reclassified according to the following rules:
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
		mat <- cbind(xy, (j-10)*2e4 + bl, v.map, v.tmf, v.cci, yr) # (j-10)*2e4 is a trick to have an unique id for each combination (coarse area, year)
		colnames(mat) <- c('x1', 'x2', 'block', 'map', 'tmf', 'cci', 'year')
		mat
	})
	tmp <- do.call(rbind, tmp)
	if(is.null(tmp)) return(list(NULL, NULL, NULL, NULL, NULL))

	## In order to calculate time since disturbance, a separate dataset containing only TMF data is created
	tmf <- lapply(10:20, function(j){
		if(sum(na) == length(bl)) return(NULL)

		v.tmf <- readValues(get(paste0('tmf', j)), row = bss$row[i], nrows = bss$nrow[i], mat = TRUE)
		# Original TMF: 1 = Undisturbed, 2 = Disturbed, 3 = Deforested, 4 = Regrowth, 6 = OtherCover
		# Modified TMF: 1 = Undisturbed, 2 = Regrowth, 3 = Disturbed, 4 = Deforested, 5 = OtherCover, 8 = NA
		v.tmf[v.tmf == 2] <- 7
		v.tmf[v.tmf == 5] <- 8
		v.tmf[v.tmf == 3] <- 9
		v.tmf[v.tmf == 4] <- 2
		v.tmf[v.tmf == 7] <- 3
		v.tmf[v.tmf == 9] <- 4
		v.tmf[v.tmf == 6] <- 5

		v.tmf
	})
	tmf <- do.call(cbind, tmf)
	t.age <- lapply(1:ncol(tmf), function(j){
		if(j == 1) return(rep(0, nrow(tmf)))
		else return(age(tmf[,1:j]))
	})
	t.age <- do.call(c, t.age)

	tmax <- lapply(10:20, function(j){
		v.tem = data.frame(t0 = readValues(get(paste0('tem0', j)), row = bss$row[i], nrows = bss$nrow[i]),
			t1 = readValues(get(paste0('tem1', j)), row = bss$row[i], nrows = bss$nrow[i]),
			t2 = readValues(get(paste0('tem2', j)), row = bss$row[i], nrows = bss$nrow[i]),
			t3 = readValues(get(paste0('tem3', j)), row = bss$row[i], nrows = bss$nrow[i]))
		v.tem
	})
	tmax <- do.call(rbind, tmax)
	L <- 1 + 0 * tmax

	mcwd <- lapply(10:20, function(j){
		v.mcw <- data.frame(t0 = readValues(get(paste0('mcw0', j)), row = bss$row[i], nrows = bss$nrow[i]),
			t1 = readValues(get(paste0('mcw1', j)), row = bss$row[i], nrows = bss$nrow[i]),
			t2 = readValues(get(paste0('mcw2', j)), row = bss$row[i], nrows = bss$nrow[i]),
			t3 = readValues(get(paste0('mcw3', j)), row = bss$row[i], nrows = bss$nrow[i]))
		v.mcw <- log(1 + abs(v.mcw))
	})
	mcwd <- do.call(rbind, mcwd)
	
	## The pixels that are not in MapBiomas are removed, as well as those NA pixels in TMF
	tmp[tmp[,'map'] == 0 | tmp[,'tmf'] %in% c(0, 8), 'block'] <- NA
	na <- getRowIndex_NA(as.matrix(cbind(tmp, tmax, mcwd)))
	if(length(na) == nrow(tmp)) return(list(NULL, NULL, NULL, NULL, NULL))

	## This part is a little bit tricky. The coarse-level average is calculated based on the AGC of forest pixels + non-forest pixels
	### Since we assume that AGC(non-forest) = 0, we have to count how many non-forest pixels we have in order to add to the calculation
	### I call this "extra zeros", reliable pixels that are removed from the model, but not from the average calculation
	s0 <- data.frame(s0 = 0, block = tmp[,'block'])
	s0[which(tmp[,'tmf'] %in% c(0, 8)), 1] <- 1
	s0 <- na.omit(s0)
	s0 <- cbind(sort(unique(s0[,2])), rowsum(s0[,1], group = s0[,2]))
	s0 <- s0[which(s0[,2] > 0),]

	## I convert the dataframe to a list (for practical reasons)
	if(length(na) > 0){
		tmp <- list(map = tmp[,'map'][-na], tmf = tmp[,'tmf'][-na], cci = tmp[,'cci'][-na],
			x1 = tmp[,'x1'][-na], x2 = tmp[,'x2'][-na], block = tmp[,'block'][-na], age = t.age[-na],
			tmax = as.matrix(tmax[-na,]), mcwd = as.matrix(mcwd[-na,]), L = as.matrix(L[-na,]))
	} else {
		tmp <- list(map = tmp[,'map'], tmf = tmp[,'tmf'], cci = tmp[,'cci'],
			x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'], age = t.age,
			tmax = as.matrix(tmax), mcwd = as.matrix(mcwd), L = as.matrix(L))
	}

	## We save the ordering of groups in the aggregation function, which allows us identify the coarse-scale pixel containing each fine-scale pixel
	groups <- sort(unique(tmp$block))

	## Then, we evaluate the smoothers at the pixel values extracted, and create the model matrix at the fine scale ('tmpgam')
	### Based on this model matrix, we derive the matrices at the coarse scale, that allow us to recover the coefficients while approximating the disaggregation constraints
	tmpgam.X1 <- mgcv:::Predict.matrix3(sm1[[1]], tmp)$X
	tmpgam.X1 <- tmpgam.X1 %*% qrQ1[,2:ncol(tmpgam.X1)]
	i_st <- 2
	i_en <- i_st + ncol(tmpgam.X1) - 1
	z. <- tmpgam.X1 %*% b[[1]][i_st:i_en]

	tmpgam.X2 <- mgcv:::Predict.matrix3(sm2[[1]], tmp)$X
	tmpgam.X2 <- tmpgam.X2 %*% qrQ2[,2:ncol(tmpgam.X2)]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X2) - 1
	z. <- z. + tmpgam.X2 %*% b[[1]][i_st:i_en]

	tmpgam.X3 <- FastPredictMat3(sm3[[1]], tmp, qrQ3)
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X3) - 1
	z. <- z. + tmpgam.X3 %*% b[[1]][i_st:i_en]

	tmpgam.X4 <- mgcv:::Predict.matrix3(sm4[[1]], tmp)$X
	tmpgam.X4 <- tmpgam.X4 %*% qrQ4[,2:ncol(tmpgam.X4)]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X4) - 1
	z. <- z. + tmpgam.X4 %*% b[[1]][i_st:i_en]

	tmpgam.age <- mgcv:::Predict.matrix3(sm5[[1]], tmp)$X
	tmpgam.age <- tmpgam.age %*% qrQ5[,2:ncol(tmpgam.age)]

	tmpgam.X5 <- tmpgam.age
	na <- which(tmp$tmf != 1) # Undisturbed
	if(length(na) > 0) tmpgam.X5[na,] <- 0
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X5) - 1
	z. <- z. + tmpgam.X5 %*% b[[1]][i_st:i_en]

	tmpgam.X6 <- tmpgam.age
	na <- which(tmp$tmf != 2) # Regrowth
	if(length(na) > 0) tmpgam.X6[na,] <- 0
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X6) - 1
	z. <- z. + tmpgam.X6 %*% b[[1]][i_st:i_en]

	tmpgam.X7 <- tmpgam.age
	na <- which(tmp$tmf != 3) # Degraded
	if(length(na) > 0) tmpgam.X7[na,] <- 0
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X7) - 1
	z. <- z. + tmpgam.X7 %*% b[[1]][i_st:i_en]

	z. <- z. + b[[1]][1]

	## From here on, we derive the matrices at the coarse scale. It is necessary to read the 2022 disaggregation paper to understand what they mean.
	### Please refer to its Section 2.2
	gz. <- g(z.)
	v.out <- rowsum(gz., group = tmp$block)
	tmp_Agz. <- v.out

	#MMt
	dg2 <- 0 * gz. + 1
	v.out <- rowsum(dg2, group = tmp$block)
	tmp_MMt <- v.out

	#MX
	tmpgam.X1 <- as.numeric(dg(z.)) * tmpgam.X1
	tmpgam.X2 <- as.numeric(dg(z.)) * tmpgam.X2
	tmpgam.X3 <- as.numeric(dg(z.)) * tmpgam.X3
	tmpgam.X4 <- as.numeric(dg(z.)) * tmpgam.X4
	tmpgam.X5 <- as.numeric(dg(z.)) * tmpgam.X5
	tmpgam.X6 <- as.numeric(dg(z.)) * tmpgam.X6
	tmpgam.X7 <- as.numeric(dg(z.)) * tmpgam.X7
	v.out <- cbind(rowsum(dg(z.), group = tmp$block),
		rowsum(tmpgam.X1, group = tmp$block),
		rowsum(tmpgam.X2, group = tmp$block),
		rowsum(tmpgam.X3, group = tmp$block),
		rowsum(tmpgam.X4, group = tmp$block),
		rowsum(tmpgam.X5, group = tmp$block),
		rowsum(tmpgam.X6, group = tmp$block),
		rowsum(tmpgam.X7, group = tmp$block))
	tmp_AGX <- v.out
	rm(tmpgam.X1, tmpgam.X2, tmpgam.X3, tmpgam.X4, tmpgam.X5, tmpgam.X6, tmpgam.X7); gc()

	return(list(tmp_Agz.[,1], tmp_AGX, tmp_MMt[,1], s0, groups))
}

## Here, the cluster is created and the variables, exported
clusterEvalQ(cl, rm(list=ls()))
clusterExport(cl, list('bss', 'sm1', 'sm2', 'sm3', 'sm4', 'sm5', 'df.Y', 'ag.m', 'g', 'dg', 'cellFromRow'), envir = nenv)
clusterApply(cl, seq_along(cl), function(i) sink(paste0('output/worker_', i)))
clusterEvalQ(cl, {
	library(terra)
	library(mgcv)
	library(Rcpp)
	library(numDeriv)
	source('utils/PredictMat_new.R')
	library(disag.v8)

	qrQ1 <- qr.Q(attr(sm1[[1]], 'qrc'), complete = TRUE)
	qrQ2 <- qr.Q(attr(sm2[[1]], 'qrc'), complete = TRUE)
	qrQ3 <- qr.Q(attr(sm3[[1]], 'qrc'), complete = TRUE)
	qrQ4 <- qr.Q(attr(sm4[[1]], 'qrc'), complete = TRUE)
	qrQ5 <- qr.Q(attr(sm5[[1]], 'qrc'), complete = TRUE)
	
	return(TRUE)
})

## While we don't calculate the analytical gradient, we evaluate it numerically
grfitCpp <- function(x, ...) grad(fitCpp, x, ..., method = 'simple')
clusterExport(cl, list('grfitCpp'), envir = nenv)

## This is the function to coordinate the generation of model matrices in parallel
out.fn <- function(b){
	## All target variables are started with zero values, and then updated after the result of each raster block is received
	b <- list(matrix(b, ncol = 1))
	Agz. <- matrix(0, nrow = dim(df.Y)[1], ncol = 2)
	Agz.[,1] <- df.Y[,1]
	X_ <- matrix(0, nrow = dim(df.Y)[1], ncol = bd)
	diagMMt <- matrix(0, nrow = dim(df.Y)[1], ncol = 2)
	diagMMt[,1] <- df.Y[,1]

	## Sends the function to the workers, and updates the target variables
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i, b = b), tag = i)
	}
	for (i in 1:d){
		k <- recvOneData(cl)

		ni <- nodes + i
		if (ni <= d) {
			sendCall(cl[[k$node]], ag.m, list(i = ni, b = b), tag = ni)
		}
		k <- k$value$value

		if(i == 1 | i %% round(d/5, 0) == 0) print(paste0('received: ', i, ' - ag.m - ', Sys.time()))

		if(!is.null(k[[5]])){
			groups <- k[[5]]

			v.out <- k[[1]]
			v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
			Agz.[,2] <- Agz.[,2] + v.out

			v.out <- k[[2]]
			v.out <- rbind(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1),]
			X_ <- X_ + as.matrix(v.out)

			v.out <- k[[3]]
			v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
			diagMMt[,2] <- diagMMt[,2] + v.out

			# Again, for the extra zeros
			v.out <- k[[4]]
			v.out <- rbind(0, v.out)[match(df.Y[,1], c(NA, v.out[,1]), nomatch = 1),]
			diagMMt[,2] <- diagMMt[,2] + v.out[,2]

			rm(v.out, k)
		}
	}
	print(paste0('received: all - ag.m - ', Sys.time()))
	b <- b[[1]]

	## After all coarse-scale model matrices are received, we have to derive yet more auxiliary matrices
	### Again, it is necessary to refer to the paper's Section 2.2.3 

	## Any linear operation between the fine-scale and the coarse-scale data is supported (sum, mean, weighted sum etc.)
	### In the AGC case, the following block of lines is necessary because the coarse-scale values are the average of the fine-scale ones
	n <- diagMMt[,2]
	diagMMt[,2] <- diagMMt[,2] * (1/n)^2
	Agz.[,2] <- Agz.[,2] * (1/n)
	X_ <- X_ * (1/n)

	## Calculates the auxiliary matrices
	Y_ <- df.Y[,2] - Agz.[,2]
	ql <- diagMMt[,2] > 0
	if(exists('ll.t2')){
		diagL <- sqrt(1/ll.t2)
	} else {
		diagL <- sqrt(1/diagMMt[,2])
		ll.t2 <- diagMMt[,2]
	}
	y. <- diagL * Y_
	X. <- diagL * X_
	X.qr <- qr(X., LAPACK = TRUE)
	Qty. <- qr.qty(X.qr, y.)
	R <- qr.R(X.qr, complete = TRUE)[,sort.list(X.qr$pivot)]

	r <- crossprod(y.) - crossprod(Qty.)
	N <- sum(ql)
	er <- crossprod(Y_)
	print(er) ## This will print the quadratic error between the true coarse values and the predicted ones. Ideally this decreases at each iteration

	saveRDS(list(N, n), 'output/1Nn.rds')
	saveRDS(list(Y_, X_, n), 'output/1gen_mtc-mu.rds')
	return(list(Qty., R, r, N, b, sum(log(ll.t2))))
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
set.seed(1)
if(file.exists('output/1mtcs_sol.rds')){
	b0 <- readRDS('output/1mtcs_sol.rds')[[5]]
	st <- readRDS('output/1mtcs_sol.rds')[[7]]
	np <- length(st)
} else {
	set.seed(1)
	b0 <- rnorm(bd, 0, 0.01); b0[1] = 4.5
	np <- 10
	st <- rep(0, np)
}
source('2fit_ML2.R', local = nenv)

