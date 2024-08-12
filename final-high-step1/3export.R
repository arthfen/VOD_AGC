################### Script: 3export.R ##################################

## This script is pretty much a copy of '2fit_ML.R', but with a modification:
### Instead of using the fine-scale model matrices to generate the coarse-scale model matrices, we multiply them by the optimized coefficients to get our predictions.
### This will generate the predictions for the 'year' specified below. It will then save several files to "output/tifs<year>"
### Therefore, after running this script, it is necessary to run the postprocessing script "output/rds2tif_year.R" to convert these files to a .tif. After that, "output/tifs<year>" can be deleted.
### The variable 'extrapol' extrapolates on non-observed L-VOD pixels if TRUE, or not if FALSE

## Since this follows '2fit_ML.R', I omit most comments below, keeping only the new ones.

rm(list=ls())
library(terra)
library(mgcv)
library(snow)
library(parallel)
library(Rcpp)

uncond <- FALSE
uncond_approx <- FALSE

year <- 2011
nodes <- 2
extrapol <- TRUE
nr <- nrow(rast('../large_domain_allpixels/input/1blocks_high.tif'))
blks <- 3000

source('utils/PredictMat_new.R')
library(disag.v8)
g <- function(x) exp(x)
dg <- function(x) exp(x)

cl <- makeCluster(nodes, type = 'SOCK')
clusterExport(cl = cl, list('extrapol'))
clusterEvalQ(cl, {
	library(terra)
	library(mgcv)
	library(Rcpp)

	source('utils/PredictMat_new.R')
	library(disag.v8)
	return(TRUE)
})

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

mtcs = readRDS('output/1mtcs.rds')
it = max((1:50)[sapply(1:50, function(i) !any(is.na(mtcs[[i]])))])

## Loads the two results from the estimation procedure
st = mtcs[[it]][[7]]
b0 = mtcs[[it]][[5]]

cts <- readRDS('output/1gen_mtc-mu.rds')
ll.t2 <- as.numeric(readRDS('../final-high-step4/output/1mtcs_sol-var.rds')[[1]])
Y_ <- cts[[1]]
X_ <- cts[[2]]
n <- cts[[3]]
diagL <- sqrt(1/ll.t2)
y. <- diagL * Y_
X. <- diagL * X_
X.qr <- qr(X., LAPACK = TRUE)
Qty. <- qr.qty(X.qr, y.)
R <- qr.R(X.qr, complete = TRUE)[,sort.list(X.qr$pivot)]

rho <- exp(st)
Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4 + rho[6]*S5 + rho[7]*S6 + rho[8]*S7 + rho[9]*S8 + rho[10]*S9
V <- crossprod(R)/rho[1] + Sl
br <- b0

set.seed(1)
eigenV <- eigen(V)
sqrtinv <- sqrt(1/eigenV$values)
hat.cov <- crossprod(sqrtinv * t(eigenV$vectors))

## This allows the possibility of getting unconditional standard deviations. Currently not in use.
if(uncond){
	if(uncond_approx){
		# First-order correction for unconditional s.d. by Kass and Steffey (1989), described in Wood (2018)
		library(numDeriv)
		J <- jacobian(hbCpp_exp, exp(st), Qty = Qty, R = R, RtR = crossprod(R), r = r, N = N, b0 = mtcs[[it-1]][[5]],
			S1 = S1, S2 = S2, S3 = S3, S4 = S4, S5 = S5, S6 = S6, S7 = S7, S8 = S8, S9 = S9)
		Vp <- hessian(fitCpp_exp, exp(st), Qty = Qty, R = R, RtR = crossprod(R), r = r, N = N, b0 = mtcs[[it-1]][[5]],
			S1 = S1, S2 = S2, S3 = S3, S4 = S4, S5 = S5, S6 = S6, S7 = S7, S8 = S8, S9 = S9, it = it)
		eigenV.p <- eigen(Vp)
		sqrtinv.p <- sqrt(1/eigenV.p$values)
		hat.cov.p <- crossprod(sqrtinv.p * t(eigenV.p$vectors))
		ext <- J %*% hat.cov.p %*% t(J)
		
		hat.cov <- hat.cov + ext
	}
}

Cv <- chol(hat.cov, pivot = TRUE)
Cv <- Cv[,order(attr(Cv, 'pivot'))]
n.rep = 500 # Number of simulations, must be increased for more precise confidence intervals
nb <- dim(hat.cov)[1]
br <- t(Cv) %*% matrix(rnorm(n.rep*nb), nb, n.rep)
br <- as.numeric(b0) + br

bss <- bssize(nr, blks)

ag.m <- function(i, b){
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

	cci <- rast('../large_domain_allpixels/input/1CCI_biomass_2010.tif')
	block <- rast('../large_domain_allpixels/input/1blocks_high.tif')
	map <- rast('../large_domain_allpixels/input/1Mapbiomas_2010.tif')
	tmf <- rast(paste0('../large_domain_allpixels/input/1TMF_extent_', year, '.tif'))

	tem0 <- rast(paste0('../large_domain_allpixels/input/1Tmax_', year-0, '.tif'))[[1]]
	tem1 <- rast(paste0('../large_domain_allpixels/input/1Tmax_', year-1, '.tif'))[[1]]
	tem2 <- rast(paste0('../large_domain_allpixels/input/1Tmax_', year-2, '.tif'))[[1]]
	tem3 <- rast(paste0('../large_domain_allpixels/input/1Tmax_', year-3, '.tif'))[[1]]
	mcw0 <- rast(paste0('../large_domain_allpixels/input/1MCWD_', year-0, '.tif'))
	mcw1 <- rast(paste0('../large_domain_allpixels/input/1MCWD_', year-1, '.tif'))
	mcw2 <- rast(paste0('../large_domain_allpixels/input/1MCWD_', year-2, '.tif'))
	mcw3 <- rast(paste0('../large_domain_allpixels/input/1MCWD_', year-3, '.tif'))

	readStart(cci); readStart(map); readStart(tmf); readStart(block)
	readStart(tem0); readStart(tem1); readStart(tem2); readStart(tem3)
	readStart(mcw0); readStart(mcw1); readStart(mcw2); readStart(mcw3)

	## The block areas
	bl <- readValues(block, row = bss$row[i], nrows = bss$nrow[i])
	bl <- matrix(bl, ncol = 1)
	v.map = readValues(map, row = bss$row[i], nrows = bss$nrow[i])
	v.map[v.map == 0] <- NA
	if(extrapol) bl <- matrix(0 * v.map, ncol = 1)
	na = is.na(bl)
	if(sum(na) == nrow(bl)) return(list(NULL, NULL, NULL, NULL))

	xy <- xyFromCell(block, cellFromRow(block, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1)))
	v.cci <- readValues(cci, row = bss$row[i], nrows = bss$nrow[i])

	## Reminder - TMF:
	# Original TMF: 1 = Undisturbed, 2 = Disturbed, 3 = Deforested, 4 = Regrowth, 6 = OtherCover
	# Modified TMF: 1 = Undisturbed, 2 = Regrowth, 3 = Disturbed, 4 = Deforested, 5 = OtherCover, 8 = NA
	v.tmf <- readValues(tmf, row = bss$row[i], nrows = bss$nrow[i])
	v.tmf[v.tmf == 2] <- 7
	v.tmf[v.tmf == 5] <- 8
	v.tmf[v.tmf == 3] <- 9
	v.tmf[v.tmf == 4] <- 2
	v.tmf[v.tmf == 7] <- 3
	v.tmf[v.tmf == 9] <- 4
	v.tmf[v.tmf == 6] <- 5

	tmax = data.frame(t0 = readValues(tem0, row = bss$row[i], nrows = bss$nrow[i]),
		t1 = readValues(tem1, row = bss$row[i], nrows = bss$nrow[i]),
		t2 = readValues(tem2, row = bss$row[i], nrows = bss$nrow[i]),
		t3 = readValues(tem3, row = bss$row[i], nrows = bss$nrow[i]))
	mcwd = data.frame(t0 = readValues(mcw0, row = bss$row[i], nrows = bss$nrow[i]),
		t1 = readValues(mcw1, row = bss$row[i], nrows = bss$nrow[i]),
		t2 = readValues(mcw2, row = bss$row[i], nrows = bss$nrow[i]),
		t3 = readValues(mcw3, row = bss$row[i], nrows = bss$nrow[i]))
	mcwd <- log(1 + abs(mcwd))
	L <- 1 + tmax * 0

	outval <- NA * v.tmf

	y <- as.numeric(substr(year, 3, 4))
	t.age <- lapply(10:y, function(j){
		if(sum(na) == length(bl)) return(NULL)

		assign(paste0('tmf', j), rast(paste0('../large_domain_allpixels/input/1TMF_extent_20', j, '.tif')))
		readStart(get(paste0('tmf', j)))
		v.tmf <- readValues(get(paste0('tmf', j)), row = bss$row[i], nrows = bss$nrow[i], mat = TRUE)
		readStop(get(paste0('tmf', j)))
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

		v.tmf
	})
	t.age <- do.call(cbind, t.age)
	t.age <- lapply(1:ncol(t.age), function(j){
		if(j == 1) return(rep(0, nrow(t.age)))
		else return(age(t.age[,1:j]))
	})
	t.age <- t.age[[length(t.age)]]

	tmp <- cbind(xy, bl, v.map, v.tmf, v.cci)
	colnames(tmp) <- c('x1', 'x2', 'block', 'map', 'tmf', 'cci')
	tmp[tmp[,'tmf'] %in% c(0, 8), 'tmf'] <- NA # Remove non-land and water bodies

	readStop(cci); readStop(map); readStop(tmf); readStop(block)
	readStop(tem0); readStop(tem1); readStop(tem2); readStop(tem3)
	readStop(mcw0); readStop(mcw1); readStop(mcw2); readStop(mcw3)

	na <- getRowIndex_NA(as.matrix(cbind(tmp, tmax, mcwd)))
	if(length(na) == nrow(tmp)) return(list(NULL, NULL, NULL, NULL))

	if(length(na) > 0){
		tmp <- list(map = tmp[,'map'][-na], tmf = tmp[,'tmf'][-na], cci = tmp[,'cci'][-na],
			x1 = tmp[,'x1'][-na], x2 = tmp[,'x2'][-na], block = tmp[,'block'][-na], age = t.age[-na],
			tmax = as.matrix(tmax[-na,]), mcwd = as.matrix(mcwd[-na,]), L = as.matrix(L[-na,]))
	} else {
		tmp <- list(map = tmp[,'map'], tmf = tmp[,'tmf'], cci = tmp[,'cci'],
			x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'], age = t.age,
			tmax = as.matrix(tmax), mcwd = as.matrix(mcwd), L = as.matrix(L))
	}

	# The ordering of groups in the aggregation function
	groups <- sort(unique(tmp$block))

	tmpgam.X1 <- mgcv:::Predict.matrix3(sm1[[1]], tmp)$X
	tmpgam.X1 <- tmpgam.X1 %*% qrQ1[,2:ncol(tmpgam.X1)]

	tmpgam.X2 <- mgcv:::Predict.matrix3(sm2[[1]], tmp)$X
	tmpgam.X2 <- tmpgam.X2 %*% qrQ2[,2:ncol(tmpgam.X2)]

	tmpgam.X3 <- FastPredictMat3(sm3[[1]], tmp, qrQ3)

	tmpgam.X4 <- mgcv:::Predict.matrix3(sm4[[1]], tmp)$X
	tmpgam.X4 <- tmpgam.X4 %*% qrQ4[,2:ncol(tmpgam.X4)]

	tmpgam.age <- mgcv:::Predict.matrix3(sm5[[1]], tmp)$X
	tmpgam.age <- tmpgam.age %*% qrQ5[,2:ncol(tmpgam.age)]

	tmpgam.X5 <- tmpgam.age
	wk <- which(tmp$tmf != 1) # Undisturbed
	if(length(wk) > 0) tmpgam.X5[wk,] <- 0

	tmpgam.X6 <- tmpgam.age
	wk <- which(tmp$tmf != 2) # Regrowth
	if(length(wk) > 0) tmpgam.X6[wk,] <- 0

	tmpgam.X7 <- tmpgam.age
	wk <- which(tmp$tmf != 3) # Disturbed
	if(length(wk) > 0) tmpgam.X7[wk,] <- 0

	tmpgam <- cbind(1, tmpgam.X1, tmpgam.X2, tmpgam.X3, tmpgam.X4, tmpgam.X5, tmpgam.X6, tmpgam.X7)
	rm(tmpgam.X1, tmpgam.X2, tmpgam.X3, tmpgam.X4, tmpgam.X5, tmpgam.X6, tmpgam.X7)

	outval <- list(outval, outval)
	z.med = tpmm_vec_noblock(tmpgam, b0)
	z.sd = tpmm_vec_q_noblock(tmpgam, b, type = 4)

	if(length(na) > 0){
		outval[[1]][-na] <- z.med
		outval[[2]][-na] <- z.sd
	} else {
		outval[[1]] <- z.med
		outval[[2]] <- z.sd
	}
	rm(tmp)

	outval[[1]] <- exp(outval[[1]])
	outval[[2]] <- outval[[2]]
	# Reminder: what we have is y ~ lognormal(logmean = log(outval[[1]]), logsd = outval[[2]])

	return(outval)
}

clusterExport(cl=cl, list('sm1', 'sm2', 'sm3', 'sm4', 'sm5', 'bss', 'cellFromRow', 'ag.m', 'br', 'b0', 'g', 'year', 'br'))
clusterEvalQ(cl, {
	qrQ1 <- qr.Q(attr(sm1[[1]], 'qrc'), complete = TRUE)
	qrQ2 <- qr.Q(attr(sm2[[1]], 'qrc'), complete = TRUE)
	qrQ3 <- qr.Q(attr(sm3[[1]], 'qrc'), complete = TRUE)
	qrQ4 <- qr.Q(attr(sm4[[1]], 'qrc'), complete = TRUE)
	qrQ5 <- qr.Q(attr(sm5[[1]], 'qrc'), complete = TRUE)

	return(TRUE)
})

exporta <- function(){
	l = 0
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i+l, b = br), tag = i+l)
	}

	for (i in 1:bss$n) {
		d <- recvOneData(cl)
		if (!d$value$success) {
			saveRDS(d, 'erro.rds')
			cat('erro no numero: ', d$value$tag, '\n'); flush.console();
			stop('cluster error')
		}

		ni <- nodes + i + l
		if (ni <= bss$n){
			sendCall(cl[[d$node]], ag.m, list(i = ni, b = br), tag = ni)
		}
		b <- d$value$tag
		print(paste0('recebido: ', b))
		saveRDS(d$value$value, paste0('output/tifs', year, '/v_', b))
		rm(d)
	}
}
exporta()

