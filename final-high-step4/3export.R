################### Script: 4visualize.R ##################################

## This script is pretty much a copy of '2fit_ML.R', but with a modification:
### Instead of using the fine-scale model matrices to generate the coarse-scale model matrices, we multiply them by the optimized coefficients to get our predictions.
### This will generate the predictions for the 'year' specified below. It will then save several files to "output/tifs<year>"
### Therefore, after running this script, it is necessary to run the postprocessing script "output/rds2tif_year.R" to convert these files to a .tif. After that, "output/tifs<year>" can be deleted.
### The variable 'extrapol' extrapolates on non-observed L-VOD pixels if TRUE, or not if FALSE

rm(list=ls())
library(snow)
library(terra)
library(mgcv)
library(INLA)
library(disag.v8)

year <- 2010
nodes <- 2
extrapol <- TRUE
nr <- nrow(rast('../large_domain_allpixels/input/1blocks_high.tif'))
blks <- 3000

#auxiliary functions
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

sm = readRDS('../final-high-step1/output/1arable_sm1.rds')
sm_var1 = readRDS('mesh.rds')
sm_var2 <- sm[[4]]
qrQv2 <- qr.Q(attr(sm_var2[[1]], 'qrc'), complete = TRUE)

st = readRDS('../final-high-step1/output/1mtcs_sol.rds')[[7]]
rho = exp(st)
v = readRDS('output/1mtcs_sol-var.rds')
v = v[[3]]

Nn <- readRDS('../final-high-step1/output/1Nn.rds')
rtNn <- Nn[[1]]/sum(Nn[[2]])

bss <- bssize(nr, blks)

ag.m <- function(i, v){
	block <- rast('../large_domain_allpixels/input/1blocks_high.tif')
	map <- rast('../large_domain_allpixels/input/1Mapbiomas_2010.tif')
	tmf <- rast(paste0('../large_domain_allpixels/input/1TMF_extent_', year, '.tif'))

	readStart(block); readStart(map); readStart(tmf)

	bl <- readValues(block, row = bss$row[i], nrows = bss$nrow[i])
	bl <- matrix(bl, ncol = 1)
	v.map = readValues(map, row = bss$row[i], nrows = bss$nrow[i])
	v.map[v.map == 0] <- NA
	if(extrapol) bl <- matrix(0 * v.map, ncol = 1)

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

	outval <- NA * v.tmf
	xy <- xyFromCell(block, cellFromRow(block, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1)))
	tmp <- cbind(xy, bl, v.map, v.tmf)
	colnames(tmp) <- c('x1', 'x2', 'block', 'map', 'tmf')
	tmp[tmp[,'tmf'] %in% c(0, 8), 'tmf'] <- NA # Remove non-land and water bodies

	readStop(block); readStop(map); readStop(tmf)
	na <- getRowIndex_NA(as.matrix(tmp))
	if(length(na) == nrow(tmp)) return(list(NULL, NULL, NULL, NULL))

	if(length(na) > 0){
		tmp <- list(map = tmp[,'map'][-na], tmf = tmp[,'tmf'][-na],
			x1 = tmp[,'x1'][-na], x2 = tmp[,'x2'][-na], block = tmp[,'block'][-na])
	} else {
		tmp <- list(map = tmp[,'map'], tmf = tmp[,'tmf'],
			x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'])
	}

	tmpgam_var1 <- inla.spde.make.A(sm_var1, cbind(tmp$x1, tmp$x2))
	tmpgam_var2 <- mgcv:::Predict.matrix3(sm_var2[[1]], tmp)$X
	tmpgam_var2 <- tmpgam_var2 %*% qrQv2[,2:ncol(tmpgam_var2)]
	X_var <- cbind(tmpgam_var1, tmpgam_var2)
	y_var <- as.numeric(exp(X_var %*% v))
	diagAVA <- y_var * rho[1]

	# This was supposed to be a bias-correction procedure using the lower bounds derived, but it doesn't work if diag(V) != 1, so it must be corrected later
	# I am keeping this line just for compatibility with the scripts used for further analysis
	diagAVA <- diagAVA * rtNn

	if(length(na) > 0){
		outval[-na] <- diagAVA
	} else {
		outval <- diagAVA
	}
	
	return(outval)
}

cl <- makeCluster(nodes, type = 'SOCK')
clusterExport(cl = cl, list('extrapol', 'cellFromRow', 'sm_var1', 'sm_var2', 'qrQv2', 'rho', 'bss', 'ag.m', 'year', 'v', 'rtNn'))
clusterEvalQ(cl, {
	library(terra)
	library(mgcv)
	library(INLA)
	library(disag.v8)

	return(TRUE)
})

exporta <- function(){
	l = 0
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i+l, v = v), tag = i+l)
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
			sendCall(cl[[d$node]], ag.m, list(i = ni, v = v), tag = ni)
		}
		b <- d$value$tag
		print(paste0('recebido: ', b))
		saveRDS(d$value$value, paste0('output/tifs', year, '/v_', b))
		rm(d)
	}
}
exporta()

