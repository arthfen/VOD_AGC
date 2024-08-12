################### Script: 2fit_ML.R ##################################
## This script creates the reference smoothers for the mean parameter, based on a low-resolution sample of the data.

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
## In this step, we just define some necessary functions and loads the coarse variable

set.seed(1)
nr <- nrow(rast('../large_domain_allpixels/input/1blocks_high.tif'))
blks <- nr

## The bssize and cellFromRow functions (see 0minmax.R)
### 'coordinates' returns the xy coordinates of grid cells
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
coordinates <- function(object) {
	xyFromCell(object, cell=1:ncell(object))
}

#####################################################
##### Step 2: creating the reference smoothers ######
#####################################################
## In this step, we define the basis functions of the multi and univariate smoothers that will be used
### During estimation step (Step 4, 2fit_ML2.R), we will find the coefficients that multiply these basis functions in order to generate the disaggregation curves

## Loads the raster
cci <- rast('../large_domain_allpixels/input/1CCI_biomass_2010.tif')
map <- rast('../large_domain_allpixels/input/1Mapbiomas_2010.tif')
tem <- rast('../large_domain_allpixels/input/1Tmax_2010.tif')[[1]]
tmf <- rast('../large_domain_allpixels/input/1TMF_extent_2010.tif')
block <- rast('../large_domain_allpixels/input/1blocks_high.tif')
mcw <- rast('../large_domain_allpixels/input/1MCWD_2010.tif')

## Creates a dataframe 'tmp' containing the cell values of a lower resolution version of the rasters
### The smoother construction changes with the variable range and quantiles, so to have an approximation of the distribution of the original raster, we use a lower resolution
### The original raster (at the fine scale) can not be used because it contains too many pixels
tmpr <- spatSample(block, 1e4, method = 'regular', as.raster = TRUE)
tmp <- data.frame(coordinates(tmpr[[1]]), readValues(tmpr, mat = TRUE))
colnames(tmp) <- c('x1', 'x2', 'block')

ndf <- spatSample(c(cci, map, tem, tmf, mcw), 1e4, method = 'regular', as.raster = TRUE)
ndf <- data.frame(readValues(ndf, mat = TRUE))
colnames(ndf) <- c('cci', 'mapbiomas', 'tmax', 'tmf', 'mcwd') # format: tmf_<pixel value>

## Reminder - TMF:
# Original TMF: 1 = Undisturbed, 2 = Disturbed, 3 = Deforested, 4 = Regrowth, 6 = OtherCover
# Modified TMF: 1 = Undisturbed, 2 = Regrowth, 3 = Disturbed, 4 = Deforested, 5 = OtherCover, 8 = NA
ndf[ndf$tmf == 2, 'tmf'] <- 7
ndf[ndf$tmf == 5, 'tmf'] <- 8
ndf[ndf$tmf == 3, 'tmf'] <- 9
ndf[ndf$tmf == 4, 'tmf'] <- 2
ndf[ndf$tmf == 7, 'tmf'] <- 3
ndf[ndf$tmf == 9, 'tmf'] <- 4
ndf[ndf$tmf == 6, 'tmf'] <- 5

na <- (ndf$mapbiomas == 0 | ndf$tmf %in% c(0, 8))

mcwd = data.frame(t0 = log(1 + abs(ndf$mcwd)))
mcwd$t3 <- mcwd$t2 <- mcwd$t1 <- mcwd$t0
tmax = data.frame(t0 = ndf$tmax)
tmax$t3 <- tmax$t2 <- tmax$t1 <- tmax$t0
L <- 1 + mcwd * 0
tmp <- data.frame(x1 = tmp$x1, x2 = tmp$x2, block = tmp$block,
	tmf = ndf$tmf,
	map = ndf$mapbiomas,
	cci = ndf$cci,
	age = sample(0:10, nrow(ndf), replace = TRUE))
tmp$tmax <- as.matrix(tmax)
tmp$mcwd <- as.matrix(mcwd)
tmp$L = as.matrix(L)
tmp <- tmp[!na,]
tmp <- na.omit(tmp)

## Here, we build the reference uni and multivariate smoothers to be used. We ask them to absorb the centering constraint because we add an intercept term later.
### The values set in file 1call.R are being used to define the dimensions.
sm1 <- smoothCon(te(x1, x2, bs = c(h.sp, h.sp), m = h.m, k = c(h.b0, h.b1)), data = tmp, absorb.cons = TRUE); sm1[[1]]$X <- NULL
sm2 <- smoothCon(s(cci, bs = h.sp, m = h.m, k = h.b2), data = tmp, absorb.cons = TRUE); sm2[[1]]$X <- NULL
sm3 <- smoothCon(te(tmax, mcwd, by = L, bs = c(h.sp, h.sp), m = h.m, k = c(h.b3, h.b4)), data = tmp, absorb.cons = TRUE); sm3[[1]]$X <- NULL
sm4 <- smoothCon(s(tmf, bs = h.sp, m = h.m, k = h.b5), data = tmp, absorb.cons = TRUE); sm4[[1]]$X <- NULL
sm5 <- smoothCon(s(age, bs = h.sp, m = h.m, k = h.b6), data = tmp, absorb.cons = TRUE); sm5[[1]]$X <- NULL # age: undisturbed (sm5), regrowth (sm6), degraded (sm7)

saveRDS(list(sm1, sm2, sm3, sm4, sm5), 'output/1arable_sm1.rds')
sm1[[1]]$S <- sm2[[1]]$S <- sm3[[1]]$S <- sm4[[1]]$S <- sm5[[1]]$S <- NULL
for(j in 1:length(sm1[[1]]$margin)) sm1[[1]]$margin[[j]]$S <- sm1[[1]]$margin[[j]]$X <- NULL
rm(tmp)

