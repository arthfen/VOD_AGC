library(raster)

args <- commandArgs(trailingOnly=TRUE)
year <- args[1]

blocks <- 3000

bssize <- function(nr, nb) {
	size <- ceiling(nr/nb)
	nb <- ceiling(nr / size)

	row <- (0:(nb-1))*size + 1
	nrows <- rep(size, length(row))
	dif <- nb * size - nr 
	nrows[length(nrows)] = nrows[length(nrows)] - dif
	return(list(row=row, nrows=nrows, n=nb))
}

n17 <- raster('../../large_domain_allpixels/input/1CCI_biomass_2010.tif')

bss <- bssize(nrow(n17), blocks)

out_mean <- raster(n17)
out_mean <- raster::writeStart(out_mean, filename = paste0('3out_mean_', year), datatype = 'FLT4S', format='GTiff', options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "PREDICTOR=2"))
out_sd <- raster(n17)
out_sd <- raster::writeStart(out_sd, filename = paste0('3out_sd_', year), datatype = 'FLT4S', format='GTiff', options = "COMPRESS=LZW")

for (i in 1:bss$n) {
	k <- readRDS(paste0('tifs', year, '/v_', i))
	print(paste0('received: ', i, ' - ag.m - ', Sys.time()))

	if(is.null(k[[1]])){
		nkp <- ncol(n17) * bss$nrows[i]
		k <- list(rep(NA, nkp), rep(NA, nkp))
	}
	out_mean <- raster::writeValues(out_mean, k[[1]], bss$row[i])
	out_sd <- raster::writeValues(out_sd, k[[2]], bss$row[i])
}
out_mean <- raster::writeStop(out_mean)
out_sd <- raster::writeStop(out_sd)

