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

out_var <- raster(n17)
out_var <- raster::writeStart(out_var, filename = paste0('3out_var_', year), datatype = 'FLT4S', format='GTiff', options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "PREDICTOR=2"))

for (i in 1:bss$n) {
	k <- readRDS(paste0('tifs', year, '/v_', i))
	print(paste0('received: ', i, ' - ag.m - ', Sys.time()))

	if(class(k) == 'list'){
		nkp <- ncol(n17) * bss$nrows[i]
		k <- rep(NA, nkp)
	}
	out_var <- raster::writeValues(out_var, k, bss$row[i])
}
out_var <- raster::writeStop(out_var)

