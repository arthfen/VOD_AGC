################### Script: 4visualize.R ##################################

## In this script, we generate the curves extracted by the estimation procedure
### I didn't add almost any detailed comments because it basically replicates the same ideas described before
#### Exception: lines 77-89

rm(list=ls())
library(raster)
library(mgcv)

mask <- raster('../large_domain_allpixels/input/1Mapbiomas_2010.tif'); mask = sampleRegular(mask, 1e5, asRaster = TRUE, useGDAL = TRUE)
xy <- coordinates(mask)
nk <- nrow(mask)

g <- function(x) exp(x)
dg <- function(x) exp(x)

## Builds a reference smoother, but stores the basis indexes on the variable 'dims'
sm = readRDS('output/1arable_sm1.rds')
sm1 = sm[[1]]; sm2 = sm[[2]]; sm3 = sm[[3]]; sm4 = sm[[4]]; sm5 = sm[[5]]
bd <- sapply(1:4, function(i) get(paste0('sm', i))[[1]]$df)
bd <- c(bd, rep(sm5[[1]]$df, 3))
dims <- 1 + cumsum(bd)
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

mtcs = readRDS('output/1mtcs_sol.rds')
st = mtcs[[7]]
b0=mtcs[[5]]
Qty=mtcs[[1]]
R=mtcs[[2]]
r=mtcs[[3]]
N=mtcs[[4]]

rho <- exp(st)
Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4 + rho[6]*S5 + rho[7]*S6 + rho[8]*S7 + rho[9]*S8 + rho[10]*S9
V <- crossprod(R)/rho[1] + Sl
eigenV <- eigen(V)
sqrtinv <- sqrt(1/eigenV$values)
hat.cov <- crossprod(sqrtinv * t(eigenV$vectors))
hat.beta <- b0

# generate a fake model matrix
cci <- raster('../large_domain_allpixels/input/1CCI_biomass_2010.tif'); cci = sampleRegular(cci, 1e5, asRaster = TRUE, useGDAL = TRUE)
tem <- raster('../large_domain_allpixels/input/1Tmax_2010.tif'); tem = sampleRegular(tem, 1e5, asRaster = TRUE, useGDAL = TRUE)
mcw <- raster('../large_domain_allpixels/input/1MCWD_2010.tif'); mcw = sampleRegular(mcw, 1e5, asRaster = TRUE, useGDAL = TRUE)

# We calculate age in order to be able to plot the 'rug'
## Reminder - TMF:
# Original TMF: 1 = Undisturbed, 2 = Disturbed, 3 = Deforested, 4 = Regrowth, 6 = OtherCover
# Modified TMF: 1 = Undisturbed, 2 = Regrowth, 3 = Disturbed, 4 = Deforested, 5 = OtherCover, 8 = NA
l <- lapply(10:20, function(i){
	tmf <- raster(paste0('../large_domain_allpixels/input/1TMF_extent_20', i, '.tif')); tmf = sampleRegular(tmf, 1e5, asRaster = TRUE, useGDAL = TRUE)
	tmf <- getValues(tmf)
	tmf[tmf == 2] <- 7
	tmf[tmf == 5] <- 8
	tmf[tmf == 3] <- 9
	tmf[tmf == 4] <- 2
	tmf[tmf == 7] <- 3
	tmf[tmf == 9] <- 4
	tmf[tmf == 6] <- 5
	return(tmf)
})
l <- do.call(cbind, l)
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
l.age <- lapply(1:ncol(l), function(j){
	if(j == 1) return(rep(0, nrow(l)))
	else return(age(l[,1:j]))
})
l.age <- do.call(c, l.age)
l = c(l)

set.seed(1)
idk = sample(1:length(l), length(mcw[]))
l = l[idk]
l.age = l.age[idk]

mcwd <- data.frame(t0 = log(1 + abs(mcw[])))
mcwd$t3 <- mcwd$t2 <- mcwd$t1 <- mcwd$t0
tmax <- data.frame(t0 = tem[])
tmax$t3 <- tmax$t2 <- tmax$t1 <- tmax$t0
L <- 1 + mcwd * 0
for(i in 2:ncol(L)) L[,i] <- 0
tmp <- data.frame(x1 = xy[,1],
	x2 = xy[,2],
	cci = cci[],
	tmf = l,
	mask = mask[],
	age = l.age)
tmp$tmax <- as.matrix(tmax)
tmp$mcwd <- as.matrix(mcwd)
tmp$L <- as.matrix(L)
tmp$mask[tmp$mask == 0] <- NA

na <- which(is.na(tmp), arr.ind = TRUE)[,1]
tmp <- na.omit(tmp)
tmpgam.X1 <- PredictMat(sm1[[1]], tmp)
tmpgam.X2 <- PredictMat(sm2[[1]], tmp)
tmpgam.X3 <- PredictMat(sm3[[1]], tmp)
tmpgam.X4 <- PredictMat(sm4[[1]], tmp)
tmpgam.age <- PredictMat(sm5[[1]], tmp)

X <- cbind(1, tmpgam.X1, tmpgam.X2, tmpgam.X3, tmpgam.X4, tmpgam.age, tmpgam.age, tmpgam.age)
rm(tmpgam.X1, tmpgam.X2, tmpgam.X3, tmpgam.X4, tmpgam.age)

## Here, we sample from the posterior covariance matrix
### In order to understand this, please see the disaggregation paper's Equation 10.
### This generate multiple realization of the estimated curves, considering their estimated uncertainty
### This is what allow us to generate the standard deviation and quantiles of model predictions
set.seed(1)
Cv <- chol(hat.cov, pivot = TRUE)
Cv <- Cv[,order(attr(Cv, 'pivot'))]
n.rep = 2500
nb <- dim(hat.cov)[1]
br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep)
br <- as.numeric(hat.beta) + br

## From here on, we just export the curves to a PDF format
# s(x1, x2)
out <- X[,2:dims[1]] %*% br[2:dims[1],]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
spk.mean <- raster(mask)
v <- rep(NA, length(spk.mean[]))
v[-na] <- out.mean
spk.mean[] <- v
spk.sd <- raster(mask)
v <- rep(NA, length(spk.mean[]))
v[-na] <- out.sd
spk.sd[] <- v
writeRaster(spk.mean, paste0('1rout_s(x1,x2)-mean.tif'), overwrite = TRUE)
writeRaster(spk.sd, paste0('1rout_s(x1,x2)-sd.tif'), overwrite = TRUE)

# Remove non-mapped land cover classes
X <- X[tmp$tmf != 0,]
tmp <- tmp[tmp$tmf != 0,]

#s(x3)
tmp$x3 <- tmp$cci
out <- X[,(dims[1]+1):(dims[2])] %*% br[(dims[1]+1):(dims[2]),]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
om <- out.mean[order(tmp$x3)]
os <- out.sd[order(tmp$x3)]
d <- tmp$x3[order(tmp$x3)]
qs <- quantile(d, c(0.05, 0.95, 0.99))
om <- om[which(d < qs[3])]
os <- os[which(d < qs[3])]
d <- d[which(d < qs[3])]
pdf(paste0('s2-AGC.pdf'))
plot(om ~ d, type = 'l', cex = 1.5,
	ylim = c(min(om - os)*1, max(om + os)*1),
	xlim = c(min(d), max(d)), main = 's(x3) estimated - CCI',
	xlab = 'AGB (t/ha)',
	ylab = 's(AGB)')
polygon(c(d, rev(d)), c(om - os, rev(om + os)), col = rgb(80, 80, 80, 50, maxColorValue = 255), border = NA)
abline(h = 0, lty = 'dashed', lwd = 1.75, col = 'black')
abline(v = qs[1], lty = 'dotted', lwd = 1.35, col = 'black')
abline(v = qs[2], lty = 'dotted', lwd = 1.35, col = 'black')
lines(om ~ d, col = 'red', lwd = 2.5)
#histogram instead of rug - https://github.com/densitymodelling/alternate-rugplots
#rug(jitter(sample(tmp$x3, min(250, length(tmp$x3))), amount = 0.01), lwd = 0.8, col = adjustcolor("grey80", alpha.f=0.4))
oo <- par()
u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v <- c(v[1], v[2], v[3], v[3] + (v[4]+v[3])*0.1)
oo <- par(fig=v, new=TRUE, mar=c(0,0,0,0))
hist(tmp$x3, axes=FALSE, xlab="", ylab="", main="", breaks=100, col = adjustcolor("grey30", alpha.f=0.5), border=NA)
par(oo)
dev.off()

#s(x6)
Xk <- X # backup
tmpk <- tmp # backup
#
X <- X[tmp$tmf != 8,]
tmp <- tmp[tmp$tmf != 8,]
X <- X[!duplicated(tmp$tmf),]
tmp <- tmp[!duplicated(tmp$tmf),]
#
tmp$x6 <- tmp$tmf
out <- X[,(dims[3]+1):(dims[4])] %*% br[(dims[3]+1):(dims[4]),]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
p <- tmp$x6
om <- out.mean[order(p)]
os <- out.sd[order(p)]
d <- p[order(p)]
pdf('s4-TMF.pdf')
plot(om ~ d,
	ylim=range(c(om - os, om + os)),
	pch=19, xlab='x6', ylab='s(x6)',
	main = 's(x6) estimated - TMF', col = 'red', xaxt = 'n')
arrows(d, om - os, d, om + os, length = 0.05, angle = 90, code = 3)
abline(h = 0, lty = 'dashed', lwd = 1.75, col = 'black')
axis(1, at=1:5, labels = c('Undisturbed', 'Regrowth', 'Disturbed', 'Deforested', 'OtherCover'))
#
X <- Xk # restore backup
tmp <- tmpk # restore backup
dev.off()

#s(x7_1)
Xk <- X # backup
tmpk <- tmp # backup
#
X <- X[tmp$tmf == 1,]
tmp <- tmp[tmp$tmf == 1,]
#
tmp$x7 <- tmp$age
out <- X[,(dims[4]+1):(dims[5])] %*% br[(dims[4]+1):(dims[5]),]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
om <- out.mean[order(tmp$x7)]
os <- out.sd[order(tmp$x7)]
d <- tmp$x7[order(tmp$x7)]
pdf(paste0('s5-Age.pdf'))
plot(om ~ d, type = 'l', cex = 1.5,
	ylim = c(min(om - os)*1, max(om + os)*1),
	xlim = c(min(d), max(d)), main = 's(x7) Estimated - Age 1',
	xlab = 'age 1',
	ylab = 's(age)')
polygon(c(d, rev(d)), c(om - os, rev(om + os)), col = rgb(80, 80, 80, 50, maxColorValue = 255), border = NA)
abline(h = 0, lty = 'dashed', lwd = 1.75, col = 'black')
lines(om ~ d, col = 'red', lwd = 2.5)
#histogram instead of rug - https://github.com/densitymodelling/alternate-rugplots
#rug(jitter(sample(tmp$x7, min(250, length(tmp$x7))), amount = 0.1), lwd = 0.8, col = adjustcolor("grey80", alpha.f=0.4))
oo <- par()
u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v <- c(v[1], v[2], v[3], v[3] + (v[4]+v[3])*0.1)
oo <- par(fig=v, new=TRUE, mar=c(0,0,0,0))
hist(tmp$x7, axes=FALSE, xlab="", ylab="", main="", breaks=100, col = adjustcolor("grey30", alpha.f=0.5), border=NA)
par(oo)
#
X <- Xk # restore backup
tmp <- tmpk # restore backup
dev.off()

#s(x7_2)
Xk <- X # backup
tmpk <- tmp # backup
#
X <- X[tmp$tmf == 2,]
tmp <- tmp[tmp$tmf == 2,]
#
tmp$x7 <- tmp$age
out <- X[,(dims[5]+1):(dims[6])] %*% br[(dims[5]+1):(dims[6]),]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
om <- out.mean[order(tmp$x7)]
os <- out.sd[order(tmp$x7)]
d <- tmp$x7[order(tmp$x7)]
pdf(paste0('s6-Age.pdf'))
plot(om ~ d, type = 'l', cex = 1.5,
	ylim = c(min(om - os)*1, max(om + os)*1),
	xlim = c(min(d), max(d)), main = 's(x7) Estimated - Age 1',
	xlab = 'age 1',
	ylab = 's(age)')
polygon(c(d, rev(d)), c(om - os, rev(om + os)), col = rgb(80, 80, 80, 50, maxColorValue = 255), border = NA)
abline(h = 0, lty = 'dashed', lwd = 1.75, col = 'black')
lines(om ~ d, col = 'red', lwd = 2.5)
#histogram instead of rug - https://github.com/densitymodelling/alternate-rugplots
#rug(jitter(sample(tmp$x7, min(250, length(tmp$x7))), amount = 0.1), lwd = 0.8, col = adjustcolor("grey80", alpha.f=0.4))
oo <- par()
u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v <- c(v[1], v[2], v[3], v[3] + (v[4]+v[3])*0.1)
oo <- par(fig=v, new=TRUE, mar=c(0,0,0,0))
hist(tmp$x7, axes=FALSE, xlab="", ylab="", main="", breaks=100, col = adjustcolor("grey30", alpha.f=0.5), border=NA)
par(oo)
#
X <- Xk # restore backup
tmp <- tmpk # restore backup
dev.off()

#s(x7_3)
Xk <- X # backup
tmpk <- tmp # backup
#
X <- X[tmp$tmf == 3,]
tmp <- tmp[tmp$tmf == 3,]
#
tmp$x7 <- tmp$age
out <- X[,(dims[6]+1):(dims[7])] %*% br[(dims[6]+1):(dims[7]),]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
om <- out.mean[order(tmp$x7)]
os <- out.sd[order(tmp$x7)]
d <- tmp$x7[order(tmp$x7)]
pdf(paste0('s7-Age.pdf'))
plot(om ~ d, type = 'l', cex = 1.5,
	ylim = c(min(om - os)*1, max(om + os)*1),
	xlim = c(min(d), max(d)), main = 's(x7) Estimated - Age 1',
	xlab = 'age 1',
	ylab = 's(age)')
polygon(c(d, rev(d)), c(om - os, rev(om + os)), col = rgb(80, 80, 80, 50, maxColorValue = 255), border = NA)
abline(h = 0, lty = 'dashed', lwd = 1.75, col = 'black')
lines(om ~ d, col = 'red', lwd = 2.5)
#histogram instead of rug - https://github.com/densitymodelling/alternate-rugplots
#rug(jitter(sample(tmp$x7, min(250, length(tmp$x7))), amount = 0.1), lwd = 0.8, col = adjustcolor("grey80", alpha.f=0.4))
oo <- par()
u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v <- c(v[1], v[2], v[3], v[3] + (v[4]+v[3])*0.1)
oo <- par(fig=v, new=TRUE, mar=c(0,0,0,0))
hist(tmp$x7, axes=FALSE, xlab="", ylab="", main="", breaks=100, col = adjustcolor("grey30", alpha.f = 0.5), border=NA)
par(oo)
#
X <- Xk # restore backup
tmp <- tmpk # restore backup
dev.off()


# tensor product of tmax and mcwd
library(alphahull)
db <- data.frame(x = c(tmp$mcwd), y = c(tmp$tmax))
db <- db[!duplicated(db),]
set.seed(0)
db <- db[sample(1:nrow(db), 500),]
cxh <- ahull(x = db$x, y = db$y, alpha = 3)
#pdf('Rplots1.pdf'); plot(cxh); dev.off()

set.seed(1)
v1 <- seq(min(tmp$mcwd), max(tmp$mcwd), length.out = 300) # x
v2 <- seq(min(tmp$tmax), max(tmp$tmax), length.out = 300) # y

Xk <- matrix(data = 0, nrow = length(v1), ncol = length(v2))
qrQ3 <- qr.Q(attr(sm3[[1]], 'qrc'), complete = TRUE)
for (i in 1:length(v1)) {
	xk <- data.frame(mcwd = v1[i], tmax = v2)
	Xs <- lapply(1:2, function(j){
		Predict.matrix2(sm3[[1]]$margin[[j]], xk)
	})
	Xs <- tensor.prod.model.matrix(Xs)
	Xs <- Xs %*% qrQ3[,2:ncol(Xs)]

	incxh <- 1 * inahull(cxh, as.matrix(xk))
	incxh[incxh == 0] <- NA
	
#	Xk[i,] <- incxh * apply(Xs %*% br[(dims[2]+1):(dims[3]),], 1, mean)
	Xk[i,] <- incxh * Xs %*% hat.beta[(dims[2]+1):(dims[3]),]
}
pdf('s3-tmax-mcwd.pdf')
#plot(x = db$x, y = db$y, pch = 16, bg = NA, col = adjustcolor("grey50", alpha.f = 0.4), cex = 0.6)
filled.contour(x = v1, y = v2, z = Xk, color = terrain.colors, zlim = c(-0.2, 0.2), xlab = 'log[abs( Maximum cumulative water deficit )]', ylab = 'Temperature (K)',
  plot.axes = {
    axis(1)
    axis(2)
    points(db$x, db$y, pch = 16, col = adjustcolor("grey50", alpha.f = 0.4), cex = 0.6, bg = NA)
  }
)
dev.off()

