################### Script: 2fit_ML2.R ##################################
## This script repeatedly iterates the model, searching for better solutions to the disaggregation problem. This is the estimation procedure for the variance parameter.
## The estimation of the variance parameter is not so simple, as the exact Hessian can not be easily obtained, and the approximations are not good.
## In order to solve this issue, a BFGS-scheme with a fixed (small) alpha is performed. This means that several iterations are needed, but they are not as slow as those for the mean parameter.

#####################################################
######## Step 4: iterative model estimation #########
#####################################################
## We define the log-likelihood functions with a term for the variance appended.
## Such a term was previously disregarded during the estimation of the mean parameter, because it was held constant.

## 'fit_at_min' is the log-likelihood function
## 'gtfit_at_min' is the gradient of the log-likelihood function wrt. the smoothing parameters ('rho'), calculated numerically for simplicity
## 'hv_at_min' outputs the new vector of parameters
fit_at_min <- function(rho, ref, v0, Sv1, Sv2, Sv3, Sv4, g0_var, H0_var, alpha){
	rho <- exp(rho)
	if(min(rho) < 1e-10 | max(rho) > 1e10) return(Inf)
	t2 <- 1/(4 * pi * rho[1] * exp(st)[1])
	Slv = t2 * (rho[1]*rho[1]*Sv1 + 2*rho[1]*Sv2 + Sv3) + rho[2] * Sv4

	g0v = g0_var + Slv %*% v0;
	H0v = 1/alpha * (H0_var + Slv)
	tryCatch({
		hatv = v0 - solve(H0v, g0v)
	}, error = function(e){
		return(Inf)
	})
	
	logLik = ref + 1/2 * t(v0) %*% Slv %*% v0 + 1/2 * t(g0v) %*% (hatv - v0) + 1/4 * t(hatv - v0) %*% H0v %*% (hatv - v0)
	logLik
}
grfit_at_min <- function(x, ...) grad(fit_at_min, x, ..., method = 'simple')
hv_at_min <- function(rho, ref, v0, Sv1, Sv2, Sv3, Sv4, g0_var, H0_var, alpha){
	rho <- exp(rho)
	t2 <- 1/(4 * pi * rho[1] * exp(st)[1])
	Slv = t2 * (rho[1]*rho[1]*Sv1 + 2*rho[1]*Sv2 + Sv3) + rho[2] * Sv4

	g0v = g0_var + Slv %*% v0;
	H0v = 1/alpha * (H0_var + Slv)
	hatv = v0 - solve(H0v, g0v)

	attr(hatv, 'H') <- H0v
	hatv
}

## We set the stop criteria
### 'n.dif' equals the norm between the vector of coefficients in two successive iterations
### low 'n.dif' values indicate that the model has converged
it0 <- it
while(it < maxit){
	it <- it + 1

	## At each iteration, we print its number and calculate the necessary matrices (using 'out.fn')
	### At each step, everything is saved to help debug in case of problems
	print(paste0('iteration #', it))
	ll.t2 <- out.fn(b=b0, v=v0)

	mtcs[[it]] <- ll.t2

	print(paste0('starting optimization: ', Sys.time()))	

	Y_ <- ll.t1[[1]]
	rho <- exp(st)
	Sl <- rho[1] * S1 + rho[2] * S2 + rho[3] * S3 + rho[4] * S4 + rho[5] * S5 + rho[6] * S6 + rho[7] * S7 + rho[8] * S8 + rho[9] * S9

	## We calculate the gradient of the log-likelihood wrt. the vector v0
	## All terms below are exact, but one term varies during the optimization procedure and must be appended later (see L.21)
	nk <- ncol(ll.t2[[2]])
        invAVA <- as.numeric(1/ll.t2[[1]])
        Vinv <- solve(crossprod(ll.t1[[2]] * sqrt(invAVA))/rho[1] + Sl, tol = 1e-30)
        g0_var <- as.numeric(t(Y_) %*% (- invAVA^2 * ll.t2[[2]] * Y_)/rho[1]) # gradient of the squared error
        g0_var <- g0_var + colSums(ll.t2[[2]] * invAVA) # gradient of the log determinant (likelihood)
        clusterExport(cl, list('rho', 'invAVA', 'Vinv', 'll.t2', 'll.t1'), envir = nenv)
        g0_var <- g0_var + parSapply(cl, 1:nk, function(j){ # gradient of the log determinant (Laplace)
                sum(t(Vinv) * (- t(ll.t1[[2]]) %*% (c(invAVA^2 * ll.t2[[2]][,j]/rho[1]) * ll.t1[[2]]) ))
        })
        g0_var <- g0_var * 1/2 # because we actually return 0.5 * logLik in fitCpp

	## We calculate an approximation to the Hessian matrix, using a BFGS-like scheme with alpha fixed
	## Again, one term varies during the optimization procedure and must be appended later (see L.22)
        if(it == (it0 + 1)){
                alpha <- 1
                dg <- 1/(0.1 + sqrt(as.numeric(crossprod(g0_var))))
                Hk <- diag(0 * g0_var + 1/dg) # non-inverted Hessian
        } else {
                yk <- (g0_var - g0_var0)
                Hk <- Hk + yk %*% t(yk)/as.numeric(crossprod(yk, sk)) - Hk %*% sk %*% t(sk) %*% t(Hk)/as.numeric(t(sk) %*% Hk %*% sk)
        }

	## We calculate the current log-likelihood value, disregarding the appended variance terms
	diagMMt <- ll.t2[[1]]
	ql <- diagMMt > 0
	diagL <- 0 * diagMMt
	diagL[ql] <- sqrt(1/diagMMt[ql])
	diagL <- as.numeric(diagL)
	y. <- diagL * ll.t1[[1]]
	X. <- diagL * ll.t1[[2]]
	X.qr <- qr(X., LAPACK = TRUE)
	Qty. <- qr.qty(X.qr, y.)
	R <- qr.R(X.qr, complete = TRUE)[,sort.list(X.qr$pivot)]
	r <- crossprod(y.) - crossprod(Qty.)
	sldiagA <- sum(log(diagMMt))
	ref <- ll_at_min(st, Qty=Qty., R=R, RtR=crossprod(R), r=r, N=length(ll.t1[[1]]), b0=b0,
		S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, sldiagA=sldiagA)
	rm(diagMMt, diagL, y., X., X.qr, Qty., R, r, sldiagA, invAVA, Vinv); gc()
	print(ref)

	## We optimize the log-likelihood wrt. the smoothing parameters ('rho'), now with the appended variance terms
	l <- lapply(1:5, function(i){
		print(i)
                set.seed(i)
		jj = rnorm(np_var)

                optim(jj, fn=fit_at_min, gr=grfit_at_min, ref=ref, v0=v0, Sv1=Sv1, Sv2=Sv2, Sv3=Sv3, Sv4=Sv4,
                	g0_var=g0_var, H0_var=Hk, alpha=alpha, method = 'BFGS')
        })
	vs <- sapply(l, function(i) i$value)
	print(vs)
	l <- l[[which.min(vs)]]

	vnew <- hv_at_min(l$par, ref=ref, v0=v0, Sv1=Sv1, Sv2=Sv2, Sv3=Sv3, Sv4=Sv4, g0_var=g0_var, H0_var=Hk, alpha=alpha)
	Hk <- attr(vnew, 'H'); attr(vnew, 'H') <- NULL
	sk <- vnew - v0

	### If the likelihood increases compared to the previous iteration, then the model is no longer improving and we can stop.
	ml[it] <- as.numeric(l$value)
	if(it > 2 & !res){
		print(ml[it])
		if(ml[it] > ml[it-1]){
			print('end')
			break
		}
	} else {
		res = FALSE
	}

	mtcs[[it]][[4]] <- l$par
	saveRDS(mtcs, 'output/1mtcs-var.rds')
	if(it > 1){
		saveRDS(mtcs[[it]], 'output/1mtcs_sol-var.rds')
		mtcs[[it-1]] <- NA
		gc()
	}

	## Then, we update the variables before starting the next iteration, and print some relevant information
	n.dif <- crossprod(v0 - vnew)
	v0 <- vnew
	g0_var0 <- g0_var
	st.var <- l$par

	print(paste0('n.dif: ', round(n.dif, 4)))
	print(paste0('rhos: ', paste0(round(st.var, 4), collapse = ' ')))
	print(paste0('ll: ', round(ml[it], 4)))

	rm(ll.t2); gc()
}
