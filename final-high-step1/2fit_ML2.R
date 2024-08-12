################### Script: 2fit_ML2.R ##################################
## This script repeatedly iterates the model, searching for better solutions to the disaggregation problem. This is the estimation procedure for the mean parameter.

#####################################################
######## Step 4: iterative model estimation #########
#####################################################
## In this step, we used all previously generated results to estimate the coefficients of the disaggregation model

## We set the stop criteria: number of iterations
while(it < maxit){
	it <- it + 1

	## At each iteration, we print its number and calculate the necessary matrices (using 'out.fn')
	### At each step, everything is saved to help debug in case of problems
	print(paste0('iteration #', it))
	ll.t1 <- out.fn(b=b0)

	mtcs[[it]] <- ll.t1

	## Then, we optimize the log-likelihood function considering all the matrices generated
	### In order to improve the chance of finding the true minima, we do it multiple times
	### Here I decided to do it 10 times due to the high number of coarse pixels, but this may vary in each case
	### We then take the result with the lowest log-likelihood
	print(paste0('starting optimization: ', Sys.time()))	
	l <- lapply(1:10, function(i){
		print(i)
                set.seed(i)
		jj = rnorm(np)
		jj[1] = 10 # Attempt to generate feasible starting values

                optim(jj, fn=fitCpp, gr=grfitCpp, Qty=ll.t1[[1]], R=ll.t1[[2]],
			RtR = crossprod(ll.t1[[2]]), r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
			S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9,
			sldiagA=ll.t1[[6]], it=it, method='BFGS')
        })
	vs <- sapply(l, function(i) i$value)
	print(vs)
	l <- l[[which.min(vs)]]

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

	mtcs[[it]][[7]] <- l$par
	saveRDS(mtcs, 'output/1mtcs.rds')
	if(it > 1){
		saveRDS(mtcs[[it]], 'output/1mtcs_sol.rds')
		mtcs[[it-1]] <- NA
		gc()
	}

	## After model optimization, we can calculate:
	### 'bnew': the new updated vector of coefficients
	### 'gg': the norm of the gradient at the minimum. Ideally, this should be close to zero, indicating a successful optimization
	bnew <- hbCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]], RtR = crossprod(ll.t1[[2]]), r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
		S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9)
	gg <- crossprod(grfitCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]], RtR = crossprod(ll.t1[[2]]), r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
		S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, sldiagA=ll.t1[[6]], it=it))

	## Then, we update the variables before starting the next iteration, and print some relevant information
	n.dif <- crossprod(b0 - bnew)
	st <- l$par
	b0 <- bnew
	print(paste0('n.dif: ', round(n.dif, 4)))
	print(paste0('norm gr @ minimum: ', gg))
	print(paste0('ll: ', round(ml[it], 4)))
	print(paste0('rhos: ', paste0(round(st, 4), collapse = ' ')))
	print(paste0('s2: ', round(exp(st[1]), 8)))
	
	rm(ll.t1); gc()
}
