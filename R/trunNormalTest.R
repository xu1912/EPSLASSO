# This function tests the significance of all the predictors given the solution of beta, sigma from EPSLASSO.
#
# @param A Matrix of predictors. Required.
# @param Y Trait values. Required.
# @param c1 Right censored point. Required.
# @param c2 Left censored point. Required.
# @param beta_e A sequence of beta_i values estimated from epsLassoSolver. Required.
# @param sigma Sigma estimated from epsLassoSolver. Required.
# @param m_w Methods used to estimate W matrix. Default is "lso" for LASSO solution using glmnet. Another method is "dzg" for Danzig-type estimator.
# @param paral Parallel computing indicator. Default is FALSE, not using parallel.
# @param paral_n Number of cores that are used for parallel computing. Default is NULL. When paral is TRUE, default is the number of system available cores - 1.
# @param resol The refining step when m_w="dzg". Default is 1.3. A large resol results in faster convergence speed, but rough solutions.
# @param tol The convergence threshold for refining when m_w="dzg". Default is 0.001.
# @param maxTry The maximum refining steps when m_w="dzg". Default is 10.
# @param verbose Print debugging info or not.
# @import glmnet
#' @import MASS foreach
# 


trunAllTest_parallel=function(A, Y, c1, c2, beta_e, sigma, m_w="lso", paral=FALSE, paral_n=NULL, resol=1.3, tol=1e-3, maxTry=10, verbose = TRUE){

	p <- ncol(A);
	nX=nrow(A);
	W <- matrix(0, p, p-1);

	if(paral==TRUE){
			if(is.na(paral_n)){
					cl <- parallel::makeCluster(parallel::detectCores() - 1)
			}else{
					cl=parallel::makeCluster(paral_n)
			}
			doParallel::registerDoParallel(cl)

			par_W=foreach (i = 1:p, .combine=rbind, .export=c('linprogPD'), .packages=c('MASS','glmnet')) %dopar% {

					## Null model
					beta_i=beta_e
					beta_i[i]=0
					Ax=A %*% beta_i
					Axy=Ax-Y;
					a1<-(c1-Ax)/sigma
					a2<-(c2-Ax)/sigma
					f1<-dnorm(a1)
					f2<-dnorm(a2)
					F1<-pnorm(-a1)
					F2<-pnorm(a2)
					F12<-F2 + F1
					l_null=nX/2*log(2*pi*sigma^2) + sum(Axy^2)/(2*sigma^2) + sum(log(F12))

					mv=(f2-f1)/F12
					vv=(a2*f2-a1*f1)/F12 + mv^2
					SM=((t(A)%*%A%*%beta_i - t(A)%*%Y)/sigma^2-t(A)%*%mv/sigma)/nX
					HM=t(A)%*%diag(as.numeric(1-vv))%*%A/sigma^2/nX

					if (m_w=="lso"){
						la  = rep(0,nX)                           # Gradient w.r.t parameter of interest
						lb  = matrix(0,nrow = nX, ncol=p-1)
						for(k in 1:nX){
								la[k]=(Y[k]-A[k,]%*%beta_i)*A[k,i]/sigma^2-mv[k]*A[k,i]/sigma
								lb[k,]=(Y[k]-A[k,]%*%beta_i)*A[k,-i]/sigma^2-mv[k]*A[k,-i]/sigma

						}
					
						cv.fit = glmnet(lb,la,standardize = FALSE,intercept = FALSE)
						fit     = cv.fit$glmnet.fit
						tmp     = which(fit$lambda == cv.fit$lambda.min)
						if (sum(fit$beta[,tmp]!=0)==0){
							W[i,] = rep(0,p-1)
						} else {
							W[i,]    = fit$beta[,tmp]
						}
					}else{
						#Danzig-type estimator for W
						tm=sqrt(log(p)/nX)
						
						#W[i,] = Wisolver(HM, i, lambda=tm,search=search,resol=resol)
						wi=tryCatch(HM[i,-i]%*%solve(HM[-i,-i]), error=function(e) {HM[i,-i]%*%ginv(HM[-i,-i])})
						W[i,] <- linprogPD(wi, HM[-i,-i], HM[i,-i], tm)

						tm.stop <- 0;
						try.no <- 1;
						while ((tm.stop != 1)&&(try.no<maxTry)){
							last.wi <- W[i,]
					  		try.no <- try.no+1
					  		tm <- tm*resol;

							wi_e <- tryCatch(linprogPD(wi, HM[-i,-i], HM[i,-i], tm), warning=function(w) w)
							if (is(wi_e,"warning")){
									W[i,] <- last.wi;
						    		tm.stop <- 1;
							}else{
								W[i,]=wi_e
							}
							diff.wi=sum((last.wi-W[i,])^2)
							if(diff.wi<tol){
								tm.stop <- 1;
							}
						}
					}
					
					if (i == 1){
							var   = max(HM[i,i] - W[i,]%*%HM[c((i+1):p),i])
					} else if (i == p){
								var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1)),i])
							} else {
								var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1),(i+1):p),i])
							}

					## Decorrelated Score
					S = SM[i,] - W[i,] %*% SM[-i,]
					pval_score=2*pnorm(-abs(sqrt(nX)*S/sqrt(max(var,0.1))))

					#Partial likelihood test
					#alternative model
					beta_a=beta_e
					Ax=A %*% beta_e
					Axy=Ax-Y;
					a1<-(c1-Ax)/sigma
					a2<-(c2-Ax)/sigma
					f1<-dnorm(a1)
					f2<-dnorm(a2)
					F1<-pnorm(-a1)
					F2<-pnorm(a2)
					F12<-F2 + F1
					mv=(f2-f1)/F12
					vv=(a2*f2-a1*f1)/F12 + mv^2
					SM=((t(A)%*%A%*%beta_e - t(A)%*%Y)/sigma^2-t(A)%*%mv/sigma)
					S = SM[i,] - W[i,] %*% SM[-i,]
					HM=t(A)%*%diag(as.numeric(1-vv))%*%A/sigma^2
					HM_i=tryCatch(solve(HM[i,i]-W[i,] %*% HM[-i,i]), error=function(e) {ginv(HM[i,i]-W[i,] %*% HM[-i,i])})

					# Wald test
					beta_a[i]=beta_e[i]-HM_i*S
					tmp=(nX)*beta_a[i]^2*(max(var,1e-8))
					pval_wald  = 1-pchisq(tmp,1)

					beta_a[-i]=beta_e[-i]-beta_a[i]*W[i,]
					Ax=A %*% beta_a
					Axy=Ax-Y;
					a1<-(c1-Ax)/sigma
					a2<-(c2-Ax)/sigma
					f1<-dnorm(a1)
					f2<-dnorm(a2)
					F1<-pnorm(-a1)
					F2<-pnorm(a2)
					F12<-F2 + F1
					l_altn=nX/2*log(2*pi*sigma^2) + sum(Axy^2)/(2*sigma^2) + sum(log(F12))

					tmp   = -2*(l_altn-l_null)
					pval_plr  = 1-pchisq(tmp,1)

					return( data.frame(pval_score, pval_wald, pval_plr) )
			}
			parallel::stopCluster(cl)
			return(par_W)
	}else{

			pval_score=c()
			pval_wald=c()
			pval_plr=c()
			xperc = 0;
			xp = round(p/10);
			tm=sqrt(log(p)/nX)
			for (i in 1:p) {
					if ((i %% xp)==0){
							xperc = xperc+10;
							if (verbose) {
									print(paste(xperc,"% done",sep=""));
							}
					}
										
					## Null model
					beta_i=beta_e
					beta_i[i]=0
					Ax=A %*% beta_i
					Axy=Ax-Y;
					a1<-(c1-Ax)/sigma
					a2<-(c2-Ax)/sigma
					f1<-dnorm(a1)
					f2<-dnorm(a2)
					F1<-pnorm(-a1)
					F2<-pnorm(a2)
					F12<-F2 + F1
					l_null=nX/2*log(2*pi*sigma^2) + sum(Axy^2)/(2*sigma^2) + sum(log(F12))

					mv=(f2-f1)/F12
					vv=(a2*f2-a1*f1)/F12 + mv^2
					SM=((t(A)%*%A%*%beta_i - t(A)%*%Y)/sigma^2-t(A)%*%mv/sigma)/nX
					HM=t(A)%*%diag(as.numeric(1-vv))%*%A/sigma^2/nX

					if(m_w=="lso"){ ## Use Lasso-type estimator for W solved by glmnet
						la  = rep(0,nX)                           # Gradient w.r.t parameter of interest
						lb  = matrix(0,nrow = nX, ncol=p-1)
						for(k in 1:nX){
								la[k]=(Y[k]-A[k,]%*%beta_i)*A[k,i]/sigma^2-mv[k]*A[k,i]/sigma
								lb[k,]=(Y[k]-A[k,]%*%beta_i)*A[k,-i]/sigma^2-mv[k]*A[k,-i]/sigma

						}
					
						cv.fit = glmnet(lb,la,standardize = FALSE,intercept = FALSE)
						fit     = cv.fit$glmnet.fit
						tmp     = which(fit$lambda == cv.fit$lambda.min)
						if (sum(fit$beta[,tmp]!=0)==0){
							W[i,] = rep(0,p-1)
						} else {
							W[i,]    = fit$beta[,tmp]
						}
					
					}else{
						#Danzig-type estimator for W
						#W[i,] = Wisolver(HM, i, lambda=tm,search=search,resol=resol)
						wi=tryCatch(HM[i,-i]%*%solve(HM[-i,-i]), error=function(e) {HM[i,-i]%*%ginv(HM[-i,-i])})
						W[i,] <- linprogPD(wi, HM[-i,-i], HM[i,-i], tm)

						tm.stop <- 0;
						try.no <- 1;
						while ((tm.stop != 1)&&(try.no<maxTry)){
							last.wi <- W[i,]
					  		try.no <- try.no+1
					  		tm <- tm*resol;

							wi_e <- tryCatch(linprogPD(wi, HM[-i,-i], HM[i,-i], tm), warning=function(w) w)
							if (is(wi_e,"warning")){
									W[i,] <- last.wi;
						    		tm.stop <- 1;
							}else{
								W[i,]=wi_e
							}
							diff.wi=sum((last.wi-W[i,])^2)
							if(diff.wi<tol){
								tm.stop <- 1;
							}
						}
					}

					if (i == 1){
							var   = max(HM[i,i] - W[i,]%*%HM[c((i+1):p),i])
					} else if (i == p){
								var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1)),i])
							} else {
								var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1),(i+1):p),i])
							}

					## Decorrelated Score
					S = SM[i,] - W[i,] %*% SM[-i,]
					pval_score[i]=2*pnorm(-abs(sqrt(nX)*S/sqrt(max(var,0.1))))

					#Partial likelihood test
					#alternative model
					beta_a=beta_e
					Ax=A %*% beta_e
					Axy=Ax-Y;
					a1<-(c1-Ax)/sigma
					a2<-(c2-Ax)/sigma
					f1<-dnorm(a1)
					f2<-dnorm(a2)
					F1<-pnorm(-a1)
					F2<-pnorm(a2)
					F12<-F2 + F1
					mv=(f2-f1)/F12
					vv=(a2*f2-a1*f1)/F12 + mv^2
					SM=((t(A)%*%A%*%beta_e - t(A)%*%Y)/sigma^2-t(A)%*%mv/sigma)
					S = SM[i,] - W[i,] %*% SM[-i,]
					HM=t(A)%*%diag(as.numeric(1-vv))%*%A/sigma^2
					HM_i=tryCatch(solve(HM[i,i]-W[i,] %*% HM[-i,i]), error=function(e) {ginv(HM[i,i]-W[i,] %*% HM[-i,i])})

					# Wald test
					beta_a[i]=beta_e[i]-HM_i*S
					tmp=(nX)*beta_a[i]^2*(max(var,1e-8))
					pval_wald[i]  = 1-pchisq(tmp,1)

					beta_a[-i]=beta_e[-i]-beta_a[i]*W[i,]
					Ax=A %*% beta_a
					Axy=Ax-Y;
					a1<-(c1-Ax)/sigma
					a2<-(c2-Ax)/sigma
					f1<-dnorm(a1)
					f2<-dnorm(a2)
					F1<-pnorm(-a1)
					F2<-pnorm(a2)
					F12<-F2 + F1
					l_altn=nX/2*log(2*pi*sigma^2) + sum(Axy^2)/(2*sigma^2) + sum(log(F12))

					tmp   = -2*(l_altn-l_null)
					pval_plr[i]  = 1-pchisq(tmp,1)
					
			}
			return( data.frame(pval_score, pval_wald, pval_plr) )
	}

}
