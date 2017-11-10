#' Test for high-dimensional sparse regression for extreme phenotype sampling data.
#'
#' This function tests the effects of predictors modelled jointly for high-dimensional extreme phenotype sampling data.
#' @param X Matrix of predictors. Required.
#' @param Y Trait values. Required.
#' @param c1 Right censored point. Required.
#' @param c2 Left censored point. Required.
#' @param lam0 A sequence of lambda values. Default is the sequence used in GLMNET.
#' @param m_w Methods used to estimate W matrix. Default is "lso" for LASSO solution using glmnet. Another method is "dzg" for Danzig-type estimator.
#' @param scal Scale matrix X or not. Default is TRUE.
#' @param paral Parallel computing indicator. Default is FALSE, not using parallel.
#' @param paral_n Number of cores that are used for parallel computing. Default is NULL. When paral is TRUE, default is the number of system available cores - 1.
#' @param resol The refining step when m_w="dzg". Default is 1.3. A large resol results in faster convergence speed, but rough solutions.
#' @param tol The convergence threshold for refining when m_w="dzg". Default is 0.001.
#' @param maxTry The maximum refining steps when m_w="dzg". Default is 10.
#' @param verbose Print debugging info or not.
#' @import glmnet foreach
#' @importFrom methods is
#' @importFrom stats dnorm lm model.frame model.matrix pchisq pnorm qnorm uniroot
#' @export
#' @examples
#' library(mvtnorm)
#' sd=1
#' n=100
#' p1=0.2
#' p2=0.2
#' p=100
#' nc=10
#' eff=0.5
#'
#' beta_eff=c(rep(eff,nc),rep(0,p-nc))
#' 
#' cov_m=diag(p)
#' X_b=rmvnorm(n,mean=c(rep(1,p/2),rep(2,p/2)), sigma=cov_m)
#' Y_b=X_b%*%beta_eff+rnorm(n,0,sd)
#' sample_threshold_low=ceiling(n*p1)
#' sample_threshold_up=ceiling(n*p2)
#'
#' ind_sample = order(Y_b)[c(1:sample_threshold_low,(n-sample_threshold_up+1):n)]
#' Y_p = Y_b[ind_sample]
#' c1 = min(Y_b[ind_sample[-(1:sample_threshold_low)]])
#' c2 = max(Y_b[ind_sample[1:sample_threshold_low]])
#' X_p = X_b[ind_sample,]
#'
#' res=epsLasso(X_p,Y_p,c1,c2)
#' res
#'

epsLasso=function(X, Y, c1, c2, lam0=NULL, m_w="lso", scal=TRUE, paral=FALSE, paral_n=NULL, resol=1.3, tol=1e-3, maxTry=10, verbose = TRUE){

	try(if (missing(X) || missing(Y) || missing(c1) || missing(c2) ) stop('\n Inputs: X, Y, c1 and c2 should be specified!\n'))
	
	if(scal){
		A=scale(X);
	}else{
		A=X
	}
	
	y_mu=mean(Y);
	Y=Y-y_mu;
	c1=c1-y_mu;
	c2=c2-y_mu;

	nX=nrow(A)
	pX=ncol(A)

	if(is.null(lam0)){
		lamMax=max(abs(t(A)%*%Y))
		if(nX<pX){
			lamMin=0.01*lamMax
		}else{
			lamMin=0.0001*lamMax
		}
		lamR=exp((log(lamMax)-log(lamMin))/99)
		lam0=c(lamMax,lamMax*lamR^(-seq(1:99)))
		lam0=lam0[order(lam0)];
	}
	##initial guess from glmnet
	if(F){
		rlasso=cv.glmnet(A,Y,intercept=F,standardize = F)
        rlasso=glmnet(A,Y,lambda=rlasso$lambda.min,intercept=F,standardize = F)
        f_idx_glmnet=rlasso$beta@i+1
		cep_fit=lm_eps(Y~A[,f_idx_glmnet]-1,c1,c2)
	}
	##initial guess from square-root lasso
	if(T){
		lambda_e <- sqrt(qnorm(1-(0.1/pX))/nX);
		outLas <- flare::slim(A,Y,lambda=c(lambda_e),method="lq",q=2,verbose=FALSE);
		beta_e=outLas$beta
		f_idx_las=which(beta_e!=0)
		if(length(f_idx_las)<nX){
			cep_fit=lm_eps(Y~A[,f_idx_las]-1,c1,c2)
			sigma_i=max(sqrt(sum((Y-A%*%outLas$beta)^2)/nX),cep_fit$sigma)
		}else{
			sigma_i=sqrt(sum((Y-A%*%outLas$beta)^2)/nX)
		}
	}
	
	##tuning lambda -- cross validation
	if(F){	
		nfolds=10
			foldid = sample(rep(seq(nfolds), length = nX))

		if(paral==T){
			if(is.na(paral_n)){
				cl <- parallel::makeCluster(parallel::detectCores() - 1)
			}else{
				cl=parallel::makeCluster(paral_n)
			}
			doParallel::registerDoParallel(cl)
			cv_eps=foreach(i = 1:length(lam0), .combine=c, .export=c('epsLassoSolver','epsLeastR','linprogPD'), .packages='MASS') %dopar% {
						mse=c()
						sigma=cep_fit$sigma
						outlist = as.list(seq(nfolds))
						for (j in seq(nfolds)) {
								which = foldid == j
								if (is.matrix(Y))
										Y_sub = Y[!which, ]
								else Y_sub = Y[!which]
								outlist[[j]] = epsLassoSolver(A[!which, , drop = FALSE], Y_sub, c1, c2, lambda = lam0[i], sigma)
								Ax=A[which,] %*% outlist[[j]]$beta
								Axy=Ax-Y[which];
								mse[j]=sum(Axy^2)
						}
				return(mean(mse))
				}
			parallel::stopCluster(cl)

		}else{
			cv_eps=c()
			for(i in 1:length(lam0)){
					mse=c()
					sigma=cep_fit$sigma
					outlist = as.list(seq(nfolds))
					for (j in seq(nfolds)) {
							which = foldid == j
							if (is.matrix(Y))
									Y_sub = Y[!which, ]
							else Y_sub = Y[!which]
							outlist[[j]] = epsLassoSolver(A[!which, , drop = FALSE], Y_sub, c1, c2, lambda = lam0[i], sigma)
							Ax=A[which,] %*% outlist[[j]]$beta
							Axy=Ax-Y[which];
							mse[j]=sum(Axy^2)
					}
					cv_eps[i]=mean(mse)
			}
		}

		lam_e=lam0[which(cv_eps==min(cv_eps))[1]]
			res_eps=epsLassoSolver(A,Y,c1,c2,lambda=lam_e,sigma=cep_fit$sigma)
			beta_e=res_eps$beta
			sigma_e=res_eps$sigma
			lambda_e=res_eps$lambda
	}
	

	if(paral==T){
		if(is.na(paral_n)){
			cl <- parallel::makeCluster(parallel::detectCores() - 1)
		}else{
			cl=parallel::makeCluster(paral_n)
		}
		doParallel::registerDoParallel(cl)
		ebic=foreach(i = 1:length(lam0), .combine=c, .export=c('epsLassoSolver','epsLeastR','linprogPD'), .packages='MASS') %dopar% {
			#sigma=cep_fit$sigma
			out_eps=epsLassoSolver(A,Y,c1,c2,lambda=lam0[i],sigma_i)
			Ax=A %*% out_eps$beta
			Axy=Ax-Y;
			a1<-(c1-Ax)/out_eps$sigma
			a2<-(c2-Ax)/out_eps$sigma
			F1<-pnorm(-a1)
			F2<-pnorm(a2)
			F12<-F2 + F1
			K=sum(abs(out_eps$beta)>0)
			ebic_i=2*(nX/2*log(2*pi*out_eps$sigma^2) + sum(Axy^2)/(2*out_eps$sigma^2) + sum(log(F12)))+K*log(nX)
			return(ebic_i)
		}
		parallel::stopCluster(cl)
		lam_e=lam0[which(ebic==min(ebic))[1]]
		res_eps=epsLassoSolver(A,Y,c1,c2,lambda=lam_e,sigma=sigma_i)
		beta_e=res_eps$beta
		sigma_e=res_eps$sigma
		lambda_e=res_eps$lambda
		if(sum(beta_e!=0)>0){
			eps_fit=lm_eps(Y~A[,which(beta_e!=0)]-1,c1,c2)
			beta_e[which(beta_e!=0)]=eps_fit$coef
			sigma_e=eps_fit$sigma
		}
    }else{
		out_eps=list()
		ebic=c()
		for(i in 1:length(lam0)){
				sigma=cep_fit$sigma
				out_eps[[i]]=epsLassoSolver(A,Y,c1,c2,lambda=lam0[i],sigma)
				Ax=A %*% out_eps[[i]]$beta
				Axy=Ax-Y;
				a1<-(c1-Ax)/out_eps[[i]]$sigma
				a2<-(c2-Ax)/out_eps[[i]]$sigma
				F1<-pnorm(-a1)
				F2<-pnorm(a2)
				F12<-F2 + F1
				K=sum(abs(out_eps[[i]]$beta)>0)
				ebic[i]=2*(nX/2*log(2*pi*out_eps[[i]]$sigma^2) + sum(Axy^2)/(2*out_eps[[i]]$sigma^2) + sum(log(F12)))+K*(pX^(1/3))
		}
		beta_e=out_eps[[which(ebic==min(ebic))[1]]]$beta
		sigma_e=out_eps[[which(ebic==min(ebic))[1]]]$sigma
		lambda_e=out_eps[[which(ebic==min(ebic))[1]]]$lambda
		if(sum(beta_e!=0)>0){
			eps_fit=lm_eps(Y~A[,which(beta_e!=0)]-1,c1,c2)
			beta_e[which(beta_e!=0)]=eps_fit$coef
			sigma_e=eps_fit$sigma
		}
	}

	#if(method=="score"){
		p_val=trunAllTest_parallel(A, Y, c1, c2, beta_e, sigma_e, m_w=m_w, paral=paral, paral_n=paral_n, resol=resol, tol=tol, maxTry=maxTry, verbose = verbose)
	#}else{
	#	p_val=trunPLRTest_parallel(A, Y, c1, c2, beta_e, sigma_e, lambda_e, paral=paral, paral_n=paral_n,search=F)
	#}
	return(list(pvals=p_val,sigma=sigma_e))
	
}
