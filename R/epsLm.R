#' Linear regression for extreme phenotype sampling data under low-dimensional situation (sample size n > p number of predictors).
#'
#' This function solves the Maximum Likelihood Estimate of the low-dimensional linear model for extreme phenotype sampling data using Newton-Raphson (NR) procedure. This function is prepared based on functions from the R package CEPSKAT. 
#'
#' @param formula Regression model to be fit. Required.
#' @param c1 Right censored point. Required.
#' @param c2 Left censored point. Required.
#' @param delta Convergence threshold for NR procedure. Default is 0.001.
#' @param MAXITERNUM Maximum iteration number for NR procedure. Default is 1000.
#' @param data The dataframe stores data for the formula. Default is NULL.
#' @param verbose Print debugging info or not. Default is FALSE.
#' @export
#' @examples
#' n=100
#' p1=0.2
#' p2=0.2
#' X=rnorm(n)
#' Y=1+0.5*X+rnorm(n)
#'
#' Y_eps=Y[order(Y)[c(1:(n*p1),(n-n*p2+1):n)]]
#' X_eps=X[order(Y)[c(1:(n*p1),(n-n*p2+1):n)]]
#' c1=Y_eps[n*p1+1]
#' c2=Y_eps[n*p1]
#' 
#' res=lm_eps(Y_eps~X_eps, c1, c2)
#' res
#'

lm_eps <-
function(formula, c1, c2, delta=0.001, MAXITERNUM=1000, data=NULL, verbose=FALSE){
	preerror = 0
	X<-model.matrix(formula,data=data)
	Y<-model.frame(formula, data=data)[,1]
	out.lm<-lm(formula, data=data)
	sigma.ols =mean(abs(out.lm$residuals))
	c = (c1-c2)/2 
	gamma0 = c(out.lm$coefficients,sigma.ols)
	n = nrow(X)
	K = ncol(X)
	iternum = 0
	sigma = gamma0[K+1]
	alpha0 = gamma0[-(K+1)]

	while(T){
		Xa = (X%*%alpha0)[,1]
		J = matrix(0,nrow=K,ncol=K)
		a1<-(c1-Xa)/sigma
		a2<-(c2-Xa)/sigma
		f1<-dnorm(a1)
		f2<-dnorm(a2)
		F1<-pnorm(-a1)
		F2<-pnorm(a2)
		Y_Xa<-(Y-Xa)
		f12<-f1-f2
		F12<-F2 + F1
		mv=(f2-f1)/F12
		vv=(a2*f2-a1*f1)/F12 + mv^2
		J[1:K, 1:K]<-t(X) %*%diag(as.numeric(-1+vv))%*%X / sigma^2
		Jinv = solve(J)
		V = matrix(0,nrow=K,ncol=1)
		V[1:K,1] = colSums((X/sigma)*as.vector(Y_Xa/sigma-f12/F12))
	    sigma_next = tryCatch(uniroot(function(sigmaT) sum(-1/sigmaT+(Y_Xa)^2/sigmaT^3+(((c2-Xa)/sigmaT^2)*f2-((c1-Xa)/sigmaT^2)*f1)/F12),c(sigma/500,500*sigma))$root,  error=function(e){if(verbose){print("Warning")};sigma_next =sigma})
		alphaNext = alpha0-Jinv%*%V
		iternum = iternum + 1
		curerror = sum(abs(alphaNext-alpha0))+abs(sigma-sigma_next)
		if(iternum > MAXITERNUM && curerror > preerror){
			print("Newton-Raphson diverging.")
			break
		}
		preerror = curerror  
		if(curerror < delta){
			break
		}
		alpha0 = alphaNext
		sigma=sigma_next
	}
	MLEs = c(alphaNext,sigma)
	MLEs = data.frame(MLEs)
	rownames = rep(0,K+1)
	for(i in 1:K){
        rownames[i] = paste("alpha",i-1,sep="")
	}
	#rownames[1] = "intercept"
	rownames[K+1] = "sigma"
	row.names(MLEs) = rownames
	alpha0 = alphaNext
	pval=c()
	for(i in 1:length(alpha0)){
		alpha_t=alpha0
		alpha_t[i]=0
		Ax=X %*% alpha_t
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
		SM=((t(X)%*%X%*%alpha_t - t(X)%*%Y)/sigma^2-t(X)%*%mv/sigma)/n
		HM=t(X)%*%diag(as.numeric(1-vv))%*%X/sigma^2/n
		pval[i]=2*pnorm(-abs(sqrt(n)*SM[i,]/sqrt(HM[i,i])))
	}    
  return(list(coef=alpha0, sigma=sigma, pvals=pval))
}


