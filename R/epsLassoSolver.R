## This function solves EPSLASSO with an initial guess of beta and sigma. 
## 
## @param A Matrix of predictors. Required.
## @param Y Trait values. Required.
## @param c1 Right censored point. Required.
## @param c2 Left censored point. Required.
## @param lambda Penalized parameter (lambda>0) controlling the sparsity.
## @param sigma Initial estimate of sigma.
## @param maxIter Maximum iteration number, default 1000.
## @param verbose Print debugging info or not.
## 


epsLassoSolver=function(A, Y, c1, c2, lambda, sigma, maxIter=1000, verbose=FALSE){
	
	pX=ncol(A)
	nX=nrow(A)
	x=matrix(0,pX,1);
	sigma=as.numeric(sigma)
	esp_beta=1
	esp_var=1
	fit=epsLeastR(A,y=Y,z=lambda,c1,c2,sigma)
        x=fit$x
	step=0
	#sigma_tr=c()
	while ( (esp_beta>0.001 || esp_var>0.001) & step<=maxIter ){
		step=step+1
		sigma_old=sigma
		x_old=x
		Ax = A%*%x
		a1<-(c1-Ax)/sigma
		a2<-(c2-Ax)/sigma
		f1<-dnorm(a1)
		f2<-dnorm(a2)
		F1<-pnorm(-a1)
		F2<-pnorm(a2)
		f12<-f1-f2
		Y_Ax<-(Y-Ax)
		F12<-F2 + F1
		#sigma = tryCatch(uniroot(function(sigmaT) sum(-1/sigmaT+(Y_Ax)^2/sigmaT^3+(((c2-Ax)/sigmaT^2)*f2-((c1-Ax)/sigmaT^2)*f1)/F12),c(sigma/500,500*sigma))$root, 
		#		error=function(e){if(verbose){print("Warning")};sigma=sigma_old})
		bb=-sum(((c2-Ax)*f2-(c1-Ax)*f1)/F12)
		cc=-sum(Y_Ax^2)
		sigma=(-bb+sqrt(bb^2-4*cc*nX))/(2*nX)
		#sigma_tr[step]=sigma
		if(sigma==sigma_old){
			break;
		}
		fit=epsLeastR(A,y=Y,z=lambda,c1,c2,sigma)
                if(fit$ill_tag==1){
                        sigma=sigma_old
                        x=x_old
                        break
                }else{
                        x=fit$x
                }

		esp_beta=norm(x-x_old)
		esp_var=abs(sigma-sigma_old)
		if(step>1 && esp_beta>pre_esp_beta && esp_var>pre_esp_var){
			if(verbose){
				print("Diverging.")
			}
			sigma=sigma_old
			x=x_old
			break;
		}
		pre_esp_beta = esp_beta
		pre_esp_var = esp_var
	}
	out = list(lambda = lambda, sigma=sigma,beta=x)
}
