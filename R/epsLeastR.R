#' This function solves the sparse penalized Maximum Likelihood Estimate of beta given an estimate of sigma and regularization parameter. This function is prepared based on functions from the SLEP LASSO Matlab package developed by Jun Liu, and Jieping Ye. 
#' @param A Matrix of predictors. Required.
#' @param y Response vector. Required.
#' @param z L_1 norm regularization parameter (z >=0). Required.
#' @param c1 Right censored point. Required.
#' @param c2 Left censored point. Required.
#' @param sigma Estimate of sigma. Required.
#' @param tol Convergence threshold for approximation procedure. Default is 0.001.
#' @param maxIter Maximum iteration number for approximation procedure. Default is 1000.
#' @export


epsLeastR=function(A, y, z, c1, c2, sigma, maxIter=1000, tol=1e-3){

# Verify the number of input parameters
try(if (missing(A) || missing(y) || missing(z) || missing(c1) || missing(c2) || missing(sigma)) stop('\n Inputs: A, y, z and sigma should be specified!\n'))
sigma=as.numeric(sigma)
ill_tag=0
# Get the size of the matrix A
m=nrow(A);
n=ncol(A);

# Verify the length of y
try(if (length(y) != m) stop('\n Check the length of y!\n'));

# Verify the value of z
try(if (z<0) stop('\n z should be nonnegative!\n'));

## Detailed initialization
## Starting point initialization

# compute AT y
ATy=t(A)%*%y;

# process the regularization parameter
# L1 norm regularization
# z here is the scaling factor lying in [0,inf)
try(if (z<0) stop('\n z should be <0'))
lambda=z;


# initialize a starting point
x=ATy;  # if .x0 is not specified, we use ratio*ATy, where ratio is a positive value.


# compute A x
Ax=A%*%x;

x_norm=sum(abs(x));
x_2norm=sum(x*x);
if (x_norm>=1e-6){
    #ratio=initFactor(x_norm, Ax, y, lambda,'LeastR', rsL2, x_2norm);
	ratio=  as.numeric((t(Ax)%*%y - lambda %*% x_norm) / (t(Ax)%*%Ax));
    x=ratio*x;
	Ax=ratio*Ax;
}

## The main program

## The Armijo Goldstein line search scheme + accelearted gradient descent
bFlag=0; # this flag tests whether the gradient step only changes a little

L=1;
# We assume that the maximum eigenvalue of A'A is over 1

# assign xp with x, and Axp with Ax
xp=x;
Axp=Ax;
xxp=matrix(0,nrow=n,ncol=1);

# alphap and alpha are used for computing the weight in forming search point
alphap=0;
alpha=1;
funVal=c();
preerror=0;
for (iterStep in 1:maxIter){
    # --------------------------- step 1 ---------------------------
    # compute search point s based on xp and x (with beta)
    beta=(alphap-1)/alpha;
	s=x + beta* xxp;
	x_old=x;

    # --------------------------- step 2 ---------------------------
    # line search for L and compute the new approximate solution x

    # compute the gradient (g) at s
    As=Ax + beta* (Ax-Axp);
	
	# compute mg
	# c1 - up threshold; c2 - bottom threshold
	a1<-(c1-As)/as.numeric(sigma)
	a2<-(c2-As)/as.numeric(sigma)
	f1<-dnorm(a1)
	f2<-dnorm(a2)
	F1<-pnorm(-a1)
	F2<-pnorm(a2)
	F12<-F2 + F1	
	mg=(f2-f1)/F12
	
    # compute AT As
	ATAs=t(A)%*%As;

    # obtain the gradient g
    g=(ATAs - ATy)/sigma^2 - t(A)%*%mg/sigma ;

    # copy x and Ax to xp and Axp
    xp=x;
	Axp=Ax;

    while (1){
        # let s walk in a step in the antigradient of s to get v
        # and then do the l1-norm regularized projection
        v=s-g/L;

        # L1-norm regularized projection
        x=sign(v)*pmax(abs(v)-lambda / L,0);

        v=x-s;  # the difference between the new approximate solution x and the search point s

        # compute A x
        Ax=A %*% x;

        Av=Ax - As;
        r_sum=sum(v^2);
		l_sum=sum(Av^2);

        if (is.na(r_sum)||is.na(l_sum)){
            bFlag=1; # this shows that, the gradient step makes little improvement
			x=x_old;
            break;
        }
            
        if (r_sum <=1e-20){
            bFlag=1; # this shows that, the gradient step makes little improvement
            break;
        }

        # the condition is ||Av||_2^2 <= (L - rsL2) * ||v||_2^2
        if(l_sum <= r_sum * L){
            break;
		} else{
            L=max(2*L, l_sum/r_sum);
            #print(L);
            }
    }
        
    # --------------------------- step 3 ---------------------------
    # update alpha and alphap, and check whether converge
    alphap=alpha;
	alpha= (1+ sqrt(4*alpha*alpha +1))/2;

    xxp=x-xp;
	Axy=Ax-y;
	a1<-(c1-Ax)/sigma
	a2<-(c2-Ax)/sigma
	F1<-pnorm(-a1)
	F2<-pnorm(a2)
	F12<-F2 + F1
    funVal[iterStep] = m/2*log(2*pi*sigma^2) + sum(Axy^2)/(2*sigma^2) + sum(log(F12)) + sum(abs(x)) * lambda;
        
    if (bFlag){
        #print('\n The program terminates as the gradient step changes the solution very small.');
        break;
    }

    if (iterStep>=2){
		if(funVal[iterStep]- funVal[iterStep-1]>1){
			x=xp
			#print("Cannot converge!")
			if(iterStep==2){
                ill_tag=1
            }
			break;
		}
        if (abs( funVal[iterStep] - funVal[iterStep-1] ) <= tol){
            break;}
    }
	
}
return(list(x=x,ill_tag=ill_tag))
}
