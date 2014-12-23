# Reza Asad	

# This clears all the varaibles previously
# stored.
rm(list = setdiff(ls(), lsf.str()))

# Get the Data
data(Boston, package = 'MASS');

# Fit regression with no covariets and
# full covariets.Improve the fits using
# Mallow's Cp forward and backward feature
# selection.
null.reg <- lm(medv~1, data = Boston);
full.reg <- lm(medv~., data = Boston);
sih <- summary(full.reg)$sigma;
null.forward <- step(null.reg,scope = list(lower = null.reg, upper = full.reg)
	, data = Boston)
full.backward <- step(full.reg,scope = list(lower = null.reg, upper = full.reg)
	, data = Boston)
# They result in the same model and the AIC was 1585.76.

# 10 Fold CV for Prediction Error 
k = 10;
n = nrow(Boston)
N = 100;
ii = rep(1:k, each = floor(n/k))
if(length(ii)<n){
	ii = c(ii, rep(1:(n-length(ii))));
}
	
cv.err = rep(0, length = N)
cvF.err = rep(0, length = N)
cvN.err = rep(0, length = N)

for( j in 1:N){
	ii = sample(ii);
	for( i in 1:k){
		lm.cv <- step(lm(medv~., data = Boston[ii != i, ]),
			direction = 'backward', scale = sih, trace = 0);
		lm.pr <- predict(lm.cv, newdata = Boston[ii == i,]);
		lmF.cv <- lm(medv~., data = Boston, subset = (ii != i));
		lmF.pr <- predict(lmF.cv, newdata = Boston[ii == i,]);
		lmN.cv <- lm(medv~1, data = Boston, subset = (ii != i));
		lmN.pr <- predict(lmN.cv, newdata = Boston[ii == i,]);
		cvN.err[j] <- cvN.err[j] + sum((lmN.pr-Boston[ii == i,]$medv)^2);
		cvF.err[j] <- cvF.err[j] + sum((lmF.pr-Boston[ii == i,]$medv)^2);
		cv.err[j] <- cv.err[j] + sum((lm.pr-Boston[ii == i,]$medv)^2);
	}
}	
print(mean(cv.err)/n)
print(mean(cvF.err)/n)
print(mean(cvN.err)/n)

# Stagewise Selection
#library(lars)
#b <- lars(x = as.matrix(Boston)[,-14], y = Boston[,14], 
#	type = 'forward.stagewise')
#cvlas <- cv.lars(as.matrix(Boston)[,-14], y = Boston[,14], type="lasso")
#cvlas
#frac <- cvlas$fraction[which.min(cvlas$cv)]
#frac
#las.coef <- predict.lars(b, type="coefficients", mode="fraction", s=1)

# Regression Tree, Another Way to Select Variables
library(tree)
regT <- tree(medv~., data = Boston)

# 10 Fold CV for prediction error
k = 10;
n = nrow(Boston)
N = 100;
ii = rep(1:k, each = floor(n/k))
if(length(ii)<n){
	ii = c(ii, rep(1:(n-length(ii))));
}
cvT.err = rep(0, N)
for( j in 1:N){
	ii = sample(ii);	
	for( i in 1:k){
		regT <- tree(medv~., data = Boston[ii != i,])
		regT.pr <- predict(regT, newdata = Boston[ii == i,]);
		cvT.err[j] <- cvT.err[j] + sum((regT.pr-Boston[ii == i,]$medv)^2);
	}
}	
print(mean(cvT.err)/n)

# Box Plot Comparing Prediction Errors
boxplot(data.frame(cv.err/n, cvT.err/n), ylab='10-fold CV',
        names=c('Cp', 'RegTree'), ylim=c(16, 28), col='gray80')

# Ridge Regression Comparision
library(MASS)

# 5000 Values of Lambda Between 0 and 1000
lam = seq(0,1000, length = 5000)

# Ridge Regression using the lambdas.
lmRidg <- lm.ridge(medv~., data = Boston, lambda = lam)

# Plot Change of Coefficients as a Function of Lambdas
#matplot(lam, coef(lmRidg)[,-1], lwd = 4,lty =1, type = 'l',
#	xlab = expression(lambda), ylab = expression(beta[lambda]))

# Compute Degress of Freedom for Each Lambda
df.lam <- rep(0, length(lam));
x = as.matrix(Boston[,-14])
x = scale(x, scale = TRUE)
for(i in 1:length(lam)){
	xhat <- x%*%solve(t(x)%*%x+ lam[i]*diag(ncol(x)),t(x))
	df.lam[i] <- sum(diag(xhat))
}
df <- seq(1:(ncol(x)));
edf.lam <- rep(0, length = length(df));
# Find the Lambda Corresponding to Each Degree of Freedom  
for(i in 1:length(df)){
	edf.lam[i] <- lam[which.min(abs(df.lam - df[i]))]
}

# 10 Fold CV to Pick the Best Lambda
k <- 10;
n <- nrow(x);
N <- 100;

ii <- rep(1:k, each = floor(n/k));
if(length(ii)<n){
	ii = c(ii, rep(1:k, length = n - length(ii)))
}
cvR.err = matrix(0, N,length(edf.lam))
for(j in 1:N){
	ii = sample(ii)
	for(i in 1:k){
		for(m in 1:length(edf.lam)){
			lmRidg <- lm.ridge(medv~., data = Boston[ii!= i,],
				 lam = edf.lam[m])
			lmR.pr <- as.vector(as.matrix(Boston[ii == i,-14])%*% coef(lmRidg)[-1]+coef(lmRidg)[1]);
			cvR.err[j,m] <- cvR.err[j,m] + sum((Boston[ii == i,]$medv - lmR.pr)^2);  
		}
	}
}
print(colMeans(cvR.err/n))

# Box Plot For comparing Ridge Regression and Cp
boxplot(data.frame(cvR.err/n, cv.err/n), ylab = '10-fold CV',
	names = c(paste('edf-', df, sep = ''), 'Cp'), ylim = c(23, 32), col = 'gray80')




