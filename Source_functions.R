###################################################################
# date : Oct 13th
# this file provides the source functions for debiasing method
# under hidden confounding effects.
###################################################################
library(glmnet) 
library(cate) 
library(MASS)
library(purrr)
library(paran)
set.seed(1)
#source("cate_src.R") # uncomment and use open source file if library "cate" installation failed

### Data generating function
generate_data = function(n=300, p=600, K=3, s=1, rho=1, mod="logistic", loading="sparse", sigma_e=1)
{
  ## FUNCTION: generate the covariates and response according to the respective 
  ##           generalized linear model with unmeasured confounders considered.
  ##
  ## INPUTS: n - the sample size of the data
  ##         p - the dimension of covariates
  ##         K - the dimension of the unmeasured confounders
  ##         s - number of nonzero elements in covariate coefficients
  ##         rho - diagonal entry of error variance matrix for factor model X = W^T U + E
  ##         mod - linear or logistic (can be generated to more GLMs)
  ##         loading - sparse or dense (when loading is sparse, p should be divisible by 3)
  ##         sigma_e - error variance for linear model of y
  ##
  ## OUTPUTS: X - covariates
  ##          Y - response
  ##          U - unmeasured confounders
  
  ## Set true U, beta and W 
  U <- matrix(rnorm(n*K, 0, 1), nrow = n, ncol = K)   
  beta <- rep(1,K)
  subp <- p/3
  if(loading == "sparse"){
    i3_p1 <- matrix(data=c(rep(0.5,subp), rep(0,subp),rep(0,subp)), 
                    nrow=subp, ncol=K, byrow=F) 
    i3_p2 <- matrix(data=c(rep(0,subp), rep(1,subp),rep(0,subp)), 
                    nrow=subp, ncol=K, byrow=F)
    i3_p3 <- matrix(data=c(rep(0,subp), rep(0,subp),rep(1.5,subp)), 
                    nrow=subp, ncol=K, byrow=F)
    W <- cbind(t(i3_p1),t(i3_p2),t(i3_p3)) 
  }
  if(loading == "dense"){
    W <- matrix(runif(K*p, min = 0, max = 1), nrow = K, ncol = p)
  }
  gamma <- c(0, rep(1, s), rep(0, p-s-1)) 
  
  ## Generate X = W^T U + E, E ~ N( 0, rho * I_p)
  sigma_E <- diag(rho,p) 
  E <- mvrnorm(n, rep(0,p), sigma_E)
  X <- U %*% W + E

  ## Generate Y 
  if(mod == "linear")
  {
    e <- rnorm(n, 0, sigma_e)
    Y <- U%*% beta + X %*% gamma + e
  }
  
  if(mod == "logistic")
  {
    xb <- U%*% beta + X %*% gamma
    prob <- 1/(1 + exp(-xb))
    Y <- rbinom(n=n, size=1, prob=prob)
  }
  
  return(list("X"=X, "Y"=Y))
}

### point and interval estimation function
estimate_coefficient = function(Y, X, mod="logistic")
{
  ## FUNCTION: estimate the coefficient of interest using the proposed methodology
  ##           and take the hidden confounding effects into account
  ##
  ## INPUTS: Y - response
  ##         X - covariates
  ##         mod - linear or logistic (can be generated to more GLMs)
  ##
  ## OUTPUTS: K_hat - estimated number of confounders
  ##          theta - debiased estimator for parameter of interest
  ##          upper - upper bound for estimated confidence interval
  ##          lower - lower bound for estimated confidence interval
  ##          length - the length for the estimated confidence interval
  
  D <- X[,1]  # interested covariates
  Q <- X[,-1] # nuisance covariates
  
  ## Parallel Analysis to estimate unmeasured confounders dimension
  PA <- paran(X, iterations = 500, centile = 0)
  K_hat <- PA$Retained

  ## MLE for W and Sigma_E by EM algorithm under working identifiability condition
  fa <- factor.analysis(X, K_hat, method = "ml")  
  Wupdate.t <- fa$Gamma
  Sgm_inv <- solve(diag(fa$Sigma))
  orthg <- t(Wupdate.t) %*% Sgm_inv %*% Wupdate.t/p
  V <- eigen(orthg)$vectors
  Wupdate <- t(Wupdate.t %*% V)
  
  ## GLS Estimator for U_hat 
  Uupdate <- t(solve(Wupdate %*% Sgm_inv %*% t(Wupdate)) %*% Wupdate %*% Sgm_inv %*% t(X)) 
  M <- cbind(Uupdate,Q)
  Z <- cbind(Uupdate,X)
  
  if(mod == "linear"){
    # Initial estimator  
    cv.fit1 <- cv.glmnet(Z, Y, standardize=FALSE, family="gaussian", 
                         penalty.factor=c(rep(0,K_hat), rep(1,p)))
    coef_hat <- coef(cv.fit1, s="lambda.1se", complete=TRUE)[-1]
    beta_hat <- coef_hat[1:K_hat]
    theta_hat <- coef_hat[K_hat+1]
    gamma_hat <- coef_hat[-1:-(K_hat+1)]
    
    # Debiasing step
    nzero_bar <- 0.5  # empirically control for the sparsity level 
    bprime <- Z %*% coef_hat
    bprime2 <- rep(1,p)
    
    cv.fit2 <- cv.glmnet(M, D, type.measure="mse", alpha=1, family="gaussian")
    nzero_level <- cv.fit2$nzero/p                 # nonzero proportion at each lambda 
    lmd <- cv.fit2$lambda[nzero_level < nzero_bar] # filter lambda
    cvsd <- cv.fit2$cvsd[nzero_level < nzero_bar]  # filter sd
    cvm <- cv.fit2$cvm[nzero_level < nzero_bar]    # filter mse
    three_se <- 3*cvsd[which.min(cvm)]             
    
    # largest lambda such that the error is within three standard error of the minimum
    lambda.3se <- max(lmd[abs(cvm-min(cvm)) < three_se]) 
    fit2 <- glmnet(M, D, lambda = lambda.3se)
    w <- coef(fit2)[-1]
    I_hat <- D %*% D/n - w %*% t(M) %*% D /n 
    score <- -t(D - M %*% w) %*% (Y - Z %*% coef_hat)/n
  }
  
  if(mod == "logistic"){
    # Initial estimator
    cv.fit1 <- cv.glmnet(Z, Y, standardize=FALSE, alpha=1, family="binomial", 
                        penalty.factor=c(rep(0,K_hat), rep(1,p)))
    coef_hat <- coef(cv.fit1, s="lambda.1se", complete=TRUE)[-1]
    beta_hat <- coef_hat[1:K_hat]
    theta_hat <- coef_hat[K_hat+1]
    gamma_hat <- coef_hat[-1:-(K_hat+1)]
    
    # Debiasing step
    nzero_bar <- 0.5 # empirically control for the sparsity level 
    bprime <- exp(Z %*% coef_hat)/(1+exp(Z %*% coef_hat))
    bprime2 <- exp(Z %*% coef_hat)/(1+exp(Z %*% coef_hat))^2
    cv.fit2 <- cv.glmnet(M, D, type.measure="mse", alpha=1, family="gaussian", 
                        weights = bprime2)
    
    nzero_level <- cv.fit2$nzero/p                  # nonzero proportion at each lambda 
    lmd <- cv.fit2$lambda[nzero_level < nzero_bar]  # filter lambda
    cvsd <- cv.fit2$cvsd[nzero_level < nzero_bar]   # filter sd
    cvm <- cv.fit2$cvm[nzero_level < nzero_bar]     # filter mse
    three_se <- 3*cvsd[which.min(cvm)] 
    
    # largest lambda such that the error is within three standard error of the minimum
    lambda.3se <- max(lmd[abs(cvm-min(cvm)) < three_se])
    fit2 <- glmnet(M, D, weights = bprime2, lambda = lambda.3se)
    w <-  coef(fit2)[-1]
    I_hat <- mean(bprime2 * D * (D - M %*% w))
    score <- -mean((Y - bprime)*(D - M %*% w))
    
  }

  
  theta_debiased <- theta_hat - score/I_hat
  ci_upper <- theta_debiased + 1.96 * I_hat^(-1/2)/sqrt(n) 
  ci_lower <- theta_debiased - 1.96 * I_hat^(-1/2)/sqrt(n) 
  ci_length <- 2* 1.96 * I_hat^(-1/2)/sqrt(n) 

  
  return(list("K_hat"=K_hat, "theta"=theta_debiased,
              "upper"=ci_upper, "lower"=ci_lower,
              "length"=ci_length))
  
}

