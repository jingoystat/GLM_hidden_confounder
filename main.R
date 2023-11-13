###################################################################
# date : Oct 13th
# this file is the main file for debiasing method simulation
# under hidden confounding effects.
###################################################################
library(glmnet) 
library(cate) 
library(MASS)
library(purrr)
library(paran)
set.seed(1)
#source("cate_src.R") # uncomment this open source file in case library "cate" installation failed
source("Source_functions.R") 


p <- 300        # X_i dimension 
n <- 100        # sample size  
nsim <- 300     # replications

ci_length <- ci_upper <- ci_lower <- 
  coverage_full <- theta_debiased <- rep(NA, nsim)


for(sim in 1:nsim){
  simu_data <- generate_data(n=n, p=p)
  debiase_results <- estimate_coefficient(simu_data$Y, simu_data$X)
  
  theta_debiased[sim] <- debiase_results$theta
  ci_upper[sim] <- debiase_results$upper 
  ci_lower[sim] <- debiase_results$lower
  coverage_full[sim] <- 1*(ci_upper[sim]>0)*(ci_lower[sim]<0)
  ci_length[sim] <- debiase_results$length
}

coverage_full
(coverage <- mean(coverage_full))
(ci_length_mean <- mean(ci_length))
