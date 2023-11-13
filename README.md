# High-Dimensional Inference for Generalized Linear Models with Hidden Confounding

## Description 

Statistical inferences for high-dimensional regression models have been extensively studied for their wide applications ranging from genomics, neuroscience, to economics. However, in practice, there are often potential unmeasured confounders associated with both the response and covariates, which can lead to invalidity of standard debiasing methods. This paper focuses on a generalized linear regression framework with hidden confounding and proposes a debiasing approach to address this high-dimensional problem, by adjusting for the effects induced by the unmeasured confounders. We establish consistency and asymptotic normality for the proposed debiased estimator. The finite sample performance of the proposed method is demonstrated through extensive numerical studies and an application to a genetic data set.

## Usage

The files include an example of n=100 observations with covariate dimension p=300. The default model is that the response y given covariates and unmeasured confounders follow a logistic regression. We estimate the unmeasured confounders and then perform debias step to obtain point and interval estimation results for coefficient of interest. 

``main.R`` includes the main part of this example. It computes the coverage of the estimated interval as well as the average confidence interval width over 300 replications.
Source_functions.R includes the data generating function and the point and interval estimation function.
cate_src.R is a file includes high dimensional factor analysis function. Just in case library(cate) does not work on your side, this open source file serves as a back-up to use the functions we need to adjust for unmeasured confounders.
