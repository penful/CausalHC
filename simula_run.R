#### 
# Example Code run simulation as an example of the simulations proposed on the 
# paper titled: 
# "A causal latent transition model with multivariate outcomes and unobserved heterogeneity: 
# Application  to human capital development" by 
# - F. Bartolucci (University of Perugia, IT)
# - F. Pennoni (University of Milano-Bicocca, IT)
# - G. Vittadini (University of Milano-Bicocca, IT)
####

# load packages
rm(list=ls())
require(LMest)
require(mvtnorm)
source("sim_data.R")

## set up
# to change for more real values as: # n value 1000,2000
n = 10    # example
k = 2        # value 2,3
rhoc = 1     # choice 1,2
de0c = 1     # choice 1,2
de1c = 1     # choice 1,2
tol = 10^-10 # tolerance for convergence
nit = 10 # to change for more real values as: #nit = 1000   # number of replicates

## other parameters
if(rhoc==1) rho = 0.25
if(rhoc==2) rho = 0.75
if(de0c==1) de0 = -1
if(de0c==2) de0 = 0
if(de1c==1) de1 = rep(0,3)
if(de1c==2) de1 = rep(0.5,3)
Si = rho+(1-rho)*diag(2)
if(k==2){
  Mu = matrix(c(-1,-2,1,2),2)
  Be0alpha = c(0,0)
  Be0beta = c(0,1)
  Be1 = matrix(c(0,0,0,1,0,-1),3)
  Ga0alpha = matrix(0,2,2)
  Ga0beta = matrix(c(0,-1,1,0),2)
  Ga1 = array(0,c(3,2,2))
  Ga1[,1,2] = c(0,1,-1)
  Ga1[,2,1] = c(0,-1,1)
}
if(k==3){
  Mu = matrix(c(-1.5,-3,0,0,1.5,3),2)
  Be0alpha = c(0,1,0.25)
  Be0beta = c(0,2,1.25)
  Be1 = matrix(c(0,0,0,1,0,-1,1,0,-1),3)
  Ga0alpha = matrix(c(0,-1,0,0,0,0,0,-1,0),3)
  Ga0beta = matrix(c(0,-2,-1,1,0,-1,1,0,0),3)
  Ga1 = array(0,c(3,3,3))
  Ga1[,1,2] = Ga1[,1,3] = Ga1[,2,3] = c(0,1,-1)
  Ga1[,2,1] = Ga1[,3,1] = Ga1[,3,2] = c(0,-1,1)
}
file_name_tmp = paste("sim_tmp_n",n,"_k",k,"_rhoc",rhoc,"_de0c",de0c,"_de1c",de1c,".RData",sep = "")
file_name = paste("sim_n",n,"_k",k,"_rhoc",rhoc,"_de0c",de0c,"_de1c",de1c,".RData",sep = "")

## simulate
est0v = est1v = vector("list",nit)
for(it in 1:nit){
  print(it)
# simualate covariates
  Xs = rmvnorm(n,rep(0,3),0.5+0.5*diag(3))
  X = Xs; X[,3] = 2*(Xs[,3]>=0)-1
# simulate responses
  out = sim_data(X,Mu,Si,Be0alpha,Be0beta,Be1,Ga0alpha,Ga0beta,Ga1,de0,de1)
  Y1 = out$Y1; Y2 = out$Y2; Z = out$Z

# estimate models
  YY = array(NA, c(n, 2, 2))
  YY[,1,] = Y1
  YY[,2,] = Y2
  X1 = X
  X22 = array(NA, c(n, 1, 4))
  X22[,1,1:3]= X
  X22[,1,4] = Z

# estimate model without covariates and with covariates
  est0 = try(est_lm_cov_latent_cont(YY,X1 = NULL,X2 = X22[,1,4], k = k, out_se = TRUE, tol=tol))
  est1 = try(est_lm_cov_latent_cont(YY,X1 = X1,X2 = X22[,1,1:4], k = k, out_se = TRUE, tol=tol))
  if(!inherits(est0,"try-error")) est0v[[it]] = est0
  if(!inherits(est1,"try-error")) est1v[[it]] = est1
  
# save results
  save.image(file_name_tmp)
}
save.image(file_name)