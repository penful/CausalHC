#### 
# Example Code to generate dasta for an example of the simulations proposed on the 
# paper titled: 
# "A causal latent transition model with multivariate outcomes and unobserved heterogeneity: 
# Application  to human capital development" by 
# - F. Bartolucci (University of Perugia, IT)
# - F. Pennoni (University of Milano-Bicocca, IT)
# - G. Vittadini (University of Milano-Bicocca, IT)
####


sim_data <- function(X,Mu,Si,Be0alpha,Be0beta,Be1,Ga0alpha,Ga0beta,Ga1,de0,de1){

#simulate data 
#
# INPUT:
# X        = matrix of covariates
# Mu       = matrix of means of the responses (r x k)
# Si       = variance-covarinace matrix of the responses (r x r)
# Be0alpha = intercept for latent state at 1st time without treatment
# Be0beta  = intercept for latent state at 1st time with treatment
# Be1      = regression coefficients for latent state at 1st time
# Ga0alpha = interpect for transition to latent state at 2nd time without treatment
# Ga0beta  = interpect for transition to latent state at 2nd time without treatment
# Ga1      = regression coefficients for transition to latent state at 2nd
# de0      = intercept in the logistic model for the treatment
# de1      = regression coefficients in the logistic model for the treatment
#
# OUTPUT:
# Y1       = response at the 1st time occasion
# Y2       = response at the 2nd time occasion
# Z        = treatment
# H1       = latent state at the 1st time occasion
# H2       = latent state at the 2sn time occasion

# preliminaries
  n = nrow(X)
  r = nrow(Mu)
  k = ncol(Mu)

# simulate latent state and responses at the 1st time occasion
  PH1 = exp(rep(1,n)%o%Be0alpha+X%*%Be1)
  PH1 = (1/rowSums(PH1))*PH1
  H1 = rep(0,n)
  for(i in 1:n) H1[i] = which(rmultinom(1,1,PH1[i,])==1)
  Y1 = matrix(0,n,r)
  for(h in 1:k){
    ind = which(H1==h)
    lind = length(ind)
    if(lind>0) Y1[ind,] = rmvnorm(lind,Mu[,h],Si)
  }
  
# simulate treatment
  pZ = exp(de0+X%*%de1); pZ = pZ/(1+pZ)
  Z = 1*(runif(n)<=pZ)

# simulate latent state and responses at the 2nd time occasion
  PH2 = matrix(0,n,k)
  for(i in 1:n){
    if(Z[i]==0) PH2[i,] = exp(Ga0alpha[H1[i],]+X[i,]%*%Ga1[,H1[i],])
    else PH2[i,] = exp(Ga0beta[H1[i],]+X[i,]%*%Ga1[,H1[i],])
  }
  PH2 = (1/rowSums(PH2))*PH2
  H2 = rep(0,n)
  for(i in 1:n) H2[i] = which(rmultinom(1,1,PH2[i,])==1)
  Y2 = matrix(0,n,r)
  for(h in 1:k){
    ind = which(H2==h)
    lind = length(ind)
    if(lind>0) Y2[ind,] = rmvnorm(lind,Mu[,h],Si)
  }

# output
  out = list(Y1=Y1,Y2=Y2,Z=Z,H1=H1,H2=H2)
  
}