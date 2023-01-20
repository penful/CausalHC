#### 
# Example Code to estimate the models proposed  in the paper titled 
# "A causal latent transition model with multivariate outcomes and unobserved heterogeneity: 
# Application  to human capital development" by 
# - F. Bartolucci (University of Perugia, IT)
# - F. Pennoni (University of Milano-Bicocca, IT)
# - G. Vittadini (University of Milano-Bicocca, IT)
####

#### load pseudo data as an example  ####
rm(list=ls())
if(!("LMest"%in%installed.packages())) install.packages("LMest")
if(!("mvtnorm"%in%installed.packages())) install.packages("mvtnorm")
#
require(LMest)
require(mvtnorm)
#
load("data.Rdata")

# id: unit identifier
# time: time period (1, first period, 2, second period)
# TT: treatment
# RAI: school motivations
# ec2: parents escs index (index of economic, social and cultural status)
# q:  quality of class relations
# s: external support
# b: well-being at school
# m: discomfort at school
# b1: bulling acting
# b2: bullying right away
# gen: gender
# npita: Italian nationality of the father
# pscore: employment status of the father (continous)
# pscoremiss: dummy variable for missing employment status of the parents 
# Y0IT: achievement score in Italian in t=0
# Y0M: achievement score in Maths in t=0
# Y1: achievement score in Italian in t=1
# Y2: achievement score in Maths in t=1
#
####

### Code employed to estimate the proposed models ###
#### Causal Latent Transition (CLT) Model ####
set.seed(143)
out3Y <- lmestCont(responsesFormula = Y1 +Y2 ~ NULL,
                   latentFormula = ~  
                     RAI + ec2 + q + s + b + 
                     m + b1 + b2 + gen + npita + pscore + +  
                     pscoremiss |
                     TT +  ec2 + gen + npita +  pscore +  
                     pscoremiss  +  Y0IT + Y0M,
                     index=c("id","time"),
                     k = 2, 
                     data = data, 
                     out_se = TRUE, output = TRUE,
                     tol=10^-12)

#### Results of the CLT model ####
summary(out3Y)

# Estimated conditional means 
out3Y$Mu

# Estimated variance-covariance matrix
out3Y$Si

# Estimated initial probabilities 
apply(out3Y$Piv,2,mean)

# Estimated averaged transition matrix
apply(out3Y$PI[, , , 2], c(1, 2), mean)

# Results on the initial probabilities
# Estimated regression coefficients for the covariates (Be) 
# and standard errors (seBe)
TabBeYY <-  cbind(out3Y$Be,out3Y$seBe)
colnames(TabBeYY)<-c("Be", "seBe")
round(TabBeYY,3)

# Results on the transition probabilities
# Estimated regression coefficients for the covariates (Ga) 
TabGa1Y <- round(cbind(out3Y$Ga,out3Y$seGa),3)
colnames(TabGa1Y) <- c("Ga1", "seGa1", "Ga2", "seGa2")

# Estimated treatment effect on the probability to transit from the first 
# to the second latent state (Ga1) and related standard error (seGa1)
# Estimated treatment effect on the probability to transit from the second 
# to the third latent state (Ga2) and related standard error (seGa2)
round(TabGa1Y[2:2,],1)
#

# Density plots

# contour plot of the joint density
plot(out3Y, "density")

# Contour plot of the marginal densities
plot(out3Y,what="density",components=c(1,2))


#### Difference in Difference estimations ####

#### DiD models  ####

## Model 1 (all students) ITALIAN score 
#
DID2IT5 <- lm(Y1 ~ TT + RAI + ec2 + q + s + b + 
               m + b1 + b2 + gen + npita + pscore +
               pscoremiss,
               data=data[data$time==1,])
summary(DID2IT5)

## DiD models according to the students scores ##
# selection of the first time period
dd <- data[data$time==1,]
# selection of students  performing below the median value for Italian
ind1 <- which(dd$Y0IT <= median(dd$Y0IT))

# DiD model for the Italian score of  best performing students
# (Table 8, top panel, third column)
#
DID2IT5b <- lm(Y1 ~ TT + RAI +  ec2 + q + s + b + 
                 m + b1 + b2 + gen + npita + pscore +
                 pscoremiss ,
                 data=dd[-ind1,])
summary(DID2IT5b)

# DiD model for the Italian score of worst performing students
DID2IT5w <- lm(Y1 ~ TT + RAI +   ec2 + q + s + b+ 
                 m + b1 + b2+ gen + npita + pscore +
                 pscoremiss,
                 data=dd[ind1,])
summary(DID2IT5w)

# DiD model for the Math score of all students
#
DID2M5 <- lm(Y2 ~ TT + RAI +  ec2 + q + s + b+ 
               m + b1 + b2 + gen + npita + pscore +
               pscoremiss,
               data=data[data$time==1,])
summary(DID2M5)

# Selection of students performing below the median value for maths
ind2 <- which(dd$Y0M<=median(dd$Y0M))
#
# DiD model for the Math score of  best performing students
# (Table 9, top panel, second column)
#
DID2M5b <- lm(Y2 ~ TT + RAI + ec2 + q + s + b + 
                m + b1 + b2 + gen + npita + pscore +
                pscoremiss,
                data=dd[-ind2,])
summary(DID2M5b)

# DiD model for the Math score of  worst performing students
# (Table 9, top panel, third column)
#
DID2M5w <- lm(Y2 ~ TT + RAI + ec2+ q+ s+ b+ 
                m + b1 + b2 + gen + npita + pscore +
                pscoremiss ,
                data=dd[ind2,])
summary(DID2M5w)

#### DiD EQUATION  #### 
# with covariates for the second time occasion
#
# calculate differences
data$Yd1 <- data$Y1-data$Y0IT
data$Yd2 <- data$Y2-data$Y0M

# DiD EQUATION (15) Model 1 (all) Italian
# (Table 8, bottom panel)
#
DID2IT <- lm(Yd1 ~ TT + ec2 + gen +       
              npita +  pscore + pscoremiss  + Y0IT,
              data=data[data$time==2,])
summary(DID2IT)
 
# DiD EQUATION (15) Model 2 (best) Italian
# (Table 8, bottom panel)
#
DL2 <- data[-ind1,]
DL2$Yd1 <- DL2$Y1-DL2$Y0IT
DL2$Yd2 <- DL2$Y2-DL2$Y0M
#
DID2ITb <- lm(Yd1 ~ TT + ec2 + gen + npita +  pscore +  
               pscoremiss  + Y0IT,
               data=DL2[DL2$time==2,])
summary(DID2ITb)


# DiD EQUATION (15) Model 3 (worst) Italian
# (Table 8, bottom panel)
#
dd <- data[data$time==1,]
ind1 <- which(dd$Y0IT<= median(dd$Y0IT))
DL1 <- data[ind1,]

DL1$Yd1 <- DL1$Y1-DL1$Y0IT
DID2ITw <- lm(Yd1 ~ TT + ec2 + gen + npita +  pscore +  
              pscoremiss  + Y0IT,
              data=DL1[DL1$time==2,])
summary(DID2ITw)



# DiD EQUATION (15) Model 1 (all) Maths
# (Table 9, bottom panel)
#
DID2M <- lm(Yd2 ~ TT + ec2 + gen + npita +  pscore +  
           pscoremiss  +  Y0M,
           data=data[data$time==2,])

summary(DID2M)

# DiD EQUATION (15) Model 2 (best) Maths
# (Table 9, bottom panel)
dd <- data[data$time==1,]
ind3 <- which(dd$Y0M<=median(dd$Y0M))
DL3 <- data[ind3,]

DL4 <- data[-ind3,]
DID2Mb <- lm(Yd2 ~ TT + ec2 + gen + npita +  pscore +  
             pscoremiss  +  Y0M,
             data=DL4[DL4$time==2,])
summary(DID2Mb)

# DiD EQUATION (15) Model 3 (worse) Maths
# (Table 9, bottom panel)

dd <- data[data$time==1,]
ind3 <- which(dd$Y0M<=median(dd$Y0M))
DL3 <- data[ind3,]

DL3$Yd2 <- DL3$Y2-DL3$Y0M
DID2Mw <- lm(Yd2 ~ TT + ec2 + gen + npita +  pscore +  
            pscoremiss  + Y0M,
            data=DL3[DL3$time==2,])
summary(DID2Mw)

# Double robust DID model 
# Italian score for all studends
library(DRDID) 
data$TT <-as.numeric(data$TT)
data$time <-as.numeric(data$time)
outIT1 <- drdid(yname = "Y1", tname = "time", idname = "id", 
               dname = "TT",
               xformla = ~ ec2 + gen + npita + pscore +  
               pscoremiss  + Y0IT,
               data = data, panel = TRUE)
summary(outIT1)


