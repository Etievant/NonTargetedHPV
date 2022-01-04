## Implementation of the simulation for the randomized setting when the assumption of no treatment effect on the non-targeted types is violated (simulation of WebSection E.3.1.1).
# We only consider n = 5,000 and 10,000.
# Part 1, for beta2j = nu_1.
# Code for beta2j = nu_2 and nu_3 are given in Simul_Randomized_TreatmentEffectY2_Part2_revision.R and Simul_Randomized_TreatmentEffectY2_Part3_revision.R

# We only consider n = 5,000 and 10,000.
# Part 1, for beta2j = nu.
# Code for beta2j = - nu is given in Simul_Randomized_TreatmentEffectY2_Part2_revision.R

## Variables:
# Y1^(k): infection with HPV type k, targeted by the vaccine. We have Y1i^(k) = 1 if the ith person is infected with this type, and 0 otherwise.
# Y1c = Y1^(16) + Y1^(18) - Y1^(16)Y1^(18): binary variable indicating an infection with HPV type 16, type 18, or both. This is the primary outcome we consider.
# Y2^(j): infection with HPV type j, not targeted by the vaccine. We have Y2i^(j) = 1 if the ith person is infected with this a type, and 0 otherwise.
# Y2s = sum_j Y2^(j): categorical variable indicating the total number of infections by HPV types not targeted by the vaccine. This is the secondary outcome we consider.
# T: treatment, that is HPV vaccination. We have Ti = 1 if the ith person has been vaccinated, and 0 otherwise.
# A: covariate e.g. sexual activity.
# W: covariate e.g. age and geographical region.

## Comments:
# We assume that treatment has been randomized.
# We assume that the vaccine actually effects the non targeted HPV types.
# We consider 20 HPV types that are not targeted by the vaccine.
# The objective is to assess the performance of of method Aug, compared to UnAug, under violation of the assumption of no effect of the vaccine on the non targeted HPV types.

source("EstimationFunctions.R")
load("Parameters.RData") # file with certain of the fixed parameters
library(parallel)

## All the possible scenarios/configurations:

# Incidence of HPV types targeted by the vaccine
P_Y1  = rbind(c(0.14, 0.07), c(0.05, 0.05), c(0.032, 0.015)) 

# Different categories for A, that will induce different levels of correlation between Y1 and Y2
Aa    = rbind(c(0, 1, 2.5), c(0, 1, 2), c(0, 0.75, 2))       
# A is assumed to be categorical, with 3 categories. We assume that someone who is not sexually active can't be infected by any of the HPV types (first category is zero)

# Sample sizes
N     = c(5 * 10^3, 10^4)      

k = 1
PARAM = matrix(nrow = length(P_Y1[,1])*length(Aa[,1])*length(N), ncol = length(P_Y1[1,]) + length(Aa[1,]) + 1) 
for(l in 1:length(P_Y1[,1])){
  for(m in 1:length(Aa[,1])){
    for(i in 1:length(N)){
      PARAM[k,] = c(P_Y1[l,], Aa[m,], N[i])
      k = k+1
    }}}

## Fixed parameters:
# Incidence of the 20 non-targeted HPV types
pY2             = c(0.07, 0.03, 0.0145, 0.055, 0.115, 0.04, 0.02, 0.055, 0.065, 0.175, 0.19, 0.13, 0.095, 0.12, 0.09, 0.07, 0.14, 0.07, 0.085, 0.12)
l               = length(pY2)

# Treatment effect on each of the HPV types
beta_1.16       = - 0.73 # Effect on Y1^(16). Value taken from the Costa Rica Vaccine Trial data
beta_1.18       = - 0.86 # Effect on Y1^(18).
beta_1          = c(beta_1.16, beta_1.18)
beta_2          = c(0.01, -0.03, -0.32, -0.17, 0.04, 0.07, 0.15, 0.09, 0.18, 0.08, -0.06, 0.12, -0.02, -0.07, -0.16, -0.01, -0.03, -0.18, 0.04, -0.13) # We assume that the vaccine has an effect on HPV viruses that are not targeted by the vaccine.
# The beta2j are taken from the CVT (with incident + prevalent infections); for each non targeted type we used the estimated effect of the vaccine.

# W = (WAge, WRegion)
# We consider 13 age groups for WAge
wAge            = seq(15, 21, length.out = 13) # Different age groups, 15, ..., 21 years old
# For less strata we could use 15-18 and 18-21 age groups

# We assume that an individual can come from Region 1, Region 2 or Region 3
wRegion         = c(0,1,2)
k               = length(wRegion) # number of strata, without taking into account the age groups
PwRegion        = rep(1/k, k) # P(Wregion = wregion) = 1/3 for each region

# P(Wage = wage) for each wage in wAge. They will be "normalized" so that Sum_wage P(Wage = wage) = 1.
# We actually use slightly different probabilities for each of the k strata
PwAge           = param$PwAge
PwAge           = PwAge / rowSums(PwAge)

# Effect of W on each of the HPV types
q_WAge    = 1 / 10 * c(0.01, 0.1) # Effect of WAge on the targeted types
q_WRegion = rbind(c(-1.45, 0.06), c(0, -0.26), c(0.2, 0.5))       # Effect of WRegion on the targeted types         
s_WAge    = c(0.0035, 0.0026, 0.0071, 0.0156, 0.0004, 0.0090, 0.0073, 0.0078, 0.0054, 0.0015, 0.0082, 0.0036, 0.0085, 0.0052, 
                    0.0077, 0.0120, 0.0112, 0.0143, 0.0011, 0.0019) # Effect of WAge on each of the non targeted HPV types.
s_WRegion = c(-0.2504, -0.1048, -0.0994, -0.3612, -0.1164, -0.2218, -0.2030, -0.0325, 0.1126,  0.3296, -0.1547,  0.3212, 
                    -0.2316, 0.1313,  0.5098, -0.0070, -0.1339, -0.0015,  0.3554, -0.2277) # Effect of WRegion on each of the non targeted HPV types.

# A can take 3 different values: a_low, a_medium and a_high. The associated probabilities depends on the value of W
# The following probabilities will be "normalized" so that P(A = a_low | W = w) + P(A = a_medium | W = w) 
# + P(A = a_high | W = w) = 1, in each strata of W.
Pahigh    = param$Pahigh
Pamedium  = param$Pamedium 
Palow     = param$Palow 
Ptot      = Pahigh + Pamedium + Palow
Pahigh    = Pahigh / Ptot
Pamedium  = Pamedium / Ptot
Palow     = Palow / Ptot

# Number of data sets replications
Nreplic   = 10^4

# Function to run for each scenario/configuration
Onerun = function(p){
  
  set.seed(1234)
  # Considered scenario
  pY1             = PARAM[p, 1:2]
  a               = PARAM[p, 3:5]
  a_low           = a[1]
  a_medium        = a[2]
  a_high          = a[3]
  n               = PARAM[p, 6]
  
  # Generation of WRegion of each of the n individuals.
  WWRegion = sample(wRegion, n, replace = TRUE, prob = PwRegion)        
  
  # Generation of WAge of each of the n individuals.
  WWAge = rep(NA, n)
  for(i in 1:length(wRegion)){
    WWAge[which(WWRegion == wRegion[i])] = sample(wAge, sum(WWRegion == wRegion[i]), replace = TRUE, prob = PwAge[i,])
  }
  
  ## Computation of the data-set characteristics
  # Mean of A, E(A)
  meanA = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        meanA = meanA + a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        meanA = meanA + a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        meanA = meanA + a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  
  # Variance of A
  varA = - meanA^2
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        varA = varA + a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        varA = varA + a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        varA = varA + a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  
  # So that P(Y1 = 1) = pY1, for each targeted type
  meanY1_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        pt = 1 / 2
        
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
      }}}
  
  # So that P(Y2 = 1) = pY2, for each non-targeted type
  meanY2_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_2) * pt + 1 - pt )
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_2) * pt + 1 - pt )
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_2) * pt + 1 - pt )
      }}}
  
  # Computation of the intercepts, so that E(Y1) = pY1 and E(Y2) = pY2 for each HPV type.
  mu_1            = log(pY1 / meanY1_part1) # so that E(Y1) = pY1, for each targeted type.
  mu_2            = log(pY2 / meanY2_part1) # so that E(Y2) = pY2, for each non-targeted type.
  Mu_1            = matrix(rep(mu_1,n), nrow = n, byrow = TRUE)
  Mu_2            = matrix(rep(mu_2,n), nrow = n, byrow = TRUE)
  
  # Mean of the composite primary outcome, E(Y1c)
  meanY1c = sum(pY1)
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        pt = 1 / 2
        
        meanY1c = meanY1c - exp(sum(mu_1)) * exp(sum(wAge[j] %*% t(q_WAge)) + sum(t(q_WRegion[i,]))) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * (exp(sum(beta_1)) * pt + 1 - pt )
        meanY1c = meanY1c - exp(sum(mu_1)) * exp(sum(wAge[j] %*% t(q_WAge)) + sum(t(q_WRegion[i,]))) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * (exp(sum(beta_1)) * pt + 1 - pt )
        meanY1c = meanY1c - exp(sum(mu_1)) * exp(sum(wAge[j] %*% t(q_WAge)) + sum(t(q_WRegion[i,]))) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * (exp(sum(beta_1)) * pt + 1 - pt )
      }}}
  
  # Variance of the composite primary outcome
  varY1c = meanY1c * (1 - meanY1c)
  
  # Mean of the secondary outcome
  meanY2s = sum(pY2) # E(Y2s)
  
  # Variance of the secondary outcome
  covY2 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        for(t in 1:0){
          covY2 = covY2 + t(exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]  * ((1 - pt) * (1 - t) + pt * t)
          covY2 = covY2 + t(exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
          covY2 = covY2 + t(exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
        }}}}
  varY2s = meanY2s * (1 - meanY2s) + (sum(covY2 - diag(diag(covY2))))
  
  # Covariance between primary composite and secondary outcomes
  covY1cY2s_part1 = 0
  covY1cY2s_part2 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        for(t in 0:1){
          
          pt = 1 / 2
          
          covY1cY2s_part1 = covY1cY2s_part1 + t(exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% (exp(mu_1 + t %*% t(beta_1) + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
          covY1cY2s_part1 = covY1cY2s_part1 + t(exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% (exp(mu_1 + t %*% t(beta_1) + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
          covY1cY2s_part1 = covY1cY2s_part1 + t(exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% (exp(mu_1 + t %*% t(beta_1) + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
          
          covY1cY2s_part2 = covY1cY2s_part2 - (exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * exp(sum(mu_1 + t %*% t(beta_1) + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^3 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
          covY1cY2s_part2 = covY1cY2s_part2 - (exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * exp(sum(mu_1 + t %*% t(beta_1) + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^3 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
          covY1cY2s_part2 = covY1cY2s_part2 - (exp(mu_2 + t %*% t(beta_2) + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * exp(sum(mu_1 + t %*% t(beta_1) + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^3 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * ((1 - pt) * (1 - t) + pt * t)
        }}}}
  covY1cY2s = sum(covY1cY2s_part1) + sum(covY1cY2s_part2) - meanY1c * meanY2s # cov(Y1c, Y2s)
  
  # Correlation between primary and secondary outcomes
  corrY1cY2s = covY1cY2s / sqrt(varY1c * varY2s)
  
  # Correlation between primary and secondary outcomes
  corrY1cY2s = covY1cY2s / sqrt(varY1c * varY2s)
  
  # Parameter of interest
  EY1c_T1_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY1c_T1_part1 = EY1c_T1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T1_part1 = EY1c_T1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T1_part1 = EY1c_T1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  
  EY1c_T1_part2 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY1c_T1_part2 = EY1c_T1_part2 + exp(sum(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T1_part2 = EY1c_T1_part2 + exp(sum(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T1_part2 = EY1c_T1_part2 + exp(sum(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  
  EY1c_T1 = sum(EY1c_T1_part1 * exp(mu_1 + beta_1)) - EY1c_T1_part2 * exp(sum(mu_1 + beta_1))
  
  EY1c_T0_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY1c_T0_part1 = EY1c_T0_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T0_part1 = EY1c_T0_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T0_part1 = EY1c_T0_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  
  EY1c_T0_part2 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY1c_T0_part2 = EY1c_T0_part2 + exp(sum(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T0_part2 = EY1c_T0_part2 + exp(sum(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY1c_T0_part2 = EY1c_T0_part2 + exp(sum(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  
  EY1c_T0       = sum(EY1c_T0_part1 * exp(mu_1)) - EY1c_T0_part2 * exp(sum(mu_1))
  
  beta_1.true   = log(EY1c_T1 / EY1c_T0) # Parameter of interest
  
  res           = NULL
  
  for(nrep in 1:Nreplic){
    
    ## Data generation
    sample  = sample(n, n)
    
    WRegion = WWRegion[sample]
    WAge = WWAge[sample]
    
    A = rep(NA, n)
    for(i in 1:length(wRegion)){
      for(j in 1:length(wAge)){
        A[which((WRegion == wRegion[i])&(WAge == wAge[j]))] = sample(c(a_low, a_medium, a_high), size = sum((WRegion == wRegion[i])&(WAge == wAge[j])), replace = TRUE, prob = c(Palow[i,j], Pamedium[i,j], Pahigh[i,j]))
      }}
    
    T = rbinom(n, size = 1, prob = 1 / 2)
    
    WRegion_Y1 = 0
    for(i in 1:length(wRegion)){
      WRegion_Y1 = WRegion_Y1 + 1 * (WRegion == wRegion[i]) %*% t(q_WRegion[i, ])
    }
    
    p1              = exp(Mu_1 + T %*% t(beta_1) + WAge %*% t(q_WAge) + WRegion_Y1)
    PY1             = A * p1
    Y1              = matrix(rbinom(2*n, size = 1, prob = PY1), nrow = n, byrow = FALSE) # Generation of the infections with the targeted types for each of the n individuals.
    Y1c             = 1 * (rowSums(Y1) > 0) # Primary composite outcome
    
    p2              = exp(Mu_2 + T %*% t(beta_2) + WAge %*% t(s_WAge) + WRegion %*% t(s_WRegion))
    PY2             = A * p2 
    Y2              = matrix(rbinom(l*n, size = 1, prob = PY2), nrow = n, byrow = FALSE)  # Generation of the infections with the non targeted types for each of the n individuals.
    Y2s             = rowSums(Y2)   # Secondary outcome, counting the total number of infections with non-targeted types.
    
    W = interaction(WAge, WRegion, sep = "_")
    
    ## Estimation of the parameters
    
    # Method Aug_W (Aug with observed covariates only, as proposed by Zhang, Tsiatis and Davidian, 2008)
    est.Aug_W = AugmentedEE(Y1 = Y1c, T = T, Z = as.data.frame(cbind(WAge, WRegion)))
    CI.Aug_W = cbind(est.Aug_W$beta_1.hat - est.Aug_W$se.beta_1.hat * qnorm(0.975), est.Aug_W$beta_1.hat + est.Aug_W$se.beta_1.hat * qnorm(0.975))
    Aug_W = cbind("Aug_W", est.Aug_W$beta_1.hat, est.Aug_W$se.beta_1.hat, CI.Aug_W)
    
    # Method Aug (with secondary outcome only)
    est.Aug_Y2 = AugmentedEE(Y1 = Y1c, T = T, Z = as.data.frame(Y2s))
    CI.Aug_Y2 = cbind(est.Aug_Y2$beta_1.hat - est.Aug_Y2$se.beta_1.hat * qnorm(0.975), est.Aug_Y2$beta_1.hat + est.Aug_Y2$se.beta_1.hat * qnorm(0.975))
    Aug_Y2 = cbind("Aug_Y2", est.Aug_Y2$beta_1.hat, est.Aug_Y2$se.beta_1.hat, CI.Aug_Y2)
    
    # Method Aug_WY2 (with secondary outcome and observed covariates)
    est.Aug_WY2 = AugmentedEE(Y1 = Y1c, T = T, Z = as.data.frame(cbind(WAge, WRegion, Y2s)))
    CI.Aug_WY2 = cbind(est.Aug_WY2$beta_1.hat - est.Aug_WY2$se.beta_1.hat * qnorm(0.975), est.Aug_WY2$beta_1.hat + est.Aug_WY2$se.beta_1.hat * qnorm(0.975))
    Aug_WY2 = cbind("Aug_WY2", est.Aug_WY2$beta_1.hat, est.Aug_WY2$se.beta_1.hat, CI.Aug_WY2)
    
    # Method UnAug
    est.UnAug = NaiveEE(Y1 = Y1c, T = T)
    CI.UnAug = cbind(est.UnAug$beta_1.hat - est.UnAug$se.beta_1.hat * qnorm(0.975), est.UnAug$beta_1.hat + est.UnAug$se.beta_1.hat * qnorm(0.975))
    UnAug = cbind("UnAug", est.UnAug$beta_1.hat, est.UnAug$se.beta_1.hat, CI.UnAug)
    
    recap = rbind(UnAug, Aug_W, Aug_Y2, Aug_WY2)     
    colnames(recap) = c("Approach", "beta_1.hat", "se.beta_1.hat", "CI.left", "CI.right")
    
    res = rbind(res, cbind(recap, n = n, pY1 = paste(pY1, collapse = "_"), pY2 = paste(pY2, collapse = "_"), A = paste(a, collapse = "_"), beta_1.true = beta_1.true, corrY1cY2s = corrY1cY2s, varA = varA))
  }
  
  myfile  = paste0("RES_Randomized_TreatmentEffectY2_Part1-n", n, "-pY1", paste(pY1, collapse = "_"), "-A", paste(a, collapse = "_"), "-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData")
  save(res, file = myfile)
}

P = nrow(PARAM)
resultat = mclapply(1:P, Onerun, mc.cores = 18)

# Results for the randomized settings
RES = NULL
for(p in 1: nrow(PARAM)){
  
  pY1             = PARAM[p, 1:2]
  a               = PARAM[p, 3:5]
  n               = PARAM[p, 6]
  
  load(paste0("RES_Randomized_TreatmentEffectY2_Part1-n", n, "-pY1", paste(pY1, collapse = "_"), "-A", paste(a, collapse = "_"), "-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData"))
  RES = rbind(RES, res)
}

RECAP = as.data.frame(RES)
RECAP$beta_1.hat = as.numeric(RECAP$beta_1.hat)
RECAP$Approach = as.factor(RECAP$Approach)
RECAP$Approach = factor(RECAP$Approach, labels = c("Aug_W", "Aug_WY2", "Aug_Y2", "UnAug"))
RECAP$se.beta_1.hat = as.numeric(RECAP$se.beta_1.hat)
RECAP$beta_1.true = as.numeric(RECAP$beta_1.true)
RECAP$pY1 = as.factor(RECAP$pY1)
RECAP$pY1 = factor(RECAP$pY1, labels = c(expression(P[Y[1]^(16) == 1]==~0.032~","~P[Y[1]^(18) == 1]==~0.015), expression(P[Y[1]^(16) == 1]==~0.05~","~P[Y[1]^(18) == 1]==~0.05), expression(P[Y[1]^(16) == 1]==~0.14~","~P[Y[1]^(18) == 1]==~0.07) ))
RECAP$A = as.factor(RECAP$A)
RECAP$A = factor(RECAP$A, labels = c(expression(a[low]==~0~","~a[medium]==~0.75~","~a[high]==~2), expression(a[low]==~0~","~a[medium]==~1~","~a[high]==~2), expression(a[low]==~0~","~a[medium]==~1~","~a[high]==~2.5) ))
RECAP$corrY1cY2s = as.numeric(RECAP$corrY1cY2s)
RECAP$CI.left = as.numeric(RECAP$CI.left)
RECAP$CI.right = as.numeric(RECAP$CI.right)
RECAP$n = as.factor(RECAP$n)
RECAP$n = factor(RECAP$n,levels(RECAP$n)[c(2,1)])

myfile  = paste0("RECAP_Randomized_TreatmentEffectY2_Part1-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData")
save(RECAP, file = myfile)

library(ggplot2)
library(gtable)
library(grid)
library(xtable)
# with more readable facets names
RECAP$pY1 = factor(RECAP$pY1, labels = c("Low", "Medium", "High"))
RECAP$A = factor(RECAP$A, labels = c("Small", "Medium", "Large"))
plot = ggplot(RECAP, aes(x = n, y = beta_1.hat, color = Approach)) + geom_boxplot(coef = NULL) + geom_hline(aes(yintercept = beta_1.true)) + theme_light() + theme(plot.title = element_text(size = 11), axis.text = element_text(size = 11), axis.title = element_text(size = 11), legend.text = element_text(size = 11), strip.background = element_rect(color="black", fill="white", size = 0.5, linetype="solid"), strip.text.x = element_text(size = 11, color = "black"), strip.text.y = element_text(size = 11, color = "black")) + ylab((expression(hat(beta[1]^C) ))) + facet_grid(pY1~A, labeller = label_parsed, scales = "free_y") + scale_color_manual(values=c("#932d48","#015fc6", "#ea95ff", "#53b31f"), labels = c( expression(Aug[W]), expression(Aug[WY2]), "Aug", "UnAug"))
labelT = "Variance of the unmeasured covariate A"
labelR = "Risk of infections with targeted types 16 and 18"
# Get the ggplot grob
z = ggplotGrob(plot)
# Get the positions of the strips in the gtable: t = top, l = left, ...
posR = subset(z$layout, grepl("strip-r", name), select = t:r)
posT = subset(z$layout, grepl("strip-t", name), select = t:r)
# Add a new column to the right of current right strips, 
# and a new row on top of current top strips
width = z$widths[max(posR$r)]    # width of current right strips
height = z$heights[min(posT$t)]  # height of current top strips
z = gtable_add_cols(z, width, max(posR$r))  
z = gtable_add_rows(z, height, min(posT$t)-1)
# Construct the new strip grobs
stripR = gTree(name = "Strip_right", children = gList(
  rectGrob(gp = gpar(color = "black", fill = "white")),
  textGrob(labelR, rot = -90, gp = gpar(fontsize = 11, col = "black"))))
stripT = gTree(name = "Strip_top", children = gList(
  rectGrob(gp = gpar(color = "black", fill = "white")),
  textGrob(labelT, gp = gpar(fontsize = 11, col = "black"))))
# Position the grobs in the gtable
z = gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
z = gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
# Add small gaps between strips
z = gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
z = gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
grid.newpage()
grid.draw(z)
# Save the plot
pdf(paste0("Comparison_Randomized_TreatmentEffectY2_Part1-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".pdf"),  width = 10, height = 7)
grid.draw(z) # print it
dev.off() # Stop writing to the PDF file

# Computation of beta_2, for all the considered scenarios, as defined in Web Section D.1
BETA_2 = NULL
for(p in 1:nrow(PARAM)){
  pY1             = PARAM[p, 1:2]
  a               = PARAM[p, 3:5]
  a_low           = a[1]
  a_medium        = a[2]
  a_high          = a[3]
  n               = PARAM[p, 6]
  
  ## Computation of the data-set characteristics
  # Mean of A, E(A)
  meanA = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        meanA = meanA + a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        meanA = meanA + a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        meanA = meanA + a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  # Variance of A
  varA = - meanA^2
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        varA = varA + a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        varA = varA + a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        varA = varA + a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  # So that P(Y1 = 1) = pY1, for each targeted type
  meanY1_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        pt = 1 / 2 
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
      }}}
  # So that P(Y2 = 1) = pY2, for each non-targeted type
  meanY2_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_2) * pt + 1 - pt )
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_2) * pt + 1 - pt )
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_2) * pt + 1 - pt )
      }}}
  # Computation of the intercepts, so that E(Y1) = pY1 and E(Y2) = pY2 for each HPV type.
  mu_1            = log(pY1 / meanY1_part1) # so that E(Y1) = pY1, for each targeted type.
  mu_2            = log(pY2 / meanY2_part1) # so that E(Y2) = pY2, for each non-targeted type.
  # True beta_2
  EY2c_T1_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY2c_T1_part1 = EY2c_T1_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T1_part1 = EY2c_T1_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T1_part1 = EY2c_T1_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  EY2c_T1_part2 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY2c_T1_part2 = EY2c_T1_part2 + exp(sum(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T1_part2 = EY2c_T1_part2 + exp(sum(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T1_part2 = EY2c_T1_part2 + exp(sum(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  EY2c_T1 = sum(EY2c_T1_part1 * exp(mu_2 + beta_2)) - EY2c_T1_part2 * exp(sum(mu_2 + beta_2))
  EY2c_T0_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY2c_T0_part1 = EY2c_T0_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T0_part1 = EY2c_T0_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T0_part1 = EY2c_T0_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  EY2c_T0_part2 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        EY2c_T0_part2 = EY2c_T0_part2 + exp(sum(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T0_part2 = EY2c_T0_part2 + exp(sum(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        EY2c_T0_part2 = EY2c_T0_part2 + exp(sum(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  EY2c_T0       = sum(EY2c_T0_part1 * exp(mu_2)) - EY2c_T0_part2 * exp(sum(mu_2))
  beta_2.true = log(EY2c_T1 / EY2c_T0)
  BETA_2 = cbind(BETA_2, beta_2.true)
}

## Details of the results
Nreplic   = 10^4
Eff = NULL
for(i in 1:nrow(PARAM)){
  
  RECAP1 = RECAP[((i-1)*(4*Nreplic) + 1):(i*(4*Nreplic)),]
  
  beta_1 = RECAP1$beta_1.true[1]
  
  # Coverage of the confidence intervals (for the methods which return estimates of the variance)
  cov.UnAug   = sum((RECAP1[which(RECAP1$Approach == "UnAug"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "UnAug"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "UnAug"),5])
  cov.Aug_W   = sum((RECAP1[which(RECAP1$Approach == "Aug_W"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Aug_W"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Aug_W"),5])
  cov.Aug_Y2  = sum((RECAP1[which(RECAP1$Approach == "Aug_Y2"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Aug_Y2"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Aug_Y2"),5])
  cov.Aug_WY2 = sum((RECAP1[which(RECAP1$Approach == "Aug_WY2"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Aug_WY2"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Aug_WY2"),5])
  
  sd.UnAug    = sd(RECAP1[which(RECAP1$Approach == "UnAug"),2])         # empirical standard deviation for hat beta1 estimated with the UnAug approach.
  sd.Aug_W    = sd(RECAP1[which(RECAP1$Approach == "Aug_W"),2])         # empirical standard deviation for hat beta1 estimated with the augmented approach (with observed covariates only).
  sd.Aug_Y2   = sd(RECAP1[which(RECAP1$Approach == "Aug_Y2"),2])        # empirical standard deviation for hat beta1 estimated with the augmented approach (with secondary outcome).
  sd.Aug_WY2  = sd(RECAP1[which(RECAP1$Approach == "Aug_WY2"),2])       # empirical standard deviation for hat beta1 estimated with the augmented approach (with secondary outcome and observed covariates).
  
  sandwich_sd.UnAug   = mean(RECAP1[which(RECAP1$Approach == "UnAug"),3]) # mean of the sandwich standard deviations for hat beta1, with the UnAug approach.
  sandwich_sd.Aug_W   = mean(RECAP1[which(RECAP1$Approach == "Aug_W"),3]) # mean of the sandwich standard deviations for hat beta1, with the augmented approach (with observed covariates only).
  sandwich_sd.Aug_Y2  = mean(RECAP1[which(RECAP1$Approach == "Aug_Y2"),3]) # mean of the sandwich standard deviations for hat beta1, with the augmented approach (with secondary outcome).
  sandwich_sd.Aug_WY2 = mean(RECAP1[which(RECAP1$Approach == "Aug_WY2"),3]) # mean of the sandwich standard deviations for hat beta1, with the augmented approach (with secondary outcome and observed covariates).
  
  mean.UnAug    = mean(RECAP1[which(RECAP1$Approach == "UnAug"),2]) # mean of hat beta1, when estimated with the UnAug approach.
  MSE.UnAug     = mean((RECAP1[which(RECAP1$Approach == "UnAug"),2] - beta_1)^2) # MSE, when beta1 is estimated with the UnAug approach.
  
  mean.Aug_W    = mean(RECAP1[which(RECAP1$Approach == "Aug_W"),2]) # mean of hat beta1, when estimated with the augmented approach (with observed covariates only).
  MSE.Aug_W     = mean((RECAP1[which(RECAP1$Approach == "Aug_W"),2] - beta_1)^2) # MSE, when beta1 is estimated with the augmented approach (with observed covariates only).
  
  mean.Aug_Y2   = mean(RECAP1[which(RECAP1$Approach == "Aug_Y2"),2]) # mean of hat beta1, when estimated with the augmented approach (with secondary outcome).
  MSE.Aug_Y2    = mean((RECAP1[which(RECAP1$Approach == "Aug_Y2"),2] - beta_1)^2) # MSE, when beta1 is estimated with the augmented approach (with secondary outcome).
  
  mean.Aug_WY2  = mean(RECAP1[which(RECAP1$Approach == "Aug_WY2"),2]) # mean of hat beta1, when estimated with the augmented approach (with secondary outcome and observed covariates).
  MSE.Aug_WY2   = mean((RECAP1[which(RECAP1$Approach == "Aug_WY2"),2] - beta_1)^2) # MSE, when beta1 is estimated with the augmented approach (with secondary outcome and observed covariates).
  
  eff_UnAug_Aug_W     = MSE.UnAug / MSE.Aug_W # ratio of MSE for the UnAug and augmented approach (with observed covariates only)
  eff_UnAug_Aug_Y2    = MSE.UnAug / MSE.Aug_Y2 # ratio of MSE for the UnAug and augmented approach (with secondary outcome)
  eff_UnAug_Aug_WY2   = MSE.UnAug / MSE.Aug_WY2 # ratio of MSE for the UnAug and augmented approach (with secondary outcome and observed covariates)
  
  Eff = rbind(Eff, c(eff_beta1_Aug_W = eff_UnAug_Aug_W, eff_beta1_Aug_Y2 = eff_UnAug_Aug_Y2, eff_beta1_Aug_WY2 = eff_UnAug_Aug_WY2, bias.UnAug = mean.UnAug - beta_1, bias.Aug_W = mean.Aug_W - beta_1, bias.Aug_Y2 = mean.Aug_Y2 - beta_1, bias.Aug_WY2 = mean.Aug_WY2 - beta_1, empir_sd.UnAug = sd.UnAug, sandwich_sd.UnAug = sandwich_sd.UnAug, empir_sd.Aug_W = sd.Aug_W, sandwich_sd.Aug_W = sandwich_sd.Aug_W, empir_sd.Aug_Y2 = sd.Aug_Y2, sandwich_sd.Aug_Y2 = sandwich_sd.Aug_Y2, empir_sd.Aug_WY2 = sd.Aug_WY2, sandwich_sd.Aug_WY2 = sandwich_sd.Aug_WY2,  CIcov.UnAug = cov.UnAug, CIcov.Aug_W = cov.Aug_W, CIcov.Aug_Y2 = cov.Aug_Y2, CIcov.Aug_WY2 = cov.Aug_WY2, n = as.character(RECAP1[1,]$n), pY1 = as.character(RECAP1[1,]$pY1), A =  as.character(RECAP1[1,]$A), beta_1 = beta_1, corr = RECAP1$corrY1cY2s[1] ) )
}
Eff = as.data.frame(Eff)
Eff = cbind(Eff, beta_2.true = c(BETA_2))
ColNames = colnames(Eff[,c(1:20, 23:25)])
Eff[ColNames] = sapply(Eff[ColNames], as.numeric)
#print(xtable(Eff, type = "latex", digits = 3), include.rownames=FALSE) # LaTeX table
beta_1          = c(beta_1.16, beta_1.18)
myfile  = paste0("Eff_Randomized_TreatmentEffectY2_Part1-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData") # save as Rdata
save(Eff, file = myfile)
Eff[ColNames] = round(Eff[ColNames], digits = 3)
Eff = Eff[which(Eff$n == 10000),] # save only the scenarios with n = 10,000
#Eff = Eff[-which(Eff$A == "a[low] == ~0 ~ \",\" ~ a[medium] == ~1 ~ \",\" ~ a[high] == ~2.5" ),]
#Eff = Eff[-which(Eff$pY1 == "P[Y[1]^(16) == 1] == ~0.05 ~ \",\" ~ P[Y[1]^(18) == 1] == ~0.1"),]
#Eff = Eff[-which(Eff$pY1 == "P[Y[1]^(16) == 1] == ~0.1 ~ \",\" ~ P[Y[1]^(18) == 1] == ~0.07"),]
write.csv(Eff, file = paste0("Eff_Randomized_TreatmentEffectY2_Part1-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".csv")) # save as csv

