## Implementation of the simulation for the non-randomized setting when the observed covariate only has a few large strata  (simulation of WebSection E.2.2).
# We only consider n = 5,000 and 10,000.

## Variables:
# Y1^(k): infection with HPV type k, targeted by the vaccine. We have Y1i^(k) = 1 if the ith person is infected with this type, and 0 otherwise.
# Y1c = Y1^(16) + Y1^(18) - Y1^(16)Y1^(18): binary variable indicating an infection with HPV type 16, type 18, or both. This is the primary outcome we consider.
# Y2^(j): infection with HPV type j, not targeted by the vaccine. We have Y2i^(j) = 1 if the ith person is infected with this a type, and 0 otherwise.
# Y2s = sum_j Y2^(j): categorical variable indicating the total number of infections by HPV types not targeted by the vaccine. This is the secondary outcome we consider.
# T: treatment, that is HPV vaccination. We have Ti = 1 if the ith person has been vaccinated, and 0 otherwise.
# A: unobserved confounder e.g. sexual activity.
# W: observed confounder e.g. age and geographical region.

## Comments:
# We assume that treatment has not been randomized.
# We assume that some confounders for the relationship between treatment and primary outcome are observed.
# We consider that the observed covariate only has a few large strata.
# We assume that other confounders have not been observed.
# We assume that the vaccine does not have an effect on the HPV viruses that are not targeted by the vaccine.
# We consider 20 HPV types that are not targeted by the vaccine.
# The objective is to assess the performance of our proposed approaches (Joint-Reg, Joint-MH and SS-Joint), compared to more naive approaches (e.g MH).

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
beta_2          = 0 # We assume that the vaccine does not have an effect on HPV viruses that are not targeted by the vaccine.

# W = (WAge, WRegion)
# We consider 2 categories for WAge
wAge            = c(mean(seq(15,18,by = 0.5)), mean(seq(18.5,21,by = 0.5))) # Different age groups: 15-18 and 18-21 age groups.

# We assume that an individual can come from Region 1, Region 2 or Region 3
wRegion         = c(0,1,2)
k               = length(wRegion) # number of strata, without taking into account the age groups
PwRegion        = rep(1/k, k) # P(Wregion = wregion) = 1/3 for each region

# P(Wage = wage) for each wage in wAge. They will be "normalized" so that Sum_wage P(Wage = wage) = 1.
# We actually use slightly different probabilities for each of the k strata
PwAge           = cbind(rowSums(param$PwAge[,1:7]), rowSums(param$PwAge[,8:13]))
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
Pahigh    = cbind(rowSums(param$Pahigh[,1:7]), rowSums(param$Pahigh[,8:13]))
Pamedium  = cbind(rowSums(param$Pamedium[,1:7]), rowSums(param$Pamedium[,8:13]))
Palow     = cbind(rowSums(param$Palow[,1:7]), rowSums(param$Palow[,8:13]))
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
  
  # Generation of the region of each of the n individuals.
  WWRegion = sample(wRegion, n, replace = TRUE, prob = PwRegion)        
  
  # Generation of the age of each of the n individuals.
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
  
  # So that P(T = 1) = 0.5
  alpha_T = - 1 / 18
  beta_T = 1.5
  gamma_T = 1
  delta_T = -0.91
  meanT = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        meanT = meanT + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]) / (1 + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p])) * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        meanT = meanT + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]) / (1 + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p])) * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        meanT = meanT + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]) / (1 + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p])) * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  
  # So that P(Y1 = 1) = pY1, for each targeted type
  meanY1_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        
        pt = 1 / 2 / meanT * exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]) / (1 + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]))
        
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
        meanY1_part1 = meanY1_part1 + exp(wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * (exp(beta_1) * pt + 1 - pt )
      }}}
  
  # So that P(Y2 = 1) = pY2, for each non-targeted type
  meanY2_part1 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] 
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        meanY2_part1 = meanY2_part1 + exp(wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p] * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
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
        
        pt = 1 / 2 / meanT * exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]) / (1 + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]))
        
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
        covY2 = covY2 + t(exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] 
        covY2 = covY2 + t(exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        covY2 = covY2 + t(exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion)) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
      }}}
  varY2s = meanY2s * (1 - meanY2s) + (sum(covY2 - diag(diag(covY2))))
  
  # Covariance between primary composite and secondary outcomes
  covY1cY2s_part1 = 0
  covY1cY2s_part2 = 0
  for(i in 1:length(wRegion)){
    for(j in 1:length(wAge)){
      for(p in 1:length(a)){
        
        pt = 1 / 2 / meanT * exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]) / (1 + exp(delta_T + alpha_T * wAge[j] + beta_T * wRegion[i] + gamma_T * a[p]))
        
        covY1cY2s_part1 = covY1cY2s_part1 + t(exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% (exp(mu_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * (t(exp(beta_1)) * pt + 1 - pt)) * a[p]^2 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i]
        covY1cY2s_part1 = covY1cY2s_part1 + t(exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% (exp(mu_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * (t(exp(beta_1)) * pt + 1 - pt)) * a[p]^2 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i]
        covY1cY2s_part1 = covY1cY2s_part1 + t(exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) %*% (exp(mu_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,])) * (t(exp(beta_1)) * pt + 1 - pt)) * a[p]^2 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i]
        
        covY1cY2s_part2 = covY1cY2s_part2 - (exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * exp(sum(mu_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^3 * (p == 1) * Palow[i,j] * PwAge[i,j] * PwRegion[i] * (exp(sum(beta_1)) * pt + 1 - pt)
        covY1cY2s_part2 = covY1cY2s_part2 - (exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * exp(sum(mu_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^3 * (p == 2) * Pamedium[i,j] * PwAge[i,j] * PwRegion[i] * (exp(sum(beta_1)) * pt + 1 - pt)
        covY1cY2s_part2 = covY1cY2s_part2 - (exp(mu_2 + wAge[j] %*% t(s_WAge) + wRegion[i] %*% t(s_WRegion))) * exp(sum(mu_1 + wAge[j] %*% t(q_WAge) + t(q_WRegion[i,]))) * a[p]^3 * (p == 3) * Pahigh[i,j] * PwAge[i,j] * PwRegion[i] * (exp(sum(beta_1)) * pt + 1 - pt)
      }}}
  covY1cY2s = sum(covY1cY2s_part1) + sum(covY1cY2s_part2) - meanY1c * meanY2s # cov(Y1c, Y2s)
  
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
    
    T = rbinom(n, size = 1, prob = 1 / 2 / meanT * exp(delta_T + alpha_T * WAge + beta_T * WRegion + gamma_T * A) / (1 + exp(delta_T + alpha_T * WAge + beta_T * WRegion + gamma_T * A)))
    
    WRegion_Y1 = 0
    for(i in 1:length(wRegion)){
      WRegion_Y1 = WRegion_Y1 + 1 * (WRegion == wRegion[i]) %*% t(q_WRegion[i, ])
    }
    
    p1              = exp(Mu_1 + T %*% t(beta_1) + WAge %*% t(q_WAge) + WRegion_Y1)
    PY1             = A * p1
    Y1              = matrix(rbinom(2*n, size = 1, prob = PY1), nrow = n, byrow = FALSE) # Generation of the infections with the targeted types for each of the n individuals.
    Y1c             = 1 * (rowSums(Y1) > 0) # Primary composite outcome
    
    p2              = exp(Mu_2 + WAge %*% t(s_WAge) + WRegion %*% t(s_WRegion))
    PY2             = A * p2 
    Y2              = matrix(rbinom(l*n, size = 1, prob = PY2), nrow = n, byrow = FALSE)  # Generation of the infections with the non targeted types for each of the n individuals.
    Y2s             = rowSums(Y2)   # Secondary outcome, counting the total number of infections with non-targeted types.
    
    W = interaction(WAge, WRegion, sep = "_")
    
    ## Estimation of the parameters
    
    est.jointReg    = JointReg(Y1 = Y1c, Y2 = Y2s, T = T, Wcat = W) # Estimation using the joint approach, when the observed covariates are included in the estimating equations.
    est.jointNC     = JointNC(Y1 = Y1c, Y2 = Y2s, T = T) # Estimation using the joint approach, when the observed covariates are neglected.
    est.jointMH     = JointStratificationMH(Y1c, Y2s, T = T, W = W) # Estimation of the treatment effect with the joint standardized approach with M-H weights.
    est.MH          = StratificationMH(Y1 = Y1c, T = T, W = W) # Estimation of the treatment effect with the usual standardization, using MH weights. Observed covariates are used, but not the non-targeted HPV types.
    est.SSjoint     = SSJoint(Y1c, Y2s, T, W) # Estimation using the stratum-specific joint approach
    
    jointMH         = cbind("JointMH", est.jointMH$beta_1.hat, est.jointMH$se.beta_1.hat,  est.jointMH$beta_1.hat - est.jointMH$se.beta_1.hat * qnorm(0.975), est.jointMH$beta_1.hat + est.jointMH$se.beta_1.hat * qnorm(0.975))
    jointNC         = cbind("JointNC", est.jointNC$beta_1.hat, est.jointNC$se.beta_1.hat, est.jointNC$beta_1.hat - est.jointNC$se.beta_1.hat * qnorm(0.975), est.jointNC$beta_1.hat + est.jointNC$se.beta_1.hat * qnorm(0.975))
    MH              = cbind("MH", est.MH$beta_1.hat, est.MH$se.beta_1.hat, est.MH$beta_1.hat - est.MH$se.beta_1.hat * qnorm(0.975), est.MH$beta_1.hat + est.MH$se.beta_1.hat * qnorm(0.975))
    jointReg        = cbind("JointReg", est.jointReg$beta_1.hat, est.jointReg$se.beta_1.hat, est.jointReg$beta_1.hat - est.jointReg$se.beta_1.hat * qnorm(0.975), est.jointReg$beta_1.hat + est.jointReg$se.beta_1.hat * qnorm(0.975))
    SSjoint         = cbind("SSJoint", est.SSjoint$beta_1.hat, est.SSjoint$se.beta_1.hat,  est.SSjoint$beta_1.hat - est.SSjoint$se.beta_1.hat * qnorm(0.975), est.SSjoint$beta_1.hat + est.SSjoint$se.beta_1.hat * qnorm(0.975))
    
    recap = rbind(jointMH, jointNC, MH, jointReg, SSjoint)
    colnames(recap) = c("Approach", "beta_1.hat", "se.beta_1.hat", "CI.left", "CI.right")
    
    res = rbind(res, cbind(recap, n = n, pY1 = paste(pY1, collapse = "_"), pY2 = paste(pY2, collapse = "_"), A = paste(a, collapse = "_"), beta_1.true = beta_1.true, corrY1cY2s = corrY1cY2s, varA = varA, beta_1.hat.modgee = est.jointReg$beta_1.hat.naive))
  }
  
  myfile  = paste0("RES_NonRandomized_LargeStrataW-n", n, "-pY1", paste(pY1, collapse = "_"), "-A", paste(a, collapse = "_"), "-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData")
  save(res, file = myfile)
}

P = nrow(PARAM)
resultat = mclapply(1:P, Onerun, mc.cores = 18)

RES = NULL
for(p in 1: nrow(PARAM)){
  pY1             = PARAM[p, 1:2]
  a               = PARAM[p, 3:5]
  n               = PARAM[p, 6]
  
  load(paste0("RES_NonRandomized_LargeStrataW-n", n, "-pY1", paste(pY1, collapse = "_"), "-A", paste(a, collapse = "_"), "-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData"))
  RES = rbind(RES, res)
}

RECAP = as.data.frame(RES)
RECAP$beta_1.hat = as.numeric(RECAP$beta_1.hat)
RECAP$Approach = as.factor(RECAP$Approach)
RECAP$Approach = factor(RECAP$Approach, labels = c("Joint-MH", "Joint-NC", "Joint-Reg", "MH", "SS-Joint"))
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

myfile  = paste0("RECAP_NonRandomized_LargeStrataW-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData")
save(RECAP, file = myfile)

library(ggplot2)
library(gtable)
library(grid)
library(xtable)

## Plot of all scenarios
RECAP$pY1 = factor(RECAP$pY1, labels = c("Low", "Medium", "High"))
RECAP$A = factor(RECAP$A, labels = c("Small", "Medium", "Large"))
RECAP$Approach = factor(RECAP$Approach,levels(RECAP$Approach)[c(3,1,5,2,4)])
plot = ggplot(RECAP, aes(x = n, y = beta_1.hat, color = Approach)) + geom_boxplot(coef = NULL) + geom_hline(aes(yintercept = beta_1.true)) + theme_light() + theme(plot.title = element_text(size = 11), axis.title = element_text(size = 11), axis.text = element_text(size = 11), legend.text = element_text(size = 11), strip.background = element_rect(color="black", fill="white", size = 0.5, linetype="solid"), strip.text.x = element_text(size = 11, color = "black"), strip.text.y = element_text(size = 11, color = "black")) + ylab((expression(hat(beta[1]^C) ))) + facet_grid(pY1~A, labeller = label_parsed, scales = "free_y") + scale_color_manual(values=c("#F8766D", "#7CAE00", "#FF65AE", "#00BFC4", "#C77CFF"))

labelT = "Variance of the unmeasured confounder A"
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
pdf(paste0("Comparison_NonRandomized_LargeStrataW-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".pdf"),  width = 10, height = 7)
grid.draw(z) # print it
dev.off() # Stop writing to the PDF file

## Details of the results
Nreplic   = 10^4
Eff = NULL
for(i in 1:nrow(PARAM)){
  
  RECAP1 = RECAP[((i-1)*(5*Nreplic) + 1):(i*(5*Nreplic)),]
  
  beta_1 = RECAP1$beta_1.true[1]
  
  # Coverage of the confidence intervals 
  cov.JointNC     = sum((RECAP1[which(RECAP1$Approach == "Joint-NC"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Joint-NC"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Joint-NC"),5])
  cov.JointMH     = sum((RECAP1[which(RECAP1$Approach == "Joint-MH"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Joint-MH"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Joint-MH"),5])
  cov.MH          = sum((RECAP1[which(RECAP1$Approach == "MH"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "MH"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "MH"),5])
  cov.JointReg    = sum((RECAP1[which(RECAP1$Approach == "Joint-Reg"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "Joint-Reg"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "Joint-Reg"),5])
  cov.SSJoint     = sum((RECAP1[which(RECAP1$Approach == "SS-Joint"),4] < beta_1)&(RECAP1[which(RECAP1$Approach == "SS-Joint"),5] > beta_1)) / length(RECAP1[which(RECAP1$Approach == "SS-Joint"),5])
  
  sd.JointNC      = sd(RECAP1[which(RECAP1$Approach == "Joint-NC"),2])      # empirical standard deviation for hat.beta1 estimated with the joint approach (neglected covariates).
  sd.JointMH      = sd(RECAP1[which(RECAP1$Approach == "Joint-MH"),2])      # empirical standard deviation for hat.beta1 estimated with the joint M-H approach.
  sd.MH           = sd(RECAP1[which(RECAP1$Approach == "MH"),2])            # empirical standard deviation for hat.beta1 estimated with the naive standardized approach.
  sd.JointReg     = sd(RECAP1[which(RECAP1$Approach == "Joint-Reg"),2])     # empirical standard deviation for hat.beta1 estimated with the joint approach, when the observed covariates are included in the estimating equations.
  sd.SSJoint      = sd(RECAP1[which(RECAP1$Approach == "SS-Joint"),2])      # empirical standard deviation for hat.beta1 estimated with the joint Stratun Specific approach.
  
  sandwich_sd.JointNC   = mean(RECAP1[which(RECAP1$Approach == "Joint-NC"),3]) # mean of the sandwich standard deviations for hat beta1, with the joint approach (neglected covariates).
  sandwich_sd.JointMH   = mean(RECAP1[which(RECAP1$Approach == "Joint-MH"),3]) # mean of the sandwich standard deviations for hat beta1, with the joint M-H approach.
  sandwich_sd.MH        = mean(RECAP1[which(RECAP1$Approach == "MH"),3]) # mean of the sandwich standard deviations for hat beta1, with the naive standardized approach.
  sandwich_sd.JointReg  = mean(RECAP1[which(RECAP1$Approach == "Joint-Reg"),3]) # mean of the sandwich standard deviations for hat beta1, with the joint approach, when the observed covariates are included in the estimating equations.
  sandwich_sd.SSJoint   = mean(RECAP1[which(RECAP1$Approach == "SS-Joint"),3]) # mean of the sandwich standard deviations for hat beta1, with SS-joint.
  
  mean.JointNC    = mean(RECAP1[which(RECAP1$Approach == "Joint-NC"),2]) # mean of hat beta1, when estimated with the joint approach (neglected covariates).
  MSE.JointNC     = mean((RECAP1[which(RECAP1$Approach == "Joint-NC"),2] - beta_1)^2) # MSE, when beta1 is estimated with the joint naive approach (neglected covariates).
  
  mean.JointMH    = mean(RECAP1[which(RECAP1$Approach == "Joint-MH"),2]) # mean of hat beta1, when estimated with the joint M-H approach.
  MSE.JointMH     = mean((RECAP1[which(RECAP1$Approach == "Joint-MH"),2] - beta_1)^2) # MSE, when beta1 is estimated with the joint M-H approach.
  
  mean.MH         = mean(RECAP1[which(RECAP1$Approach == "MH"),2]) # mean of hat beta1, when estimated with the naive standardized MH approach.
  MSE.MH          = mean((RECAP1[which(RECAP1$Approach == "MH"),2] - beta_1)^2) # MSE, when beta1 is estimated with the naive standardized MH approach.
  
  mean.JointReg   = mean(RECAP1[which(RECAP1$Approach == "Joint-Reg"),2]) # mean of hat beta1, when estimated with the joint approach, when the observed covariates are included in the estimating equations.
  MSE.JointReg    = mean((RECAP1[which(RECAP1$Approach == "Joint-Reg"),2] - beta_1)^2) # MSE, when beta1 is estimated with the joint approach, when the observed covariates are included in the estimating equations.
  
  mean.SSJoint    = mean(RECAP1[which(RECAP1$Approach == "SS-Joint"),2]) # mean of hat beta1, when estimated with the joint Stratum Specific approach.
  MSE.SSJoint     = mean((RECAP1[which(RECAP1$Approach == "SS-Joint"),2] - beta_1)^2)
  
  eff_MH_JointMH    = MSE.MH / MSE.JointMH # ratio of MSE for the naive method MH and Joint methods
  eff_MH_JointNC    = MSE.MH / MSE.JointNC
  eff_MH_JointReg   = MSE.MH / MSE.JointReg
  eff_MH_SSJoint   = MSE.MH / MSE.SSJoint
  
  Eff = rbind(Eff, c(eff_beta1_JointMH = eff_MH_JointMH, eff_beta1_JointNC = eff_MH_JointNC, eff_beta1_JointReg = eff_MH_JointReg, bias.MH = mean.MH - beta_1, bias.JointMH = mean.JointMH - beta_1,  bias.JointNC = mean.JointNC - beta_1,  bias.JointReg = mean.JointReg - beta_1, bias.SSJoint = mean.SSJoint - beta_1,  empir_sd.MH = sd.MH, sandwich_sd.MH = sandwich_sd.MH, empir_sd.JointMH = sd.JointMH, sandwich_sd.JointMH = sandwich_sd.JointMH, empir_sd.JointNC = sd.JointNC, sandwich_sd.JointNC = sandwich_sd.JointNC, empir_sd.JointReg = sd.JointReg, sandwich_sd.JointReg = sandwich_sd.JointReg,  empir_sd.SSJoint = sd.SSJoint, sandwich_sd.SSJoint = sandwich_sd.SSJoint, CIcov.MH = cov.MH, CIcov.JointMH = cov.JointMH, CIcov.JointNC = cov.JointNC,  CIcov.JointReg = cov.JointReg,  CIcov.SSJoint = cov.SSJoint, n = as.character(RECAP1[1,]$n), pY1 = as.character(RECAP1[1,]$pY1), A =  as.character(RECAP1[1,]$A), beta_1 = beta_1, corr = RECAP1$corrY1cY2s[1] ) )
}

Eff = as.data.frame(Eff)
ColNames = colnames(Eff[,c(1:24, 27:28)])
Eff[ColNames] = sapply(Eff[ColNames], as.numeric)
# print(xtable(Eff, type = "latex", digits = 3), include.rownames=FALSE) # LaTeX table
beta_1          = c(beta_1.16, beta_1.18)
myfile  = paste0("Eff_NonRandomized_LargeStrataW-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".RData")  # save as Rdata
save(Eff, file = myfile)
Eff[ColNames] = round(Eff[ColNames], digits = 3)
Eff = Eff[which(Eff$n == 10000),] # save only the scenarios with n = 10,000
write.csv(Eff, file = paste0("Eff_NonRandomized_LargeStrataW-beta1", paste(round(beta_1, digits = 3), collapse = "_"), ".csv")) # save as csv


