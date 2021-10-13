library(gee)

## Functions to estimate the treatment effect and its variance, with different methods. 

## 1- NaiveEE: estimation of the treatment effect with the naive, or un-augmented, approach.

## 2- AugmentedEE: estimation of the treatment effect on the primary outcome of interest with method Aug, following Zhang, Tsiatis
# and Davidian (Biometrics, 2008), when the primary outcome is binary and the secondary outcome is categorical. The primary outcome
# estimating equation is augmented with a function of the secondary outcome, and two logit models are used for E(Y1 | Y2, T = 1) 
# and  E(Y1 | Y2, T = 0).

## 3- JointNC: estimation of the treatment effect on the primary outcome of interest with method JointNC. Relies on the joint 
# estimation of the treatment effect on primary and secondary outcomes. Information on potential observed covariates is not used. 

## 4- StratificationMH: estimation of the treatment effect with usual stratification on observed covariates, using Mantel-Haenszel
# weights. 

## 5- JointStratificationMH: estimation of the treatment effect on the primary outcome of interest with method JointMH. Relies on
# the joint estimation of the treatment effect on primary and secondary outcomes, using stratification on observed covariates with
# Mantel-Haenszel weights. 

## 6- JointReg: estimation of the treatment effect on the primary outcome of interest with method JointReg. Relies on joint 
# estimation of the treatment effect on primary  and secondary outcomes, using regression models where observed covariates are
# included. Covariates may be categorical, continuous, or a combination of both. For continuous  covariates, quadratic polynomial
# functions are used.

## 7- JointReg_CatCov: estimation of the treatment effect on the primary outcome of interest with method JointReg, when observed
# covariates are categorical. It relies on joint estimation of the treatment on primary and secondary outcomes, using regression
# models where the observed categorical covariates are included. 

## 8- JointReg_ContCov: estimation of the treatment effect on the primary outcome of interest with method JointReg, when observed
# covariates are continuous. It relies on joint estimation of the treatment on primary and secondary outcomes, using regression
# models where the observed continuous covariates are included (using quadratic polynomial functions for the continuous covariates).

## 9- SSJoint: estimation of the treatment effect on the primary outcome of interest with method JointNC, applied stratum by 
# stratum of an observed categorical covariate (with only a few strata). The final estimate is computed as a weighted average of 
# the stratum-specific estimates. 


## NaiveEE:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.

## Values:
# mu_1.hat: estimate of the intercept.
# se.mu_1.hat: sandwich estimate for the standard deviation of mu_1.hat.
# beta_1.hat: estimate of the treatment effect.
# se.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.

NaiveEE = function(Y1 = Y1, T = T){
  n                 = length(Y1) # sample size

  e.mu_1.hat        = sum((1 - T) * Y1) / sum(1-T) # estimation of the intercept.
  mu_1.hat          = log(e.mu_1.hat)
  e.beta_1.hat      = sum(T * Y1) / sum(T) / e.mu_1.hat # estimation of the treatment effect.
  beta_1.hat        = log(e.beta_1.hat)
  
  # Sandwich estimation of the variance
  p1.hat            = e.mu_1.hat * (e.beta_1.hat * T + 1 - T)
  U1.hat            = rbind((Y1 - p1.hat) / (1 - p1.hat), T * (Y1 - p1.hat) / (1 - p1.hat))
  drond_U1.hat      = matrix(NA, nrow = 2, ncol = 2)
  drond_U1.hat[1,]  = rbind(sum(p1.hat * (Y1 - 1) / (1 - p1.hat)^2), sum(T * p1.hat * (Y1 - 1) / (1 - p1.hat)^2)) / n
  drond_U1.hat[2,]  = rbind(sum(T * p1.hat * (Y1 - 1) / (1 - p1.hat)^2), sum(T^2 * p1.hat * (Y1 - 1) / (1 - p1.hat)^2)) / n
  
  if(e.beta_1.hat == 0){
    print("No treated cases")
    beta_1.hat      = NA
    Var_theta1.hat  = matrix(NA, nrow = 2, ncol = 2)
  }else{
    Var_theta1.hat  = 1 / n * solve(drond_U1.hat) %*% (U1.hat %*% t(U1.hat) / n) %*% t(solve(drond_U1.hat)) # estimated variance, with the sandwich formula.
  }
  
  return(list(mu_1.hat = mu_1.hat, se.mu_1.hat = sqrt(diag(Var_theta1.hat))[1], beta_1.hat = beta_1.hat, se.beta_1.hat = sqrt(diag(Var_theta1.hat))[2]))
}


## AugmentedEE:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: secondary outcome. Two logit links are used to relate Y1 to Y2, for individuals with T = 0 and T = 1, respectively.

## Values:
# mu_1.hat: estimate of the intercept.
# se.mu_1.hat: sandwich estimate for the standard deviation of mu_1.hat.
# beta_1.hat: estimate of the treatment effect.
# se.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.

AugmentedEE = function(Y1 = Y1, T = T, Y2 = Y2){
  n                     = length(Y1) # sample size
  Pi_T1                 = sum(T) / n # P(T = 1)
  T0.individuals        = which(T == 0) # Untreated individuals
  T1.individuals        = which(T == 1) # Treated individuals
  dataset               = as.data.frame(cbind(T = T, Y1 = Y1, Y2 = Y2))
  dataset.T0            = dataset[T0.individuals,]
  dataset.T1            = dataset[T1.individuals,]
  
  # Two logit models are used to predict E(Y1 | Y2, T = 0) and E(Y1 | Y2, T = 1)
  model.T0              = glm(as.factor(Y1) ~ Y2, family = "binomial", data = dataset.T0)
  model.T1              = glm(as.factor(Y1) ~ Y2, family = "binomial", data = dataset.T1)
  meanY1_Y2T0           = predict.glm(model.T0, newdata = as.data.frame(Y2), type = "response")
  meanY1_Y2T1           = predict.glm(model.T1, newdata = as.data.frame(Y2), type = "response")

  # Parameter estimation
  e.mu_1.hat            = sum((1-T) * Y1 + (T - Pi_T1) * meanY1_Y2T0) / (sum(1 - T) + sum(Pi_T1 - T)) # estimation of the intercept for Y1.
  mu_1.hat              = log(e.mu_1.hat)
  e.beta_1.hat          = sum(T * Y1 - (T - Pi_T1) * meanY1_Y2T1) /  (sum(T) - sum(Pi_T1 - T))  / e.mu_1.hat # estimation of the treatment effect for Y1.
  beta_1.hat            = log(e.beta_1.hat)
  
  # Sandwich estimation of the variance
  p1.hat                = e.mu_1.hat * (e.beta_1.hat * T + 1 - T)
  p11.hat               = e.mu_1.hat * e.beta_1.hat
  p10.hat               = e.mu_1.hat
  U1.hat                = rbind((Y1 - p1.hat) / (1 - p1.hat), T * (Y1 - p1.hat) / (1 - p1.hat))
  a_1.hat               = cbind( (meanY1_Y2T1 - p11.hat) / (1 - p11.hat) , (meanY1_Y2T1 - p11.hat) / (1 - p11.hat))
  a_0.hat               = cbind( (meanY1_Y2T0 - p10.hat) / (1 - p10.hat) , 0)
  augmenting_term.hat   = rbind((T - Pi_T1), (T - Pi_T1)) * (t(a_1.hat) - t(a_0.hat))
  drond_U1.hat          = matrix(NA, nrow = 2, ncol = 2)
  drond_U1.hat[1,]      = rbind(sum(p1.hat * (Y1 - 1) / (1 - p1.hat)^2), sum(T * p1.hat * (Y1 - 1) / (1 - p1.hat)^2)) / n
  drond_U1.hat[2,]      = rbind(sum(T * p1.hat * (Y1 - 1) / (1 - p1.hat)^2), sum(T^2 * p1.hat * (Y1 - 1) / (1 - p1.hat)^2)) / n
  Var_theta_1.hat       = 1 / n * solve(drond_U1.hat) %*% (( (U1.hat - augmenting_term.hat) %*% t(U1.hat - augmenting_term.hat) )/ n) %*% t(solve(drond_U1.hat))
  
  return(list(mu_1.hat = mu_1.hat, se.mu_1.hat = sqrt(diag(Var_theta_1.hat))[1], beta_1.hat = beta_1.hat, se.beta_1.hat = sqrt(diag(Var_theta_1.hat))[2]))
}


## JointNC:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: secondary outcome.

## Values:
# beta_1.hat: estimate of the treatment effect.
# se.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.

JointNC = function(Y1 = Y1, Y2 = Y2, T = T){
  n                     = length(Y1) # sample size
  
  # Parameter estimation
  e.beta_1.tilde.hat    = sum(T * Y1) / sum(T) * sum(1 - T) / sum((1 - T) * Y1) # estimation of the treatment effect on Y1.
  beta_1.tilde.hat      = log(e.beta_1.tilde.hat) 
  e.beta_2.tilde.hat    = sum(T * Y2) / sum(T) * sum(1 - T) / sum((1 - T) * Y2) # estimation of the treatment effect on Y2.
  beta_2.tilde.hat      = log(e.beta_2.tilde.hat) 
  beta_1.hat            = beta_1.tilde.hat - beta_2.tilde.hat # de-biased treatment effect on Y1.
  
  # Sandwich estimation of the variance
  U.hat1                = rbind((T * Y1) / sum(T) - e.beta_1.tilde.hat * ((1 - T) * Y1) / sum(1-T), (T * Y2) / sum(T) - e.beta_2.tilde.hat * ((1 - T) * Y2) / sum(1-T)  )
  drond_U.hat1          = matrix(NA, nrow = 2, ncol = 2)
  drond_U.hat1[1,]      = rbind(sum(- e.beta_1.tilde.hat * sum((1 - T) * Y1) / sum(1-T)), 0) / n
  drond_U.hat1[2,]      = rbind(0, - e.beta_2.tilde.hat * sum((1 - T) * Y2) / sum(1-T)) / n
  if(is.na(e.beta_1.tilde.hat)){
    print("No untreated cases")
    beta_1.hat          = NA
    Var_theta.hat1      = matrix(NA, nrow = 2, ncol = 2)
  }else{
    if(e.beta_1.tilde.hat == 0){
      print("No treated cases")
      beta_1.hat        = NA
      Var_theta.hat1    = matrix(NA, nrow = 2, ncol = 2)
    }else{
      if(e.beta_1.tilde.hat == Inf){
        print("No untreated cases")
        beta_1.hat      = NA
        Var_theta.hat1  = matrix(NA, nrow = 2, ncol = 2)
      }else{
        Var_theta.hat1  = 1 / n * solve(drond_U.hat1) %*% (U.hat1 %*% t(U.hat1) / n) %*% t(solve(drond_U.hat1))
      }}}
  
  return(list(beta_1.hat = beta_1.hat, se.beta_1.hat = sqrt(Var_theta.hat1[1,1] + Var_theta.hat1[2,2] - 2*Var_theta.hat1[2,1])))
}


## StratificationMH:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary.
# W: covariate (stratification variable) used for adjustment.

## Values:
# beta_1.hat: final estimate of the treatment effect.
# se.beta_1.hat: estimate of the standard deviation of beta_1.hat.

StratificationMH = function(Y1 = Y1, T = T, W = NULL){
  if(is.null(W)){
    NaiveEE(Y1 = Y1, T = T)
  }else{
    
    n               = length(Y1)
    W               = as.factor(W)
    levelsW         = levels(W)
    K               = length(levelsW) # number of strata
    
    # Parameter estimation
    res             = NULL
    for(l in levelsW){
      nl            = sum(W == l) # number of individuals in the stratum
      nl1           = sum((W == l)&(T == 1)) # number of individuals who received the treatment in the stratum
      nl0           = sum((W == l)&(T == 0)) # number of individuals who did not receive the treatment in the stratum
      pl            = nl0 * nl1 / nl
      pl1           = nl0 / nl 
      pl0           = nl1 / nl
      Num_Y1        = sum(Y1[(W == l)&(T == 1)]) # number of treated cases in the stratum
      Denom_Y1      = sum(Y1[(W == l)&(T == 0)]) # number of un-treated cases in the stratum
      res           = rbind(res, cbind(Num_Y1 = Num_Y1, Denom_Y1 = Denom_Y1, pl = pl, pl1 = pl1, pl0 = pl0))
    }
    res             = as.data.frame(res)
    if(sum(is.na(res$pl))){
      res           = res[-which(is.na(res$pl)), ] 
    }
    if(sum(res$pl == 0)){
      res           = res[-which(res$pl == 0), ] 
    }
    beta_1.hat      = log(sum(res$Num_Y1 * res$pl1) / sum(res$Denom_Y1 * res$pl0)) # estimated treatment effect on the primary outcome via stratification on W with MH weights
    
    # Sandwich estimation of the variance
    U.hat           = res$pl * (res$Num_Y1 - exp(beta_1.hat) * res$Denom_Y1)
    bar_U.hat       = mean(U.hat)
    drond_U.hat     = - exp(beta_1.hat) * sum(res$Denom_Y1 * res$pl) / K
    Var_theta.hat   = 1 / (K) * solve(drond_U.hat) %*% (sum((U.hat)^2) / (K - 1))  %*% t(solve(drond_U.hat))
    
    return(list(beta_1.hat = beta_1.hat, se.beta_1.hat = sqrt(Var_theta.hat)))
  }
}


## JointStratificationMH:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: secondary outcome.
# W: categorical covariates used for adjustment.

## Values:
# beta_1.hat: final estimate of the treatment effect.
# se.beta_1.hat: estimate of the standard deviation of beta_1.hat.

JointStratificationMH = function(Y1 = Y1, Y2 = Y2, T = T, W = NULL){
  if(is.null(W)){
    # in the absence of measured covariates, perform the joint approach with no covariates (JointNC)
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    n                 = length(Y1)
    W                 = as.factor(W)
    levelsW           = levels(W)
    K                 = length(levelsW)
    
    # Parameter estimation
    res               = NULL
    for(l in levelsW){
      nl              = sum(W == l) # number of individuals in the stratum
      nl1             = sum((W == l)&(T == 1)) # number of individuals who received the treatment in the stratum
      nl0             = sum((W == l)&(T == 0)) # number of individuals who received the treatment in the stratum
      pl              = nl0 * nl1 / nl
      pl1             = nl0 / nl 
      pl0             = nl1 / nl
      Num_Y1          = sum(Y1[(W == l)&(T == 1)]) # number of treated cases in the stratum, for the primary outcome
      Denom_Y1        = sum(Y1[(W == l)&(T == 0)]) # number of untreated cases in the stratum, for the primary outcome
      Num_Y2          = sum(Y2[(W == l)&(T == 1)]) # number of treated cases in the stratum, for the secondary outcome
      Denom_Y2        = sum(Y2[(W == l)&(T == 0)]) # number of treated cases in the stratum, for the secondary outcome
      res             = rbind(res, cbind(Num_Y1 = Num_Y1, Denom_Y1 = Denom_Y1, Num_Y2 = Num_Y2, Denom_Y2 = Denom_Y2, pl = pl, pl1 = pl1, pl0 = pl0))
    }
    res               = as.data.frame(res)
    if(sum(is.na(res$pl))){
      res             = res[-which(is.na(res$pl)), ] 
    }
    if(sum(res$pl == 0)){
      res             = res[-which(res$pl == 0), ] 
    }
    beta_1.tilde.hat  = log(sum(res$Num_Y1 * res$pl1) / sum(res$Denom_Y1 * res$pl0)) # estimated treatment effect on the primary outcome via stratification on W with MH weights
    beta_2.tilde.hat  = log(sum(res$Num_Y2 * res$pl1) / sum(res$Denom_Y2 * res$pl0)) # estimated treatment effect on the secondary outcome via stratification on W with MH weights
    beta_1.hat        = beta_1.tilde.hat - beta_2.tilde.hat # using the estimated non-null treatment effect on the secondary outcome to reduce confounding bias in the estimated treatment effect on the primary outcome
    
    # Sandwich variance estimation
    U.hat             = rbind(res$pl * (res$Num_Y1 - exp(beta_1.tilde.hat) * res$Denom_Y1), res$pl * (res$Num_Y2 - exp(beta_2.tilde.hat) * res$Denom_Y2))
    bar_U.hat         = rowMeans(U.hat)
    drond_U.hat       = matrix(NA, nrow = 2, ncol = 2)
    drond_U.hat[1,]   = rbind(- exp(beta_1.tilde.hat) * sum(res$Denom_Y1 * res$pl), 0) / K
    drond_U.hat[2,]   = rbind(0, - exp(beta_2.tilde.hat) * sum(res$Denom_Y2 * res$pl)) / K
    Var_theta.hat     = 1 / (K) * solve(drond_U.hat) %*% ((U.hat - bar_U.hat) %*% t(U.hat - bar_U.hat) / (K - 1)) %*% t(solve(drond_U.hat))
    
    return(list(beta_1.hat = beta_1.hat, se.beta_1.hat = sqrt(Var_theta.hat[1,1] + Var_theta.hat[2,2] - 2*Var_theta.hat[1,2]), bar_U.hat = bar_U.hat ))
  }
}


## JointReg:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: secondary outcome, assumed to be categorical. 
# Wcat: categorical covariates used for adjustment.
# Wcont: continuous covariates used for adjustment. We assume Wcont affects Y1 and Y2 through linear quadratic functions.

## Values:
# beta_1.hat: estimate of the treatment effect. ("de-biased" treatment effect on Y1).
# se.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.
# beta_1.hat.naive: naive estimate of the treatment effect, using Y1 only. ("naive" treatment effect on Y1).
# beta_2.hat.naive: estimate of the treatment effect in the Y2 model. (treatment effect on Y2).

JointReg = function(Y1 = Y1, Y2 = Y2, T = T, Wcont = NULL, Wcat = NULL){
  n                     = length(Y1)
  
  if((is.null(Wcat))|(is.null(Wcont))){
    if((is.null(Wcat))&(is.null(Wcont))){
      JointNC(Y1 = Y1, Y2 = Y2, T = T)
    }else{
      if((is.null(Wcat))){
        JointReg_ContCov(Y1 = Y1, Y2 = Y2, T = T, Wcont = Wcont)
      }else{
        JointReg_CatCov(Y1 = Y1, Y2 = Y2, T = T, Wcat = Wcat)
      }}}
  else{
    Wcat                = as.matrix(Wcat)
    Wcont               = as.matrix(Wcont)
    Wcont2              = Wcont^2 # using a quadradic term for the continuous covariates
    levelsW             = levels(as.factor(Wcat)) # building dummy variables for the categories of W
    m                   = length(levelsW)
    WW                  = matrix(0, nrow = n, ncol = m)
    for(i in 1:m){
      WW[which(Wcat == levelsW[i]),i] = 1
    }
    if(sum((colSums(WW) == 0))){
      levelsWW          = levelsW[-which(colSums(WW) == 0)]
      WW                = WW[,-which(colSums(WW) == 0)]
      m                 = ncol(WW)
    }
    WW                  = WW[,-m]
    WW                  = as.matrix(WW)
    
    # regression models for the primary and secondary outcome
    # alternatively, one could have solved the associated estimating equations
    mod1.gee            = glm(Y1 ~ T + Wcont + Wcont2 + WW, family = binomial(link = "log"))
    mod2.gee            = glm(Y2 ~ T + Wcont + Wcont2 + WW, family = poisson)
    
    # parameters estimates
    mu_1.hat            = mod1.gee$coefficients[1]
    beta_1.hat          = mod1.gee$coefficients[2]
    alpha_1cont.hat     = c(mod1.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    alpha_1cat.hat      = c(mod1.gee$coefficients[(3 + ncol(Wcont) * 2):(length(mod1.gee$coefficients))])
    mu_2.hat            = mod2.gee$coefficients[1]
    beta_2.hat          = mod2.gee$coefficients[2]
    alpha_2cont.hat     = c(mod2.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    alpha_2cat.hat      = c(mod2.gee$coefficients[(3 + ncol(Wcont) * 2):(length(mod2.gee$coefficients))])
    
    # Variance of the estimates, using the estimating equations
    p1.hat              = exp(mu_1.hat + beta_1.hat * T + cbind(Wcont,Wcont2) %*% alpha_1cont.hat + WW %*% alpha_1cat.hat) 
    p2.hat              = exp(mu_2.hat + beta_2.hat * T + cbind(Wcont,Wcont2) %*% alpha_2cont.hat + WW %*% alpha_2cat.hat) 
    
    ee1                       = (Y1 - p1.hat) / (1 - p1.hat)             # First estimating equation for Y1
    ee1T                      = T * (Y1 - p1.hat) / (1 - p1.hat)         # "T" estimating equation for Y1
    ee1Wcont                  = cbind(Wcont, Wcont2) * matrix(rep((Y1 - p1.hat) / (1 - p1.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F)         # "Wcont" estimating equation for Y1
    ee1wcat                   = as.matrix((Y1 - p1.hat) / (1 - p1.hat))  # Used for the "Wcat" estimating equation for Y1
    deriv.ee1TWcont           = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the "T" estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to T.
    deriv.ee1Wcontwcat        = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # Used for the derivative of the "Wcont" estimating equation (for Y1) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to "Wcont".
    deriv.ee1Twcat            = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # Used for the derivative of the "T" estimating equation (for Y1) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to T.  Will be also used for the derivative of the "Wcat" estimating equation wr to "Wcat".
    deriv.ee1wcat             = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # Used for the derivative of the first estimating equation (for Y1) wr to "Wcat", or the derivative of the "Wcat" estimating equations wr to the intercept.
    deriv.ee1wcont            = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the first estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to the intercept. Will be also used for the derivative of the "Wcont" estimating equation wr to "Wcont".
    deriv.ee1T                = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # Derivative of the first estimating equation (for Y1) wr to T, or the derivative of the "T" estimating equations wr to the intercept.
    deriv.ee1                 = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # Derivative of the first estimating equation (for Y1) wr to the intercept.
    
    ee2                       = (Y2 - p2.hat)                            # First estimating equation for Y2
    ee2T                      = T * (Y2 - p2.hat)                        # "T" estimating equation for Y2
    ee2Wcont                  = cbind(Wcont, Wcont2) * matrix(rep((Y2 - p2.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F) # "Wcont" estimating equation for Y2
    ee2wcat                   = as.matrix((Y2 - p2.hat))  # Used for the "Wcat" estimating equation for Y2                
    deriv.ee2TWcont           = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * T * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the "T" estimating equation (for Y2) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to T.
    deriv.ee2Wcontwcat        = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) # Used for the derivative of the "Wcont" estimating equation (for Y2) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to "Wcont".
    deriv.ee2Twcat            = as.matrix(p2.hat * T * (- 1))    # Used for the derivative of the "T" estimating equation (for Y2) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to T.  Will be also used for the derivative of the "Wcat" estimating equation wr to "Wcat".         
    deriv.ee2wcat             = as.matrix(p2.hat * (- 1))    # Used for the derivative of the first estimating equation (for Y2) wr to "Wcat", or the derivative of the "Wcat" estimating equations wr to the intercept.              
    deriv.ee2wcont            = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the first estimating equation (for Y2) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to the intercept. Will be also used for the derivative of the "Wcont" estimating equation wr to "Wcont".
    deriv.ee2T                = as.matrix(p2.hat * T * (- 1))  # Derivative of the first estimating equation (for Y2) wr to T, or the derivative of the "T" estimating equations wr to the intercept.
    deriv.ee2                 = as.matrix(p2.hat * (- 1))  # Derivative of the first estimating equation (for Y2) wr to the intercept.
    
    ee1Wcat                   = WW    # "Wcat" estimating equation for Y1
    ee2Wcat                   = WW    # "Wcat" estimating equation for Y2
    deriv.ee1Wcat             = WW    # Derivative of the T estimating equation (for Y1) wr to Wcat, or derivative of Wcat estimating equations wr to T coefficient.
    deriv.ee2Wcat             = WW    # Derivative of the T estimating equation (for Y2) wr to Wcat, or derivative of the Wcat estimating equations wr to T coefficient.
    deriv.ee1TWcat            = WW    # Derivative of the first estimating equation (for Y1) wr to Wcat, or derivative of the Wcat estimating equations wr to the intercept.
    deriv.ee2TWcat            = WW    # Derivative of the first estimating equation (for Y2) wr to Wcat, or derivative of the Wcat estimating equations wr to the intercept.
    deriv.ee1WcontWcat        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(WW))    # Derivative of the Wcat estimating equation (for Y1) wr to Wcont, or derivative of the Wcont estimating equations wr to Wcat.
    deriv.ee2WcontWcat        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(WW))    # Derivative of the Wcat estimating equation (for Y2) wr to Wcont, or derivative of the Wcont estimating equations wr to Wcat.
    for(i in 1:(m-1)){
      ee1Wcat[,i]             = WW[,i] * ee1wcat 
      ee2Wcat[,i]             = WW[,i] * ee2wcat 
      deriv.ee1Wcat[,i]       = WW[,i] * deriv.ee1wcat 
      deriv.ee2Wcat[,i]       = WW[,i] * deriv.ee2wcat 
      deriv.ee1TWcat[,i]      = WW[,i] * deriv.ee1Twcat 
      deriv.ee2TWcat[,i]      = WW[,i] * deriv.ee2Twcat 
      deriv.ee1WcontWcat[,i]  = colSums(WW[,i] * deriv.ee1Wcontwcat)
      deriv.ee2WcontWcat[,i]  = colSums(WW[,i] * deriv.ee2Wcontwcat)
    }
    
    deriv.ee1Wcont        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # Derivative of the "Wcont" estimating equation for Y1
    deriv.ee2Wcont        = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # Derivative of the "Wcont" estimating equation for Y2
    for(i in 1:(ncol(Wcont) * 2)){
      deriv.ee1Wcont[,i]  = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee1wcont)
      deriv.ee2Wcont[,i]  = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee2wcont)
    }
    U.hat                 = cbind(ee1, ee1T, ee1Wcont, ee1Wcat, ee2, ee2T, ee2Wcont, ee2Wcat) # estimating function evaluated at theta.hat for each individual
    U.hat.mean            = matrix(rep(colMeans(U.hat), n), nrow = n, byrow = TRUE)
    
    # Derivative of the estimating equation evaluated in theta.hat
    # used for the sandwich variance
    drond_U.hat                                                                                             = matrix(NA, nrow = ncol(U.hat), ncol = ncol(U.hat))
    drond_U.hat[1,]                                                                                         = c(sum(deriv.ee1), sum(deriv.ee1T), colSums(deriv.ee1wcont), colSums(deriv.ee1Wcat), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[2,]                                                                                         = c(sum(deriv.ee1T), sum(deriv.ee1T), colSums(deriv.ee1TWcont), colSums(deriv.ee1TWcat), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[3:(2 + ncol(Wcont) * 2),]                                                                   = cbind(colSums(deriv.ee1wcont), colSums(deriv.ee1TWcont), deriv.ee1Wcont, deriv.ee1WcontWcat, matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F)) / n 
    UU1                                                                                                     = matrix(0, nrow = (m-1), ncol = 2*(2 +  ncol(Wcont) * 2 + m-1))
    UU1[,1]                                                                                                 = colSums(deriv.ee1Wcat) / n
    UU1[,2]                                                                                                 = colSums(deriv.ee1TWcat) / n
    UU1[,3:(2 + ncol(Wcont) * 2)]                                                                           = t(deriv.ee1WcontWcat) / n
    UU1[,(3 + ncol(Wcont) * 2):(2 +  ncol(Wcont) * 2 + m-1)]                                                = diag(colSums(deriv.ee1Wcat), nrow = length((3 + ncol(Wcont) * 2):(2 +  ncol(Wcont) * 2 + m-1))) / n
    drond_U.hat[(3 + ncol(Wcont) * 2):(2 +  ncol(Wcont) * 2 + m-1),]                                        = UU1
    
    drond_U.hat[1 + ncol(U.hat) / 2,]                                                                       = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2), sum(deriv.ee2T), colSums(deriv.ee2wcont), colSums(deriv.ee2Wcat)) / n
    drond_U.hat[2 + ncol(U.hat) / 2,]                                                                       = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2T), sum(deriv.ee2T), colSums(deriv.ee2TWcont), colSums(deriv.ee2TWcat)) / n
    drond_U.hat[(3 + ncol(U.hat) / 2):(2 + ncol(Wcont) * 2 + ncol(U.hat) / 2),]                             = cbind(matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F), colSums(deriv.ee2wcont), colSums(deriv.ee2TWcont), deriv.ee2Wcont, deriv.ee2WcontWcat) / n 
    UU2                                                                                                     = matrix(0, nrow = (m-1), ncol = 2*(2 +  ncol(Wcont) * 2 + m-1))
    UU2[,ncol(U.hat) / 2 + 1]                                                                               = colSums(deriv.ee2Wcat) / n
    UU2[,ncol(U.hat) / 2 + 2]                                                                               = colSums(deriv.ee2TWcat) / n
    UU2[,(3 + ncol(U.hat) / 2):(2+ ncol(U.hat) / 2 + ncol(Wcont) * 2)]                                      = t(deriv.ee2WcontWcat) / n
    UU2[,(3 + ncol(U.hat) / 2 + ncol(Wcont) * 2):(2+ ncol(U.hat) / 2 +  ncol(Wcont) * 2 + m-1)]             = diag(colSums(deriv.ee2Wcat), nrow = length((3 + ncol(U.hat) / 2 + ncol(Wcont) * 2):(2+ ncol(U.hat) / 2 +  ncol(Wcont) * 2 + m-1))) / n
    drond_U.hat[(3 + ncol(U.hat) / 2 + ncol(Wcont) * 2):(2 + ncol(U.hat) / 2 +  ncol(Wcont) * 2 + m-1),]    = UU2
    
    # Sandwich variance
    Var_theta.hat                 = 1 / n * solve(drond_U.hat) %*% (t(U.hat - U.hat.mean) %*% (U.hat - U.hat.mean) / n) %*% t(solve(drond_U.hat)) 
  
    return(list(beta_1.hat = beta_1.hat - beta_2.hat, se.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)] - 2 * Var_theta.hat[(2 + ncol(U.hat) / 2),2]), beta_1.hat.naive = beta_1.hat, beta_2.hat.naive = beta_2.hat))
  }
}


## JointReg_ContCov:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: (categorical) secondary outcome. 
# Wcont: continuous covariates used for adjustment.

## Values:
# beta_1.hat: estimate of the treatment effect. ("de-biased" treatment effect on Y1).
# se.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.
# beta_1.hat.naive: naive estimate of the treatment effect, using Y1 only. ("naive" treatment effect on Y1).
# beta_2.hat.naive: estimate of the treatment effect in the Y2 model. (treatment effect on Y2).

JointReg_ContCov = function(Y1 = Y1, Y2 = Y2, T = T, Wcont = NULL){
  n                   = length(Y1)
  
  if(is.null(Wcont)){
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    
    Wcont             = as.matrix(Wcont)
    Wcont2            = Wcont^2
    
    # Solving the estimating equations to obtain theta.hat = c(mu_1.hat, beta_1.hat, alpha_1.hat, mu_2.hat, beta_2.hat, alpha_2.hat)
    mod1.gee          = glm(Y1 ~ T + Wcont + Wcont2 , family = binomial(link = "log"))
    mod2.gee          = glm(Y2 ~ T + Wcont + Wcont2, family = poisson)
    
    mu_1.hat          = mod1.gee$coefficients[1]
    beta_1.hat        = mod1.gee$coefficients[2]
    alpha_1cont.hat   = c(mod1.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    
    mu_2.hat          = mod2.gee$coefficients[1]
    beta_2.hat        = mod2.gee$coefficients[2]
    alpha_2cont.hat   = c(mod2.gee$coefficients[3:(2 + ncol(Wcont) * 2)])
    
    # Variance of the estimates
    p1.hat          = exp(mu_1.hat + beta_1.hat * T + cbind(Wcont,Wcont2) %*% alpha_1cont.hat) 
    p2.hat          = exp(mu_2.hat + beta_2.hat * T + cbind(Wcont,Wcont2) %*% alpha_2cont.hat) 
    
    ee1                     = (Y1 - p1.hat) / (1 - p1.hat)             # First estimating equation for Y1
    ee1T                    = T * (Y1 - p1.hat) / (1 - p1.hat)         # "T" estimating equation for Y1
    ee1Wcont                = cbind(Wcont, Wcont2) * matrix(rep((Y1 - p1.hat) / (1 - p1.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F)         # "Wcont" estimating equation for Y1
    deriv.ee1TWcont         = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the "T" estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to T.
    deriv.ee1wcont          = cbind(Wcont, Wcont2) * matrix(rep(p1.hat * (Y1 - 1) / (1 - p1.hat)^2, ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the first estimating equation (for Y1) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to the intercept. Will be also used for the derivative of the "Wcont" estimating equation wr to "Wcont".
    deriv.ee1T              = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2)  # Derivative of the first estimating equation (for Y1) wr to T, or the derivative of the "T" estimating equations wr to the intercept.
    deriv.ee1               = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # Derivative of the first estimating equation (for Y1) wr to the intercept.
  
    ee2                     = (Y2 - p2.hat)                            # First estimating equation for Y2
    ee2T                    = T * (Y2 - p2.hat)                        # "T" estimating equation for Y2
    ee2Wcont                = cbind(Wcont, Wcont2) * matrix(rep((Y2 - p2.hat), ncol(Wcont) * 2 ), nrow = n, byrow = F) # "Wcont" estimating equation for Y2
    deriv.ee2TWcont         = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * T * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the "T" estimating equation (for Y2) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to T.
    deriv.ee2wcont          = cbind(Wcont, Wcont2) * matrix(rep(p2.hat * (- 1), ncol(Wcont) * 2), nrow = n, byrow = F) # Derivative of the first estimating equation (for Y2) wr to Wcont, or the derivative of the "Wcont" estimating equations wr to the intercept. Will be also used for the derivative of the "Wcont" estimating equation wr to "Wcont".
    deriv.ee2T              = as.matrix(p2.hat * T * (- 1)) # Derivative of the first estimating equation (for Y2) wr to T, or the derivative of the "T" estimating equations wr to the intercept.
    deriv.ee2               = as.matrix(p2.hat * (- 1))     # Derivative of the first estimating equation (for Y2) wr to the intercept.
    
    deriv.ee1Wcont = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # Derivative of the "Wcont" estimating equation for Y1
    deriv.ee2Wcont = matrix(NA, nrow = ncol(Wcont) * 2, ncol = ncol(Wcont) * 2) # Derivative of the "Wcont" estimating equation for Y2
    for(i in 1:(ncol(Wcont) * 2)){
      deriv.ee1Wcont[,i] = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee1wcont)
      deriv.ee2Wcont[,i] = colSums(cbind(Wcont,Wcont2)[,i] * deriv.ee2wcont)
    }
    
    U.hat             = cbind(ee1, ee1T, ee1Wcont, ee2, ee2T, ee2Wcont) # estimating function evaluated in theta.hat for each individual
    U.hat.mean        = matrix(rep(colMeans(U.hat), n), nrow = n, byrow = TRUE)
    
    # derivative of the estimating equation evaluated in theta.hat
    # used for the sandwich variance
    drond_U.hat       = matrix(NA, nrow = ncol(U.hat), ncol = ncol(U.hat))
    drond_U.hat[1,]   = c(sum(deriv.ee1), sum(deriv.ee1T), colSums(deriv.ee1wcont), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[2,]   = c(sum(deriv.ee1T), sum(deriv.ee1T), colSums(deriv.ee1TWcont), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[3:(2 + ncol(Wcont) * 2),]                                           = cbind(colSums(deriv.ee1wcont), colSums(deriv.ee1TWcont), deriv.ee1Wcont,  matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F)) / n 
    drond_U.hat[1 + ncol(U.hat) / 2,]                                               = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2), sum(deriv.ee2T), colSums(deriv.ee2wcont)) / n
    drond_U.hat[2 + ncol(U.hat) / 2,]                                               = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2T), sum(deriv.ee2T), colSums(deriv.ee2TWcont)) / n
    drond_U.hat[(3 + ncol(U.hat) / 2):(2 + ncol(Wcont) * 2 + ncol(U.hat) / 2),]     = cbind(matrix(rep(0, ncol(Wcont) * 2 * ncol(U.hat) / 2), nrow = ncol(Wcont) * 2, byrow = F), colSums(deriv.ee2wcont), colSums(deriv.ee2TWcont), deriv.ee2Wcont) / n 
 
    # Sandwich variance
    Var_theta.hat     = 1 / n * solve(drond_U.hat) %*% (t(U.hat - U.hat.mean) %*% (U.hat - U.hat.mean) / n) %*% t(solve(drond_U.hat)) 
    
    return(list(beta_1.hat = beta_1.hat - beta_2.hat, se.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)] - 2 * Var_theta.hat[(2 + ncol(U.hat) / 2),2]), beta_1.hat.naive = beta_1.hat, beta_2.hat.naive = beta_2.hat))
  }
}


## JointReg_CatCov:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: (categorical) secondary outcome. 
# Wcat: categorical covariates used for adjustment.

## Values:
# beta_1.hat: estimate of the treatment effect. ("de-biased" treatment effect on Y1).
# se.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.
# beta_1.hat.naive: naive estimate of the treatment effect, using Y1 only. ("naive" treatment effect on Y1).
# beta_2.hat.naive: estimate of the treatment effect in the Y2 model. (treatment effect on Y2).

JointReg_CatCov = function(Y1 = Y1, Y2 = Y2, T = T, Wcat = NULL){
  
  n = length(Y1)
  
  if(is.null(Wcat)){
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    
    Wcat          = as.matrix(Wcat)
    
    levelsW       = levels(as.factor(Wcat)) # building dummy variables for the categories of W
    m             = length(levelsW)
    WW            = matrix(0, nrow = n, ncol = m)
    for(i in 1:m){
      WW[which(Wcat == levelsW[i]),i] = 1
    }
    if(sum((colSums(WW) == 0))){
      levelsWW    = levelsW[-which(colSums(WW) == 0)]
      WW          = WW[,-which(colSums(WW) == 0)]
      m           = ncol(WW)
    }
    WW            = WW[,-m]
    WW            = as.matrix(WW)
    
    mod1.gee      = glm(Y1 ~ T + WW, family = binomial(link = "log"))
    mod2.gee      = glm(Y2 ~ T + WW, family = poisson)
    
    mu_1.hat      = mod1.gee$coefficients[1]
    beta_1.hat    = mod1.gee$coefficients[2]
    alpha_1cat.hat= c(mod1.gee$coefficients[(3):(length(mod1.gee$coefficients))])
    
    mu_2.hat      = mod2.gee$coefficients[1]
    beta_2.hat    = mod2.gee$coefficients[2]
    alpha_2cat.hat= c(mod2.gee$coefficients[(3):(length(mod2.gee$coefficients))])
    
    # Variance of the estimates
    p1.hat        = exp(mu_1.hat + beta_1.hat * T + WW %*% alpha_1cat.hat) 
    p2.hat        = exp(mu_2.hat + beta_2.hat * T + WW %*% alpha_2cat.hat) 
    
    ee1                     = (Y1 - p1.hat) / (1 - p1.hat)             # First estimating equation for Y1
    ee1T                    = T * (Y1 - p1.hat) / (1 - p1.hat)         # "T" estimating equation for Y1
    ee1wcat                 = as.matrix((Y1 - p1.hat) / (1 - p1.hat))  # Used for the "Wcat" estimating equation for Y1
    deriv.ee1Twcat          = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # Used for the derivative of the "T" estimating equation (for Y1) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to T.  Will be also used for the derivative of the "Wcat" estimating equation wr to "Wcat".
    deriv.ee1wcat           = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # Used for the derivative of the first estimating equation (for Y1) wr to "Wcat", or the derivative of the "Wcat" estimating equations wr to the intercept.
    deriv.ee1T              = as.matrix(p1.hat * T * (Y1 - 1) / (1 - p1.hat)^2) # Derivative of the first estimating equation (for Y1) wr to T, or the derivative of the "T" estimating equations wr to the intercept.
    deriv.ee1               = as.matrix(p1.hat * (Y1 - 1) / (1 - p1.hat)^2) # Derivative of the first estimating equation (for Y1) wr to the intercept.
    
    ee2                     = (Y2 - p2.hat)                            # First estimating equation for Y2
    ee2T                    = T * (Y2 - p2.hat)                        # "T" estimating equation for Y2
    ee2wcat                 = as.matrix((Y2 - p2.hat))                 # Used for the "Wcat" estimating equation for Y2
    deriv.ee2Twcat          = as.matrix(p2.hat * T * (- 1))            # Used for the derivative of the "T" estimating equation (for Y2) wr to Wcat, or the derivative of the "Wcat" estimating equations wr to T.  Will be also used for the derivative of the "Wcat" estimating equation wr to "Wcat".  
    deriv.ee2wcat           = as.matrix(p2.hat * (- 1))                # Used for the derivative of the first estimating equation (for Y2) wr to "Wcat", or the derivative of the "Wcat" estimating equations wr to the intercept.
    deriv.ee2T              = as.matrix(p2.hat * T * (- 1))            # Derivative of the first estimating equation (for Y2) wr to T, or the derivative of the "T" estimating equations wr to the intercept.
    deriv.ee2               = as.matrix(p2.hat * (- 1))                # Derivative of the first estimating equation (for Y2) wr to the intercept.
    
    ee1Wcat                 = WW    # "Wcat" estimating equation for Y1
    ee2Wcat                 = WW    # "Wcat" estimating equation for Y2
    deriv.ee1Wcat           = WW    # Derivative of the T estimating equation (for Y1) wr to Wcat, or derivative of Wcat estimating equations wr to T coefficient.
    deriv.ee2Wcat           = WW    # Derivative of the T estimating equation (for Y2) wr to Wcat, or derivative of the Wcat estimating equations wr to T coefficient.
    deriv.ee1TWcat          = WW    # Derivative of the first estimating equation (for Y1) wr to Wcat, or derivative of the Wcat estimating equations wr to the intercept.
    deriv.ee2TWcat          = WW    # Derivative of the first estimating equation (for Y2) wr to Wcat, or derivative of the Wcat estimating equations wr to the intercept.
    for(i in 1:(m-1)){
      ee1Wcat[,i]           = WW[,i] * ee1wcat 
      ee2Wcat[,i]           = WW[,i] * ee2wcat 
      deriv.ee1Wcat[,i]     = WW[,i] * deriv.ee1wcat 
      deriv.ee2Wcat[,i]     = WW[,i] * deriv.ee2wcat 
      deriv.ee1TWcat[,i]    = WW[,i] * deriv.ee1Twcat 
      deriv.ee2TWcat[,i]    = WW[,i] * deriv.ee2Twcat 
    }
    U.hat             = cbind(ee1, ee1T, ee1Wcat, ee2, ee2T, ee2Wcat) # estimating function evaluated in theta.hat for each individual
    U.hat.mean        = matrix(rep(colMeans(U.hat), n), nrow = n, byrow = TRUE)

    # derivative of the estimating equation evaluated in theta.hat
    # used for the sandwich variance
    drond_U.hat       = matrix(NA, nrow = ncol(U.hat), ncol = ncol(U.hat))
    drond_U.hat[1,]   = c(sum(deriv.ee1), sum(deriv.ee1T), colSums(deriv.ee1Wcat), rep(0, ncol(U.hat) / 2)) / n
    drond_U.hat[2,]   = c(sum(deriv.ee1T), sum(deriv.ee1T), colSums(deriv.ee1TWcat), rep(0, ncol(U.hat) / 2)) / n
    UU1               = matrix(0, nrow = (m-1), ncol = 2*(2 +  + m-1))
    UU1[,1]           = colSums(deriv.ee1Wcat) / n
    UU1[,2]           = colSums(deriv.ee1TWcat) / n
    UU1[,(3 ):(2  + m-1)]                                               = diag(colSums(deriv.ee1Wcat), nrow = length((3 ):(2  + m-1)), ncol = length((3 ):(2  + m-1))) / n
    drond_U.hat[(3 ):(2 + m-1),]                                        = UU1
    
    drond_U.hat[1 + ncol(U.hat) / 2,]                                   = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2), sum(deriv.ee2T), colSums(deriv.ee2Wcat)) / n
    drond_U.hat[2 + ncol(U.hat) / 2,]                                   = c(rep(0, ncol(U.hat) / 2), sum(deriv.ee2T), sum(deriv.ee2T),  colSums(deriv.ee2TWcat)) / n
    UU2                                                                 = matrix(0, nrow = (m-1), ncol = 2*(2 + m-1))
    UU2[,ncol(U.hat) / 2 + 1]                                           = colSums(deriv.ee2Wcat) / n
    UU2[,ncol(U.hat) / 2 + 2]                                           = colSums(deriv.ee2TWcat) / n
    UU2[,(3 + ncol(U.hat) / 2 ):(2+ ncol(U.hat) / 2  + m-1)]            = diag(colSums(deriv.ee2Wcat), nrow = length((3 + ncol(U.hat) / 2 ):(2+ ncol(U.hat) / 2  + m-1))) / n
    drond_U.hat[(3 + ncol(U.hat) / 2 ):(2 + ncol(U.hat) / 2 + m-1),]    = UU2
  
    # Sandwich variance
    Var_theta.hat     = 1 / n * solve(drond_U.hat) %*% (t(U.hat - U.hat.mean) %*% (U.hat - U.hat.mean) / n) %*% t(solve(drond_U.hat)) 
    
    return(list(beta_1.hat = beta_1.hat - beta_2.hat, se.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[(2 + ncol(U.hat) / 2),(2 + ncol(U.hat) / 2)] - 2 * Var_theta.hat[(2 + ncol(U.hat) / 2),2]), beta_1.hat.naive = beta_1.hat, beta_2.hat.naive = beta_2.hat))
  }
}


## SSJoint:

## Arguments:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: secondary outcome.
# W: categorical covariates (= stratification variable with only a few strata) used for adjustment.

## Values:
# beta_1.hat: final estimate of the treatment effect.
# se.beta_1.hat: estimate of the standard deviation of beta_1.hat.

SSJoint = function(Y1 = Y1, Y2 = Y2, T = T, W = NULL){
  if(is.null(W)){
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{

    n               = length(Y1)
    W               = as.factor(W)
    levelsW         = levels(W)
    
    res             = NULL
    for(l in levelsW){
      est           = JointNC(Y1 = Y1[which(W == l)], Y2 = Y2[which(W == l)], T = T[which(W == l)])
      res           = rbind(res, cbind(beta_1.hat = est$beta_1.hat, se.beta_1.hat = est$se.beta_1.hat))
    }
    res             = as.data.frame(res)
    if(sum(is.na(res$beta_1.hat))){
      res           = res[-which(is.na(res$beta_1.hat)), ] 
    }
    pl              = 1 / res$se.beta_1.hat^2 # weights inversely proportional to the variance of the estimate in each stratum
    sum_pl          = sum(pl)
    beta_1.hat      = sum(res$beta_1.hat * pl) / sum_pl # parameter estimation
    se.beta_1.hat   = sqrt(sum((res$se.beta_1.hat * pl)^2)) / sum_pl # variance estimation
    
    return(list(beta_1.hat = beta_1.hat, se.beta_1.hat = se.beta_1.hat))
  }
}
