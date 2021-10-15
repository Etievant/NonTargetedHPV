# NonTargetedHPV

##  Increasing efficiency and reducing bias when assessing HPV vaccination efficacy by using non-targeted HPV strains

Code for replicating the simulations in “Increasing efficiency and reducing bias when assessing HPV vaccination efficacy by using non-targeted HPV strains”, by Etievant, Sampson and Gail.

### Packages required 

```
gee, ggplot2, grid, gtable, parallel and xtable.
```

### Scripts

* Script `Simul_Randomized.R` allows to replicate the simulations proposed by Etievant, Sampson and Gail in Section 3.2.

* Scripts `Simul_NonRandomized_Part1.R` and `Simul_NonRandomized_Part1.R` allow to replicate the simulations proposed in Section 3.3. 

* Scripts `Simul_Randomized_TreatmentEffectY2_Part1.R` and `Simul_Randomized_TreatmentEffectY2_Part2.R` allows to replicate the simulations proposed in Web Appendix E.2.1.1. 

* Scripts `Simul_NonRandomized_TreatmentEffectY2_Part1.R` and `Simul_NonRandomized_TreatmentEffectY2_Part2.R` allows to replicate the simulations proposed in Web Appendix E.2.1.2.

* Script `Simul_NonRandomized_LargeStrataW.R` allows to replicate the simulations proposed in Web Appendix E.2.2. 

* Script `Simul_NonRandomized_LargeStrataW_SingleTargetedType.R` allows to replicate the simulations proposed in Web Appendix E.2.3. 

Each script also relies on functions provided in `EstimationFunctions.R` and certain parameters provided in `Parameters.RData`.


### Functions provided in `EstimationFunctions.R`

* **NaiveEE** - estimation of the treatment effect with the ``unaugmented'' approach, UnAug.

* **AugmentedEE** - estimation of the treatment effect on the primary outcome of interest with method Aug, following Zhang, Tsiatis and Davidian (Biometrics, 2008), when the primary outcome is binary and the secondary outcome is categorical. The primary outcome estimating equation is augmented with a function of the secondary outcome, and two logit models are used  for E(Y1 | Y2, T = 1) and  E(Y1 | Y2, T = 0).

* **JointNC** - estimation of the treatment effect on the primary outcome of interest with method JointNC. Relies on the joint estimation of the treatment effect on primary and secondary outcomes. Information on potential observed covariates is not used. 

* **StratificationMH** - estimation of the treatment effect with usual stratification on observed covariates, using Mantel-Haenszel weights. 

* **JointStratificationMH** - estimation of the treatment effect on the primary outcome of interest with method JointMH. Relies on the joint estimation of the treatment effect on primary and secondary outcomes, using stratification on observed covariates with Mantel-Haenszel weights. 

* **JointReg** - estimation of the treatment effect on the primary outcome of interest with method JointReg. Relies on joint estimation of the treatment effect on primary  and secondary outcomes, using regression models where observed covariates are included. Covariates may be categorical, continuous, or a combination of both. For continuous  covariates, quadratic polynomial functions are used.

* **JointReg_CatCov** - estimation of the treatment effect on the primary outcome of interest with method JointReg, when observed covariates are categorical. It relies on joint estimation of the treatment on primary and secondary outcomes, using regression models where the observed categorical covariates are included. 

* **JointReg_ContCov** - estimation of the treatment effect on the primary outcome of interest with method JointReg, when observed covariates are continuous. It relies on joint estimation of the treatment on primary and secondary outcomes, using regression models where the observed continuous covariates are included (using quadratic polynomial functions for the continuous covariates).

* **SSJoint** - estimation of the treatment effect on the primary outcome of interest with method JointNC, applied stratum by stratum of an observed categorical covariate (with only a few strata). The final estimate is computed as a weighted average of the stratum-specific estimates. 

### Instructions to run each script

* Save the chosen script(s) and files `EstimationFunctions.R` and `Parameters.RData` in the same directory.

* Open and run the whole script(s).

* The results of the simulations are saved automatically in figures and tables. For example, when running script `Simul_Randomized.R`, file Comparison_Randomized-beta1-0.73_-0.86.pdf and table Eff_Randomized-beta1-0.73_-0.86.csv will give the same results as displayed in Figure 2 and Table 2 of Section 3.2.



