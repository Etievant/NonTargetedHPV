# NonTargetedHPV

Replication of the simulation studies in "Increasing efficiency and reducing bias when assessing HPV vaccination efficacy by using non-targeted HPV strains", by Etievant, Sampson and Gail (Biometrics, 2023). The methods described in the article are implemented in the file `EstimationFunctions.R`.

### Required packages 

```
gee, ggplot2, grid, gtable, parallel and xtable.
```

### Scripts

* Script `Simul_Randomized.R` replicates the simulations proposed by Etievant, Sampson and Gail in Section 3.2.

* Scripts `Simul_NonRandomized_Part1.R` and `Simul_NonRandomized_Part1.R` replicate the simulations proposed in Section 3.3. 

* Scripts `Simul_Randomized_TreatmentEffectY2_Part1.R` and `Simul_Randomized_TreatmentEffectY2_Part2.R` and `Simul_Randomized_TreatmentEffectY2_Part3.R` replicate the simulations proposed in Web Section E.3.1.1. 

* Scripts `Simul_NonRandomized_TreatmentEffectY2_Part1.R`, `Simul_NonRandomized_TreatmentEffectY2_Part2.R` and `Simul_NonRandomized_TreatmentEffectY2_Part3.R` replicate the simulations proposed in Web Section E.3.1.2.

* Script `Simul_NonRandomized_LargeStrataW.R` replicates the simulations proposed in Web Section E.3.2. 

* Script `Simul_NonRandomized_LargeStrataW_SingleTargetedType.R` replicates the simulations proposed in Web Section E.3.3. 

Each script relies on functions provided in `EstimationFunctions.R` and certain parameters provided in `Parameters.RData`.


### Instructions to run each script

* Save the chosen script(s) and files `EstimationFunctions.R` and `Parameters.RData` in the same directory.

* Open and run the whole script(s).

* The results of the simulations are saved in figures and tables. For example, when running script `Simul_Randomized.R`, file Eff_Randomized-beta1-0.73-0.86.csv will give the same results as displayed in Table 2 of Section 3.2.


### Functions provided in `EstimationFunctions.R`

* **NaiveEE** - estimation of the treatment effect with the naive, or un-augmented, approach (UnAug).

* **AugmentedEE** - estimation of the treatment effect with augmentation to gain efficiency. The primary outcome estimating equation is augmented with a function of Z, and two logit models are used for E(Y1 | Z, T = 1) and  E(Y1 | Z, T = 0). One can use data on baseline covariates for the augmentation (method Aug_W, as proposed by Zhang, Tsiatis and Davidian in Biometrics, 2008), data on a secondary outcome unaffected by the treatment (method Aug), or data on both (method Aug_WY2).

* **JointNC** - estimation of the treatment effect on the primary outcome of interest with method JointNC. Relies on the joint estimation of the treatment effect on primary and secondary outcomes. Information on potential observed covariates is not used. 

* **StratificationMH** - estimation of the treatment effect with usual stratification on observed covariates, using Mantel-Haenszel weights. 

* **JointStratificationMH** - estimation of the treatment effect on the primary outcome of interest with method JointMH. Relies on the joint estimation of the treatment effect on primary and secondary outcomes, using stratification on observed covariates with Mantel-Haenszel weights. 

* **JointReg** - estimation of the treatment effect on the primary outcome of interest with method JointReg. Relies on joint estimation of the treatment effect on primary  and secondary outcomes, using regression models where observed covariates are included. Covariates may be categorical, continuous, or a combination of both. For continuous  covariates, quadratic polynomial functions are used.

* **JointReg_CatCov** - estimation of the treatment effect on the primary outcome of interest with method JointReg, when observed covariates are categorical. It relies on joint estimation of the treatment on primary and secondary outcomes, using regression models where the observed categorical covariates are included. 

* **JointReg_ContCov** - estimation of the treatment effect on the primary outcome of interest with method JointReg, when observed covariates are continuous. It relies on joint estimation of the treatment on primary and secondary outcomes, using regression models where the observed continuous covariates are included (using quadratic polynomial functions for the continuous covariates).

* **SSJoint** - estimation of the treatment effect on the primary outcome of interest with method JointNC, applied stratum by stratum of an observed categorical covariate (with only a few strata). The final estimate is computed as a weighted average of the stratum-specific estimates. 


