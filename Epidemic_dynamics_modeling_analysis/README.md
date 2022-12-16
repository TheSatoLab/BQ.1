# Epidemic dynamics modeling analysis

## Summary
Motivated by the model established by [Obermeyer et al](https://www.science.org/doi/10.1126/science.abm1208), we developed a method to model the relationship between viral epidemic dynamics and S substitutions. This model can simultaneously estimate i) the effect of each S substitution on effective reproduction number (Re) and ii) the relative Re of a viral group represented by each S haplotype. The key concept of the model used in this study is the same as the one in Obermeyer et al. However, our method is independent of the predefined viral classification such as PANGO lineage but based on the viral classification according to the profile of S substitutions. Therefore, our method can link the effect of S substitutions to viral epidemic dynamics in a more direct manner. Also, in our method, a Markov Chain Monte Carlo (MCMC) method is used for parameter estimation instead of variational inference, an approximation method.

We constructed a Bayesian hierarchal model, which represents the epidemic dynamics of each S haplotype according to growth rate parameters for each S haplotype, which is represented by a linear combination of the effect of S substitutions. Arrays in the model index over one or more indices: L = 254 viral lineages (i.e., S haplotypes) l; S = 107 substitutions/substitution clusters s; and T = 229 days t. The model is:

$$ \sigma_1\sim StudentT^+(5,0,10) $$

$$ f_m \sim Laplace(0,10) $$

$$ \beta_l\sim StudentT(5,\sum_{m}{f_mX_{lm}},\sigma_1) $$

$$ y_{.t}\sim Multinomial(\sum_{l} y_{lt},softmax(\alpha.+\beta_.t)) $$


The count of viral lineage $l$ at time $t$, $y_{lt}$, is modeled as a hierarchal Multinomial logistic regression with intercept $\alpha_l$ and slope $\beta_l$ parameters for lineage $l$. The slope (or viral lineage growth) parameter $\beta_l$ is generated from Student’s t distribution with five degrees of freedom, the mean value represented by $f_mX_{lm}$, and standard deviation, $\sigma_1$. $f_mX_{lm}$ denotes the linear combination of the effect of each substitution, where $f_m$ and $X_{lm}$ are the effect of substitution m and the profile of substitution $m$ in lineage $l$ (i.e., the substitution profile matrix constructed in the above paragraph), respectively. As a prior of $f_m$, the Laplace distribution with the mean 0 and the standard deviation 10 was set. In other words, we estimated the parameter $f_m$ in the framework of Bayesian least absolute shrinkage and selection operator (LASSO). As a prior of $\sigma_1$, a half Student’s t distribution with the mean 0 and the standard deviation 10 was set. For the other parameters, non-informative priors were set.
The relative Re of each viral lineage, $r_l$, was calculated according to the slope parameter $\beta_l$ as:

$$ r_l=exp\left(\gamma\beta_l\right) $$

where $\gamma$ is the average viral generation time (2.1 days) (http://sonorouschocolate.com/covid19/index.php?title=Estimating_Generation_Time_Of_Omicron). Similarly, the effect size of substitution m on the relative Re, $F_l$, was calculated according to the coefficient $f_l$ as:

$$ F_l=exp\left(\gamma f_l\right) $$

Parameter estimation was performed via the MCMC approach implemented in CmdStan v2.30.1 (https://mc-stan.org) with CmdStanr v0.5.3 (https://mc-stan.org/cmdstanr/). Four independent MCMC chains were run with 500 and 2,000 steps in the warmup and sampling iterations, respectively. We confirmed that all estimated parameters showed <1.01 R-hat convergence diagnostic values and >200 effective sampling size values, indicating that the MCMC runs were successfully convergent.

## Contents:
* **input/lineage\_frequency.txt:** A count matrix where the row is day and the column is viral lineage (e.g., S haplotype)
* **input/S\_substitution\_profile.txt:** A matrix where the row is viral lineage (e.g., S haplotype) and the column is substitution
* **script/multinomial\_mut\_regression.stan:** The Bayesian hierarchal model, implemented in Stan
* **script/multinomial\_mut\_regression.R:** An R code to run the stan model file

* **output/mutation_effect.txt:** Estimated effect of each mutation on Re
* **output/S_haplotype_Re.txt:** Estimated Re of each S haplotype

## Usage:
R --vanilla --slave < script/multinomial_mut_regression.R

## System requirements (R libraries)
* **tidyverse** 1.3.1
* **cmdstanr** 0.5.3




