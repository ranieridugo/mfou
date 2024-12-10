f_ stands for (collection of) functions, p_ stand for program

f_covariance: covariance matrix of the mfBm and mfOU
f_covariance_derivatives: analytical gradient of the cross-covariance function of the mfOU used in the MDE estimation
f_estimator_wxy23: estimation routine for univariate marginal parameters proposed by Wang, Xiao, and Yu (2023), JoE
f_estimator_dgp24.R: estimators of rho and eta starting from given values for the univariate parameters (Dugo Giorgio Pigato 2024a)
f_estimator_mde.R: functions performing the Minimum Distance Estimation (MDE) advocated in Dugo, Giorgio, and Pigato 2024b
f_spillovers.R: forecast error variance decomposition in the mfOU and mfBm models and spillover indices according to Diebold and Yilmaz (2012), IJF

Dependencies: pracma, calculus, ggplo2

--------------------------------------------------------------------------------------------------------------

p_ stands for program
p_empirics.R: replicate the results in the empirical section of Dugo Giorgio and Pigato 2024b using the functions defined in f_

Dependencies: tidyverse

--------------------------------------------------------------------------------------------------------------

Correspondence to: dugoranieri@gmail.com