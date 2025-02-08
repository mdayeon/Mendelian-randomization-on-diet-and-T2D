# Mendelian-randomization-on-diet-and-T2D

Mendelian randomization approach to assess causality in diet and type 2 diabetes mellitus association.
3 specific dietary habits as exposures
T2D as the outcome of interest

The MR approach consists of univariable and multivariable including two-step mediation analyses. 
  For univariable MR, TwoSampleMR R package was used to replicate the 3 dietary habits reported to be associated with cardiometabolic health.
  For multivariable MR (MVMR), MVMR along with TwoSampleMR packages was used to evaluate those 3 dietary associations in joint with potential mediating traits as secondary exposures.
  For two-step MR mediation, RMediation R package was used to estimate indirect (mediated) effect of a given potential mediating trait, calculate confidence intervals and thereby confirm/validate the MVMR results.
  
To conduct these analyses the following codes/scripts and tables of genetic instruments were used.

  1. UVMR_diet_T2D_TwoSampleMR.R

  2. MVMR_diet_T2D_using_both_2samplemr_mvmr_with_proxy.R
     
  3. Supplementary Table 9 for the table of genetic instruments for UVMR and MVMR analyses
