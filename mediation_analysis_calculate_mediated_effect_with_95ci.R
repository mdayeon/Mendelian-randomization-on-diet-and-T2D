##To compute 95% CI for indirect effect (mediated effect) of an exposure (x) on outcome (y)

###Rmediation / Prodclin approach###
#https://www.public.asu.edu/~horourke/Research_in_Prevention_Laboratory_at_Arizona_State_University/PRODCLIN.html
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3233842/

library(Rmediation)

#product of coefficients approach to estimate the indirect effect (mediated effect)
#step 1 = path a : exposure (x) --> mediator (m)
#step 2 = path b : mediator (m) --(x)--> outcome (y)


#In RMediation, the function "medci" produces CIs for the product of two normal variables and mediated effects. 
#mu.x = beta estimate for path a 
#my.y = beta estimate for path b 
#se.x = se for path a
#se.y = se for path b
#alpha = significance level for CI (e.g., 90% CI, alpha = 0.1, 95% CI, alpha = 0.05)

#rho = the correlation between two variables with the default value of 0. 

#type = specified a method to estimate CI 
  #¡°prodclin¡± (default, the PRODCLIN program), 
    #Prodclin is a program that uses the distribution of the product of two normally distributed variables to compute asymmetric confidence intervals for the mediated effect. 
  #¡°dop¡± (the RDOP program) - distribution of products which is same as prodclin. dop is the best statistical approach to construct CI for mediated effect.
  #¡°MC¡± (the Monte Carlo approach), 
  #¡°asymp¡± (the AND method)
  #¡°all¡± (using all four methods).

#for the primary exposure-outcome association of interest (e.g., cheese-T2D, BMI = mediator)
medci(mu.x = beta estimate for step1/path a, 
      mu.y = beta estimate for spte2/path b, 
      se.x = standard error for step1/ path a, 
      se.y = standard error for step2/ path b, 
      rho = 0, #default
      alpha = 0.1, #default significance level for CI
      type = "prodclin" #specified method for CI estimation)



#function "qprodnormal" computes the quantile for the distribution of product. 
#p = probability for quantiles
qprodnormal(p = 0.975, mu.x = 0.2, mu.y = 0.4, se.x = 1, se.y = 1, type = "all")


