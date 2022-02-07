## Henry, W. (2022). effectsize0 (Version 1.0) [Source code].
## ------- Inspired from: Goulet-Pelletier and Cousineau (2018); Fitts (2020)
rm(list = ls()) # clear workspace

effectsize0 = function(x, y, design, method) {
  coverage = 0.95 # coverage
  sd = if (design == "between"){ #sd
    sqrt((var(x)*(length(x)-1)+var(y)*(length(y)-1))/(length(x)+length(y)-2)) # (weighted) pooled standard deviation for Cohen's d_s/d_p
  } else if (design == "within_av") {
    (sd(x)+ sd(y))/2 # (unweighted) average standard deviation for Cohen's d_av
  } else {
    r = cor(x,y) # correlation coeffecient
    cov = cov(x, y) # covariance
    sqrt(var(x)+var(y)-2*r*sd(x)*sd(y)) # standard deviation for Cohen's d_z/d_D and Cohen's d_rm/d_Dc; equivalent to the standard deviation of the differnce scores
  }
  
  d_biased = (mean(x)-mean(y))/sd # Cohen's d = d biased
  if (design == "within_rm") {
    r = cor(x,y) # correlation coefficient
    cov = cov(x, y) # covariance
    d_biased = d_biased*sqrt(2*(1-r)) # Cohen's d_rm
  } else {
    d_biased = d_biased # Cohen's d (_s/_av/_z)
  }
  
  nu = if (design == "between") { # df
    length(x)+length(y)-2
  } else {
    n=length(x)
    n-1 
              # Williams' note: Goulet-Pelletier and Cousineau (2018) used 2(n-1); but Fitts (2020) argued that the appropriate df is function of Rho, that is, lies between 2(n-1) and n-1 as a function of the increasing value of Rho
              # but compared to 2(n-1) using n-1 offers an appropriate coverage rate at multiple effect sizes, Rho values and sample sizes (Fitts, 2020)
  }     
  
  J = gamma(0.5*nu)/(sqrt(nu/2)*gamma(0.5*(nu-1))) # Hedges' correction factor (with a gamma function) to unbiased d (Hedges, 1981) 
  d_unbiased = d_biased*J # Cohen's d_unbiased/ Hedges' g, g*, Hedges'h (see Nakagawa & Cuthill, 2007)
  
  ## compute se 
  library(psych)
  ms = c(length(x), length(y))
  hm =  harmonic.mean(ms); hm # Harmonic mean = 2*length(x)*length(y)/(length(x)+length(y))
  
  var_effs = if (design == "between") { # variance around Cohen's d
    nu/(nu-2)*(2/hm)*(1+d_biased^2*(hm/2))-d_biased^2/J^2
  } else if (design == "within_av" && design == "within_rm") {
    r = cor(x,y) # correlation coefficient
    cov = cov(x, y) # covariance
    (nu/(nu-2))*(2*(1-r)/n)*(1+d_biased^2*(n/(2*(1-r))))-d_biased^2/J^2 
  } else {
    r = cor(x,y) # correlation coefficient
    cov = cov(x, y) # covariance
    (nu/(nu-2))*(1/n)*(1+d_biased^2*n)-d_biased^2/J^2
  }
  variance_g = var_effs* J^2 # variance around Hedges'g
  se_d = sqrt(var_effs) # se
  se_g = sqrt(variance_g) # se
  
  ## Compute 95% CI from a noncentral t distribution (method == "noncentral") and a central t distribution (method == "central")
  # Estimate a noncentrality paramter lambda
  lambda_d = if (design == "between") {
    d_biased*sqrt(hm/2)
  } else if (design == "within_av" && design == "within_rm") {
    r = cor(x, y)
    d_biased*sqrt(n/(2*(1-r)))
  } else {
    d_biased*sqrt(n)
  }
  lambda_g = if (design == "between") {
    d_unbiased*sqrt(hm/2)
  } else if (design == "within_av" && design == "within_rm") {
    d_unbiased*sqrt(n/(2*(1-r)))
  } else {
    d_unbiased*sqrt(n)
  } 
  if (method == "noncentral") {
    # Confidence interval around lambda
    tll_d = qt(0.5-coverage/2, df=nu, ncp=lambda_d)
    tul_d = qt(0.5+coverage/2, df=nu, ncp=lambda_d)
    tll_g = qt(0.5-coverage/2, df=nu, ncp=lambda_g) 
    tul_g = qt(0.5+coverage/2, df=nu, ncp=lambda_g)
    # Cohen's d 95%CI
    dll = tll_d/lambda_d*d_biased
    dul = tul_d/lambda_d*d_biased
    # Hedges'g 95%CI
    gll = tll_g/lambda_g*d_unbiased
    gul = tul_g/lambda_g*d_unbiased
  } else if (method == "central") { # not sure whether the following central method is what Goulet-Pelletier & Cousineau (2018) refer about
    tll = qt(0.5-coverage/2, df=nu)
    tul = qt(0.5+coverage/2, df=nu)
    # Cohen's d 95%CI
    dll = d_biased+se_d*tll
    dul = d_biased+se_d*tul
    # Hedges'g 95%CI
    gll = d_unbiased+se_d*tll
    gul = d_unbiased+se_d*tul
  } else {
    # Cohen's d 95%CI
    dll = d_biased-1.96*se_d
    dul = d_biased+1.96*se_d
    # Hedges'g 95%CI
    gll = d_unbiased-1.96*se_d
    gul = d_unbiased+1.96*se_d
  }
  #----------------------------------------------------------------------------------------------------------------
  if (design == "between") {
      cat("Cohen's d_s = ", d_biased, "\n", coverage*100, "% CI = [", dll, dul, "]\n") 
      cat("Hedges'g_s = ", d_unbiased, "\n", coverage*100, "% CI = [", gll, gul, "]\n")
      cat("Variance around d_s = ", var_effs,"\n")
      cat("Variance around g_s = ", variance_g,"\n")
    } else if (design == "within_av") {
      cat("variance-covariance matrix", "{",var(x),"
                            ",cov(x,y),var(y), "}\n")
      cat("Cohen's d_av = ", d_biased, "\n", coverage*100, "% CI = [", dll, dul, "]\n")
      cat("Hedges'g_av = ", d_unbiased, "\n", coverage*100, "% CI = [", gll, gul, "]\n")
      cat("Variance around d_av = ", var_effs,"\n")
      cat("Variance around g_av = ", variance_g,"\n")
    } else if (design == "within_rm") {
      cat("variance-covariance matrix", "{",var(x),"
                            ",cov(x,y),var(y), "}\n")
      cat("Cohen's d_rm = ", d_biased, "\n", coverage*100, "% CI = [", dll, dul, "]\n")
      cat("Hedges'g_rm = ", d_unbiased, "\n", coverage*100, "% CI = [", gll, gul, "]\n")
      cat("Variance around d_rm = ", var_effs,"\n")
      cat("Variance around g_rm = ", variance_g,"\n")
    } else if (design == "within_z") {
      cat("variance-covariance matrix", "{",var(x),"
                            ",cov(x,y),var(y), "}\n")
      cat("Cohen's d_z = ", d_biased, "\n", coverage*100, "% CI = [", dll, dul, "]\n")
      cat("Hedges'g_z = ", d_unbiased, "\n", coverage*100, "% CI = [", gll, gul, "]\n")
      cat("Variance around d_z = ", var_effs,"\n")
      cat("Variance around g_z = ", variance_g,"\n")
    }
}


#-----------------------------------------------------#
#######################################################
#-----------------------------------------------------#

## design : choose which effect size to calculate according to your desgin
# between --> Cohen's d_s/Hedges'g_s
# within_av --> Cohen's d_av/Hedges'g_av
# within_rm --> Cohen's d_rm/Hedges'g_rm
# within_z --> Cohen's d_z/Hedges'g_z

## method : choose from which method 95%CI should be estimated
# noncentral
# central
# classic
## effectsize0(x=x, y=y, design ="between" or "within_av" or "within_rm", or "within_z", method = "noncentral" or "central" or "classic")
