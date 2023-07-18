## Henry, W. (2023). effectsize2.
## ------- Advanced version taking into account some works to have the most accurate CI
# see Fitts (2021); Cousineau and Goulet-Pelletier (2020, 2021)

rm(list = ls()) # clear workspace

effectsize2 <- function(x, y, design = c("between", "within"), coverage) {
  library(DPQ)
  library(psych)
  library(sadists)
  library(MBESS)
  if (design == "between"){ 
    sd <- sqrt((var(x)*(length(x)-1)+var(y)*(length(y)-1))/(length(x)+length(y)-2)) # (weighted) pooled standard deviation for Cohen's d_s/d_p
   } else if (design == "within") {
    r <- cor(x,y)
    sd <- sqrt(var(x)+var(y)-2*r*sd(x)*sd(y)) # standard deviation for Cohen's d_z/d_D and Cohen's d_rm/d_Dc; equivalent to the standard deviation of the diffrence scores
  }
  d_biased <- (mean(x)-mean(y))/sd # Cohen's d = d biased
  if (design == "between") { # df
    nu <- length(x)+length(y)-2
  } else {
    n <- length(x)
    nu <- length(x) - 1
  }     
  J <- exp (lgamma(nu/2) - log(sqrt(nu/2)) - lgamma((nu-1)/2))  # from Robert Calin-Jageman as cited in many reports (e.g., Cousineau & Goulet-Pelletier, 2020; 2021)
  d_unbiased <- d_biased*J # Cohen's d_unbiased/ Hedges' g, g*, Hedges'h (see Nakagawa & Cuthill, 2007)
  ## compute se 
  ms <- c(length(x), length(y))
  hm <-  harmonic.mean(ms); hm # Harmonic mean = 2*length(x)*length(y)/(length(x)+length(y))
  if (design == "between") { # variance around Cohen's d
    var_effs <- nu/(nu-2)*(2/hm)*(1+d_biased^2*(hm/2))-d_biased^2/J^2
  } else if (design == "within") {
    var_effs <- (nu/(nu-2))*(1/n)*(1+d_biased^2*n)-d_biased^2/J^2
  }
  variance_g <- var_effs*(J^2) # variance around Hedges'g
  ## Compute 95% CIs  
  # Estimate a non centrality parameter lambda
  if (design == "between") {
    lambda_d <- d_biased*sqrt(hm/2)
  } else if (design == "within") {
    lambda_d <- d_biased*sqrt(n)
  }
  if (design == "between") {
    lambda_g <- d_unbiased*sqrt(hm/2)
  } else if (design == "within") {
    lambda_g <- d_unbiased*sqrt(n)
  }
  if (design == "between") {
    dll <- qlambdap(1/2-coverage/2, df = nu, t = lambda_d)
    dul <- qlambdap(1/2+coverage/2, df = nu, t = lambda_d)
    gll <- qlambdap(1/2-coverage/2, df = nu, t = lambda_g)
    gul <- qlambdap(1/2+coverage/2, df = nu, t = lambda_g)
  } else if (design == "within") {
    # Confidence interval around lambda --> Steiger and Fouladi (1997)
    tCI <- conf.limits.nct(lambda_d, nu, conf.level = coverage)
    tll_d <- tCI$Lower.Limit
    tul_d <- tCI$Upper.Limit
    # Cohen's d 95%CI
    dll <- tll_d/lambda_d*d_biased
    dul <- tul_d/lambda_d*d_biased
    # Confidence interval around lambda --> Hedges and Olkin (1985)
    tll_g <- qtAppr(0.5-coverage/2, df=nu, ncp=lambda_g) 
    tul_g <- qtAppr(0.5+coverage/2, df=nu, ncp=lambda_g)
    # Hedges'g 95%CI
    gll <- tll_g/lambda_g*d_unbiased
    gul <- tul_g/lambda_g*d_unbiased
  }
  #----------------------------------------------------------------------------------------------------------------
  if (design == "between") {
    cat("Cohen's d_p = ", d_biased, "\n", coverage*100, "%CI = [", dll, dul, "]\n") 
    cat("Hedges'g_p = ", d_unbiased, "\n", coverage*100, "%CI = [", gll, gul, "]\n")
    cat("Variance around d_p = ", var_effs,"\n")
    cat("Variance around g_p = ", variance_g,"\n")
  } else if (design == "within") {
    cat("Correlation coefficient r =",cor(x, y),"\n")
    cat("Cohen's d_D = ", d_biased, "\n", coverage*100, "%CI = [", dll, dul, "]\n")
    cat("Hedges'g_D = ", d_unbiased, "\n", coverage*100, "%CI = [", gll, gul, "]\n")
    cat("Variance around d_D = ", var_effs,"\n")
    cat("Variance around g_D = ", variance_g,"\n")
  }
}


#-----------------------------------------------------#
#######################################################
#-----------------------------------------------------#

## design : choose which effect size need to be calculated according to your desgin
# between --> Cohen's d_p/Hedges'g_p
# within_z --> Cohen's d_D/Hedges'g_D

## effectsize2(x=x, y=y, design = c("between", "within", coverage = c(0.95, 0.99, 0.995))
