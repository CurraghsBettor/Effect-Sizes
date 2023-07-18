# Effect Sizes
The R function `effectsize0` allows to compute four pairs of effect sizes, that is, Cohen’s d_s and Hedges’ g_s for between-subject design, Cohen’s d_av and Hedges’ g_av, Cohen’s d_rm and Hedges’ g_rm, and Cohen’s d_z and Hedges’ g_z for within-subject design(s). CIs can be computed according to different methods, for instance, from a noncentral t-distribution. Furthermore, CIs around each effect size can be computed with different coverage rates. The variance around each effect size is provided when displaying the output. This R script is highly inspired by Goulet-Pelletier and Cousineau (2018) and Fitts (2020). This R code is currently under modification in order to also estimate Cohen’s d_1 effect sizes (one-sample design), as well as corresponding confidence intervals.

This code regarding the estimation of Cohen's d_av and d_rm is mostly false.


Fitts, D. A. (2020). Commentary on “A review of effect sizes and their confidence intervals, Part I: The Cohen’s d family”: The degrees of freedom for paired samples      designs. The Quantitative Methods for Psychology, 16(4), 250-261.

Goulet-Pelletier, J. C., & Cousineau, D. (2018). A review of effect sizes and their confidence intervals, Part I: The Cohen’sd family. The Quantitative Methods for Psychology, 14(4), 242-265.

I propose `effectsize2` 
