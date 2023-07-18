# Effect Sizes
The R function `effectsize0` allows to compute four pairs of effect sizes, that is, Cohen’s $d_{s}$ and Hedges’ $g_{s}$ for between-subjects design, Cohen’s $d_{av}$ and Hedges’ $g_{av}$, Cohen’s $d_{rm}$ and Hedges’ $g_{rm}$, and Cohen’s $d_{z}$ and Hedges’ $g_{z}$ for within-subjects design(s). CIs can be computed according to different methods, for instance, from a noncentral t-distribution. Furthermore, CIs around each effect size can be computed with different coverage rates. The variance around each effect size is provided when displaying the output. This R script is highly inspired by Goulet-Pelletier and Cousineau (2018) and Fitts (2020). This R code is currently under modification in order to also estimate Cohen’s d_1 effect sizes (one-sample design), as well as corresponding confidence intervals.

This code regarding the estimation of Cohen's d_av and d_rm is mostly false.


Fitts, D. A. (2020). Commentary on “A review of effect sizes and their confidence intervals, Part I: The Cohen’s d family”: The degrees of freedom for paired samples designs. The Quantitative Methods for Psychology, 16(4), 250-261.

Goulet-Pelletier, J. C., & Cousineau, D. (2018). A review of effect sizes and their confidence intervals, Part I: The Cohen’sd family. The Quantitative Methods for Psychology, 14(4), 242-265.

I propose `effectsize2`, a function allowing to obtain Cohen's $d_{p}$ ($d_{s}$) and Hedges'  $g_{p}$ ($g_{s}$) for a between-subjects design. The protocol to estimate confidence intervals relies on the noncentral $\Lambda\$' distribution (Lecoutre, 1999, 2007) that was found to be an exact method (Cousineau and Goulet-Pelletier, 2020; see also Cousineau and Goulet-Pelletier, 2021). The function also allows to estimate Cohen's $d_{D}$ ($d_{z}$) and Hedges' $g_{D}$  ($g_{z}$ ) for a within-subjects design. The protocol to estimate Cohen's d_D confidence intervals relies on the method proposed by Steiger and Fouladi (1997), whereas, the protocol to estimate Hedges' g_D confidence intervals relies on the method proposed by Hedges and Olkin (1985). the choices were in accordance with Fitts'(2021) findings. 
