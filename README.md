# Bayesian meta-analysis with systematic error

## R

- `get_eumaeus_data.R` can be used to extract Eumaeus results to apply this method to (requires database key). Those data should be saved into the `data` directory (git-ignored). 
- `get_est_functions.R` contains functions for applying stan models to historical controls data
- `compare_likelihoods.R` compares the timing/results for those models, which use different likelihoods

## stan

- `metanalysis-likelihood.stan` is the stan model that implements the method for an arbitrary likelihood. Requires a set of negative control estimates and standard errors per site and a approximation of the profile likelihood of the effect of interest for each site. The prior distributions are specified here but the parameters in the R code. This model assumes that the likelihood was evaluated along a regular series of points (now forced to be in the R code).
- `metanalysis-poisson.stan` uses the poisson likelihood (e.g., for historical controls data) to do the same thing. Requires the event counts/person-time for the exposed and comparison group for both the effect of interest and the negative controls (see R code).
- `metanalysis-normal.stan` uses the normal likelihood for both the effect of interest and the negative controls. Requires point estimates and standard errors for each.

## data

- Where the Eumaeus results get saved.

## renv

- R packages

## documents

## results
