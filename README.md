# Bayesian meta-analysis with systematic error

## code

- `get_eumaeus_data.R` can be used to extract Eumaeus results to apply this method to (requires database key). Those data should be saved into the `data` directory (git-ignored). 
- `approx-eumaeus.R` applies stan model described below to a set of Eumaeus results. 

## stan

- `NCs-multiple-sites-priors-real-effect-likelihood.stan` is the stan model that implements the method. Requires a set of negative control estimates and standard errors per site and a approximation of the profile likelihood of the effect of interest for each site. The prior distributions are specified here but the parameters in the R code. This model assumes that the likelihood was evaluated along a regular series of points.
- `NCs-multiple-sites-priors-real-effect-likelihood-irregular.stan` is a stan model that does the same thing, but if the set of points is not regular. This is realllllly slow and should be avoided (or sped up).
- The rest is a series of me building more and more layers into the model and can be mostly ignored.

## data

- Where the Eumaeus results get saved.

## renv

- R packages

## documents