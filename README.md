
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CGLS

<!-- badges: start -->
<!-- badges: end -->

# Phylogenetic bayesian mixed models using INLA

# (Last update: 06/03/2024)

While species are known to be non-independent due to shared evolutionary
history, this is rarely considered in community-level analysees.
Additionally, communities are not statisitcally indepedent if there are
species overlap (which is often the case).

The analysis accounts for non-independence between species within
communities using GLS. The main function would be “get_comm_pair_r”,
which requries two inputs: a site-species matrix, and a
variance-covariance matrix based on phylogenetic tree.

Compositional matrices in presence/absence, number of individuals,
percent cover, biomass etc can be used. Matrices based on scoring
systems with uneven intervals can also be used (e.g., DAFOR), but this
can lead to interpretation issues.

Users can calculate the variance-covariance matrix based on any
phylogenetic model (e.g. vcv in package ape is based on Brownian
Motion).

If the resulting variance-covariance matrix is not positive definite,
nearPD from the Matrix package will be applied to find the nearest
positive definite matrix. Without a positive definite matrix PGLS will
fail to run.

Please make sure both the species compositional matrix and phylogenetic
matrix have identical species name (e.g., space denoted as ” “, not”\_”
or “.”).

“get_comm_pair_r” will return a variance-covariance matrix at the
community level, which can be used in gls (from nlme).

likelihood_lambda.R were extracted from Revell et al. (2010, Methods in
Ecology & Evolution). Similar with PGLS, the function optimizes the
phylogenetic covariance matrix by minimizing log-likelihood. By default
the lambda value is restricted to 0-1.

Other functions are used for the simulations / empirical analyses and
producing figures

Note that we are using Bayesian GLMM rather than Frequentist GLMM to
improve the speed.

# Known problems

Might need to use inla.rerun to enhance model stability

It is possible that communities are indeed independent sample even if
there are species overlap. In those cases the optimization of lambda
will produce unstable results. A fix is to 1) re-run optimization
several times to see if results are largely stable and 2) compare wAIC
between models with optimized lambda and with independence assumption.
If the model with optimized lambda doesn’t perform better, users can use
the model assuming communities are independent.

Note that even with the same phylogenetic matrix, INLA will give
different results across runs. So re-running the model will give
different slightly different results (if phylogeny is important). If
phylogeny is unimportant, INLA will give very unstable results, as the
wAIC of all lambda values are driven by random fluctutation across runs.
