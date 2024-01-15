# CGLS - Community Generalized Least Square 
# (Last update: 15/01/2024)

While species are known to be non-independent due to shared evolutionary history, this is rarely considered in community-level analysees. Additionally, communities are not statisitcally indepedent if there are species overlap (which is often the case).

The analysis accounts for non-independence between species within communities using GLS. The main function would be "get_comm_pair_r", which requries two inputs: a site-species matrix, and a variance-covariance matrix based on phylogenetic tree. 

Compositional matrices in presence/absence, number of individuals, percent cover, biomass etc can be used. Matrices based on scoring systems with uneven intervals can also be used (e.g., DAFOR), but this can lead to interpretation issues.

Users can calculate the variance-covariance matrix based on any phylogenetic model (e.g. vcv in package ape is based on Brownian Motion). 

If the resulting variance-covariance matrix is not positive definite, nearPD from the Matrix package will be applied to find the nearest positive definite matrix. Without a positive definite matrix PGLS will fail to run.

Please make sure both the species compositional matrix and phylogenetic matrix have identical species name (e.g., space denoted as " ", not "_" or ".").

"get_comm_pair_r" will return a variance-covariance matrix at the community level, which can be used in gls (from nlme)

Other functions are used for the simulations / empirical analyses.
