# Simulate data and fit 2-species occupancy model à la Rota et al. (2016) w/ Jags and Nimble

We consider a two-species static occupancy model à la [Rota et al. (2016)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12587). We simulate data from this model, and fit the model to these data using `Jags` and `Nimble`.
The equations below do no display well, you'd rather have a look to the [PDF file](https://github.com/oliviergimenez/bayes2speciesoccupancy/blob/master/simul_rota.pdf). With the [Rmd file](https://github.com/oliviergimenez/bayes2speciesoccupancy/blob/master/simul_rota.Rmd), you can run the code in RStudio and reproduce the results.

Rota and colleagues used `Stan` and we refer to the supplementary material of their paper for further details. Fidino and colleagues extended the Rota's model to multiple seasons, therefore allowing the estimation of local colonization and extinction (dynamic models), and they used `Jags', check out their paper [here](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13117) and the code [there](https://github.com/mfidino/dcom).

