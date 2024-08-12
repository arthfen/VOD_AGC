# VOD_AGC
This repository contains Supplementary Data for the manuscript "Human influence on Amazon's aboveground carbon dynamics intensified over the last decade".

## Content
* The R package developed to perform declaring the disaggregation model's log-likelihood function (folder *disag.v8*)
* The R scripts used to estimate the disaggregation model parameters

## Instructions
1. The R package can be instaled with:
```sh
cd final-high-step1/disag.v8
R CMD INSTALL .
```

2. To run the model, you must run the scripts of *final-high-step1* in order.
A job script *0job.sh* is an example of how the model was ran in the LSCE's cluster Obelix.
In that case, you could submit the job as:
```sh
cd final-high-step1
qsub 0job.sh
```
If the file *0job.sh* is used, then the variance estimatino (i.e., folder *final-high-step4*) will follow automatically.

## Important remarks
The corresponding Zenodo repositories contains files 3out_mean_{year}.tif and 3out_sd_{year}.tif.

**These files should not be interpreted the mean and standard deviation of the model predictions**.

Instead, they were obtained as follows:
1. After model estimation, a sample of size 500 was taken from the posterior distribution of model parameters.
2. For each pixel, the model matrix multiplied these samples and the mean μ, and standard deviation, σ, were calculated.
3. File 3out_mean_{year}.tif corresponds to exp(μ), and file 3out_sd_{year}.tif, to σ.

Therefore, the expected value (i.e., the expected AGC in MgC/ha) in each pixel follows a [log-normal distribution]([https://twitter.com/your_username](https://en.wikipedia.org/wiki/Log-normal_distribution)) with parameters μ and σ.
The file 3out_mean_{year}.tif, used for calculations in the paper, correspond to the median of such a distribution, but the uncertainty is also provided in case users are interested in other quantities (i.e., confidence intervals, for example).
