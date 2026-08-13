Data and code to recreate the analyses in 'Non-stationary population dynamics in Pacific salmon through an era of widespread ecological change' by Greenberg et al.

Data:
two datasets are included here, representing the compiled spawner-recruit time-series and metadata for these time-series - in the 'data' folder. Note that the compilation is completed, and may be updated in the future, in this [repository](https://github.com/Pacific-salmon-assess/dfo_salmon_compilation)

Code:
'model_fits.R' provides code to fit Stan models (see 'stan' folder in 'code' - note for cmdstanr these will need to be placed in your respective cmdstan path)
'intext_results_reporting.R' manipulates the posterior estimates to recreate the reported Results. The model fits used in the paper are included in 'outputs/model fits'
'main figures.R' includes code to recreate figures

Appendix
'index.Rmd' in 'docs' provides rmarkdown code for the Supplemental Material/Appendix. An HTML version is available to view [here](https://pacific-salmon-assess.github.io/nonstationary_dynamics/)
