# BayesianSSD_attrition

This repository contains the R code necessary to perform a simulation study about Bayesian power with attrition in longitudinal trials  using the `BayesSSD` package. Also, it contains a tutorial on how to use this package. The package can be dowloaded from `https://github.com/ulrichlosener/BayesSSD`. 

The code for the simulation study from the paper accepted for publication in *Advances in Methods and Practices in Psychological Science* can be found in the file "simulation_missing.R". Figures 1 and 4 can be reproduced using the code in "survival_hazard_plots.R". The tutorial Markdown file is called "empirical_example.Rmd".

Different Attrition patterns are taken into account by means of survival functions. To get an impression of how these survival functions can model drop-out, see out the corresponding ShinyApp on https://losener.shinyapps.io/attrition-survival/.

For questions and comments, please contact me under u.c.losener1@uu.nl. 
