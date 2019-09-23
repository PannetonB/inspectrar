# inspectrar

This R package bundles __InSpectoR__, a GUI based script for optical spectroscopy data  
visualisation and editing. It also contains tools for developing basic chemometrics models.  
Data sets are provided to get you started. A few functions that may be of general 
interest are also provided.

To install in R, use the _**devtools**_ package and the following command:  

__devtools::install_github("pannetonb/inspectrar", build_opts = c("--no-resave-data", "--no-manual"))__

The default build_opts has a "--no-build-vignettes". Removing forces building the vignettes.  

For building the vignettes, the following packages are required:  _**kable**_, _**kableExtra**_ and _**bookdown**_.
Also, make sure a recent version of _**pandoc**_ is installed. This may be installed using the _install.pandoc()_ 
function for the _**installr**_ package.

After installation, you may have to restart **R** or **RStudio**.