# inspectrar

**NOTE: *inspectrar* is no longer support**. The package *Gwidgets2* no longer works with the *RGtk2* toolkit that was used to build the graphical interface. The supported *tcltk* tool kit does not work as support for interacting with graphic objects is quite different when compared to *GTk2*. (as of April 2022). A Shiny app version is under construction and may be released in the coming months. If you need to use the package, it runs smoothly on **R version 4.0.3** with versions 1.0-8, 1.0-7 of the ***gWidgets2*** and ***gWidgets2RGtk2*** packages respectively.

This R package bundles **InSpectoR**, a GUI based script for optical spectroscopy data  
visualisation and editing. It also contains tools for developing basic chemometrics models.  
Data sets are provided to get you started. A few functions that may be of general interest are also provided.

To install in R, use the ***devtools*** package and the following command:

First for building the vignettes, the following packages are required: ***knitr***, ***kableExtra*** and ***bookdown***. Also, make sure a recent version of ***pandoc*** is installed. This may be installed using the *install.pandoc()* function for the ***installr*** package. Finally, make sure you have the latest version of package ***rlang*** installed before installing from github.

Then, you can install ***inspectrar*** using the following function from the ***devtools*** package:

**devtools::install_github("pannetonb/inspectrar", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)**

The default build_opts has a "--no-build-vignettes". Removing forces building the vignettes.

After installation, you may have to restart **R** or **RStudio**.
