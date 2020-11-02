#### compile_pkg.R ####
# This file is meant for recompilation during development.
# For regular installation use devtools::install_github("user/pkgname")
# Use build_vignettes = TRUE in install_github to get html vignettes

library(roxygen2)
library(devtools)

rm(list = ls()) # cleans current environment
gc() # releases memory removed above

# pkgbuild::compile_dll() #run if DLL error occurs or if new C functions don't appear
# create NAMESPACE and .Rd files for documentation. Compile C files if any.
devtools::document()

# install preliminary package for function testing and documentation inspection
# setting build_vignettes = FALSE will considerably speed up installation.
# if built with vignettes, use browseVignettes("pkgname") to view.
devtools::install(build_vignettes = FALSE)

# devtools::build(binary = TRUE)
