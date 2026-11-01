# Precompute the vignettes: the resampling takes too long to run when the
# package is checked, so each vignette is knitted here from its .Rmd.orig
# source and the result (with its figures) is what the package ships.
# Run from the package root with the package installed:
#   Rscript vignettes/precompute.R
old <- setwd("vignettes")
for (f in c("lmmresample", "time-courses", "multiple-comparisons")) {
  knitr::knit(paste0(f, ".Rmd.orig"), paste0(f, ".Rmd"))
}
setwd(old)
