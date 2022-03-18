## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = getwd())


## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("mamouzgar/hsslda", build_vignettes = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("mamouzgar/hsslda", build_vignettes = TRUE)

## ---- warning=FALSE, include=FALSE--------------------------------------------
library(magrittr)

## ---- warning=FALSE, include=TRUE, results = 'hide', message = FALSE----------
library(hsslda)

## ---- include = FALSE, results = 'hide', echo= FALSE--------------------------
# load("/Users/mamouzgar/phd-projects/Rpackage/hsslda/data/TcellHartmann2020_sampleData.rda")

## ----  eval=FALSE-------------------------------------------------------------
#  data(TcellHartmann2020_sampleData)
#  

## -----------------------------------------------------------------------------
colnames(TcellHartmann2020_sampleData)
head(TcellHartmann2020_sampleData,3)

