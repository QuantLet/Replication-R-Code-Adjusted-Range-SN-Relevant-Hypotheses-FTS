#!/usr/bin/env Rscript
needed <- c('digest','matrixStats')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install first: ',paste(missing,collapse=', '))

source('R/empirical_analysis.R')
run_empirical_analysis()
