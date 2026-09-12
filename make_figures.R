#!/usr/bin/env Rscript
needed <- c('digest','matrixStats','ggplot2','patchwork')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install first: ',paste(missing,collapse=', '))

source('R/figures.R')
make_figures()
