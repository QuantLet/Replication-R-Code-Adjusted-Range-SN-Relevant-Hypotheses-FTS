#!/usr/bin/env Rscript
needed <- c('digest','matrixStats','qrng','spacefillr')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install first: ',paste(missing,collapse=', '))

source('R/primary_simulations.R')
run_primary_simulations()
run_loading_simulations()
run_mean_simulations()

source('R/size_adjusted_power.R')
run_size_adjusted_power()

source('R/variance_diagnostics.R')
run_variance_diagnostics()
