# Calculate the difference of the partial sum of the covariance operators of
# two functional data objects
D.mn = function(lambda, fdata1, fdata2) {
  if(lambda == 0) {return("error, lamba = 0 not possible in D.mn")} else {
    n = length(fdata1$coefs[1,])
    m = length(fdata2$coefs[1,])

    factor1 <- sqrt((floor(lambda*n)-1)/(n-1))
    factor2 <- sqrt((floor(lambda*m)-1)/(m-1))

    cov1 <- var.fd(factor1*fdata1[1:(floor(lambda*n))])
    cov2 <- var.fd(factor2*fdata2[1:(floor(lambda*m))])

    cov1$coefs <- cov1$coefs - cov2$coefs

    return(cov1)
  }
}


# Calculate the normalizer and the statistic for two functional data samples
#   Input:
#     fdata1, fdata2: Functional data object
#
#   Output:
#     Vector where the first entry is the normalizer and the second entry is
#     the statistic
V_n = function(fdata1, fdata2) {

  n = length(fdata1$coefs[1,])

  nPoints = 20
  points = 1:nPoints/nPoints
  intValues = c()

  if(fdata1$basis$type == "fourier") {
    for(k in 1:nPoints) {
      covDiff <- D.mn(points[k], fdata1, fdata2)
      intValues[k] <- sum(covDiff$coefs^2)
    }
  } else {stop("Currently only possible for the fourier basis.")}

  integrand = (intValues-(floor(points*n)/n)^2*intValues[nPoints])^2
  res = sqrt(1/(nPoints-1)*sum(integrand) )
  statistic = intValues[nPoints]
  return(c(res, statistic))
}


H_n = function(fdata1, fdata2) {
  
  n = length(fdata1$coefs[1,])
  
  nPoints = 20
  points = 1:nPoints/nPoints
  intValues = c()
  
  if(fdata1$basis$type == "fourier") {
    for(k in 1:nPoints) {
      covDiff <- D.mn(points[k], fdata1, fdata2)
      intValues[k] <- sum(covDiff$coefs^2)
    }
  } else {stop("Currently only possible for the fourier basis.")}
  
  integrand = (intValues-(floor(points*n)/n)^2*intValues[nPoints])
  res = max(integrand) - min(integrand)
  statistic = intValues[nPoints]
  return(c(res, statistic))
}
# Create a list of lists which contains all the data for the simulation
# Input:
#   nsim: Number which specifies the number of simulation runs, i.e. in this
#         case the number of samples that are generated
#   dataInfo: A list containing the info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c(fIID", "fMA1", "nonGaussian") determining the data type
genDataList = function(nsim, dataInfo) {

  dataInfo2 = dataInfo
  dataInfo2$n = dataInfo$m
  dataInfo1 = dataInfo
  dataInfo1$factor = 1
  datalist = list()
  for(i in 1:nsim) {
    L = list()
    switch(dataInfo$type,
      fIID = {
        # generate independent functional data
        L$fdata1 = fIID(dataInfo1)
        L$fdata2 = fIID(dataInfo2)
      },
      fMA1 = {
        #generate fMA1 data
        L$fdata1 = fMA1(dataInfo1)
        L$fdata2 = fMA1(dataInfo2)
      },
      nonGaussian = {
        L$fdata1 = nonGaussian(dataInfo1)
        L$fdata2 = nonGaussian(dataInfo2)
      }
    )
    datalist[[i]] = L
  }
  return(datalist)
}


# Calculate a (nsim x 2)-matrix where the i-th row contains the normalizer (in
# the first column) and the statistic (in the second column) of the i-th
# simulation run
# Input:
#   simInfo: A list containing general info about the simulation:
#     $nsim: Number which specifies the number of simulation runs (here the
#            number of rows of the output matrix)
#     $cores: Number of cores which are used for parallel computations
#   dataInfo: A list containing the Info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c(fIID", "fMA1", "nonGaussian") determining the data type
calcStatistics = function(simInfo, dataInfo,datalist) {

  nsim = simInfo$nsim
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("V_n","D.mn"))
  clusterEvalQ(cl, library(fda))

  statisticList <-  parallel::parLapply(cl,1:nsim, function(i) V_n(fdata1 = datalist[[i]]$fdata1,
                                                                   fdata2 = datalist[[i]]$fdata2))
  stopCluster(cl)
  statistics = matrix(ncol = 2, nrow = nsim)
  for(i in 1:nsim){statistics[i,] = statisticList[[i]]}

  return(statistics)
}

calcStatistics_hong = function(simInfo, dataInfo,datalist) {
  
  nsim = simInfo$nsim
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("H_n","D.mn"))
  clusterEvalQ(cl, library(fda))
  
  statisticList <-  parallel::parLapply(cl,1:nsim, function(i) H_n(fdata1 = datalist[[i]]$fdata1,
                                                                   fdata2 = datalist[[i]]$fdata2))
  stopCluster(cl)
  statistics = matrix(ncol = 2, nrow = nsim)
  for(i in 1:nsim){statistics[i,] = statisticList[[i]]}
  
  return(statistics)
}
# Specify a vector of numbers which are used as factors for the second
# sample. For each number, the "calcStatistics" function is used and the output
# matrix is written in a csv file in a folder specified by "path"
# Input:
#   path: String containing the path where the csv files should be saved
#   aVector: Vector containing the different factors for the expectation
#   simInfo: A list containing general info about the simulation:
#     $nsim: Number which specifies the number of simulation runs
#     $cores: Number of cores which are used for parallel computations
#   dataInfo: A list containing the Info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c(fIID", "fMA1", "nonGaussian") determining the data type
writeStatistics <- function(path, factorVector, datInfo, simInfo) {
  dataInfo <- datInfo
  for(i in 1:length(factorVector)) {
    dataInfo$factor = factorVector[i]
    nsim <-  simInfo$nsim
    datalist <- genDataList(nsim, dataInfo)
    res <- calcStatistics(simInfo = simInfo, dataInfo = dataInfo,datalist)
    file <- file.path(path, paste0("shao_", i, ".csv"))
    write.table(res, file = file, row.names = FALSE, col.names = FALSE)
  }
}

writeStatistics_hong <- function(path, factorVector, datInfo, simInfo) {
  dataInfo <- datInfo
  for(i in 1:length(factorVector)) {
    dataInfo$factor = factorVector[i]
    nsim <-  simInfo$nsim
    datalist <- genDataList(nsim, dataInfo)
    res <- calcStatistics_hong(simInfo = simInfo, dataInfo = dataInfo,datalist)
    file <- file.path(path, paste0("hong_", i, ".csv"))
    write.table(res, file = file, row.names = FALSE, col.names = FALSE)
  }
}