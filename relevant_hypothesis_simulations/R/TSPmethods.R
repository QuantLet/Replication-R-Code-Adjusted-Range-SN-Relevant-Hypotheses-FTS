# Calculate the partial sum of a functional data object
# Input:
#   lambda: Number in (0,1]. The sum of the first up to the floor(lambda*n)-th
#           data point is calculated, where n is the sample size
#   fdata: Functional data object
S_n = function(lambda, fdata) {
  if(lambda == 0) {return("error, lamba = 0 not possible in S_n")} else {
    n = length(fdata$coefs[1,])
    res = 1/n * sum(fdata[1:(floor(lambda*n))])
    return(res)
  }
}

# Calculate the difference of the partial sums of two functional data samples
D_N = function(lambda, fdata1, fdata2) {
  n = length(fdata1$coefs[1,])
  m = length(fdata2$coefs[1,])

  res = S_n(lambda, fdata1)-S_n(lambda, fdata2)

  return(res)
}

H_n = function(fdata1, fdata2, simInfo) {
  n = length(fdata1$coefs[1,])
  
  nPoints = 20
  points = 1:nPoints/nPoints
  intValues = c()
  
  # if the functional data is created with a fourier basis, the integral of the
  # squared fd object is just the sum of the squared coefficients
  if(fdata1$basis$type == "fourier") {
    for(k in 1:length(points)) {
      intValues[k] = sum(D_N(points[k], fdata1, fdata2)$coefs^2)
    }
  } else {
    for(k in 1:length(points)) {
      maxSubintervals = 1000
      D.square = function(t) {eval.fd(t, D_N(points[k], fdata1, fdata2))^2}
      intValues[k] = integrate(D.square, lower = 0, upper = 1,
                               subdivisions = maxSubintervals)$value
    }
  }
  
  switch(simInfo$normalizerType,
         V_n = {
           integrand = (intValues-(floor(points*n)/n)^2*intValues[length(points)])
           normalizer =  max(integrand) - min(integrand)
         },
         # alternative normalizer proposed by the AE
         V_n.star = {
           integrand = (intValues-(floor(points*n)/n)^2*intValues[length(points)])
           normalizer =  max(integrand) - min(integrand)
         },
         V_n.2star = {
           integrand = (intValues-(floor(points*n)/n)^2*intValues[length(points)])
           normalizer =  max(integrand) - min(integrand)
         }
  )
  
  statistic = intValues[length(points)]
  return(c(normalizer, statistic))
}
# Calculate the normalizer and the statistic for two functional data samples
#   Input:
#     fdata1, fdata2: Functional data object
#     simInfo: List where only one entry is important here:
#       $normalizerType: One of c("V_n", "V_n.star", "V_n.2star")
#
#   Output:
#     Vector where the first entry is the normalizer and the second entry is
#     the statistic
V_n = function(fdata1, fdata2, simInfo) {

  n = length(fdata1$coefs[1,])

  nPoints = 20
  points = 1:nPoints/nPoints
  intValues = c()

  # if the functional data is created with a fourier basis, the integral of the
  # squared fd object is just the sum of the squared coefficients
  if(fdata1$basis$type == "fourier") {
    for(k in 1:length(points)) {
      intValues[k] = sum(D_N(points[k], fdata1, fdata2)$coefs^2)
    }
  } else {
    for(k in 1:length(points)) {
      maxSubintervals = 1000
      D.square = function(t) {eval.fd(t, D_N(points[k], fdata1, fdata2))^2}
      intValues[k] = integrate(D.square, lower = 0, upper = 1,
                               subdivisions = maxSubintervals)$value
    }
  }

  switch(simInfo$normalizerType,
    V_n = {
      integrand = (intValues-(floor(points*n)/n)^2*intValues[length(points)])^2
      normalizer = sqrt(1/(nPoints-1)*sum(integrand) )
    },
    # alternative normalizer proposed by the AE
    V_n.star = {
      integrand = abs(intValues-(floor(points*n)/n)^2*intValues[length(points)])
      normalizer = max(integrand)
    },
    V_n.2star = {
      integrand = abs(intValues-(floor(points*n)/n)^2*intValues[length(points)])
      normalizer = 1/(nPoints-1)*sum(integrand)
    }
  )

  statistic = intValues[length(points)]
  return(c(normalizer, statistic))
}


# Create a list of lists which contains all the data for the simulation
# Input:
#   nsim: Number which specifies the number of simulation runs, i.e. in this
#         case the number of samples that are generated
#   dataInfo: A list containing the info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c("BB", "fIID", "fMA1", "fMA1dependent", "nonGaussian")
#            determining the data type
genDataList = function(nsim, dataInfo) {

  dataInfo2 = dataInfo
  dataInfo2$n = dataInfo$m
  dataInfo1 = dataInfo
  dataInfo1$muInfo$a = 0
  datalist = list()
  for(i in 1:nsim) {
    L = list()
    switch(dataInfo$type,
      BB = {
        #generate independent brownian bridges
        L$fdata1 = BB(dataInfo1)
        L$fdata2 = BB(dataInfo2)
      },
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
      fMA1dependent = {
        L = dependentSamples(dataInfo)
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
#     $normalizerType: One of c("V_n", "V_n.star", "V_n.2star")
#   dataInfo: A list containing the Info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c("BB", "fIID", "fMA1", "fMA1dependent", "nonGaussian")
#            determining the data type
calcStatistics = function(simInfo, dataInfo,datalist) {
  
  nsim = simInfo$nsim
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("V_n","S_n","D_N"))
  clusterEvalQ(cl, library(fda))
  # use multiple cores for the calculation of the test statistics
  statisticList <-  parallel::parLapply(cl,1:nsim, function(i) V_n(fdata1 = datalist[[i]]$fdata1,
                                                                   fdata2 = datalist[[i]]$fdata2, simInfo = simInfo))
  stopCluster(cl)
  statistics = matrix(ncol = 2, nrow = nsim)
  for(i in 1:nsim){statistics[i,] = statisticList[[i]]}
  
  return(statistics)
}

calcStatistics_hong = function(simInfo, dataInfo,datalist) {
  
  nsim = simInfo$nsim
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("H_n","S_n","D_N"))
  clusterEvalQ(cl, library(fda))
  # use multiple cores for the calculation of the test statistics
  statisticList <-  parallel::parLapply(cl,1:nsim, function(i) H_n(fdata1 = datalist[[i]]$fdata1,
                                                                   fdata2 = datalist[[i]]$fdata2, simInfo = simInfo))
  stopCluster(cl)
  statistics = matrix(ncol = 2, nrow = nsim)
  for(i in 1:nsim){statistics[i,] = statisticList[[i]]}
  
  return(statistics)
}

# Specify a vector of numbers which are used as factors for the expectation
# function of the second sample. For each number, the "calcStatistics" function
# is used and the output matrix is written in a csv file in a folder specified
# by "path"
# Input:
#   path: String containing the path where the csv files should be saved
#   aVector: Vector containing the different factors for the expectation
#   simInfo: A list containing general info about the simulation:
#     $nsim: Number which specifies the number of simulation runs
#     $cores: Number of cores which are used for parallel computations
#     $normalizerType: One of c("V_n", "V_n.star", "V_n.2star")
#   dataInfo: A list containing the Info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c("BB", "fIID", "fMA1", "fMA1dependent", "nonGaussian")
#            determining the data type
writeStatistics <- function(path, aVector, datInfo, simInfo) {
  dataInfo <- datInfo
  for(i in 1:length(aVector)) {
    dataInfo$muInfo$a = aVector[i]
    nsim <-  simInfo$nsim
    datalist <- genDataList(nsim, dataInfo)
    res <- calcStatistics(simInfo = simInfo, dataInfo = dataInfo,datalist)
    file <- file.path(path, paste0("shao_", i, ".csv"))
    write.table(res, file = file, row.names = FALSE, col.names = FALSE)
  }
}

writeStatistics_hong <- function(path, aVector, datInfo, simInfo) {
  dataInfo <- datInfo
  for(i in 1:length(aVector)) {
    dataInfo$muInfo$a = aVector[i]
    nsim <-  simInfo$nsim
    datalist <- genDataList(nsim, dataInfo)
    res <- calcStatistics_hong(simInfo = simInfo, dataInfo = dataInfo,datalist)
    file <- file.path(path, paste0("hong_", i, ".csv"))
    write.table(res, file = file, row.names = FALSE, col.names = FALSE)
  }
}

