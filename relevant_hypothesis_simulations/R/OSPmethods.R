
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

# Calculate the normalizer and the statistic for a functional data object
#   Input:
#     fdata: Functional data object
#
#   Output:
#     Vector where the first entry is the normalizer and the second entry is
#     the statistic
V_n = function(fdata) {
  n = length(fdata$coefs[1,])
  
  nPoints = 20
  points = 1:nPoints/nPoints
  intValues = c()
  for(k in 1:length(points)) {
    if (fdata$basis$type == "fourier") {
      # for the fourier basis, the integral of the squared function is the sum
      # of the squared coeficients
      S.n <- S_n(points[k],fdata)
      intValues[k] = sum(S.n$coefs^2)
    } else {
      # calculate the integral numerically
      maxSubintervals = 1000
      S.square = function(t) {eval.fd(t, S_n(points[k], fdata))^2}
      intValues[k] = integrate(S.square, lower = 0, upper = 1,
                               subdivisions = maxSubintervals)$value
    }
  }
  integrand = (intValues-(floor(points*n)/n)^2*intValues[length(points)])^2
  res = sqrt(1/(nPoints-1)*sum(integrand) )
  statistic = intValues[length(points)]
  return(c(res, statistic))
}


H_n = function(fdata) {
  n = length(fdata$coefs[1,])
  
  nPoints = 20
  points = 1:nPoints/nPoints
  intValues = c()
  for(k in 1:length(points)) {
    if (fdata$basis$type == "fourier") {
      # for the fourier basis, the integral of the squared function is the sum
      # of the squared coeficients
      S.n <- S_n(points[k],fdata)
      intValues[k] = sum(S.n$coefs^2)
    } else {
      # calculate the integral numerically
      maxSubintervals = 1000
      S.square = function(t) {eval.fd(t, S_n(points[k], fdata))^2}
      intValues[k] = integrate(S.square, lower = 0, upper = 1,
                               subdivisions = maxSubintervals)$value
    }
  }
  integrand = intValues-(floor(points*n)/n)^2*intValues[length(points)]
  res = max(integrand) - min(integrand)
  statistic = intValues[length(points)]
  return(c(res, statistic))
}
# Create a list of lists which contains all the data for the simulation
# Input:
#   nsim: Number which specifies the number of simulation runs, i.e. in this
#         case the number of samples that are generated
#   dataInfo: A list containing the info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c("BB", "fIID", "fMA1", "fAR1") determining the data type
genDataList = function(nsim, dataInfo) {
  
  datalist = list()
  for(i in 1:nsim) {
    L = list()
    switch(dataInfo$type,
           #generate independent brownian bridges
           BB = { L$fdata = BB(dataInfo) },
           # generate independent functional data
           fIID = { L$fdata = fIID(dataInfo) },
           #generate fMA1 data
           fMA1 = { L$fdata = fMA1(dataInfo) },
           #generate fAR1 data
           fAR1 = { L$fdata = fAR1(dataInfo) }
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
#     $type: One of c("BB", "fIID", "fMA1", "fAR1") determining the data type
calcStatistics = function(simInfo, dataInfo,datalist) {
  
  nsim <-  simInfo$nsim
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("V_n","S_n"))
  clusterEvalQ(cl, library(fda))
  # use multiple cores for the calculation of the test statistics
  statisticList <-  parallel::parLapply(cl,1:nsim, function(i) V_n(fdata = datalist[[i]]$fdata))
  stopCluster(cl)
  statistics = matrix(ncol = 2, nrow = nsim)
  for(i in 1:nsim){statistics[i,] = statisticList[[i]]}
  
  return(statistics)
}

calcStatistics_hong = function(simInfo, dataInfo,datalist) {
  
  nsim <-  simInfo$nsim
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("H_n","S_n"))
  clusterEvalQ(cl, library(fda))
  # use multiple cores for the calculation of the test statistics
  statisticList <-  parallel::parLapply(cl,1:nsim, function(i) H_n(fdata = datalist[[i]]$fdata))
  stopCluster(cl)
  statistics = matrix(ncol = 2, nrow = nsim)
  for(i in 1:nsim){statistics[i,] = statisticList[[i]]}
  
  return(statistics)
}
# Specify a vector of numbers which are used as factors for the expectation
# function of the data. For each number, the "calcStatistics" function is used
# and the output matrix is written in a csv file in a folder specified by "path"
# Input:
#   path: String containing the path where the csv files should be saved
#   aVector: Vector containing the different factors for the expectation
#   simInfo: A list containing general info about the simulation:
#     $nsim: Number which specifies the number of simulation runs
#     $cores: Number of cores which are used for parallel computations
#   dataInfo: A list containing the Info for the data generation. See the
#             documentation for each data type in "genData.R" for more details.
#             Additionally:
#     $type: One of c("BB", "fIID", "fMA1", "fAR1") determining the data type
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

