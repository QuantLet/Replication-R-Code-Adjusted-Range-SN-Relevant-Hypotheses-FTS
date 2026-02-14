# Functions for calculating quantiles of the limit distribution and
# empirical rejection probabilities from the statistics which are
# calculated in advance and saved as csv files


# Calculate one realization of the limit distribution for one Brownian motion
# process
# Input:
#   fBrownian: Functional data object containing a Brownian motion
#   simInfo: List which is here only needed for the two sample problem:
#     $normalizerType: One of c("V_n", "V_n.star", "V_n.2star")
W = function(fBrownian, simInfo = NULL) {

  if(is.null(simInfo$normalizerType)) {
    normalizerType = "V_n"
  } else {normalizerType = simInfo$normalizerType}

  B.1 = eval.fd(1, fBrownian)[1,]

  nPoints = 20
  fct = function(lambda) {lambda * (eval.fd(lambda, fBrownian)-lambda*B.1)}
  points = 1:nPoints/nPoints
  values = c()
  for(i in 1:length(points)) {
    values[i] = fct(points[i])
  }
  switch(normalizerType,
    V_n = {
      res = sqrt(1/(nPoints-1)*sum(values^2))
    },
    # limit distribution when using the alternative normalizer proposed by the AE
    V_n.star = {
      res = max(abs(values))
    },
    V_n.2star = {
      res = 1/(nPoints-1)*sum(abs(values))
    }
  )
  return(B.1/res)
}

W_hong = function(fBrownian, simInfo = NULL) {
  
  if(is.null(simInfo$normalizerType)) {
    normalizerType = "V_n"
  } else {normalizerType = simInfo$normalizerType}
  
  B.1 = eval.fd(1, fBrownian)[1,]
  
  nPoints = 20
  fct = function(lambda) {lambda * (eval.fd(lambda, fBrownian)-lambda*B.1)}
  points = 1:nPoints/nPoints
  values = c()
  for(i in 1:length(points)) {
    values[i] = fct(points[i])
  }
  switch(normalizerType,
         V_n = {
           res = max(values) - min(values)
         },
         # limit distribution when using the alternative normalizer proposed by the AE
         V_n.star = {
           res = max(values) - min(values)
         },
         V_n.2star = {
           res = max(values) - min(values)
         }
  )
  return(B.1/res)
}

# Calculate quantiles of the limit distribution
# Input:
#   fBrownianSample: Functional data object containing a sample of Brownian
#                    motions
#   simInfo: List which is here only needed for the two sample problem:
#     $normalizerType: One of c("V_n", "V_n.star", "V_n.2star")
quant.fct = function(fBrownianSample, simInfo = NULL, prob = 0.95) {
  nsim = length(fBrownianSample$coefs[1,])
  qVector = c()
  for(i in 1:nsim) {
    qVector[i] = W(fBrownianSample[i], simInfo)
  }

  q = quantile(qVector , prob)
  return(q)
}

quant.fct_hong = function(fBrownianSample, simInfo = NULL, prob = 0.95) {
  nsim = length(fBrownianSample$coefs[1,])
  qVector = c()
  for(i in 1:nsim) {
    qVector[i] = W_hong(fBrownianSample[i], simInfo)
  }
  
  q = quantile(qVector , prob)
  return(q)
}




# Calculate a list where each entry contains the matrix of the statistics
# and normalizers the test decision of each simulation runs
# Input:
#   list: List where each entry contains a (nsim x 2)-matrix with the statistics
#         and normalizers (output of "readTables")
#   Delta: Real number specifying the threshold parameter for the relevant
#          hypotheses
#   q: Real number specifying the quantile
calcRejectionList = function(list, Delta, q) {
  ntables = length(list)
  rejectionList = list()
  for(i in 1:ntables) {
    normalizers = list[[i]][,1]
    statistics = list[[i]][,2]
    rejectionList[[i]] = statistics > Delta + q*normalizers
  }
  return(rejectionList)
}


# Calculate empirical rejection probabilities
# Input:
#   list: List where each entry contains a (nsim x 2)-matrix with the statistics
#         and normalizers (output of "readTables")
#   Delta: Real number specifying the threshold parameter for the relevant
#          hypotheses
#   q: Real number specifying the quantile
empRejProb = function(list, Delta, q) {

  rejectionList = calcRejectionList(list, Delta, q)
  ntables = length(rejectionList)
  empRejProbs = c()
  for(i in 1:ntables) {
    empRejProbs[i] = mean(rejectionList[[i]])
  }
  return(empRejProbs)
}


# Function for reading the statistics and normalizers from the csv files
# Input:
#   path: String which specifies the path to the csv files
#   ntables: Number of csv files
readTables = function(path, ntables) {
  statisticsList = list()
  for(i in 1:ntables) {
    statisticsList[[i]] = read.table(file.path(path, paste0("shao_", i, ".csv")))
  }
  return(statisticsList)
}


readTables_hong = function(path, ntables) {
  statisticsList = list()
  for(i in 1:ntables) {
    statisticsList[[i]] = read.table(file.path(path, paste0("hong_", i, ".csv")))
  }
  return(statisticsList)
}
