# This script contains functions that have been used for ultra-fast preselection 
# and space-time trajectory generation.

# Author: Dennis van der Meer
# Email: denniswillemvandermeer@gmail.com

# The forecast and observation results are stored in a list with depth 4:
# Days (40) -> FCs and OBS (2) -> Horizon (9) --> Station (17)

# setwd("C:/Users/SharedComputer/Documents/dennis/UltraFastPreselection")
setwd("~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection")
rm(list = ls(all = TRUE)) 
libs <- c("doSNOW", "data.table","scoringRules",
          "seismicRoll", 'xtable',"quantreg","ggplot2","copula",
          "cowplot","viridis")
invisible(lapply(libs, library, character.only = TRUE))

# Load necessary functions
source(file.path(getwd(), "/my_functions_git.R"))

# Read Hawaii_1min.RData that has been cleaned before (mainly by removing values with
# SE < 10 and setting clear-sky index values higher than 1.75 to 1.75).
# load("C:/Users/SharedComputer/Documents/dennis/UltraFastPreselection/Hawaii_1min.RData")
load("/Volumes/G-DRIVE_USB-C/Data/Hawaii/Hawaii_1min.RData")
station.names <- colnames(Data1min[1:17])

#################################################################################
# Ultra-fast preselection forecasts
#################################################################################
Days <- seq.Date(from = as.Date("2010-07-01"), to = as.Date("2010-08-09"), by = "day") # 40 Days
K <- 9 
Results <- list()
pb <- txtProgressBar(min = 0, max = length(Days), style = 3)
j <- 0
ptm <- proc.time() 
for(d in 1:length(Days)) {              # Loop over the days
  X <- Data1min[as.Date(Data1min$Time, format = "%Y-%m-%d") == Days[d],]
  Tm <- X$Time                          # Save the time 
  X <- subset(X, select = -Time)        # Deselect time for the algorithm
  station.names <- colnames(X)
  FCs <- UltraFastFCs(K = K, dat = X, station.names = station.names, model = "QR")
  Obs <- Observations(K = 1, Time = Tm, dat = X, station.names = station.names) # Extract the observations including time
  Results[[d]] <- list(FCs,Obs)         # This is a list of days --> Res,Obs --> horizons --> locations
  rm(FCs,Obs)                           # Just to be sure
  j <- j+1
  setTxtProgressBar(pb, j)              # Update the progress bar
}
proc.time() - ptm

#######################################################################################
# Univariate prediction (use only 4 time lags)
#######################################################################################
n_l <- 4 
y.test = fcs = FH.stat.k = ME.stat.k <- list()
Days <- seq.Date(from = as.Date("2010-07-01"), to = as.Date("2010-08-09"), by = "day") # 40 Days
Days <- seq.Date(from = as.Date("2010-07-01"), to = as.Date("2010-07-02"), by = "day") # 2 Days
K <- 2
Results <- list()
pb <- txtProgressBar(min = 0, max = length(Days), style = 3)
j <- 0
m <-  350                 # Not used herein but to ensure that the forecasts of the 
nt <- 50                  # QR^{endo} and QR^{exo} models are aligned time.
taus <- seq(0.05,0.95,0.05)
ptm <- proc.time()
for(d in 1:length(Days)) {              # Loop over the days
  X <- Data1min[as.Date(Data1min$Time, format = "%Y-%m-%d") == Days[d],]
  Tm <- X$Time                          # Save the time
  X <- subset(X, select = -Time)        # Deselect time for the algorithm
  station.names <- colnames(X)
  
  for(k in 1:K){
    for(s0 in 1:length(station.names)) { # Loop over the stations
      Y             <- as.data.frame(X[,s0])
      Y[,(ncol(Y)+1):(ncol(Y)+n_l)] <- data.table::shift(Y, n = c(k:(n_l+k-1)), type = "lag") 
      colnames(Y)[1]<- c("Y")
      X.train       <- Y[1:(m+nt),]
      fit           <- quantreg::rq(formula = Y~., data = X.train, tau = taus)
      X.test        <- Y[(m+nt):nrow(Y),2:ncol(Y)]
      
      y.test[[s0]]  <- Y[(m+nt):nrow(X),1]
      fcs[[s0]]     <- predict(object = fit, newdata = X.test)
    }
    FH.stat.k[[k]] <- fcs
    ME.stat.k[[k]] <- y.test
  }
  Results[[d]] <- list(FH.stat.k,ME.stat.k) # This is a list of days --> Res,Obs --> horizons --> locations
  j <- j+1
  setTxtProgressBar(pb, j) # Update the progress bar
}
proc.time() - ptm

# my.CRPS <- DailyCRPS(lstOfResults = Results, station.names = station.names)
# p.rel <- Plot.ReliabilityDiagram(Results,station.names,K)
# p.sharp <- Plot.Sharpness(Results,station.names,K)
# 
# ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/CRPS_QR_NoPreselection.pdf", plot = my.CRPS[[2]], 
#        device = "pdf", units = "cm", height = 8.5, width = 8.5) 
# ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/Reliability_QR_NoPreselection.pdf", plot = p.rel, 
#        device = "pdf", units = "cm", height = 8.5, width = 8.5) 
# ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/Sharpness_QR_NoPreselection.pdf", plot = p.sharp, 
#        device = "pdf", units = "cm", height = 8.5, width = 8.5) 

####################################### Generate scenarios #######################################
# Load the forecast results:
# load("FC_Res_QR_Preselection.RData")
# load("FC_Res_QR_NoPreselection.RData")

# The following functions can be downloaded from Thorarinsdottir et al.:
# https://www.tandfonline.com/doi/suppl/10.1080/10618600.2014.977447?scroll=top
source("ucgs_a_977447_sm3893/CreateFiguresAndTables.R")
S <- 5  # Number of scenarios (100 in the paper)
copulas <- c("empirical","normal") # ,"clayton","t"
# copulas <- c("high","low","biased")
runTimes <- NULL
ES.t <- as.data.frame(matrix(data = NA, nrow = length(Days), ncol = 4))
colnames(ES.t) <- copulas
c <- 1
for(cop in copulas){
  print(paste("Working on",cop,"copula",sep = " "))
  res.scen <- cop_fit(S = S, K = K, Results = Results, Data1min = Data1min, copula = cop, Days = Days) # Empirical copula is fastest
  runTimes <- cbind(runTimes,res.scen[[2]])
  res.scen <- res.scen[[1]]
  # Cast the scenarios in such a format that we can use the rank histograms plot functions of Thorarinsdottir:
  res <- list()
  for(d in 1:length(Results)){
    Obs <- Results[[d]][[2]] # First day, then observations
    obs <- do.call(cbind,Obs[[1]])
    obs <- obs[,-18] # Remove the time column and start from K because of NAs introduced by lagging --> latter not necessary because I select the correct row
    lst <- list()
    for(i in K:length(res.scen[[d]])){           # Loop over the time stamps for which the scenarios were created starting from K because of NAs
      tmp <- t(flatten(res.scen[[d]][[i]]))      # Then, flatten the 3D array to a S x (ns*K) matrix, which is what the multivariate rank histogram functions expects
      tmp1 <- rbind( rep(obs[i,], K), tmp)       # Then, add a row that contains the observations (which has to be repeated of course)
      lst[[i]] <- tmp1                           
    }
    lst[sapply(lst, is.null)] <- NULL # Remove the NULL entries in lst introduced because we started from i=K
    res[[d]] <- lst
    # Calculate energy score
    tmp <- Calc.ES(scen.lst = res.scen[[d]], obs.lst = Obs, K = K) 
    ES.t[d,c] <- mean(tmp, na.rm = T)
  }
  c <- c+1
  res <- unlist(res, recursive = F)
  ARH <- avg.rhist(B = res, reps = length(res), hist_xlab="Average Rank") 
  BDR <- bd.rhist(B = res, reps = length(res), hist_xlab="Band Depth Rank")
  # Save the results
  save(ARH, file = paste("RankHistograms/ARH_",cop,"_Preselection",".RData",sep = ""))
  save(BDR, file = paste("RankHistograms/BDR_",cop,"_Preselection",".RData",sep = ""))
}
write.table(ES.t, file = "energyScore.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
colnames(runTimes) <- copulas
