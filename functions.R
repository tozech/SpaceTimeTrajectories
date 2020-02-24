# This script contains functions that have been used for ultra-fast preselection 
# and space-time trajectory generation.

# Author: Dennis van der Meer
# Email: denniswillemvandermeer@gmail.com

# The forecast and observation results are stored in a list with depth 4:
# Days (40) -> FCs and OBS (2) -> Horizon (9) --> Station (17)

#################################################################################
# Plot PICP and PINAW function
#################################################################################
PIs.tmp <- function(ens,obs,title=NULL) {
  
  stopifnot(nrow(ens) == length(obs))
  
  if(!is.null(title)) {title <- title} # Optional: add a title to the plot
  
  idx <- complete.cases(cbind(ens,obs))     # Find rows that have missing values
  ens <- ens[idx,]                          # and remove rows that contain these
  obs <- obs[idx]                           # also in the observations.
  ens <- t(apply(ens, 1, sort))             # Make sure there's no quantile crossing (requires sum(is.na)=0)
  
  R <- 1 # max(obs) - min(obs)
  PINAW <- matrix(0, nrow = floor(ncol(ens)/2), ncol = 1)
  for(c in 1:floor(ncol(ens)/2)) {
    PINAW[c,] <- (1/R) * mean(ens[,(ncol(ens)-c+1)]-ens[,c])
  }
  
  PICP <- matrix(0, nrow = floor(ncol(ens)/2), ncol = 1)
  for(c in 1:floor(ncol(ens)/2)) {
    PICP[c,] <- mean(ifelse(ens[,c] <= obs & obs <= ens[,(ncol(ens)-c+1)], 1, 0))
  }
  
  taus <- seq(0,1,(1/(length(PICP)+1)))
  taus <- sort(taus[2:(length(taus)-1)], decreasing = TRUE)
  Res <- data.frame(PINAW = PINAW, PICP = PICP, eta = taus)
  
  p.PINAW <- ggplot(Res, aes(x=eta,y=PINAW)) +
    geom_line(size=1) +
    ggtitle(title) +
    theme(axis.text=element_text(size=24),
          axis.title=element_text(size=24)) +
    theme(plot.title = element_text(hjust = 0.5, size = 24))
  
  p.PICP  <- ggplot(Res, aes(x=eta,y=PICP)) + geom_line(size=1) +
    geom_abline(slope = 1, intercept = 0, linetype=3,size=1) +
    ggtitle(title) +
    theme(axis.text=element_text(size=24),
          axis.title=element_text(size=24)) +
    theme(plot.title = element_text(hjust = 0.5, size = 24))
  
  out <- list(Res,p.PINAW,p.PICP)
  names(out) <- c("Results", "PINAW", "PICP")
  out
}
#################################################################################
# Energy score function (according to scoringRules)
#################################################################################
# requires a list of scenarios per time instance and an observation
# list per time instant
Calc.ES <- function(scen.lst,obs.lst,K){ # K for the maximum forecast horizon because we hop over it
  obs <- do.call(cbind,obs.lst[[1]])
  if(ncol(obs) > 17){obs <- obs[,-18]} # Preselection observations contain time vector
  ES.t <- matrix(data = NA, nrow = nrow(obs), ncol = 1)
  for(i in 1:length(scen.lst)){

    dat <- matrix( scen.lst[[i]], nrow = nrow(scen.lst[[i]]) )
    ES.t[i,1] <- scoringRules::es_sample(obs[i,],dat)

  }
  return(ES.t)
}
#################################################################################
# Function to generate scenarios
#################################################################################
# Function to generate scenarios with the correlated random numbers:
# For each time instant, the function stores a 3D array (ns x S x K) in a list.
# The inputs are therefore: a list of lists with the probabilistic forecasts where the elements of the
# first list are forecast horizons and that of the nested list are ns matrices of (T x taus).

scen.gen <- function(fcst.lst, S, K, ns, cop_model){ # generate scenarios
  # fcst.lst  = list with forecasts
  # S         = number of scenarios
  # K         = number of forecast horizons
  # ns        = number of stations
  # cop_model = the copula model
  
  test.length <- nrow(fcst.lst[[1]][[1]])
  my.ar <- array(data = NA, dim = c(ns,S,K))  # Store the scenarios in a 3D array
  lst <- list()                               # and store the scenarios in a list per time instant
  cl <- makeCluster(8)
  registerDoSNOW(cl)

  lst <- foreach(t = 1:test.length, .packages = c("copula")) %dopar% { # , .options.snow = opts
    set.seed(123)
    # # Some conditional statements to generate Fig. 5:
    # if(cop_model=="high"){
    #   rnd <- runif(S)
    #   ex <- replicate(ns*K,rnd)
    # } else if(cop_model=="biased"){
    #   rnd <- runif(S, min = 0, max = 0.25)
    #   ex <- replicate(ns*K,rnd)
    # } else if(cop_model=="low"){
    #   ex <- matrix(data = runif(S*ns*K), nrow = S) # Sample from a completely random distribution
    # } else {
    #   ex <- rCopula(S, cop_model) # Sample from the copula at each time step
    # }
    ex <- rCopula(S, cop_model) # Sample from the copula at each time step
    
    for(k in 1:K) {             # Loop over the forecast horizons
      for(s0 in 1:ns) {         # Loop over the stations
        for(s in 1:S) {         # Loop over the scenarios
          if(t <= K){t <- K+1}  # skip t <= K because those are NA
          my.ar[s0,s,k] <- approx(x = seq(0.05,0.95,0.05), y = sort(fcst.lst[[k]][[s0]][t,]),      
                                  xout = ex[s,s0+(k-1)*ns], yleft = min(fcst.lst[[k]][[s0]][t,]),  
                                  yright = max(fcst.lst[[k]][[s0]][t,]))$y
        }
      }
    }
    lst[[t]] <- my.ar
  }
  stopCluster(cl)
  return(lst)
}
#################################################################################
# Function to flatten 3d arrays
#################################################################################
# This basically pastes the scenarios into rows and stations * forecast horizons.
# That is what the multivariate rank histogram functions expect.
flatten <- function(ar){ # Function to create a DF of an 3D array
  dyn.GMM <- as.data.frame(matrix(0, nrow = nrow(ar)*length(ar)/(nrow(ar)*ncol(ar)), ncol = ncol(ar)))
  t <- nrow(ar)
  for(i in 1:(length(ar)/(nrow(ar)*ncol(ar)))) {  # These are all the days to loop over
    dyn.GMM[(t*i-t+1):(t*i),] <- ar[,,i]
  }
  return(dyn.GMM)
}
#################################################################################
# Functions required by the ultra-fast preselection algorithm
#################################################################################
zNorm <- function(x)
{
  x <- (x - mean(x)) / (sd(x)*sqrt((length(x)-1)/length(x))) 
  return(x)
}

mass2 <- function(x, query)
{
  require("seismicRoll")
  n <- length(x)
  m <- length(query)
  Q <- zNorm(query)
  #moving standard deviation
  S <- seismicRoll::roll_sd(x, m, align = "left")[1:(n-m+1)]
  S <- S*sqrt((m-1)/m) #get population sd
  Q <- rev(Q)
  Q[(m+1):n] <- 0
  dots <- Re(fft(fft(x)*fft(Q), inverse = TRUE)/n)
  distProfile <- sqrt(2*abs(m-dots[m:n]/S)) #abs counters the (possible) small negative value, due to 0 subtract 0, at the query index
  return(distProfile)
}

#################################################################################
# Function to parallel computing of the probabilistic forecasts per day "dat" (X)
#################################################################################

UltraFastFCs <- function(K, dat, station.names, model="QR"){

  zNorm <- function(x)
  {
    x <- (x - mean(x)) / (sd(x)*sqrt((length(x)-1)/length(x))) 
    return(x)
  }
  
  mass2 <- function(x, query)
  {
    require("seismicRoll")
    n <- length(x)
    m <- length(query)
    Q <- zNorm(query)
    #moving standard deviation
    S <- seismicRoll::roll_sd(x, m, align = "left")[1:(n-m+1)]
    S <- S*sqrt((m-1)/m) #get population sd
    Q <- rev(Q)
    Q[(m+1):n] <- 0
    dots <- Re(fft(fft(x)*fft(Q), inverse = TRUE)/n)
    distProfile <- sqrt(2*abs(m-dots[m:n]/S)) #abs counters the (possible) small negative value, due to 0 subtract 0, at the query index
    return(distProfile)
  }
  
  X <- dat
  m <-  350                 # query length
  ns <- ncol(X)             # number of stations
  nt <- 50                  # number of time lags
  start <- nt + 1           # 
  end <- nrow(X) - m + 1    # 
  taus <- seq(0.05,0.95,0.05)
  
  meas = pred = FH.k = ME.k <- list()
  
  cl <- makeCluster(6)
  registerDoSNOW(cl)
  FH.k <- foreach(k = 1:K, .packages = c("quantreg","data.table","glmnet","quantregForest")) %dopar% {

    X.design <- data.table::shift(x = X, n = k, fill = NA, type = "lag")     # Lag the entire matrix with 1
    X.design <- as.data.frame(do.call(what = cbind, X.design))          # This creates a list of vectors so cbind them again
    colnames(X.design) <- station.names

    for(s0 in 1:ns){ # Loop over the stations 
      
      # Start the procedure
      Y = X[,s0]  # Take the response variable from the original DF           
      fcst <- matrix(NA, nrow = length(Y), ncol = length(taus))
      
      for(i in (start+k):(end-1)){ # Loop over the training periods
        
        ##################################################################################
        # The main preselection routine
        ##################################################################################
        query <- Y[i:(i+m-1)]             # Take a slice of the response variable (query)
        data <- X.design[(i-nt):(i+m-2),] # Take a slice of the X including nt time lags
        x <- as.numeric(unlist(data))     # Convert matrix into vector
        dp <- mass2(x = x, query = query) # Distance profile
        #masking
        joints <- 1:(ns-1) * nrow(data)   # Indices of the joining points
        excZoneStart <- joints - m + 2
        excZoneEnd <- joints
        for(l in 1:(ns-1))
        {
          dp[excZoneStart[l]:excZoneEnd[l]] <- Inf
        }
        
        tmp <- order(dp)[1:20]            # Select 20 predictor out of this procedure

        station <- ceiling(tmp/nrow(data))# station index for the top 20 predictors
        lag <- nt - tmp%%nrow(data) + 1   # lag for the top 20 predictors
        ##################################################################################
        # This ends the preselection routine
        ##################################################################################
        
        ##################################################################################
        # Forecast with the preselection
        ##################################################################################
        predictor = predictor.new <- NULL
        for(l in 1:length(station))
        {
          predictor <- cbind(predictor, X.design[i:(i+m-1) - lag[l], station[l]])
          predictor.new <- append(predictor.new, X.design[i + m - lag[l], station[l]])
        }
        # Default model is QR, experiments with other models is possible
        if(model == "QR"){
          # QR:
          dat <- as.data.frame(predictor)
          nms <- colnames(dat)
          dat <- cbind(query,dat) 
          fit <- quantreg::rq(formula = query ~ ., tau = taus, data = dat)
          predictor.new <- as.data.frame(matrix(predictor.new,nrow = 1))
          colnames(predictor.new) <- nms
          fcst[(i+m),] <- predict(fit, predictor.new)
        } else if(model =="QRF"){
          dat <- as.data.frame(predictor)
          nms <- colnames(dat)
          fit <- quantregForest(y = query, x = dat, nthreads = 4, ntree = 1000)
          predictor.new <- as.data.frame(matrix(predictor.new,nrow = 1))
          colnames(predictor.new) <- nms
          fcst[(i+m),] <- predict(object = fit, newdata = predictor.new, what = taus)
        } else {
          # QR Lasso:
          dat <- as.data.frame(predictor)
          nms <- colnames(dat)
          dat <- cbind(query,dat)
          # The following bit is to ensure that there is no singular design matrix
          # that is sometimes caused by the lambda found in the cv procedure.
          set.seed(123)
          cv <- cv.glmnet(predictor, query, nfolds = 3)
          fit <- try(quantreg::rq(formula = query ~ ., tau = taus, data = dat,
                    method = "lasso", lambda = cv$lambda.min), silent = TRUE)
          
          attempt <- 0
          while(class(fit) == "try-error") {
            set.seed(attempt)
            cv <- cv.glmnet(predictor, query, nfolds = 3)
            fit <- try(quantreg::rq(formula = query ~ ., tau = taus, data = dat,
                          method = "lasso", lambda = cv$lambda.min), silent = TRUE)
            attempt <- attempt + 1
          }

          predictor.new <- as.data.frame(matrix(predictor.new,nrow = 1))
          colnames(predictor.new) <- nms
          fcst[(i+m),] <- predict(fit, predictor.new)
        }
      }
      
      pred[[s0]] <- fcst[(start+m):length(Y),]
      
      rm(fcst,Y)
    }
    
    FH.k[[k]] <- pred
    
  }
  
  stopCluster(cl)
  
  return(FH.k)
  
}
#################################################################################
# Function to extract the measurements of X (because that won't work in parallel)
#################################################################################
Observations <- function(K, Time, dat, station.names) {
  X <- dat
  X <- dat
  m <-  350                 # query length
  ns <- ncol(X)             # number of stations
  nt <- 50                  # number of time lags
  start <- nt + 1           # 
  end <- nrow(X) - m + 1    # 
  taus <- seq(0.05,0.95,0.05)
  
  meas = ME.k <- list()
  for(k in 1:K){ # Loop over the forecast horizons
    
    X.design <- shift(x = X, n = k, fill = NA, type = "lag")     # Lag the entire matrix with 1
    X.design <- as.data.frame(do.call(what = cbind, X.design))          # This creates a list of vectors so cbind them again
    colnames(X.design) <- station.names
    
    for(s0 in 1:ns){ # Loop over the stations 
      
      Y = X[,s0]  # Take the response variable from the original DF           
      
      meas[[s0]] <- Y[(start+m):length(Y)] # * Ioh[(start+m):length(Y)]     # Let's keep it to CSI for now.
      
      rm(Y)
    }
    meas$Time <- Tm[(start+m):length(Tm)]
    ME.k[[k]] <- meas
    
  }
  return(ME.k)
}
#################################################################################
# Calculate the uncertainty of a time series (according to Lauret et al. 2019)
#################################################################################
UNCERTAINTY <- function(obs){
  tmp <- numeric()
  for(i in 1:length(obs)){
    tmp1 <- numeric()
    for(j in 1:i){
      tmp1[j] <- abs(obs[i] - obs[j])
    }
    tmp1 <- sum(tmp1)
    tmp[i] <- tmp1
  }
  return(sum(tmp)/length(obs)^2)
}
#################################################################################
# Function to plot daily averaged CRPS in a raster (Fig. 8 in the manuscript)
#################################################################################
DailyCRPS <- function(lstOfResults, station.names){
  lst_crps <- list()
  lst_resolution <- list()
  lst_reliability <- list()
  pb <- txtProgressBar(min = 0, max = length(lstOfResults), style = 3)
  j <- 0
  for(d in 1:length(lstOfResults)){ # Loop over the days
    Daily <- lstOfResults[[d]]
    FCs <- Daily[[1]]
    OBS <- Daily[[2]]
    
    tmp_crps <- matrix(data = NA, nrow = 9, ncol = 17)
    tmp_resolution <- matrix(data = NA, nrow = 9, ncol = 17)
    tmp_reliability <- matrix(data = NA, nrow = 9, ncol = 17)
    for(k in 1:length(FCs)){ # Loop over the forecast horizons
      for(s0 in 1:length(station.names)){ # Loop over the stations
        res <- verification::crpsDecomposition(OBS[[1]][[s0]][(k+1):length(OBS[[1]][[s0]])], FCs[[k]][[s0]][(k+1):nrow(FCs[[1]][[s0]]),])
        uncertainty <- UNCERTAINTY(OBS[[1]][[s0]][(k+1):length(OBS[[1]][[s0]])]) 
        tmp_resolution[k, s0] <- uncertainty - res$CRPSpot
        tmp_crps[k,s0] <- res$CRPS
        tmp_reliability[k, s0] <- res$Reli
      }
    }
    lst_crps[[d]] <- tmp_crps
    lst_resolution[[d]] <- tmp_resolution
    lst_reliability[[d]] <- tmp_reliability
    j <- j+1
    setTxtProgressBar(pb, j) # Update the progress bar
  }
  CRPS <- Reduce("+", lst_crps) / length(lst_crps) # Calculate the average over the days
  Resolution <- Reduce("+", lst_resolution) / length(lst_resolution) # Calculate the average over the days
  Reliability <- Reduce("+", lst_reliability) / length(lst_reliability) # Calculate the average over the days
  
  CRPS <- as.data.frame(CRPS); Resolution <- as.data.frame(Resolution); Reliability <- as.data.frame(Reliability)
  colnames(CRPS) <- station.names; colnames(Resolution) <- station.names; colnames(Reliability) <- station.names
  CRPS$Horizon <- seq(1,nrow(CRPS),1); Resolution$Horizon <- seq(1,nrow(Resolution),1); Reliability$Horizon <- seq(1,nrow(Reliability),1); 
  mydf <- melt(CRPS, id.vars="Horizon"); resolution_melt <- melt(Resolution, id.vars="Horizon"); reliability_melt <- melt(Reliability, id.vars="Horizon"); 
  colnames(mydf)[3] <- "CRPS"; colnames(resolution_melt)[3] <- "Resolution"; colnames(reliability_melt)[3] <- "Reliability"
  # Prepare for the plot
  plot.size = 8; line.size = 0.1; point.size = 0.6
  # Divide the scale into quantiles
  # Adapted from https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
  quantiles <- quantile(mydf$CRPS, probs = seq(0,1,length.out = 6), na.rm = T)
  labels <- c()
  brks <- seq(0.07,0.17,0.01) # breaks maximum and minimum of data
  # round the labels for the legend
  for(idx in 1:length(brks)){
    labels <- c(labels,round(brks[idx + 1], 2))
  }
  labels <- labels[1:length(labels)-1]
  mydf$CRPS_quantiles <- cut(mydf$CRPS, 
                             breaks = brks, 
                             include.lowest = TRUE, 
                             labels = labels)
  dwnwnd <- c("AP7", "AP4", "AP3", "AP6", "DH5", "AP1", "DH2", "AP5", "DH3", "DH4", "DH1", "DH7", "DH10", "DH11", "DH9", "DH6", "DH8")
  mydf$variable <- factor(mydf$variable,levels = dwnwnd)
  brks_scale <- levels(mydf$CRPS_quantiles)
  labels_scale <- rev(brks_scale)
  
  # Plot
  plot.size = 8; line.size = 0.1; point.size = 0.6
  p <- ggplot(data = mydf, aes(x = variable, y = as.factor(Horizon))) +
    geom_raster(aes(fill=CRPS_quantiles)) +
    xlab("Station") +
    ylab("Forecast horizon (min)") +
    theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
          text = element_text(family = "Times"),
          axis.text=element_text(size=plot.size),
          axis.text.x = element_text(angle = 90),
          axis.title=element_text(size=plot.size),
          legend.position = "bottom",
          legend.justification = "center",
          legend.title = element_text(size = plot.size),
          legend.text = element_text(size = plot.size),
          legend.margin = ggplot2::margin(0,0,0,0),
          legend.box.margin = ggplot2::margin(-10,-10,0,-10)) +
    scale_fill_manual(
      values = rev(magma(10)),
      breaks = rev(brks_scale),
      name = "CRPS",
      drop = FALSE,
      labels = labels_scale,
      guide = guide_legend(
        direction = "horizontal",
        keyheight = unit(2, units = "mm"),
        keywidth = unit(70 / length(labels), units = "mm"),
        title.position = 'top',
        title.hjust = 0.5,
        label.hjust = 1,
        nrow = 1,
        byrow = T,
        reverse = T,
        label.position = "bottom"
      )
    )
  
  return(list(mydf,p,resolution_melt,reliability_melt))
}
#################################################################################
# Function to calculate daily CRPS skill using Dazhi Yang's CH-PeEn. Code is in
# part adapted from "A universal benchmarking method for probabilistic solar 
# irradiance forecasting"
#################################################################################
DailySkill <- function(Data1min, station.names, Results){
  Days <- seq.Date(from = as.Date("2010-07-01"), to = as.Date("2010-08-09"), by = "day") # Select 40 days 
  train <- Data1min[!(as.Date(Data1min$Time, format = "%Y-%m-%d") %in% Days),] # Select all other days from data
  time_fit <- paste0(lubridate::hour(train$Time), ":", lubridate::minute(train$Time))
  K <- 9 
  skillscore <- list()
  
  cl <- makeCluster(5)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(Days), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  skillscore <- foreach(d = 1:length(Days), .packages = c("lubridate","scoringRules"), .options.snow = opts) %dopar% {
    X <- Data1min[as.Date(Data1min$Time, format = "%Y-%m-%d") == Days[d],]
    Tm <- X$Time                          # Save the time
    X <- subset(X, select = -Time)        # Deselect time for the algorithm
    station.names <- colnames(X)
    Daily <- Results[[d]]     # Load the forecast results
    FCs <- Daily[[1]]         # Extract forecasts
    OBS <- Daily[[2]]         # Extract observations
    m <-  350                 # query length
    ns <- ncol(X)             # number of stations
    nt <- 50                  # number of time lags
    start <- nt + 1           
    
    meas <- NULL # Create DF with measurements and time vector from which we can construct CH-PeEn
    for(s0 in 1:ns){ # Loop over the stations 
      Y = X[,s0]  # Take the response variable from the original DF           
      meas <- cbind(meas,Y[(start+m):length(Y)])
      rm(Y)
    }
    meas <- as.data.frame(meas)
    colnames(meas) <- station.names
    meas$Time <- Tm[(start+m):length(Tm)]
    
    timestamps <- unique(paste0(lubridate::hour(meas$Time), ":", lubridate::minute(meas$Time))) # Time stamps in the test data
    time.group <- match(paste0(lubridate::hour(meas$Time), ":", lubridate::minute(meas$Time)), timestamps)
    
    # Calculate CRPS of CH-PeEn (independent of forecast horizon k)
    CRPS_chpeen <- matrix(data = NA, nrow = nrow(meas), ncol = ns) # Store CRPS of CH-PeEn
    for(s0 in 1:ns){
      samples <- NULL 
      for(i in 1:length(timestamps)){
        tmp <- train[which(time_fit==timestamps[i]),s0] # Select column based on station
        tmp <- tmp[which(!is.na(tmp)|!is.nan(tmp))] #remove na's and nan's
        tmp <- ifelse(tmp < 0, NA, tmp)
        tmp <- tmp[complete.cases(tmp)]
        if(length(tmp) != 0)
          samples[[i]] <- tmp
      }
      for(j in (K+1):nrow(meas)){
        CRPS_chpeen[j,s0] <- scoringRules::crps_sample(y = meas[j,s0], dat = samples[[time.group[j]]])
      }
    }
    # Calculate forecast CRPS and skill score. Store the skill score.
    CRPS_ss <- matrix(data = NA, nrow = K, ncol = ns) # Store daily skill score
    for(k in 1:K){ # Loop over the forecast horizons
      CRPS_fcs <- matrix(data = NA, nrow = nrow(meas), ncol = ns)
      ss <- NULL
      for(s0 in 1:ns){ # Loop over the stations
        for(j in (k+1):nrow(meas)){
          CRPS_fcs[j,s0] <- scoringRules::crps_sample(y = meas[j,s0], dat = FCs[[k]][[s0]][j,])
        }
        tmp <- 1 - mean(CRPS_fcs[,s0], na.rm = T) / mean(CRPS_chpeen[,s0], na.rm = T) # Calculate skill score for station s0
        ss <- cbind(ss,tmp)
      }
      CRPS_ss[k,] <- apply(ss,2,mean,na.rm=T) # Calculate average ss. Now it's a daily score per forecast horizon.
    }
    CRPS_ss <- as.data.frame(CRPS_ss)
    colnames(CRPS_ss) <- station.names
    CRPS_ss$Horizon <- seq(1,K,1)
    skillscore[[d]] <- CRPS_ss
  }
  close(pb)
  stopCluster(cl)
  
  ss <- Reduce("+", skillscore) / length(skillscore) # Calculate the average over the days
  mydf <- melt(ss, id.vars="Horizon")
  colnames(mydf)[3] <- "Skill"
  return(mydf)
}
#################################################################################
# Plot sharpness function
#################################################################################
Plot.Sharpness <- function(lstOfResults, station.names,K){
  # Function to plot the sharpness diagrams for all forecast horizons
  # and all locations and put them in a facet of forecast horizons.
  pb <- txtProgressBar(min = 0, max = length(lstOfResults)*9*17, style = 3)
  j <- 0
  tmp <- list()
  res <- NULL
  ns <- length(station.names)
  for(h in 1:K){ # Loop over the forecast horizons
    for(l in 1:ns){  # Loop over the stations
      for(d in 1:length(lstOfResults)){ # Loop over the days 
        Daily <- lstOfResults[[d]]
        FCs <- Daily[[1]]
        OBS <- Daily[[2]]
        tmp[[d]] <- PIs.tmp(ens = FCs[[h]][[l]], obs = OBS[[1]][[l]])[[1]] # Extract PINAW DF
        j <- j+1
        setTxtProgressBar(pb, j) # Update the progress bar
      }
      PINAW <- Reduce("+", tmp) / length(tmp) # Calculate the average over the days for 1 FH and 1 location (double check and correct)
      PINAW$Location <- station.names[l]
      PINAW$Horizon <- h
      res <- rbind(res,PINAW)
    }
  }
  plot.size = 8; line.size = 0.1; point.size = 0.6
  p <- ggplot() +
    geom_line(data = res, size = line.size, aes(x=eta,y=PINAW,colour=Location)) +
    facet_wrap(~Horizon,ncol = 3) +
    xlab("Nominal coverage level") +
    ylab("PINAW") +
    theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
          axis.text=element_text(family = "Times", size=plot.size),
          axis.title=element_text(family = "Times", size=plot.size),
          axis.text.x = element_text(angle = 90),
          legend.title = element_text(family = "Times", size=plot.size),
          legend.text = element_text(family = "Times", size=plot.size),
          legend.key.height = unit(0.75,"lines"),
          legend.margin = ggplot2::margin(0,0,0,0),
          legend.box.margin = ggplot2::margin(-5,-4,0,-10),
          strip.text = element_text(family = "Times", size=plot.size),
          panel.spacing = unit(0.05, "lines")) +
    background_grid(major = "xy")
  return(p)
}
#################################################################################
# Reliability diagram including consistency bars
#################################################################################
rel.cbars <- function(ens,obs,title=NULL) { 
  
  # Some sanity checks
  stopifnot(nrow(ens) == length(obs))
  
  if(!is.null(title)) {title <- title} # Optional: add a title to the plot
  
  idx <- complete.cases(cbind(ens,obs))     # Find rows that have missing values
  ens <- ens[idx,]                          # and remove rows that contain these
  obs <- obs[idx]                           # also in the observations.
  ens <- t(apply(ens, 1, sort))             # Make sure there's no quantile crossing (requires sum(is.na)=0)
  
  taus <- seq(0,1,(1/(ncol(ens)+1)))
  taus <- taus[2:(length(taus)-1)]
  res  <- as.data.frame(matrix(0, nrow = ncol(ens), ncol = 4))
  colnames(res) <- c("Reliability","Down","Up","Probability")
  res$Probability <- taus
  
  # Reliability:
  M_obs <- replicate(n = ncol(ens), expr = obs, simplify = "matrix") # obs vector to obs matrix.
  res$Reliability <- colSums(M_obs < ens)/nrow(ens)
  
  # Consistency bars:
  res$Up <- qbinom(p = rep(max(taus),ncol(ens)), size = replicate(n = ncol(ens), expr = nrow(ens), simplify = "vector"), prob = taus) / nrow(ens)
  res$Down <- qbinom(p = rep(min(taus),ncol(ens)), size = replicate(n = ncol(ens), expr = nrow(ens), simplify = "vector"), prob = taus) / nrow(ens)
  
  p <- ggplot() +
    geom_line(data = res, size = 1, aes(x=Probability,y=Reliability)) +
    geom_line(data = res, size = 1, aes(x=Probability,y=Down), linetype=2) +
    geom_line(data = res, size = 1, aes(x=Probability,y=Up),linetype=2) +
    geom_abline(slope = 1, intercept = 0, linetype=3, size = 1) +
    xlab("Nominal probability") +
    ylab("Observed frequency") +
    theme(axis.text=element_text(size=24),
          axis.title=element_text(size=24),
          plot.title=element_text(size=24)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  list(res,p)
}
#################################################################################
# Plot reliability function
#################################################################################
Plot.ReliabilityDiagram <- function(lstOfResults, station.names, K){
  # Function to plot the reliability diagrams for all forecast horizons
  # and all locations and put them in a facet of forecast horizons.
  pb <- txtProgressBar(min = 0, max = length(lstOfResults)*9*17, style = 3)
  j <- 0
  tmp <- list()
  res <- NULL
  ns <- length(station.names)
  for(h in 1:K){ # Loop over the forecast horizons
    for(l in 1:ns){  # Loop over the stations
      # time_instances <- NULL
      for(d in 1:length(lstOfResults)){ # Loop over the days 
        Daily <- lstOfResults[[d]]
        FCs <- Daily[[1]]
        OBS <- Daily[[2]]
        tmp[[d]] <- rel.cbars(ens = FCs[[h]][[l]], obs = OBS[[1]][[l]])[[1]] 
        # time_instances <- rbind(time_instances,nrow(FCs[[h]][[l]]))
        j <- j+1
        setTxtProgressBar(pb, j) # Update the progress bar
      }
      rMAE <- Reduce("+", tmp) / length(tmp) # Calculate the average over the days for 1 FH and 1 location (double check and correct)
      rMAE$Location <- station.names[l]
      rMAE$Horizon <- h
      res <- rbind(res,rMAE)
    }
  }
  # # Consistency bars (too narrow for this study):
  # taus <- unique(res$Probability)
  # 
  # cons_bars <- data.frame(Up=qbinom(p = rep(max(taus),length(taus)), size = replicate(n = length(taus), expr = sum(time_instances), simplify = "vector"), prob = taus) / sum(time_instances),
  #                         Down=qbinom(p = rep(min(taus),length(taus)), size = replicate(n = length(taus), expr = sum(time_instances), simplify = "vector"), prob = taus) / sum(time_instances))

  plot.size = 8; line.size = 0.1; point.size = 0.6
  p <- ggplot() +
    geom_line(data = res, size = line.size, aes(x=Probability,y=Reliability,colour=Location)) +
    # geom_line(data = cons_bars, aes(x=Probability,y=Up)) +
    # geom_line(data = cons_bars, aes(x=Probability,y=Down)) +
    geom_abline(slope = 1, intercept = 0, linetype=3, size = line.size) +
    facet_wrap(~Horizon,ncol = 3) +
    xlab("Nominal probability") +
    ylab("Observed frequency") +
    theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
          axis.text=element_text(family = "Times", size=plot.size),
          axis.text.x = element_text(angle = 90),
          axis.title=element_text(family = "Times", size=plot.size),
          legend.title = element_text(family = "Times", size=plot.size),
          legend.key.height = unit(0.75,"lines"),
          legend.margin = ggplot2::margin(0,0,0,0),
          legend.box.margin = ggplot2::margin(-5,-4,0,-10),
          legend.text = element_text(family = "Times", size=plot.size),
          strip.text = element_text(family = "Times", size = plot.size),
          panel.spacing = unit(0.05, "lines")) +
    background_grid(major = "xy")
  return(p)
}
#################################################################################
# Convert first letter of the string to upper case
#################################################################################
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
#################################################################################
# Function to fit the copulas, after which scenarios can be generated
#################################################################################
cop_fit <- function(S, K, Results, Data1min, copula, Days){
  res.scen <- list()
  pb <- txtProgressBar(min = 0, max = length(Results), style = 3)
  j <- 0
  runTimes <- NULL
  for(d in 1:length(Results)){
    FCs <- Results[[d]][[1]] # First day, then forecasts
    X <- Data1min[as.Date(Data1min$Time, format = "%Y-%m-%d") == Days[d],] # For the correlation structure
    X <- subset(X, select = -Time)        # Deselect time for the algorithm
    
    ##### To test various copulas:
    X.lagged <- NULL # Create a matrix with lagged observations of all stations
    for(l in 1:K){ # Loop over the forecast horizons
      tmp <- do.call(cbind,shift(X, n = l, type = "lag"))
      X.lagged <- cbind(X.lagged,tmp)
    }
    X.lagged <- X.lagged[1:400,] # Only the training period, not the test period (400 relates to m + nt)
    
    set.seed(123)
    ptm <- proc.time() 
    if(copula == "normal"){
      cop_model <- copula::normalCopula(dim = ncol(X.lagged), dispstr = "toep")
      cop_model@parameters <- copula::fitCopula(copula = cop_model, data = X.lagged, method = "itau")@estimate
    } else if (copula == "empirical") {
      m <- copula::pobs(X.lagged); m <- m[complete.cases(m),]
      cop_model <- copula::empCopula(m)
    } else if (copula == "clayton") {
      cop_model <- copula::claytonCopula(dim = ncol(X.lagged))
      cop_model@parameters <- copula::fitCopula(cop_model, X.lagged, method = 'itau')@estimate
    } else if (copula == "vine"){ # This is extremely slow
      cop_model <- RVineStructureSelect(data = X.lagged)
      scenarios <- RVineSim(S, cop_model)
    } else if(copula == "t"){
      cop_model <- copula::tCopula(dim = ncol(X.lagged), dispstr = "toep", df.fixed = T, df = 3)
      cop_model@parameters <- copula::fitCopula(copula = cop_model, data = X.lagged, method = "itau")@estimate
    } else if(copula == "high"){ # Correlation 1, no actual copula
      cop_model = "high"
    } else if(copula == "low"){ # Correlation 0, no actual copula
      cop_model = "low"
    } else if(copula == "biased"){ # High bias, no actual copula
      cop_model = "biased"
    }
    runTimes <- rbind(runTimes,round((proc.time() - ptm)[3],digits=3))
    ##### This ends the testing of various copulas, now 
    res.scen[[d]] <- scen.gen(fcst.lst = FCs, # Note that this is for 1 day, that's why the loop
                              S = S,
                              K = K,
                              ns = ncol(X),
                              cop_model = cop_model) 
    j <- j+1
    setTxtProgressBar(pb, j) # Update the progress bar
  }
  return(list(res.scen,runTimes))
}
