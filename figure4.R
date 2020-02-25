#######################################################################################
############################ Correlation error analysis ###############################
#######################################################################################
setwd("~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection")
load("/Volumes/G-DRIVE_USB-C/Data/Hawaii/Hawaii_1min.RData")

libs <- c("ggplot2","cowplot","copula","data.table","Cairo")
invisible(lapply(libs, library, character.only = TRUE))

Days <- seq.Date(from = as.Date("2010-07-01"), to = as.Date("2010-08-09"), by = "day") # 40 Days
lst <- list()
K <- 9
# source(file.path(getwd(), "/SupplementaryMaterials/functions.R"))
source(file.path(getwd(), "git/SpaceTimeTrajectories/functions.R"))
for(d in seq(1,40,1)){ # Loop over the days
  X <- Data1min[as.Date(Data1min$Time, format = "%Y-%m-%d") == Days[d],] # For the correlation structure
  X <- subset(X, select = -Time)        # Deselect time for the algorithm
  X.lagged <- NULL # Create a matrix with lagged observations of all stations
  for(l in 1:K){ # Loop over the forecast horizons
    tmp <- do.call(cbind,data.table::shift(X, n = l, type = "lag"))
    X.lagged <- cbind(X.lagged,tmp)
  }
  X.lagged <- X.lagged[1:400,]
  lst[[d]] <- X.lagged
}
X.lagged <- do.call(rbind,lst)
X.lagged <- X.lagged[complete.cases(X.lagged),]

# Since X.lagged is too large for the parametric copulas, randomly sample from X.lagged
set.seed(123)
X.lagged <- X.lagged[sample.int(nrow(X.lagged), 500),]

# Loop over copulas
copulas <- c("normal","t","clayton","empirical")
dat <- NULL
for(copula in copulas){
  print(paste("Working on ",copula," copula",sep = ""))
  if(copula == "normal"){
    cop_model <- copula::normalCopula(dim = ncol(X.lagged), dispstr = "toep")
    cop_model@parameters <- copula::fitCopula(copula = cop_model, data = X.lagged, method = "itau")@estimate
  } else if (copula == "empirical") {
    m <- copula::pobs(X.lagged); m <- m[complete.cases(m),]
    cop_model <- copula::empCopula(m)
  } else if (copula == "clayton") {
    cop_model <- copula::claytonCopula(dim = ncol(X.lagged))
    cop_model@parameters <- copula::fitCopula(cop_model, X.lagged, method = 'itau')@estimate
  } else if(copula == "t"){
    cop_model <- copula::tCopula(dim = ncol(X.lagged), dispstr = "toep", df.fixed = T, df = 3)
    cop_model@parameters <- copula::fitCopula(copula = cop_model, data = X.lagged, method = "itau")@estimate
  }
  # Calculate correlations 
  set.seed(123)
  cor_cop <- cor(copula::rCopula(nrow(X.lagged), cop_model))
  cor_dat <- cor(X.lagged)
  # Calculate differences and melt the matrix to plot
  cor_diff <- cor_cop - cor_dat
  melt_diffs <- melt(cor_diff)
  colnames(melt_diffs)[3] <- "Error"
  # write.table(melt_diffs, 
  #             file = paste("/Users/Dennis/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/CorrelationStudy/",copula,sep = ""),
  #             row.names = F, col.names = T)
  # # Plot spatio-temporal correlation matrix
  # p <- ggplot(data = melt_diffs, aes(x=Var1,y=Var2,fill=Error))+
  #   geom_tile() +
  #   ggtitle(paste("Spatio-temporal correlation error (",copula,"), day = all",sep = "")) +
  #   scale_fill_distiller(palette = "Spectral") +
  #   theme(axis.text=element_text(size=24),
  #         axis.title=element_text(size=24),
  #         plot.title=element_text(size=24),
  #         legend.text = element_text(size=16),
  #         legend.title = element_text(size = 18)) 
  # ggsave(filename = paste("CorrelationStudy/correlationError_",copula,".pdf",sep = ""), plot = p, 
  #        device = "pdf", units = "cm", height = 21, width = 29.7) 
  # Plot histogram of errors
  # hi <- hist(melt_diffs$Error,
  #            main = paste("Spatio-temporal correlation error (",copula,"), day = all",sep = ""),
  #            plot = F)
  dat <- rbind(dat,data.frame(Counts=melt_diffs$Error, Copula=firstup(copula)))
}
levels(dat$Copula)[levels(dat$Copula)=="T"] <- "Student-t"
plot.size = 8; line.size = 0.1; point.size = 0.6
p1 <- ggplot(data = dat, aes(x=Counts)) +
  geom_histogram(bins=100) +
  facet_wrap(~Copula, scales = "free_y") +
  theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
        axis.text=element_text(family = "Times", size=plot.size),
        axis.title=element_text(family = "Times", size=plot.size),
        plot.title=element_text(family = "Times", size=plot.size),
        panel.spacing = unit(0.05, "lines"),
        strip.text = element_text(family = "Times", size=plot.size,
                                  margin = margin(0.02,0,0.02,0, "lines"))) +
  xlab("Error (-)") +
  ylab("Count (-)") +
  background_grid(major = "xy")
print(p1)
ggsave(filename = "correlationErrorHistogram.pdf", plot = p1, 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 



