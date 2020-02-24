# This script plots the rank histograms that are stored as .RData files, first
# to study the different misspecifications in correlation and then my results.

library(ggplot2)
library(cowplot)

setwd("~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/CorrelationStudy")
RHs <- list.files(path = ".")
for(r in RHs){
  mystr <- strsplit(list.files(path = ".", pattern = r),"[_.]")[[1]]
  if(mystr[1] == "ARH"){
    load(r) # load the RH
    dat <- data.frame(Counts=ARH$counts, x=ARH$mids)
    p <- ggplot(data = dat, aes(x=x,y=Counts)) +
      geom_bar(stat = "identity") +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())
    ggsave(filename = paste(mystr[1],"_",mystr[2],".pdf",sep = ""), plot = p, 
           device = "pdf", units = "cm", height = 21, width = 29.7) # Standard: height=21,width=29.7
  } else if(mystr[1] == "BDR") {
    load(r) # load the RH
    dat <- data.frame(Counts=BDR$counts, x=BDR$mids)
    p <- ggplot(data = dat, aes(x=x,y=Counts)) +
      geom_bar(stat = "identity") +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())
    ggsave(filename = paste(mystr[1],"_",mystr[2],".pdf",sep = ""), plot = p, 
           device = "pdf", units = "cm", height = 21, width = 29.7) 
  }
}

####################################################################################################

setwd("~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/RankHistograms")
setwd("C:/Users/SharedComputer/Documents/dennis/UltraFastPreselection/RankHistograms")
RHs <- list.files(path = ".")
for(r in RHs){
  mystr <- strsplit(list.files(path = ".", pattern = r),"[_.]")[[1]]
  if(mystr[1] == "ARH" & mystr[3] != "NoPreselection"){
    load(r)
    ARH <- hist(ARH, plot = F)
    dat <- data.frame(Counts=ARH$counts, x=ARH$mids)
    p <- ggplot(data = dat, aes(x=x,y=Counts)) +
          geom_bar(stat = "identity") +
          theme(axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.line = element_blank())
    ggsave(filename = paste("Preselection",mystr[1],"_",mystr[2],".pdf",sep = ""), plot = p, 
           device = "pdf", units = "cm", height = 21, width = 29.7)
  } else if(mystr[1] == "ARH" & mystr[3] == "NoPreselection") {
    load(r)
    ARH <- hist(ARH, plot = F)
    dat <- data.frame(Counts=ARH$counts, x=ARH$mids)
    p <- ggplot(data = dat, aes(x=x,y=Counts)) +
      geom_bar(stat = "identity") +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())
    ggsave(filename = paste("NoPreselection",mystr[1],"_",mystr[2],".pdf",sep = ""), plot = p, 
           device = "pdf", units = "cm", height = 21, width = 29.7)
  } else if(mystr[1] == "BDR" & mystr[3] != "NoPreselection"){
    load(r) # load the RH
    BDR <- hist(BDR, plot = F)
    dat <- data.frame(Counts=BDR$counts, x=BDR$mids)
    p <- ggplot(data = dat, aes(x=x,y=Counts)) +
      geom_bar(stat = "identity") +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())
    ggsave(filename = paste("Preselection",mystr[1],"_",mystr[2],".pdf",sep = ""), plot = p, 
           device = "pdf", units = "cm", height = 21, width = 29.7)
  } else if(mystr[1] == "BDR" & mystr[3] == "NoPreselection") {
    load(r) # load the RH
    BDR <- hist(BDR, plot = F)
    dat <- data.frame(Counts=BDR$counts, x=BDR$mids)
    p <- ggplot(data = dat, aes(x=x,y=Counts)) +
      geom_bar(stat = "identity") +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())
    ggsave(filename = paste("NoPreselection",mystr[1],"_",mystr[2],".pdf",sep = ""), plot = p, 
           device = "pdf", units = "cm", height = 21, width = 29.7)
  }
}
