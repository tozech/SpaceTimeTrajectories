# This script postprocesses the results

# Author: Dennis van der Meer
# Email: denniswillemvandermeer@gmail.com

setwd("~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection")
libs <- c('xtable',"ggplot2","cowplot","viridis","data.table","dplyr",
          "reshape2","tidyverse","verification")
invisible(lapply(libs, library, character.only = TRUE))
load("~/Desktop/Data/Hawaii/Hawaii_1min.RData")
load("/Volumes/G-DRIVE_USB-C/Data/Hawaii/Hawaii_1min.RData")
station.names <- colnames(Data1min[1:17])
source(file.path(getwd(), "git/SpaceTimeTrajectories/functions.R"))
plot.size = 8; line.size = 0.1; point.size = 0.6
Days <- seq.Date(from = as.Date("2010-07-01"), to = as.Date("2010-08-09"), by = "day")
###################################################################################
# Generate table with CRPS, skill, reliability and resolution as a function 
# of horizon. Also plot reliability diagrams and sharpness diagrams.
###################################################################################
load("/Users/Dennis/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/FC_Res_QR_Preselection.RData")
my.CRPS_Preselection <- DailyCRPS(lstOfResults = Results, station.names = station.names)
CRPS_Preselection <- my.CRPS_Preselection[[1]]
Resolution_Preselection <- my.CRPS_Preselection[[3]]
Reliability_Preselection <- my.CRPS_Preselection[[4]]
my.Skill_Preselection <- DailySkill(Data1min = Data1min, station.names = station.names, Results = Results)
p.rel_Preselection <- Plot.ReliabilityDiagram(Results,station.names,K)
p.sharp_Preselection <- Plot.Sharpness(Results,station.names,K)
# Save figures:
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/CRPS_QR_Preselection.pdf", plot = my.CRPS_Preselection[[2]], 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/Reliability_QR_Preselection.pdf", plot = p.rel_Preselection, 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/Sharpness_QR_Preselection.pdf", plot = p.sharp_Preselection, 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 

load("/Users/Dennis/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/FC_Res_QR_NoPreselection.RData")
my.CRPS_NoPreselection <- DailyCRPS(lstOfResults = Results, station.names = station.names)
CRPS_NoPreselection <- my.CRPS_NoPreselection[[1]]
Resolution_NoPreselection <- my.CRPS_NoPreselection[[3]]
Reliability_NoPreselection <- my.CRPS_NoPreselection[[4]]
my.Skill_NoPreselection <- DailySkill(Data1min = Data1min, station.names = station.names, Results = Results)
p.rel_NoPreselection <- Plot.ReliabilityDiagram(Results,station.names,K)
p.sharp_NoPreselection <- Plot.Sharpness(Results,station.names,K)
# Save figures:
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/CRPS_QR_NoPreselection.pdf", plot = my.CRPS_NoPreselection[[2]], 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/Reliability_QR_NoPreselection.pdf", plot = p.rel_NoPreselection, 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/Sharpness_QR_NoPreselection.pdf", plot = p.sharp_NoPreselection, 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 

# Downwind station ordering
dwnwnd <- c("AP7", "AP4", "AP3", "AP6", "DH5", "AP1", "DH2", "AP5", "DH3", "DH4", "DH1", "DH7", "DH10", "DH11", "DH9", "DH6", "DH8")
# Plot CRPS as facets (Figure 8):
# Adapted from https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
tmp_Preselection <- my.CRPS_Preselection[[1]]; tmp_Preselection$Model <- "QR-exo"
tmp_NoPreselection <- my.CRPS_NoPreselection[[1]]; tmp_NoPreselection$Model <- "QR-endo"
tmp_Preselection$variable <- factor(tmp_Preselection$variable,levels = dwnwnd)
tmp_NoPreselection$variable <- factor(tmp_NoPreselection$variable,levels = dwnwnd)
mydf <- rbind(tmp_NoPreselection,tmp_Preselection)
brks_scale <- levels(mydf$CRPS_quantiles)
labels_scale <- rev(brks_scale)
plot.size = 8; line.size = 0.1; point.size = 0.6
p <- ggplot(data = mydf, aes(x = variable, y = as.factor(Horizon))) +
  geom_raster(aes(fill=CRPS_quantiles)) +
  facet_wrap(~Model, nrow = 1) +
  xlab("Station") +
  ylab("Forecast horizon (min)") +
  theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
        text = element_text(family = "Times"),
        axis.text=element_text(size=plot.size),
        axis.text.x = element_text(angle = 90),
        axis.title=element_text(size=plot.size),
        strip.text = element_text(size=plot.size),
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
      keywidth = unit(70 / length(labels_scale), units = "mm"),
      title.position = 'top',
      title.hjust = 0.5,
      label.hjust = 1,
      nrow = 1,
      byrow = T,
      reverse = T,
      label.position = "bottom"
    )
  )
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/CRPS_QR.pdf", plot = p, 
       device = cairo_pdf, units = "cm", height = 8.5, width = 17) 

# Reorder according to downwind positioning:
res_Preselection <- data.frame(Horizon=CRPS_Preselection$Horizon,
                               Station=CRPS_Preselection$variable,
                               CRPS=CRPS_Preselection$CRPS,
                               Resolution=Resolution_Preselection$Resolution,
                               Reliability=Reliability_Preselection$Reliability,
                               Skill=my.Skill_Preselection$Skill)
res_Preselection$Station <- factor(res_Preselection$Station,levels = dwnwnd)
res_NoPreselection <- data.frame(Horizon=CRPS_NoPreselection$Horizon,
                                 Station=CRPS_NoPreselection$variable,
                                 CRPS=CRPS_NoPreselection$CRPS,
                                 Resolution=Resolution_NoPreselection$Resolution,
                                 Reliability=Reliability_NoPreselection$Reliability,
                                 Skill=my.Skill_NoPreselection$Skill)
res_NoPreselection$Station <- factor(res_NoPreselection$Station,levels = dwnwnd)

# No preselection:
res_NoPres <- as.data.frame(matrix(NA, nrow = 9, ncol = 10))
colnames(res_NoPres) <- c("murMAE","sdrMAE","muCRPS","sdCRPS","muResolution","sdResolution",
                          "muReliability","sdReliability","muSkill","sdSkill")
for(i in 1:9){
  res_NoPres[i,1] <- round(mean(res_NoPreselection[res_NoPreselection$Horizon %in% i,"rMAE"]), digits = 4)
  res_NoPres[i,2] <- round(sd((res_NoPreselection[res_NoPreselection$Horizon %in% i,"rMAE"])), digits = 4)
  res_NoPres[i,3] <- round(mean(res_NoPreselection[res_NoPreselection$Horizon %in% i,"CRPS"]), digits = 4)
  res_NoPres[i,4] <- round(sd((res_NoPreselection[res_NoPreselection$Horizon %in% i,"CRPS"])), digits = 4)
  res_NoPres[i,5] <- round(mean(res_NoPreselection[res_NoPreselection$Horizon %in% i,"Resolution"]), digits = 4)
  res_NoPres[i,6] <- round(sd((res_NoPreselection[res_NoPreselection$Horizon %in% i,"Resolution"])), digits = 4)
  res_NoPres[i,7] <- round(mean(res_NoPreselection[res_NoPreselection$Horizon %in% i,"Reliability"]), digits = 4)
  res_NoPres[i,8] <- round(sd((res_NoPreselection[res_NoPreselection$Horizon %in% i,"Reliability"])), digits = 4)
  res_NoPres[i,9] <- round(mean(res_NoPreselection[res_NoPreselection$Horizon %in% i,"Skill"]), digits = 4)
  res_NoPres[i,10] <-round(sd((res_NoPreselection[res_NoPreselection$Horizon %in% i,"Skill"])), digits = 4)
}

results_NoPres <- res_NoPres %>%
  mutate_if(is.numeric, format, nsmall=1) %>%
  mutate(Model = c("QR\textsuperscript{endo}", "", "", "", "", "", "", "", "")) %>%
  mutate(Horizon = as.integer(seq(1,9,1))) %>%
  mutate(CRPS = paste(muCRPS, "$\\pm$", sdCRPS)) %>%
  mutate(Reliability = paste(muReliability, "$\\pm$", sdReliability)) %>%
  mutate(Resolution = paste(muResolution, "$\\pm$",sdResolution)) %>%
  mutate(Skill = paste(muSkill, "$\\pm$", sdSkill)) %>%
  dplyr::select(-one_of(c("murMAE", "sdrMAE", "muCRPS", "sdCRPS", "muSkill", "sdSkill","muResolution",
                          "sdResolution", "muReliability","sdReliability")))

# Preselection:
res_Pres <- as.data.frame(matrix(NA, nrow = 9, ncol = 10))
colnames(res_Pres) <- c("murMAE","sdrMAE","muCRPS","sdCRPS","muResolution","sdResolution",
                          "muReliability","sdReliability","muSkill","sdSkill")
for(i in 1:9){
  res_Pres[i,1] <- round(mean(res_Preselection[res_Preselection$Horizon %in% i,"rMAE"]), digits = 4)
  res_Pres[i,2] <- round(sd((res_Preselection[res_Preselection$Horizon %in% i,"rMAE"])), digits = 4)
  res_Pres[i,3] <- round(mean(res_Preselection[res_Preselection$Horizon %in% i,"CRPS"]), digits = 4)
  res_Pres[i,4] <- round(sd((res_Preselection[res_Preselection$Horizon %in% i,"CRPS"])), digits = 4)
  res_Pres[i,5] <- round(mean(res_Preselection[res_Preselection$Horizon %in% i,"Resolution"]), digits = 4)
  res_Pres[i,6] <- round(sd((res_Preselection[res_Preselection$Horizon %in% i,"Resolution"])), digits = 4)
  res_Pres[i,7] <- round(mean(res_Preselection[res_Preselection$Horizon %in% i,"Reliability"]), digits = 4)
  res_Pres[i,8] <- round(sd((res_Preselection[res_Preselection$Horizon %in% i,"Reliability"])), digits = 4)
  res_Pres[i,9] <- round(mean(res_Preselection[res_Preselection$Horizon %in% i,"Skill"]), digits = 4)
  res_Pres[i,10] <-round(sd((res_Preselection[res_Preselection$Horizon %in% i,"Skill"])), digits = 4)
}

results_Pres <- res_Pres %>%
  mutate_if(is.numeric, format, nsmall=1) %>%
  mutate(Model = c("QR\textsuperscript{exo}", "", "", "", "", "", "", "", "")) %>%
  mutate(Horizon = as.integer(seq(1,9,1))) %>%
  mutate(CRPS = paste(muCRPS, "$\\pm$", sdCRPS)) %>%
  mutate(Reliability = paste(muReliability, "$\\pm$", sdReliability)) %>%
  mutate(Resolution = paste(muResolution, "$\\pm$",sdResolution)) %>%
  mutate(Skill = paste(muSkill, "$\\pm$", sdSkill)) %>%
  dplyr::select(-one_of(c("murMAE", "sdrMAE", "muCRPS", "sdCRPS", "muSkill", "sdSkill","muResolution",
                          "sdResolution", "muReliability","sdReliability")))

result_table <- rbind(results_NoPres,results_Pres)
print(xtable(result_table), include.rownames=FALSE, sanitize.text.function = function(x){x})

###################################################################################
################# Generate table and plots with Energy score  #####################
###################################################################################
ES_NoPres <- read.table(file = "energyScore_NoPreselection.txt", header = T)
ES_Pres <- read.table(file = "energyScore_Preselection.txt", header = T)
colnames(ES_NoPres) <- c("Empirical","Gaussian","Clayton","Student-t")
colnames(ES_Pres) <- c("Empirical","Gaussian","Clayton","Student-t")
# Reorder:
ES_NoPres <- ES_NoPres[,c(2,4,3,1)]
ES_Pres <- ES_Pres[,c(2,4,3,1)]
# Generate tables:
tab_ES_NoPres <- data.frame(Model = "QR\textsuperscript{endo}",
                          Gaussian = paste(round(mean(ES_NoPres$Gaussian),digits = 3), "$\\pm$", round(sd(ES_NoPres$Gaussian), digits = 3)),
                          `Student-t` = paste(round(mean(ES_NoPres$`Student-t`),digits = 3), "$\\pm$", round(sd(ES_NoPres$`Student-t`), digits = 3)),
                          Clayton = paste(round(mean(ES_NoPres$Clayton),digits = 3), "$\\pm$", round(sd(ES_NoPres$Clayton), digits = 3)),
                          Empirical = paste(round(mean(ES_NoPres$Empirical),digits = 3), "$\\pm$", round(sd(ES_NoPres$Empirical), digits = 3))
                          )
tab_ES_Pres <- data.frame(Model = "QR\textsuperscript{exo}",
                          Gaussian = paste(round(mean(ES_Pres$Gaussian),digits = 3), "$\\pm$", round(sd(ES_Pres$Gaussian), digits = 3)),
                          `Student-t` = paste(round(mean(ES_Pres$`Student-t`),digits = 3), "$\\pm$", round(sd(ES_Pres$`Student-t`), digits = 3)),
                          Clayton = paste(round(mean(ES_Pres$Clayton),digits = 3), "$\\pm$", round(sd(ES_Pres$Clayton), digits = 3)),
                          Empirical = paste(round(mean(ES_Pres$Empirical),digits = 3), "$\\pm$", round(sd(ES_Pres$Empirical), digits = 3))
                          )
tab_ES <- rbind(tab_ES_NoPres,tab_ES_Pres)
print(xtable(tab_ES), include.rownames=FALSE, sanitize.text.function = function(x){x})

# Generate energy score plot:
ES_NoPres$Model <- "N"
ES_Pres$Model <- "Y"
ES_NoPres$Date <- Days
ES_Pres$Date <- Days
mydf_NoP <- melt(ES_NoPres, id.vars = c("Model","Date"))
mydf_P <- melt(ES_Pres, id.vars = c("Model","Date"))
mydf <- rbind(mydf_NoP,mydf_P)
lab1 <- c(expression(QR^{endo}), expression(QR^{exo}))
group.colors <- c(N = "red", Y = "black")

p.ES <- ggplot(data = mydf, aes(x=Date,y=value,colour=Model)) +
  facet_wrap(.~variable, nrow = 2) +
  geom_line(size=line.size,alpha=0.75) +
  theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
        text = element_text(family = "Times"),
        axis.text=element_text(family = "Times", size=plot.size),
        axis.text.x = element_text(angle = 90),
        axis.title=element_text(family = "Times", size=plot.size),
        plot.title=element_text(family = "Times", size=plot.size),
        legend.title = element_text(family = "Times", size=plot.size),
        legend.text = element_text(family = "Times", size=plot.size),
        legend.box.margin = ggplot2::margin(-5,-4,0,-10),
        legend.margin = ggplot2::margin(0,0,0,0),
        strip.text = element_text(family = "Times", size=plot.size),
        panel.spacing = unit(0.05, "lines")) +
  background_grid(major = "xy") +
  scale_color_manual(values = group.colors, labels = lab1) +
  ylab("Energy score (-)")
ggsave(filename = "~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/Paper/images/energyScorePlot.pdf", plot = p.ES, 
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5)

###################################################################################
# Generate multivariate rank histograms in Fig. 9 and Fig. 10
###################################################################################
# This assumes the multivariate rank histograms are stored as .RData (resulting
# from the final loop in main.R)
setwd("~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection/RankHistograms")
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

