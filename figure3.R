# This script is used to generate figure 3 in the manuscript

# Author: Dennis van der Meer
# Email: denniswillemvandermeer@gmail.com

setwd("~/Desktop/Drive/PhD-Thesis/My-papers/UltraFastPreselection")
libs <- c("ggplot2","cowplot","copula","Cairo")
invisible(lapply(libs, library, character.only = TRUE))

load("~/Desktop/Hawaii_1min.RData") # Load the Hawaii data
dat <- as.matrix(Data1min[strftime(Data1min$Time, format = "%Y-%m-%d") %in% "2010-07-01",1:3])
set.seed(123)
# Normal copula
cop_model <- copula::normalCopula(dim = ncol(dat), dispstr = "toep")
cop_model@parameters <- copula::fitCopula(copula = cop_model, data = dat, method = "itau")@estimate
u_norm <- copula::rCopula(nrow(dat), cop_model)

# Empirical copula
m <- copula::pobs(dat); m <- m[complete.cases(m),]
cop_model <- copula::empCopula(m)
u_emp <- copula::rCopula(nrow(dat), cop_model)

# Clayton copula
cop_model <- copula::claytonCopula(dim = ncol(dat))
cop_model@parameters <- copula::fitCopula(cop_model, dat, method = 'itau')@estimate
u_cla <- copula::rCopula(nrow(dat), cop_model)

# Student t copula
cop_model <- copula::tCopula(dim = ncol(dat), dispstr = "toep", df.fixed = T, df = 3)
cop_model@parameters <- copula::fitCopula(copula = cop_model, data = dat, method = "itau")@estimate
u_stu <- copula::rCopula(nrow(dat), cop_model)

u_norm <- data.frame(x=u_norm[,1],y=u_norm[,2],Copula="Normal")
u_emp <- data.frame(x=u_emp[,1],y=u_emp[,2],Copula="Empirical")
u_cla <- data.frame(x=u_cla[,1],y=u_cla[,2],Copula="Clayton")
u_stu <- data.frame(x=u_stu[,1],y=u_stu[,2],Copula="Student-t")

rnd <- rbind(u_norm,u_stu,u_cla,u_emp)

plot.size = 8; line.size = 0.1; point.size = 0.4

p <- ggplot(data = rnd, aes(x=x,y=y)) +
      geom_point(size=point.size, alpha=0.25) +
      geom_density_2d(size = line.size, color="black") +
      scale_x_continuous(breaks = seq(0,1,0.25), limits = c(-0.1,1.1)) +
      scale_y_continuous(breaks = seq(0,1,0.25), limits = c(-0.1,1.1)) +
      facet_wrap(.~Copula) +
      theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
            axis.text=element_text(family = "Times", size=plot.size),
            axis.title=element_text(family = "Times", size=plot.size),
            plot.title=element_text(family = "Times", size=plot.size),
            panel.spacing = unit(0.05, "lines"),
            strip.text = element_text(family = "Times", size=plot.size,
                                      margin = margin(0.02,0,0.02,0, "lines"))) +
            background_grid(major = "xy")
print(p)
ggsave(filename = "BivariateCopula.pdf", plot = p,
       device = cairo_pdf, units = "cm", height = 8.5, width = 8.5) 
