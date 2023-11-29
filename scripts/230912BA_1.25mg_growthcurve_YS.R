####################################
#growth curve function
#written by Rich Fitzjohn  in Sally Otto's lab
####################################
nderiv <- function(fit, x, eps=1e-5)
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)

spline.slope <- function(x, y, n=101, eps=1e-5)
  max(nderiv(loess(log(y) ~ x, degree=1, span=0.1),
             seq(min(x), max(x), length=n)), na.rm=TRUE)

growthrate <- function(d, ...)
  sapply(d[-1], spline.slope, x=d$t, ...)

library(here)

BA125GC_data <- read.csv(here("data_in", "GCs", "230111_125mg","230111_BA125_GC.csv"), sep = ",", header = F)

#################Naming the columns##########################
colnames(BA125GC_data) <- c("t", paste("C", seq(length=ncol(BA125GC_data)-1), sep=""))
colnames(BA125GC_data)[2:13] <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12")
colnames(BA125GC_data)[14:25] <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12")
colnames(BA125GC_data)[26:37] <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")
colnames(BA125GC_data)[38:49] <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12")
colnames(BA125GC_data)[50:61] <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12")
colnames(BA125GC_data)[62:73] <- c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")
colnames(BA125GC_data)[74:85] <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12")
colnames(BA125GC_data)[86:97] <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
##########Calculating the growth rate ############
BA125GClgr <- growthrate(BA125GC_data)

BA125GC_OD24 <- BA125GC_data[BA125GC_data$t==24, 2:length(BA125GC_data)]
BA125GC_OD24fixed <- t(BA125GC_OD24)
BA125GC_OD48 <- BA125GC_data[BA125GC_data$t==48, 2:length(BA125GC_data)]
BA125GC_OD48fixed <- t(BA125GC_OD48)
BA125GC_OD72 <- BA125GC_data[BA125GC_data$t==72, 2:length(BA125GC_data)]
BA125GC_OD72fixed <- t(BA125GC_OD72)

###################################################
BA125GCmaxOD <- apply(BA125GC_data[,2:length(BA125GC_data)], 2, function(x) max(x))

#######################################################
BA125GCfakelag <- apply(BA125GC_data[,2:length(BA125GC_data)], 2, function(x) BA125GC_data$t[min(which(x > 0.3))])
################################################################################
BA125GCwells <- read.csv(here("data_in", "GCs", "230111_125mg","230111_BA125_ANCEVOdata.csv"))

################################################################################3
BA125GCdf <- data.frame(BA125GCwells, BA125GClgr, BA125GC_OD24fixed, BA125GC_OD48fixed, BA125GC_OD72fixed, BA125GCmaxOD, BA125GCfakelag)
colnames(BA125GCdf)[8:10] <- c("OD24", "OD48", "OD72")
##############################################################################
library(dplyr)
BA125GCdf <- BA125GCdf %>%
  mutate(Strain_type = paste(Strain, Type, sep = "_")) %>%
  relocate(., Strain_type, .after = Environment)
BA125GCdf <- subset(BA125GCdf, select = -c(4,6))

#remove A18 T0 R7
BA125GCdf <- BA125GCdf[-c(28,60,92),]

BA125GC_fulldata <- BA125GCdf %>%
  group_by(Strain_type) %>%
  summarise(OD24_avg = mean(OD24, na.rm = T),
            OD48_avg = mean(OD48, na.rm = T),
            OD72_avg = mean(OD72, na.rm = T),
            OD24_median = median(OD24, na.rm = T),
            OD48_median = median(OD48, na.rm = T),
            OD72_median = median(OD72, na.rm = T),
            OD24_max = max(OD24, na.rm = T),
            OD48_max = max(OD48, na.rm = T),
            OD72_max = max(OD72, na.rm = T),
            OD24_min = min(OD24, na.rm = T),
            OD48_min = min(OD48, na.rm = T),
            OD72_min = min(OD72, na.rm = T),
            lgr_avg = mean(BA125GClgr, na.rm = T),
            fakelag_avg = mean(BA125GCfakelag, na.rm = T),
            maxOD_avg = mean(BA125GCmaxOD, na.rm = T),
            lgr_median = median(BA125GClgr, na.rm = T),
            fakelag_median = median(BA125GCfakelag, na.rm = T),
            maxOD_median = median(BA125GCmaxOD, na.rm = T),
            lgr_max = max(BA125GClgr, na.rm = T),
            fakelag_max = max(BA125GCfakelag, na.rm = T),
            maxOD_max = max(BA125GCmaxOD, na.rm = T),
            lgr_min = min(BA125GClgr, na.rm = T),
            fakelag_min = min(BA125GCfakelag, na.rm = T))

BA125GC_fulldata1 <- BA125GCdf %>%
  group_by(Strain_type, Replicate) %>%
  summarise(OD24_avg = mean(OD24, na.rm = T),
            OD48_avg = mean(OD48, na.rm = T),
            OD72_avg = mean(OD72, na.rm = T),
            OD24_median = median(OD24, na.rm = T),
            OD48_median = median(OD48, na.rm = T),
            OD72_median = median(OD72, na.rm = T),
            OD24_max = max(OD24, na.rm = T),
            OD48_max = max(OD48, na.rm = T),
            OD72_max = max(OD72, na.rm = T),
            OD24_min = min(OD24, na.rm = T),
            OD48_min = min(OD48, na.rm = T),
            OD72_min = min(OD72, na.rm = T),
            lgr_avg = mean(BA125GClgr, na.rm = T),
            fakelag_avg = mean(BA125GCfakelag, na.rm = T),
            maxOD_avg = mean(BA125GCmaxOD, na.rm = T),
            lgr_median = median(BA125GClgr, na.rm = T),
            fakelag_median = median(BA125GCfakelag, na.rm = T),
            maxOD_median = median(BA125GCmaxOD, na.rm = T),
            lgr_max = max(BA125GClgr, na.rm = T),
            fakelag_max = max(BA125GCfakelag, na.rm = T),
            maxOD_max = max(BA125GCmaxOD, na.rm = T),
            lgr_min = min(BA125GClgr, na.rm = T),
            fakelag_min = min(BA125GCfakelag, na.rm = T),
            BA125GClgr_avg = mean(BA125GClgr, na.rm = T))

BA125GCfulldata_datafixed <- BA125GC_fulldata[-c(1, 3, 5, 7), ]

library(ggforce)
library(cowplot)
library(here)
library(ggpubr)

#duplicate column "Strain_type"
BA125GC_fulldata1 <- BA125GC_fulldata1 %>%           
  mutate(Strain_type2 = Strain_type)
#separate 
library(tidyr)
BA125GC_fulldata1 <- separate(BA125GC_fulldata1, Strain_type2, into = c('Strain', 'Origin'), sep='_')

##########################plot lgr - growth rate##########################
P4ave <- ggplot(data = BA125GC_fulldata1,
                mapping = aes(x = Strain,
                              y = BA125GClgr_avg,
                              colour = Origin))+
  scale_y_continuous(expand = c(0, 0), limits= c(-0.01, 0.51))+
  geom_sina(data = BA125GC_fulldata1, mapping = aes(x = Strain_type, y = BA125GClgr_avg),
            maxwidth = 0.50,
            alpha = 0.6,
            size=4) +
  geom_point(data = BA125GC_fulldata, mapping = aes(x = Strain_type, y = lgr_median),
             color="black",
             pch="_",
             alpha=1,
             lwd = 10,
             position = position_dodge(width = 0.3)) +
  labs(y = "Growth rate (/h)", x = "") +
  scale_color_manual(values = c("grey50",
                                "purple",
                                "grey50",
                                "purple",
                                "grey50",
                                "purple",
                                "grey50",
                                "purple"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.background = element_blank(), legend.box.background = element_blank(), 
        legend.key=element_blank(),
        panel.background = element_rect(
          fill = "white",
          colour = "black",
          linewidth = 0.4),
        text=element_text(family="Times New Roman", face="bold", size = 25),
        axis.text = element_text(size = 15, hjust = 0.5)) +
  labs(title = "1.25 mg/mL BA") +
  theme(plot.title = element_text(size=20, hjust = 0.5))+
  scale_x_discrete(limits = c("A18_ancestral", "A18_evolved","A03_ancestral", "A03_evolved","A08_ancestral", "A08_evolved", "A04_ancestral", "A04_evolved"))

P4ave     


ggsave(plot = P4ave, here("figures_out", "GCs", "230111_125mg", "230912_BA125_GC_GrowthRate.jpg"), 
       width = 6.5, 
       height = 4.5)

#################################use for x-axis only#################################
write.csv(BA125GC_fulldata1, here("data_in", "GCs", "230111_125mg","BA125GC_fulldata1.csv")) 
#manually added a column (in excel) with the strain names
BA125GC_fulldata1 <- read.csv(here("data_in", "GCs", "230111_125mg", "BA125GC_fulldata1_edited.csv"))


P4axis <- ggplot(data = BA125GC_fulldata1,
                 mapping = aes(x = StrainName,
                               y = BA125GClgr_avg,
                               colour = Origin))+
  scale_y_continuous(expand = c(0, 0), limits= c(-0.01, 0.51))+
  geom_sina(data = BA125GC_fulldata1, mapping = aes(x = StrainName, y = BA125GClgr_avg),
            maxwidth = 0.50,
            alpha = 0.6,
            size=4) +
  labs(y = "Growth rate (/h)", x = "Strain") +
  scale_color_manual(values = c("grey50",
                                "purple",
                                "grey50",
                                "purple",
                                "grey50",
                                "purple",
                                "grey50",
                                "purple"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.background = element_blank(), legend.box.background = element_blank(), 
        legend.key=element_blank(),
        panel.background = element_rect(
          fill = "white",
          colour = "black",
          linewidth = 0.4)) +
  scale_x_discrete(limits = c("FH1", "GC75","P75016","P78048")) +
  theme(text=element_text(family="Times New Roman", face="bold", size = 25),
        axis.text = element_text(size = 15, hjust = 0.5)) +
  labs(title = "1.25 mg/mL BA") +
  theme(plot.title = element_text(size=20, hjust = 0.5))

P4axis

ggsave(plot = P4axis, here("figures_out", "GCs", "230111_125mg", "230912_BA125_GC_axis.jpg"), 
       width = 6.5, 
       height = 4.5)

#################################t-test#################################
t.test(BA125GClgr_avg~Origin, data = BA125GC_fulldata1 %>% filter(Strain == "A03"))
#t = -3.6643, df = 3.835, p-value = 0.02313
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.066790148 -0.008647318
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.008886032             0.046604765 

t.test(BA125GClgr_avg~Origin, data = BA125GC_fulldata1 %>% filter(Strain == "A04"))
#t = -4.0939, df = 3.0025, p-value = 0.02631
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.08981578 -0.01127214
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.005698449             0.056242407 

t.test(BA125GClgr_avg~Origin, data = BA125GC_fulldata1 %>% filter(Strain == "A08"))
#t = -5.1857, df = 3.3493, p-value = 0.01054
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.08382087 -0.02234640
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.05002657              0.10311020 

t.test(BA125GClgr_avg~Origin, data = BA125GC_fulldata1 %>% filter(Strain == "A18"))
#t = -2.2058, df = 3.0059, p-value = 0.1144
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.06290283  0.01136732
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.003233937             0.029001697