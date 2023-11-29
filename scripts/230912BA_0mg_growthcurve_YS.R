####################################
#growth curve function
#written by Rich Fitzjohn in Sally Otto's lab
####################################
nderiv <- function(fit, x, eps=1e-5)
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)

spline.slope <- function(x, y, n=101, eps=1e-5)
  max(nderiv(loess(log(y) ~ x, degree=1, span=0.1),
             seq(min(x), max(x), length=n)), na.rm=TRUE)

growthrate <- function(d, ...)
  sapply(d[-1], spline.slope, x=d$t, ...)

library(here)

BAGC_data <- read.csv(here("data_in", "GCs", "230510_0mg", "230517_BA_GC.csv"), sep = ",", header = F)

#################Naming the columns##########################
colnames(BAGC_data) <- c("t", paste("C", seq(length=ncol(BAGC_data)-1), sep=""))
colnames(BAGC_data)[2:13] <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12")
colnames(BAGC_data)[14:25] <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12")
colnames(BAGC_data)[26:37] <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")
colnames(BAGC_data)[38:49] <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12")
colnames(BAGC_data)[50:61] <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12")
colnames(BAGC_data)[62:73] <- c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")
colnames(BAGC_data)[74:85] <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12")
colnames(BAGC_data)[86:97] <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12")
##########Calculating the growth rate ############
BAGClgr <- growthrate(BAGC_data)

BAGC_OD24 <- BAGC_data[BAGC_data$t==24, 2:length(BAGC_data)]
BAGC_OD24fixed <- t(BAGC_OD24)
BAGC_OD48 <- BAGC_data[BAGC_data$t==48, 2:length(BAGC_data)]
BAGC_OD48fixed <- t(BAGC_OD48)
BAGC_OD72 <- BAGC_data[BAGC_data$t==72, 2:length(BAGC_data)]
BAGC_OD72fixed <- t(BAGC_OD72)

###################################################
BAGCmaxOD <- apply(BAGC_data[,2:length(BAGC_data)], 2, function(x) max(x))

#######################################################
BAGCfakelag <- apply(BAGC_data[,2:length(BAGC_data)], 2, function(x) BAGC_data$t[min(which(x > 0.3))])
################################################################################
BAGCwells <- read.csv(here("data_in", "GCs", "230510_0mg", "230517_BA_ANCEVOdata.csv"))

################################################################################3
BAGCdf <- data.frame(BAGCwells, BAGClgr, BAGC_OD24fixed, BAGC_OD48fixed, BAGC_OD72fixed, BAGCmaxOD, BAGCfakelag)
colnames(BAGCdf)[8:10] <- c("OD24", "OD48", "OD72")
##############################################################################
library(dplyr)
BAGCdf <- BAGCdf %>%
  mutate(Strain_type = paste(Strain, Type, sep = "_")) %>%
  relocate(., Strain_type, .after = Environment)
BAGCdf <- subset(BAGCdf, select = -c(4,6))

#remove A18 T0 R7
BAGCdf <- BAGCdf[-c(28,60,92),]

BAGC_fulldata <- BAGCdf %>%
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
            lgr_avg = mean(BAGClgr, na.rm = T),
            fakelag_avg = mean(BAGCfakelag, na.rm = T),
            maxOD_avg = mean(BAGCmaxOD, na.rm = T),
            lgr_median = median(BAGClgr, na.rm = T),
            fakelag_median = median(BAGCfakelag, na.rm = T),
            maxOD_median = median(BAGCmaxOD, na.rm = T),
            lgr_max = max(BAGClgr, na.rm = T),
            fakelag_max = max(BAGCfakelag, na.rm = T),
            maxOD_max = max(BAGCmaxOD, na.rm = T),
            lgr_min = min(BAGClgr, na.rm = T),
            fakelag_min = min(BAGCfakelag, na.rm = T))

BAGC_fulldata1 <- BAGCdf %>%
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
            lgr_avg = mean(BAGClgr, na.rm = T),
            fakelag_avg = mean(BAGCfakelag, na.rm = T),
            maxOD_avg = mean(BAGCmaxOD, na.rm = T),
            lgr_median = median(BAGClgr, na.rm = T),
            fakelag_median = median(BAGCfakelag, na.rm = T),
            maxOD_median = median(BAGCmaxOD, na.rm = T),
            lgr_max = max(BAGClgr, na.rm = T),
            fakelag_max = max(BAGCfakelag, na.rm = T),
            maxOD_max = max(BAGCmaxOD, na.rm = T),
            lgr_min = min(BAGClgr, na.rm = T),
            fakelag_min = min(BAGCfakelag, na.rm = T),
            BAGClgr_avg = mean(BAGClgr, na.rm = T))

BAGCfulldata_datafixed <- BAGC_fulldata[-c(1, 3, 5, 7),]

library(ggforce)
library(cowplot)

#duplicate column "Strain_type"
BAGC_fulldata1 <- BAGC_fulldata1 %>%           
  mutate(Strain_type2 = Strain_type)
#separate 
library(tidyr)
BAGC_fulldata1 <- separate(BAGC_fulldata1, Strain_type2, into = c('Strain', 'Origin'), sep='_')

##########################plot lgr - growth rate##########################
P4ave <- ggplot(data = BAGC_fulldata1,
                mapping = aes(x = Strain,
                              y = BAGClgr_avg,
                              colour = Origin))+
  scale_y_continuous(expand = c(0, 0), limits= c(-0.01, 0.51))+
  geom_sina(data = BAGC_fulldata1, mapping = aes(x = Strain_type, y = BAGClgr_avg),
            maxwidth = 0.50,
            alpha = 0.6,
            size=4) +
  geom_point(data = BAGC_fulldata, mapping = aes(x = Strain_type, y = lgr_median),
             colour="black",
             pch="_",
             alpha=1,
             lwd = 10,
             width = 0.5,
             position = position_dodge(width = 0.3)) +
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
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.background = element_blank(), legend.box.background = element_blank(), 
        legend.key=element_blank(),
        panel.background = element_rect(
          fill = "white",
          colour = "black",
          linewidth = 0.4),
        text=element_text(family="Times New Roman", face="bold", size = 25),
        axis.text = element_text(size = 15, hjust = 0.5)) +
  labs(title = "0 mg/mL BA") +
  theme(plot.title = element_text(size=20, hjust = 0.5)) +
  scale_x_discrete(limits = c("A18_ancestral", "A18_evolved","A03_ancestral", "A03_evolved","A08_ancestral", "A08_evolved", "A04_ancestral", "A04_evolved"))

P4ave     

ggsave(plot = P4ave, here("figures_out", "GCs", "230510_0mg", "230912_BA_GC_GrowthRate.jpg"), 
       width = 6.5, 
       height = 4.5)


#################################use for x-axis only#################################
write.csv(BAGC_fulldata1, here("data_in", "GCs", "230510_0mg", "BAGC_fulldata1.csv"))
#manually add a column (in excel) with the strain names
BAGC_fulldata1 <- read.csv(here("data_in", "GCs", "230510_0mg","BAGC_fulldata1_edited.csv"))


P4axis <- ggplot(data = BAGC_fulldata1,
                 mapping = aes(x = StrainName,
                               y = BAGClgr_avg,
                               colour = Origin))+
  scale_y_continuous(expand = c(0, 0), limits= c(-0.01, 0.51))+
  geom_sina(data = BAGC_fulldata1, mapping = aes(x = StrainName, y = BAGClgr_avg),
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
        legend.background = element_blank(), legend.box.background = element_blank(), 
        legend.key=element_blank(),
        panel.background = element_rect(
          fill = "white",
          colour = "black",
          linewidth = 0.4)) +
  scale_x_discrete(limits = c("FH1", "GC75","P75016","P78048")) +
  theme(text=element_text(family="Times New Roman", face="bold", size = 25),
        axis.text = element_text(size = 15, hjust = 0.5)) +
  labs(title = "0 mg/mL BA") +
  theme(plot.title = element_text(size=20, hjust = 0.5))

P4axis

ggsave(plot = P4axis, here("figures_out", "GCs", "230510_0mg","230912_BA_GC_axis.jpg"), 
       width = 6.5, 
       height = 4.5)

#################################t-test#################################
t.test(BAGClgr_avg~Origin, data = BAGC_fulldata1 %>% filter(Strain == "A03"))
#t = -4.9212, df = 3.0257, p-value = 0.01576
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.32698400 -0.07090809
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.2085248               0.4074708 

t.test(BAGClgr_avg~Origin, data = BAGC_fulldata1 %>% filter(Strain == "A04"))
#t = -1.8221, df = 3.0083, p-value = 0.1657
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.24197131  0.06560205
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.1825680               0.2707526 

t.test(BAGClgr_avg~Origin, data = BAGC_fulldata1 %>% filter(Strain == "A08"))
#t = -17.713, df = 4.7354, p-value = 1.645e-05
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.2685647 -0.1994828
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.1863649               0.4203887 

t.test(BAGClgr_avg~Origin, data = BAGC_fulldata1 %>% filter(Strain == "A18"))
#t = -7.9527, df = 3.7846, p-value = 0.001701
#alternative hypothesis: true difference in means between group ancestral and group evolved is not equal to 0
#95 percent confidence interval:
#  -0.2625352 -0.1243715
#sample estimates:
#  mean in group ancestral   mean in group evolved 
#0.2122811               0.4057344 