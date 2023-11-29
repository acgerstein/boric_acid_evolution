#Change "Line" names for the evolved strains to "BAEA4" or "BAEA8" in the csv
#Change "Type" to a continuous value by removing the letters "evo" or "anc" and only leaving numbers in the csv

#upload the raw data
#RAD20 from 24hr data + FoG20 from 48hr data

library(here)
library(tidyverse)
library(lmerTest)
library(ggpubr)
library(viridis)

getwd()

#strip off the "EVO/T0"
df <- read.csv(here("data_in","DDA", "fortyeighthourdataset_df_edited.csv"))
df$name <- sub("_EVO", " ", df$name)
df$name <- sub("_T0", " ", df$name)

#remove blank
df <- df[-1,]

df <- df %>% select(name, RAD20, FoG20) %>%  separate(name, c("Strain", "Replicate"), sep = "_R") 

Strain_Names <- tibble( Strain = c("A03","A04","A08","A18"), StrainName = c("GC75","P78048", "P75016", "FH1"))


df <- df %>% mutate(Replicate = gsub("([0-9])([a-zA-Z])","\\1 \\2", Replicate)) %>%
  separate(Replicate, c("Replicate", "TechnicalReplicate")) %>% 
  group_by(Strain, Replicate) %>% 
  summarise(RAD20 = mean(RAD20), 
            FoG20 = mean(FoG20)) %>%
  mutate(Strain = gsub("([a-zA-Z])([0-9])","\\1 \\2", Strain)) %>% 
  separate(Strain, c("Origin", "StrainNo"))

df$Origin[df$Origin == "A"] <- "Ancestral"
df$Origin[df$Origin == "BAEA"] <- "Evolved"
df <- df %>% mutate(A = "A") %>% unite(Strain, c(A, StrainNo), sep = "")


df <- df %>% 
  group_by(Origin, Strain, Replicate) %>% 
  summarise(RAD20 = mean(RAD20), 
            FoG20 = mean(FoG20))

#removing outlier (A18 T0 R7)
df <- df[-15,]

################################################################################################################################
#Function that will calculate the mean, median, sd, and se
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      median=median(x[[col]],na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se= sd(x[[col]])/sqrt(length(x[[col]])))}
  
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#Applying the function to RAD20 and FoG50
RAD <- data_summary(df, varname=c("RAD20"), 
                  groupnames=c("Origin", "Strain"))
FoG <- data_summary(df, varname=c("FoG20"), 
                  groupnames=c("Origin", "Strain"))

#Detach the plyr package 
detach("package:plyr", unload = TRUE)

df <- merge(df, Strain_Names)

R_p1 <- ggplot(df) +
  geom_jitter(mapping= aes(x = Origin, y=RAD20), maxwidth = 1, alpha = 0.8, size = 2.5, position = position_jitter(width = .05)) + 
  stat_summary(df,mapping = aes(x = Origin, y = RAD20), 
               fun.y = "mean", 
               geom = "crossbar",
               width = 0.3,
               size = 0.3,
               colour="red") +
  scale_y_reverse() +
  facet_wrap(~StrainName, nrow = 1) +
  xlab("Replicate") +
  ylab(expression(RAD[20]  (mm))) +
  expand_limits( y = c(30, 0))+
  theme_bw() + 
  theme(strip.background = element_rect(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 20),
        axis.text = element_text(size = 15, hjust = 0.5)) 

R_p1

T_p1 <- ggplot(df) +
  geom_jitter(mapping= aes(x = Origin, y=FoG20), maxwidth = 1,alpha = 0.7, size = 2.5, position = position_jitter(width = .05)) + 
  stat_summary(df,mapping = aes(x = Origin, y = FoG20), 
               fun.y = "mean", 
               geom = "crossbar",
               width = 0.3,
               size = 0.3,
               colour="red") +
  facet_wrap(~StrainName, nrow = 1)+ 
  expand_limits( y = c(0, 0.4))+
  xlab("Replicate") +
  ylab(expression(FoG[20])) + 
  theme_bw() + 
  theme(strip.background = element_rect(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 20),
        axis.text = element_text(size = 15, hjust = 0.5)) 
T_p1

RT_p1 <- ggarrange(R_p1, T_p1, nrow = 2, align = "hv")
RT_p1

ggsave(plot = RT_p1, here("figures_out","DDA","230912_BAE_RT1.jpg"), 
       width = 9, 
       height = 7)

################################################################################################################################
t.test(RAD20~Origin, data = df %>% filter(Strain == "A03"))
#data:  RAD20 by Origin
#t = -3.0888, df = 3.551, p-value = 0.04295
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
#  -3.9721414 -0.1111919
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#11.66667                13.70833

t.test(FoG20~Origin, data = df %>% filter(Strain == "A03"))
#data:  FoG20 by Origin
#t = -1.6757, df = 4.977, p-value = 0.1549
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
#  -0.06763075  0.01429741
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#0.1225000               0.1491667 

t.test(RAD20~Origin, data = df %>% filter(Strain == "A04"))
#data:  RAD20 by Origin
#t = -0.59062, df = 4.4118, p-value = 0.5837
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
#  -1.1526636  0.7359969
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#13.25000                13.45833 

t.test(FoG20~Origin, data = df %>% filter(Strain == "A04"))
#data:  FoG20 by Origin
#t = -0.43471, df = 5.5537, p-value = 0.6801
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
#  -0.02246753  0.01580086
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#0.1275000               0.1308333 

t.test(RAD20~Origin, data = df %>% filter(Strain == "A08"))
#data:  RAD20 by Origin
#t = 0, df = 3.3844, p-value = 1
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
#  -1.939061  1.939061
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#12.25                   12.25 

t.test(FoG20~Origin, data = df %>% filter(Strain == "A08"))
#data:  FoG20 by Origin
#t = 1.8891, df = 5.0481, p-value = 0.1169
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
#  -0.008624328  0.056957661
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#0.1541667               0.1300000 

t.test(RAD20~Origin, data = df %>% filter(Strain == "A18"))
#data:  RAD20 by Origin
#t = 2.1915, df = 4.6212, p-value = 0.0844
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
# -0.2418816  2.6307705
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#13.77778                12.58333

t.test(FoG20~Origin, data = df %>% filter(Strain == "A18"))
#data:  FoG20 by Origin
#t = -0.2582, df = 4.1468, p-value = 0.8086
#alternative hypothesis: true difference in means between group Ancestral and group Evolved is not equal to 0
#95 percent confidence interval:
# -0.01934115  0.01600782 
#sample estimates:
#  mean in group Ancestral   mean in group Evolved 
#              0.1100000               0.1116667