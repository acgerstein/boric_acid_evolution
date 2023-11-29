#Load libraries
library(tidyverse)
library(ggforce)
library(here)
library(Hmisc)

se <- function(x, ...) sqrt(var(x)/(length(x) - 1))

getwd()

#strain converstion:
#A02 - P87
#A03 - GC75
#A04 - P78048
#A08 - P75016
#A10 - P76055
#A12 - T101
#A17 - SC5314
#A18 - FH1

#order_BA <- c("A08", "A04", "A03",  "A18")
order_BA <- c("A18", "A03", "A08",  "A04")
#realnames_BA <- c("P75016", "P78048", "GC75", "FH1")
realnames_BA <- c("FH1", "GC75", "P75016", "P78048")
LA_BA <- read_csv(here("data_in", "LA", "23LA_data.csv"))

#removing outlier A18 T0 R7
library(dplyr)
LA_BA <- LA_BA %>% slice(-c(28,60,92,124,156,188,220,252,284,316,348,380,412,444,476,508,540,572,604,636,668,700,732,764)) 

#dark2 colours
Dark  <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")

#created "line" column, which is same thiug as replicate
LA_BA$line <-LA_BA$replicate

#Yana: averages across line
LA_BAave <-  LA_BA %>%
  group_by(drug, strain, line, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h))

#Yana: averages across strain
LA_BAave2 <-  LA_BA %>%
  group_by(strain, time) %>%
  summarise(OD24 = mean(OD_24h, na.rm=TRUE), OD48 = mean(OD_48h, na.rm=TRUE), OD72mix = mean(OD_72h_Mixed), OD72 = mean(OD_72h))

#Yana: if the time is zero then "darkred" and if time is not zero then "darkblue" 
LA_BAave$col <- ifelse(LA_BAave$time == 0, "darkred", "darkblue")

#Yana: this created a new column called "STL" with strain_time_line
LA_BAave$SLT <- paste(LA_BAave$strain, LA_BAave$line, LA_BAave$time, sep="_")

################################################################################################################################################################################################################################################################################################################################################
#0 mg/mL BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA0 <- subset(LA_BAave, drug==0)

BA0_t24_results <-data.frame()
BA0_t72_results <-data.frame()
for(i in unique(LA_BA0$strain)){   
  print(i)
  sub_BA0 <- subset(LA_BA0, strain==i) 
  if(length(subset(sub_BA0, time==0)$OD24) > 2){ 
    BA0test24<-t.test(subset(sub_BA0, time==0)$OD24, subset(sub_BA0, time==1)$OD24) 
    BA0test72<-t.test(subset(sub_BA0, time==0)$OD72, subset(sub_BA0, time==1)$OD72) 
    BA0_t24_results <-rbind(BA0_t24_results, c(i, round(BA0test24$estimate[2]-BA0test24$estimate[1], 3), round(BA0test24$statistic, 3), round(BA0test24$parameter, 2), round(BA0test24$p.value, 4))) 
    BA0_t72_results <-rbind(BA0_t72_results, c(i, round(BA0test72$estimate[2]-BA0test72$estimate[1], 3), round(BA0test72$statistic, 3), round(BA0test72$parameter, 2), round(BA0test72$p.value, 4))) 
  }
  else{ 
    BA0_t24_results <-rbind(BA0_t24_results, c(i, NA, NA, NA, NA)) 
    BA0_t72_results <-rbind(BA0_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA0_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA0_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

getwd()
#"C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript"

#Yana: this just saves the df as an csv file
write.csv(BA0_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA0_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA0_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA0_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230925BA0mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1)) 
for(j in order_BA) { 
  k <- k+1 
  sub_BA0<- subset(LA_BA0, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA0))), cex=2, sub_BA0$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), #cex makes the open circle points bigger
  col=sub_BA0$col, xaxt="n", yaxt="n") 
  points(jitter(rep(2, nrow(sub_BA0))), cex=2, sub_BA0$OD72, col=sub_BA0$col)
  for(i in sub_BA0$SLT){
    temp_BA0 <- subset(sub_BA0, SLT==i)
    points(c(1,2), c(temp_BA0$OD24,  temp_BA0$OD72), type="l", col=temp_BA0$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  #if(k==4) legend("bottomright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol")) #pch = the open cirlce beside anc and evol
  if(j %nin% c("A04")) text(1, max(sub_BA0$OD24)+0.2, "*", cex=3)
  if(j %nin% c("A04")) text(2, max(sub_BA0$OD72)+0.2, "*", cex=3)
  #text(0.8, 1.9, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
################################################################################################################################################################################################################################################################################################################################################
#0.2 mg/mL BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA02 <- subset(LA_BAave, drug==0.2)

BA02_t24_results <-data.frame()
BA02_t72_results <-data.frame()
for(i in unique(LA_BA02$strain)){   
  print(i)
  sub_BA02 <- subset(LA_BA02, strain==i) 
  if(length(subset(sub_BA02, time==0)$OD24) > 2){ 
    BA02test24<-t.test(subset(sub_BA02, time==0)$OD24, subset(sub_BA02, time==1)$OD24) 
    BA02test72<-t.test(subset(sub_BA02, time==0)$OD72, subset(sub_BA02, time==1)$OD72) 
    BA02_t24_results <-rbind(BA02_t24_results, c(i, round(BA02test24$estimate[2]-BA02test24$estimate[1], 3), round(BA02test24$statistic, 3), round(BA02test24$parameter, 2), round(BA02test24$p.value, 4))) #Yana: i don't get why we are binding t24_results if it's just an empty df- it's like adding a zero to a number? whats the point?
    BA02_t72_results <-rbind(BA02_t72_results, c(i, round(BA02test72$estimate[2]-BA02test72$estimate[1], 3), round(BA02test72$statistic, 3), round(BA02test72$parameter, 2), round(BA02test72$p.value, 4))) #Yana: which one is estimate[2] and which one ie estimate[1]? why are we subtracting them?
  }
  else{ 
    BA02_t24_results <-rbind(BA02_t24_results, c(i, NA, NA, NA, NA)) 
    BA02_t72_results <-rbind(BA02_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA02_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA02_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

#Yana: this just saves the df as an csv file
write.csv(BA02_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA02_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA02_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA02_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230925BA02mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1)) 
for(j in order_BA) { 
  k <- k+1 
  sub_BA02<- subset(LA_BA02, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA02))), cex=2, sub_BA02$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), col=sub_BA02$col, xaxt="n", yaxt="n") 
  points(jitter(rep(2, nrow(sub_BA02))), cex=2, sub_BA02$OD72, col=sub_BA02$col)
  for(i in sub_BA02$SLT){
    temp_BA02 <- subset(sub_BA02, SLT==i)
    points(c(1,2), c(temp_BA02$OD24,  temp_BA02$OD72), type="l", col=temp_BA02$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  #if(k==1) legend("bottomright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol")) #pch = the open cirlce beside anc and evol
  if(j %nin% c("A04")) text(1, max(sub_BA02$OD24)+0.2, "*", cex=3)
  if(j %nin% c("A10")) text(2, max(sub_BA02$OD72)+0.2, "*", cex=3)
  #text(0.8, 2.4, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
################################################################################################################################################################################################################################################################################################################################################
#0.4 mg/mL BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA04 <- subset(LA_BAave, drug==0.4)

BA04_t24_results <-data.frame()
BA04_t72_results <-data.frame()
for(i in unique(LA_BA04$strain)){   
  print(i)
  sub_BA04 <- subset(LA_BA04, strain==i) 
  if(length(subset(sub_BA04, time==0)$OD24) > 2){ 
    BA04test24<-t.test(subset(sub_BA04, time==0)$OD24, subset(sub_BA04, time==1)$OD24) 
    BA04test72<-t.test(subset(sub_BA04, time==0)$OD72, subset(sub_BA04, time==1)$OD72) 
    BA04_t24_results <-rbind(BA04_t24_results, c(i, round(BA04test24$estimate[2]-BA04test24$estimate[1], 3), round(BA04test24$statistic, 3), round(BA04test24$parameter, 2), round(BA04test24$p.value, 4))) #Yana: i don't get why we are binding t24_results if it's just an empty df- it's like adding a zero to a number? whats the point?
    BA04_t72_results <-rbind(BA04_t72_results, c(i, round(BA04test72$estimate[2]-BA04test72$estimate[1], 3), round(BA04test72$statistic, 3), round(BA04test72$parameter, 2), round(BA04test72$p.value, 4))) #Yana: which one is estimate[2] and which one ie estimate[1]? why are we subtracting them?
  }
  else{ 
    BA04_t24_results <-rbind(BA04_t24_results, c(i, NA, NA, NA, NA)) 
    BA04_t72_results <-rbind(BA04_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA04_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA04_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

#Yana: this just saves the df as an csv file
write.csv(BA04_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA04_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA04_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA04_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230922BA04mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1)) 
for(j in order_BA) { 
  k <- k+1 
  sub_BA04<- subset(LA_BA04, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA04))), cex=2, sub_BA04$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), col=sub_BA04$col, xaxt="n", yaxt="n") #Yana: WHAT IS HAPPENING
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub_BA04))), cex=2, sub_BA04$OD72, col=sub_BA04$col)
  for(i in sub_BA04$SLT){
    temp_BA04 <- subset(sub_BA04, SLT==i)
    points(c(1,2), c(temp_BA04$OD24,  temp_BA04$OD72), type="l", col=temp_BA04$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  #if(k==1) legend("bottomright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol")) #pch = the open cirlce beside anc and evol
  if(j %nin% c("A04")) text(1, max(sub_BA04$OD24)+0.2, "*", cex=3)
  if(j %nin% c("A03")) text(2, max(sub_BA04$OD72)+0.2, "*", cex=3)
  #text(0.8, 2.4, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
################################################################################################################################################################################################################################################################################################################################################
#0.78 BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA78 <- subset(LA_BAave, drug==0.78)

BA78_t24_results <-data.frame()
BA78_t72_results <-data.frame()
for(i in unique(LA_BA78$strain)){   
  print(i)
  sub_BA78 <- subset(LA_BA78, strain==i) 
  if(length(subset(sub_BA78, time==0)$OD24) > 2){ 
    BA78test24<-t.test(subset(sub_BA78, time==0)$OD24, subset(sub_BA78, time==1)$OD24) 
    BA78test72<-t.test(subset(sub_BA78, time==0)$OD72, subset(sub_BA78, time==1)$OD72) 
    BA78_t24_results <-rbind(BA78_t24_results, c(i, round(BA78test24$estimate[2]-BA78test24$estimate[1], 3), round(BA78test24$statistic, 3), round(BA78test24$parameter, 2), round(BA78test24$p.value, 4))) #Yana: i don't get why we are binding t24_results if it's just an empty df- it's like adding a zero to a number? whats the point?
    BA78_t72_results <-rbind(BA78_t72_results, c(i, round(BA78test72$estimate[2]-BA78test72$estimate[1], 3), round(BA78test72$statistic, 3), round(BA78test72$parameter, 2), round(BA78test72$p.value, 4))) #Yana: which one is estimate[2] and which one ie estimate[1]? why are we subtracting them?
  }
  else{ 
    BA78_t24_results <-rbind(BA78_t24_results, c(i, NA, NA, NA, NA)) 
    BA78_t72_results <-rbind(BA78_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA78_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA78_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

#Yana: this just saves the df as an csv file
write.csv(BA78_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA78_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA78_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA78_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230925BA78mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1)) 
for(j in order_BA) { 
  k <- k+1 
  sub_BA78<- subset(LA_BA78, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA78))), cex=2, sub_BA78$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), col=sub_BA78$col, xaxt="n", yaxt="n") #Yana: WHAT IS HAPPENING
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub_BA78))), cex=2, sub_BA78$OD72, col=sub_BA78$col)
  for(i in sub_BA78$SLT){
    temp_BA78 <- subset(sub_BA78, SLT==i)
    points(c(1,2), c(temp_BA78$OD24,  temp_BA78$OD72), type="l", col=temp_BA78$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  #if(k==1) legend("bottomright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol")) #pch = the open cirlce beside anc and evol
  if(j %nin% c("A03", "A04", "A18")) text(1, max(sub_BA78$OD24)+0.2, "*", cex=3)
  if(j %nin% c("A03","A08")) text(2, max(sub_BA78$OD72)+0.2, "*", cex=3)
  #text(0.8, 2.4, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
################################################################################################################################################################################################################################################################################################################################################
#1 mg/mL BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA1 <- subset(LA_BAave, drug==1)

BA1_t24_results <-data.frame()
BA1_t72_results <-data.frame()
for(i in unique(LA_BA1$strain)){   
  print(i)
  sub_BA1 <- subset(LA_BA1, strain==i) 
  if(length(subset(sub_BA1, time==0)$OD24) > 2){ 
    BA1test24<-t.test(subset(sub_BA1, time==0)$OD24, subset(sub_BA1, time==1)$OD24) 
    BA1test72<-t.test(subset(sub_BA1, time==0)$OD72, subset(sub_BA1, time==1)$OD72) 
    BA1_t24_results <-rbind(BA1_t24_results, c(i, round(BA1test24$estimate[2]-BA1test24$estimate[1], 3), round(BA1test24$statistic, 3), round(BA1test24$parameter, 2), round(BA1test24$p.value, 4))) #Yana: i don't get why we are binding t24_results if it's just an empty df- it's like adding a zero to a number? whats the point?
    BA1_t72_results <-rbind(BA1_t72_results, c(i, round(BA1test72$estimate[2]-BA1test72$estimate[1], 3), round(BA1test72$statistic, 3), round(BA1test72$parameter, 2), round(BA1test72$p.value, 4))) #Yana: which one is estimate[2] and which one ie estimate[1]? why are we subtracting them?
  }
  else{ 
    BA1_t24_results <-rbind(BA1_t24_results, c(i, NA, NA, NA, NA)) 
    BA1_t72_results <-rbind(BA1_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA1_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA1_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

#Yana: this just saves the df as an csv file
write.csv(BA1_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA1_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA1_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA1_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230922BA1mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1))
for(j in order_BA) { 
  k <- k+1 
  sub_BA1<- subset(LA_BA1, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA1))), cex=2, sub_BA1$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), col=sub_BA1$col, xaxt="n", yaxt="n")
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub_BA1))), cex=2, sub_BA1$OD72, col=sub_BA1$col)
  for(i in sub_BA1$SLT){
    temp_BA1 <- subset(sub_BA1, SLT==i)
    points(c(1,2), c(temp_BA1$OD24,  temp_BA1$OD72), type="l", col=temp_BA1$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  #if(k==1) legend("topright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol")) #pch = the open cirlce beside anc and evol
  if(j %nin% c("A03", "A04","A18")) text(1, max(sub_BA1$OD24)+0.2, "*", cex=3)
  if(j %nin% c("A18")) text(2, max(sub_BA1$OD72)+0.2, "*", cex=3)
  #text(0.8, 2.4, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
################################################################################################################################################################################################################################################################################################################################################
#1.25 mg/mL BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA125 <- subset(LA_BAave, drug==1.25)

BA125_t24_results <-data.frame()
BA125_t72_results <-data.frame()
for(i in unique(LA_BA125$strain)){   
  print(i)
  sub_BA125 <- subset(LA_BA125, strain==i) 
  if(length(subset(LA_BA125, time==0)$OD24) > 2){ 
    BA125test24<-t.test(subset(sub_BA125, time==0)$OD24, subset(sub_BA125, time==1)$OD24) 
    BA125test72<-t.test(subset(sub_BA125, time==0)$OD72, subset(sub_BA125, time==1)$OD72) 
    BA125_t24_results <-rbind(BA125_t24_results, c(i, round(BA125test24$estimate[2]-BA125test24$estimate[1], 3), round(BA125test24$statistic, 3), round(BA125test24$parameter, 2), round(BA125test24$p.value, 4))) #Yana: i don't get why we are binding t24_results if it's just an empty df- it's like adding a zero to a number? whats the point?
    BA125_t72_results <-rbind(BA125_t72_results, c(i, round(BA125test72$estimate[2]-BA125test72$estimate[1], 3), round(BA125test72$statistic, 3), round(BA125test72$parameter, 2), round(BA125test72$p.value, 4))) #Yana: which one is estimate[2] and which one ie estimate[1]? why are we subtracting them?
  }
  else{ 
    BA125_t24_results <-rbind(BA125_t24_results, c(i, NA, NA, NA, NA)) 
    BA125_t72_results <-rbind(BA125_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA125_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA125_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

#Yana: this just saves the df as an csv file
write.csv(BA125_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA125_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA125_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA125_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230922BA125mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1))
for(j in order_BA) { 
  k <- k+1 
  sub_BA125<- subset(LA_BA125, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA125))), cex=2, sub_BA125$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), col=sub_BA125$col, xaxt="n", yaxt="n") #Yana: WHAT IS HAPPENING
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub_BA125))), cex=2,sub_BA125$OD72, col=sub_BA125$col)
  for(i in sub_BA125$SLT){
    temp_BA125 <- subset(sub_BA125, SLT==i)
    points(c(1,2), c(temp_BA125$OD24,  temp_BA125$OD72), type="l", col=temp_BA125$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  #if(k==1) legend("topright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol")) #pch = the open cirlce beside anc and evol
  if(j %nin% c("A18", "A04", "A08")) text(1, max(sub_BA125$OD24)+0.2, "*", cex=3)
  if(j %nin% c("A04", "A18")) text(2, max(sub_BA125$OD72)+0.2, "*", cex=3)
  #text(0.8, 2.4, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
################################################################################################################################################################################################################################################################################################################################################
#1.5 mg/mL BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA15 <- subset(LA_BAave, drug==1.5)

BA15_t24_results <-data.frame()
BA15_t72_results <-data.frame()
for(i in unique(LA_BA15$strain)){   
  print(i)
  sub_BA15 <- subset(LA_BA15, strain==i) 
  if(length(subset(LA_BA15, time==0)$OD24) > 2){ 
    BA15test24<-t.test(subset(sub_BA15, time==0)$OD24, subset(sub_BA15, time==1)$OD24) 
    BA15test72<-t.test(subset(sub_BA15, time==0)$OD72, subset(sub_BA15, time==1)$OD72) 
    BA15_t24_results <-rbind(BA15_t24_results, c(i, round(BA15test24$estimate[2]-BA15test24$estimate[1], 3), round(BA15test24$statistic, 3), round(BA15test24$parameter, 2), round(BA15test24$p.value, 4))) #Yana: i don't get why we are binding t24_results if it's just an empty df- it's like adding a zero to a number? whats the point?
    BA15_t72_results <-rbind(BA15_t72_results, c(i, round(BA15test72$estimate[2]-BA15test72$estimate[1], 3), round(BA15test72$statistic, 3), round(BA15test72$parameter, 2), round(BA15test72$p.value, 4))) #Yana: which one is estimate[2] and which one ie estimate[1]? why are we subtracting them?
  }
  else{ 
    BA15_t24_results <-rbind(BA15_t24_results, c(i, NA, NA, NA, NA)) 
    BA15_t72_results <-rbind(BA15_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA15_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA15_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

#Yana: this just saves the df as an csv file
write.csv(BA15_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA15_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA15_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA15_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230922BA15mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1))
for(j in order_BA) { 
  k <- k+1 
  sub_BA15<- subset(LA_BA15, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA15))), cex=2, sub_BA15$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), col=sub_BA15$col, xaxt="n", yaxt="n") #Yana: WHAT IS HAPPENING
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub_BA15))), cex=2, sub_BA15$OD72, col=sub_BA15$col)
  for(i in sub_BA15$SLT){
    temp_BA15 <- subset(sub_BA15, SLT==i)
    points(c(1,2), c(temp_BA15$OD24,  temp_BA15$OD72), type="l", col=temp_BA15$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  #if(k==1) legend("topright", pch=21, col=c("red", "darkblue"), legend=c("Anc", "Evol")) #pch = the open cirlce beside anc and evol
  #if(j %nin% c("A18", "A03", "A08")) text(1, max(sub_BA15$OD24)+0.2, "*", cex=2)
  if(j %nin% c("A03", "A04", "A18")) text(2, max(sub_BA15$OD72)+0.2, "*", cex=3)
  #text(0.8, 2.4, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
################################################################################################################################################################################################################################################################################################################################################
#1.75 mg/mL BA
################################################################################################################################################################################################################################################################################################################################################
LA_BA175 <- subset(LA_BAave, drug==1.75)

BA175_t24_results <-data.frame()
BA175_t72_results <-data.frame()
for(i in unique(LA_BA175$strain)){   
  print(i)
  sub_BA175 <- subset(LA_BA175, strain==i) 
  if(length(subset(LA_BA175, time==0)$OD24) > 2){ 
    BA175test24<-t.test(subset(sub_BA175, time==0)$OD24, subset(sub_BA175, time==1)$OD24) 
    BA175test72<-t.test(subset(sub_BA175, time==0)$OD72, subset(sub_BA175, time==1)$OD72) 
    BA175_t24_results <-rbind(BA175_t24_results, c(i, round(BA175test24$estimate[2]-BA175test24$estimate[1], 3), round(BA175test24$statistic, 3), round(BA175test24$parameter, 2), round(BA175test24$p.value, 4))) #Yana: i don't get why we are binding t24_results if it's just an empty df- it's like adding a zero to a number? whats the point?
    BA175_t72_results <-rbind(BA175_t72_results, c(i, round(BA175test72$estimate[2]-BA175test72$estimate[1], 3), round(BA175test72$statistic, 3), round(BA175test72$parameter, 2), round(BA175test72$p.value, 4))) #Yana: which one is estimate[2] and which one ie estimate[1]? why are we subtracting them?
  }
  else{ 
    BA175_t24_results <-rbind(BA175_t24_results, c(i, NA, NA, NA, NA)) 
    BA175_t72_results <-rbind(BA175_t72_results, c(i, NA, NA, NA, NA))
  }
}

#Yana: this just names the columns cause they don't have names it's just like "NA_character_., NA_character_..1, etc."
names(BA175_t24_results) <- c("strain", "effectsize", "stat", "df", "p")
names(BA175_t72_results) <- c("strain", "effectsize", "stat", "df", "p")

#Yana: this just saves the df as an csv file
write.csv(BA175_t24_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA175_24hOD_ttest.csv", row.names=FALSE)
write.csv(BA175_t72_results, "C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/data_out/LA\\LA_BA175_72hOD_ttest.csv", row.names=FALSE)

pdf("C:/Users/oksana/Nextcloud/2022Yana/Boric_Acid/BAE_Manuscript/figures_out/LA/230925BA175mg.pdf", width=2, height=8)
k <- 0 
par(mfrow=c(4, 1), mar=c(1,0.5,1,0.5), oma=c(3.5, 3.5, 1, 1))
for(j in order_BA) { 
  k <- k+1 
  sub_BA175<- subset(LA_BA175, strain==j) 
  plot(jitter(rep(1, nrow(sub_BA175))), cex=2, sub_BA175$OD24, ylim=c(0, 2.5), xlim=c(0.8, 2.2), col=sub_BA175$col, xaxt="n", yaxt="n") #Yana: WHAT IS HAPPENING
  #points(jitter(rep(2, nrow(sub))), sub$OD48, col=sub$col)
  points(jitter(rep(2, nrow(sub_BA175))), cex=2, sub_BA175$OD72, col=sub_BA175$col)
  for(i in sub_BA175$SLT){
    temp_BA175 <- subset(sub_BA175, SLT==i)
    points(c(1,2), c(temp_BA175$OD24,  temp_BA175$OD72), type="l", col=temp_BA175$col)
  }
  if(k %% 5 ==4) axis(1, at=c(1, 2), labels=c("24", "72"), cex.axis=2)
  else axis(1, at=c(0, 0), labels=FALSE)
  
  if(k %% 4 ==1) axis(2, las=2, cex.axis = 2)
  #else axis(2, labels=FALSE)
  else axis(2, las=2, cex.axis = 2)
  
  if(k==1) legend("topright", pch=21, col=c("red", "darkblue"), legend=c(as.expression(bquote(bold("Ancestral"))), as.expression(bquote(bold("Evolved")))), cex=1.5, pt.cex = 2, bty = "n") #pch = the open cirlce beside anc and evol
  #if(j %nin% c("A18", "A03", "A08")) text(1, max(sub_BA175$OD24)+0.2, "*", cex=2)
  if(j %nin% c("A03", "A08", "A18")) text(2, max(sub_BA175$OD72)+0.2, "*", cex=3)
  #text(0.8, 2.4, realnames_BA[k], cex=1.2, pos=4) #pos values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified (x,y) coordinates.
}
#mtext("Optical density (OD)", side=2, line=2, outer=TRUE)
#mtext("Time of reading (hr)", side=1, line=2, outer=TRUE)
dev.off()
