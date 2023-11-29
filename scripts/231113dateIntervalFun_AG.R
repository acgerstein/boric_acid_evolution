library(here)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(viridis)

Transfer_df <- read_delim(here("data_in","BAE", "231113_TransferFigureDates_BAE2b.csv"), delim=",") %>% 
  select(Strain, Replicate, Concentration, Date_short, Date_long, Date)
Strain_Names <- tibble( Strain = c("A02","A03","A04","A08","A10","A12","A17","A18"), StrainName = c("P87", "GC75","P78048", "P75016", "P76055", "T101","SC5314", "FH1"))

Transfer_df <- Transfer_df[order(as.Date(Transfer_df$Date, format="%m/%d")),]
Transfer_df$Date <- as.factor(Transfer_df$Date)

Transfer_df <- Transfer_df %>%  mutate(Weeks = week(ymd(Transfer_df$Date_long)) - 41,
                                       Days = qday(ymd(Date_long))-10)

Transfer_df <- merge(Transfer_df, Strain_Names)
Transfer_df$SR <- paste(Transfer_df$Strain, Transfer_df$Replicate, sep="_")

#split up the datframe so that each "SR" = strain x rep (i.e., each isolate) becomes it's own data frame in something called a list. You can access list components using double square brackets like: list[[1]] where list[[1]] is a single dataframe and list is the collection dataframe
indivRep_ddn <- split(Transfer_df, Transfer_df$SR)
length(indivRep_ddn) # = 95 - 95 evoled isolates, each with their own dataframe

indivRep_ddn[[1]] # this is the dataframe for SR = "A02_R01"
indivRep_ddn[[12]] # this is the dataframe for SR = "A03_R01" (we are missing an A01 replicate)

#testing code - the interval function gives you the interval between one day (indivRep_ddn[[1]]$Date_long[1] = "2019-10-10") and another (indivRep_ddn[[1]]$Date_long[4] = "2019-10-16"). as.interval calculates the interval length into days (because we specified days at the end)
as.period(interval(indivRep_ddn[[1]]$Date_long[1], indivRep_ddn[[1]]$Date_long[4])) %/% days(1)

#then I use that code that I worked out to write a look where instead of hardcoding I wanted A02_R01 (indivRep_ddn[[1]]) and specific dates [1] and [4] I'm going to use a loop to iterate through all datasets (isolates) and calcuate the interval between all consecutive transfer dates within them, i.e., transfer 1-2, 2-3, 3-4 ...

#create a new list for the intervals to go into (because there will be one fewer row than the list created above - the number of intervals is the number of transfers -1)
dayInterval_ddn <- list()

# use k as the marker variable for the number of isolates (length(indivRep_ddn) = 95) and run this loop 95 times
for(k in 1:length(indivRep_ddn)){
  
  #set up vectors to hold information about the interval and drug concentration at first (conc1) and second (conc2) transfer
    dayInterval <- c()
    conc1 <- c()
    conc2 <- c()
    
    # run a loop to iterate through all consecutive transfers within a dataframe from the first transfer to the last; there is one fewer iteration than transfer since the last transfer doesn't have a parnter
    for (i in 1:nrow(indivRep_ddn[[k]])-1) {
      # calculate the interval as above
      dayInterval[i]<- as.period(interval(indivRep_ddn[[k]]$Date_long[i], indivRep_ddn[[k]]$Date_long[i+1])) %/% days(1)
      # store the drug concentration of the first transfer
      conc1[i] <- indivRep_ddn[[k]]$Concentration[i]
      # store the drug concentration of the second transfer
      conc2[i] <- indivRep_ddn[[k]]$Concentration[i+1]
    }
    # write a new dataframe for the focal SR with information about strain (using 
    # rep(indivRep_ddn[[k]]$Strain[1], nrow(indivRep_ddn[[k]])-1) to repeat the strain and replicate information X times where X is the number of transfers -1, nrow(indivRep_ddn[[k]])-1), the stored concentration information, and the calculated tranfer intervals
    df <- data.frame(strain = rep(indivRep_ddn[[k]]$Strain[1], nrow(indivRep_ddn[[k]])-1), replicate = rep(indivRep_ddn[[k]]$Replicate[1], nrow(indivRep_ddn[[k]])-1), conc1, conc2, dayInterval)
    # print(df)
    #add this dataframe into the new list that was created using 'append'
    dayInterval_ddn[[k]] <- df
  }


#Take the list just created and convert it to a dataframe
dayInterval_df <- do.call(rbind.data.frame, dayInterval_ddn)

write.csv(dayInterval_df, here("data_out","BAE", "231122_dayInterval.csv"), row.names=FALSE)

dayInterval_df_sameConc <- subset(dayInterval_df, dayInterval_df$conc1 == dayInterval_df$conc2)
dayInterval_df_diffConc <- subset(dayInterval_df, dayInterval_df$conc1 != dayInterval_df$conc2)

par(mfrow=c(2, 1), mar=c(1, 1, 2, 1), oma=c(3, 3, 1, 1))
plot(jitter(dayInterval_df_sameConc$conc1), dayInterval_df_sameConc$dayInterval, col= "black", xlim=c(0.5, 3), ylim=c(0, 8.5), xaxt="n", yaxt="n")
axis(2, las=2)
axis(1, labels=FALSE)
mtext("a) Transfer to the same concentration", side=3, adj=0.01)

plot(jitter(dayInterval_df_diffConc$conc1), dayInterval_df_diffConc$dayInterval, col= "black", xlim=c(0.5, 3), ylim=c(0, 8.5), xaxt="n", yaxt="n")
axis(2, las=2)
axis(1)
mtext("Days between transfers", outer = TRUE, side=2, line=1)
mtext("Initial BA concentration (mg/mL)", outer = FALSE, side=1, line=2)
mtext("b) Transfer to a higher concentration", side=3, adj=0.01)