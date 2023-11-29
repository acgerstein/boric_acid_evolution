library(here)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(viridis)

Transfer_df <- read_delim(here("data_in","BAE", "230726_TransferFigureDates_BAE2.csv"), delim=",") %>% 
  select(Strain, Replicate, Concentration, Date_short, Date_long, Date)
Strain_Names <- tibble( Strain = c("A02","A03","A04","A08","A10","A12","A17","A18"), StrainName = c("P87", "GC75","P78048", "P75016", "P76055", "T101","SC5314", "FH1"))

colnames(Transfer_df)[2] ="Concentration"
colnames(Transfer_df)[3] ="Replicate"

Transfer_df <- Transfer_df[order(as.Date(Transfer_df$Date, format="%m/%d")),]
Transfer_df$Date <- as.factor(Transfer_df$Date)

Transfer_df <- Transfer_df %>%  mutate(Days = yday(ymd(Date_long))-267)

Transfer_df <- merge(Transfer_df, Strain_Names)

###########################################################################
p1 <- Transfer_df %>% filter(Strain %in% c("A18", "A03", "A08", "A04")) %>% 
  ggplot(aes(x = Days, y = Concentration)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(aes(group = Replicate), alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0, 14)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5), limits = c(0, 2)) +
  xlab("Days") +
  ylab("BA Concentration (mg/mL)") +
  facet_grid(Replicate~StrainName) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines",)),
        text=element_text(family="Times New Roman", face="bold"), 
        axis.title.x = element_text(size=30), axis.title.y = element_text(size=30),
        strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

p1

ggsave(here("figures_out","BAE","230831_AB.JPG"), plot= p1, 
       height = 12,
       width = 9)


p2 <- Transfer_df %>% filter(Strain %in% c("A10", "A12", "A02", "A17")) %>% 
  ggplot(aes(x = Days, y = Concentration)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(aes(group = Replicate), alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0, 14)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5), limits = c(0, 2)) +
  xlab("Days") +
  ylab("BA Concentration (mg/mL)") +
  facet_grid(Replicate~StrainName) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines",)),
        text=element_text(family="Times New Roman", face="bold"), 
        axis.title.x = element_text(size=30), axis.title.y = element_text(size=30),
        strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
p2

ggsave(here("figures_out","BAE","230831_other.JPG"), plot= p2, 
       height = 12,
       width = 9)