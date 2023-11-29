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

###############################################################################
A18A3 <- Transfer_df %>% filter(Strain %in% c("A18", "A03")) %>% 
  ggplot(aes(x = Days, y = Concentration)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.6) +
  geom_line(aes(group = Replicate), alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 70, 7), limits = c(0, 70)) +
  scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4.25)) +
  xlab("Days") +
  ylab("BA Concentration (mg/mL)") +
  facet_grid(Replicate~StrainName) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines",)),
        text=element_text(family="Times New Roman", face="bold"),
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15))

ggsave(here("figures_out","BAE","231113_A18A03.JPG"), plot= A18A3, 
       height = 9,
       width = 6)

A8A4 <- Transfer_df %>% filter(Strain %in% c("A08", "A04")) %>% 
  ggplot(aes(x = Days, y = Concentration)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.6) +
  geom_line(aes(group = Replicate), alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 70, 7), limits = c(0, 70)) +
  scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4.25)) +
  xlab("Days") +
  ylab("BA Concentration (mg/mL)") +
  facet_grid(Replicate~StrainName) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines",)),
        text=element_text(family="Times New Roman", face="bold"),
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15))

ggsave(here("figures_out","BAE","231020_A08A04.JPG"), plot= A8A4, 
       height = 9,
       width = 6)

A10A12 <- Transfer_df %>% filter(Strain %in% c("A10", "A12")) %>% 
  ggplot(aes(x = Days, y = Concentration)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.6) +
  geom_line(aes(group = Replicate), alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 70, 7), limits = c(0, 70)) +
  scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4.25)) +
  xlab("Days") +
  ylab("BA Concentration (mg/mL)") +
  facet_grid(Replicate~StrainName) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines",)),
        text=element_text(family="Times New Roman", face="bold"),
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15))

ggsave(here("figures_out","BAE","231113_A10A12.JPG"), plot= A10A12, 
       height = 9,
       width = 6)

A02A17 <- Transfer_df %>% filter(Strain %in% c("A02", "A17")) %>% 
  ggplot(aes(x = Days, y = Concentration)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.6) +
  geom_line(aes(group = Replicate), alpha = 0.6) +
  scale_x_continuous(breaks = seq(0, 70, 7), limits = c(0, 70)) +
  scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4.25)) +
  xlab("Days") +
  ylab("BA Concentration (mg/mL)") +
  facet_grid(Replicate~StrainName) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines",)),
        text=element_text(family="Times New Roman", face="bold"),
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15))

ggsave(here("figures_out","BAE","231020_A02A17.JPG"), plot= A02A17, 
       height = 9,
       width = 6)