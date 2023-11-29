library(here)
library(tidyverse)
#upload the raw data
T_df <- read_delim(here("data_in","BAE","230831_Transfer_Data.csv"), delim=",") %>% mutate(Transfer=as.numeric(gsub("t","",Transfer)))

Strain_Names <- tibble( Strain = c("A02","A03","A04","A08","A10","A12","A17","A18"), StrainName = c("P87", "GC75","P78048", "P75016", "P76055", "T101","SC5314", "FH1"))

T_df <- merge(T_df, Strain_Names)
# How many strains survived in 3mg/mL BA 
T_df %>% filter(Concentration >= 3) %>% group_by(Strain) %>% 
  summarise(Replicate = unique(Replicate)) %>% 
  filter(!Replicate == 0)
# 9 

S <- T_df %>%filter(!Strain%in%c("A02","A04", "A08"))

SBA <- S %>%
  group_by(Concentration) %>% 
  summarise(min = min(Transfer)-0.2, max = max(Transfer)+0.2)
SBA$Concentration<-as.factor(SBA$Concentration)

S248 <- T_df %>% filter(Strain %in% c("A02","A04", "A08"))
S2BA <- S248 %>% 
  group_by(Concentration) %>% 
  summarise(min = min(Transfer)-0.2, max = max(Transfer)+0.2)
S2BA$Concentration<-as.factor(S2BA$Concentration)


col <- c("P87" = "#E41A1C", "GC75" = "#377EB8","P78048" = "#4DAF4A",
       "P75016" = "#984EA3","P76055" = "#FF7F00","T101" = "#FFFF33",
       "SC5314" = "#A65628","FH1" = "#F781BF")

colfunc <- colorRampPalette(c("blue", "red"))
colfunc(10)  

S1 <- S %>%filter(!Strain%in%c("A04","A08", "A03", "A18"))
S2 <- S %>%filter(!Strain%in%c("A10", "A12", "A17"))

p1 <- ggplot(data = S) +
  geom_point(mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName))+
  coord_cartesian(xlim = c(-0.15,27.75),ylim = c(-0.1, 12.2), clip = 'off') +
  geom_line(S1,mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName),
            size =1.2, linetype=5) +
  geom_line(S2,mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName),
            size =1.2) +
  geom_rect(SBA,mapping=aes(xmin=min,xmax=max,fill=Concentration),
            ymin=-Inf,ymax=Inf, alpha=0.15)+ 
  geom_text(SBA,mapping=aes(x=(min+max)/2,y=12,label=Concentration),
            size=4,vjust=-1.5, family="Times New Roman"
            )+
  scale_x_continuous(breaks = seq(0, 28, by = 4),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0, 12, by = 2),expand = c(0, 0))+theme_bw()+
  guides(fill = FALSE) +
  ylab("Number of Surviving Replicates") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines",)),
        axis.title.x = element_text(size=21), axis.title.y = element_text(size=21),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        text=element_text(family="Times New Roman", face="bold"))+
 scale_fill_manual(values = colfunc(10)) +
  scale_color_manual(values = col)
p1

S2481 <- S248 %>%filter(!Strain%in%c("A04", "A08"))
S2482 <- S248 %>%filter(!Strain%in%c("A02"))

p2 <- ggplot(data = S248) +
  geom_point(mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName))+
  coord_cartesian(xlim = c(-0.15,26.75),ylim = c(-0.1, 12.2), clip = 'off') +
  geom_line(S2481,mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName),size =1.2, linetype=5) +
  geom_line(S2482,mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName),size =1.2) +
  geom_rect(S2BA,mapping=aes(xmin=min,xmax=max,fill=Concentration),
            ymin=-Inf,ymax=Inf, alpha=0.2)+ 
  geom_text(S2BA,mapping=aes(x=(min+max)/2,y=12,label=Concentration),
            size=4,vjust=-1.5, family="Times New Roman")+
  scale_x_continuous(breaks = seq(0, 28, by = 4),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0, 12, by = 2),expand = c(0, 0))+
  theme_bw()+ guides(fill = FALSE)+
  ylab("Number of Surviving Replicates") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=21), axis.title.y = element_text(size=21),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.margin = margin(unit(c(20, 1, 1, 1),  "lines")), 
        text=element_text(family="Times New Roman", face="bold"))+
  scale_fill_manual(values = colfunc(10)) +
  scale_color_manual(values = col)

p2

ggsave(here("figures_out","BAE", "231121_TransferFigure_a.JPG"), plot= p1, height = 5, width = 7)
ggsave(here("figures_out", "BAE", "231121_TransferFigure_b.JPG"), plot= p2, height = 5, width = 7)

#Legend 
Legend <- ggplot(data = T_df) +
  geom_point(mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName))+
  coord_cartesian(xlim = c(-0.15,29),ylim = c(-0.1, 12.2), clip = 'off') +
  geom_line(mapping = aes(x = Transfer, y = Replicate, colour = StrainName, group = StrainName),size =1.2) +
  
  geom_text(S2BA,mapping=aes(x=(min+max)/2,y=12,label=Concentration),
            size=4,vjust=-1.5, family="Times New Roman")+
  scale_x_continuous(breaks = seq(0, 24, by = 4),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0, 12, by = 2),expand = c(0, 0))+
  theme_bw()+ guides(fill = FALSE)+
  
  theme(plot.margin = margin(unit(c(20, 1, 1, 1),  "lines")), 
        text=element_text(family="Times New Roman", size = 15))+
 scale_color_manual(values=col)

ggsave(here("figures_out","BAE","231121_JustUseForLegend.JPG"), plot= Legend,
       height = 5, width = 7)
