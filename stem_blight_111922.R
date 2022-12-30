# R script to analyze data for stem blight review
# Date: Nov. 2021

# packages
# install.packages("tidyverse")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggtext)

setwd("~/Dropbox/Botryosphaeria_manuscript/data/pathogenicity_figures")

############################################################
############ 2. figure for cultivar screening ##############
############################################################

data <- read_csv("cultivar_screening_052422.csv")
outfile <- "plot_111922"

# order cultivar based on mean lesion length across all literature
cultivarRank <- data %>%
  group_by(Cultivar) %>%
  summarise(mean = mean(LesionLength_15d, na.rm = T)) 

data_Polashock <- data %>%
  filter(Reference == "Polashock and Kramer 2006")

data_other <- data %>%
  filter(Reference != "Polashock and Kramer 2006")

# list of cultivars in Polashock paper
cultivar_Polashock <- data_Polashock$Cultivar
# list of cultivars in Smith 2009 but not included in Polashock paper
cultivar_Smith2009Unique <- (data %>% filter((Reference == "Smith 2009") & !(Cultivar %in% cultivar_Polashock)))$Cultivar
cultivar_all <- append(cultivar_Polashock, cultivar_Smith2009Unique)

# rank cultivars based on (1) lesion length in Polashock for cultivar_Polashock
# (2) lesion length in Smith 2009 for cultivar_Smith2009Unique
cultivarRank_Polashock <- data_Polashock %>%
  arrange(LesionLength_15d)

cultivarRank_Smith2009 <- data_other %>%
  filter(Cultivar %in% cultivar_Smith2009Unique) %>%
  filter(Reference == "Smith 2009") %>%
  arrange(LesionLength_15d)
cultivarRank <- append(cultivarRank_Polashock$Cultivar, cultivarRank_Smith2009$Cultivar)

data2 <- data %>%
  filter(Cultivar %in% cultivar_all) %>%
  mutate(Cultivar = factor(Cultivar, levels = cultivarRank)) %>%
  mutate(CultivarGroup = factor(ifelse(Cultivar %in% cultivar_Polashock, 
                                       "Polashock cultivar", "not Polashock cultivar"), 
                                levels = c("Polashock cultivar", "not Polashock cultivar", "NA"))) %>%
  arrange(Cultivar)
  
  
# data for plot
cultivarList <- tibble(Cultivar = factor(cultivarRank, levels = cultivarRank))
Polashock_data <- data2 %>%
  filter(Reference == "Polashock and Kramer 2006") %>%
  right_join(cultivarList, by = "Cultivar")
Smith2004_data<- data2%>%
  filter(Reference == "Smith 2004") %>%
  right_join(cultivarList, by = "Cultivar")
Smith2009_data<- data2%>%
  filter(Reference == "Smith 2009") %>%
  right_join(cultivarList, by = "Cultivar")
Creswell_data <- data2 %>%
  filter(Reference == "Creswell and Milholland 1987") %>%
  right_join(cultivarList, by = "Cultivar")

# assign different colors to species
colors <- data2 %>%
  mutate(colors = case_when(Type == "SHB" ~ "#FF6355", # red
                            Type == "RE" ~ "#FBA949", # orange
                            Type == "NHB" ~ "#8BD448", # green
                            Type == "HH" ~ "#2AA8F2", # blue
                            Type == "LB" ~ "#9C4F96", # purple
                            TRUE ~ as.character(Type)))%>%
  distinct(Cultivar, .keep_all = TRUE)

# draw plot
Polashock_plot <- ggplot(Polashock_data, aes(x = Cultivar, y = LesionLength_15d, color = CultivarGroup)) +
  scale_color_manual(values = c("blue", "gray50")) + 
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_markdown(angle = 90, color = colors$colors)) +
  ylab("Mean lesion length (mm)") + 
  ylim(0, 100)

Smith2004_plot <- ggplot(Smith2004_data, aes(x = Cultivar, y = LesionLength_15d, color = CultivarGroup)) +
  scale_color_manual(values = c("blue", "gray50")) + 
  geom_point(size = 2, shape = 8) +
  theme_bw() +
  theme(axis.text.x = element_markdown(angle = 90, color = colors$colors)) +
  ylab("Mean lesion length (mm)") + 
  ylim(0, 100)

Smith2009_plot <- ggplot(Smith2009_data, aes(x = Cultivar, y = LesionLength_15d, color = CultivarGroup)) +
  scale_color_manual(values = c("blue", "gray50")) + 
  geom_point(size = 2, shape = 8) +
  theme_bw() +
  theme(axis.text.x = element_markdown(angle = 90, color = colors$colors)) +
  ylab("Mean lesion length (mm)") + 
  ylim(0, 100)

Creswell_plot <- ggplot(Creswell_data, aes(x = Cultivar, y = LesionLength_15d, color = CultivarGroup)) +
  scale_color_manual(values = c("blue", "gray50")) + 
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_markdown(angle = 90, color = colors$colors)) +
  ylab("Mean lesion length (mm)") + 
  ylim(0, 100)

cultivarScrenningPlot <- ggarrange(Polashock_plot + rremove("legend") + rremove("xlab") + rremove("x.text") + rremove("ylab") + rremove("x.ticks"), 
                                   Creswell_plot+ rremove("legend") + rremove("xlab") + rremove("x.text") + rremove("ylab") + rremove("x.ticks"),  
                                   Smith2009_plot + rremove("legend") + rremove("xlab") + rremove("x.text") + rremove("ylab") + rremove("x.ticks"),  
                                   Smith2004_plot + rremove("legend") + rremove("ylab"), 
                                   ncol = 1, nrow = 4, heights = c(1, 1, 1, 1.8), align = 'v')

tiff(str_c(outfile, "/CultivarScreening.tiff"), width = 10, height = 8, units = "in", res = 300)
cultivarScrenningPlot
dev.off()

############################################################
############ 1. figure for pathogenicity test ##############
############################################################

p.data <- read_csv("pathogenicity_111922.csv") %>%
  drop_na(MeanLesionLength) %>%
  mutate(Reference = as.factor(Reference), Cultivar = as.factor(Cultivar), Species = as.factor(Species)) %>%
  mutate(Source = str_c(Reference, Method, sep = "-")) %>%
  filter(Method != "Detached stem/hardwood") %>%
  mutate(Species = factor(Species, levels = c("P. scoparia", "C. luteo-olivacea", "D. rudis", "B. dothidea", "N. arbuti", 
                                               "L. theobromae", "N. australe", "N. luteum", "N. parvum", "N. ribis")))

Patho.all.plot <- ggplot(p.data, aes(x = Method, y = MeanLesionLength, shape = Reference, color = Cultivar)) +
  geom_point() +
  facet_wrap(~Species, nrow = 1) + 
  scale_x_discrete(labels = c("Attached stem" = "Attached", "Detached stem" = "Detached")) +
  #scale_shape_manual(values = c(8, 19)) +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text.x = element_text(angle = 90, face = "italic")) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Mean lesion length (mm)")

tiff(str_c(outfile, "/Pathogenicity.all.tiff"), width = 10, height = 5, units = "in", res = 300)
Patho.all.plot
dev.off()

# compare detached with attached assay
m.data <- read_csv("pathogenicity_method_052522.csv") 

Pathogenicity.method.plot <- ggplot(m.data, aes(x = MeanLesionLengthDetached, y = MeanLesionLengthAttached, color = Species)) +
  geom_point() + 
  ylab("Mean lesion length from attached assay (mm)") +
  xlab("Mean lesion length from detached assay (mm)") +
  geom_errorbarh(aes(xmin=MeanLesionLengthDetached - sdDetached,
                     xmax=MeanLesionLengthDetached + sdDetached),
                 height=0.2)+
  geom_errorbar(aes(ymin=MeanLesionLengthAttached - sdAttached,
                    ymax=MeanLesionLengthAttached + sdAttached),
                width=0.1) +
  scale_y_continuous(limits = c(0, 80)) +
  scale_x_continuous(limits = c(0, 50))

tiff(str_c(outfile, "/Pathogenicity.method.tiff"), width = 10, height = 5, units = "in", res = 300)
Pathogenicity.method.plot
dev.off()


model <- lm(MeanLesionLengthDetached ~ MeanLesionLengthAttached, data = m.data)
summary(model)
summary(model)$r.squared
# 0.4893955
