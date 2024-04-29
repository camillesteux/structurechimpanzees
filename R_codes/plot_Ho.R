#!/usr/bin/env Rscript

#.libPaths("/media/camille/Donnees/R/4.2")
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(stringr)
library(scales)
library(forcats)
library(reshape)

setwd("~/stats/figures/")

# Import and format simulated Ho -----------------------------------------------

simu_het <- read.table("~/obs_het_Ranalysis.csv", sep = ",", header = TRUE)

simu_het <- simu_het[str_detect(simu_het$ind, "1"),] 
simu_het$subsp <- gsub('[[:digit:]]+', '', simu_het$ind)
simu_het$model <- as.factor(paste(str_split_i(simu_het$model, coll("_"), 4), str_split_i(simu_het$model, coll("_"), 5), str_split_i(simu_het$model, coll("_"), 6), sep="_"))
simu_het$subsp[simu_het$subsp=="schwein"] <- "schweinfurthii"
simu_het$subsp[simu_het$subsp=="troglo"] <- "troglodytes"
simu_het$emp <- "Simulated"
simu_het$empdim <- "Simu."

simu_het$common <- NA
simu_het[simu_het$subsp=="verus", ]$common <- "Western"
simu_het[simu_het$subsp=="ellioti", ]$common <- "Nigeria-Cameroon"
simu_het[simu_het$subsp=="troglodytes", ]$common <- "Central"
simu_het[simu_het$subsp=="schweinfurthii", ]$common <- "Eastern"
simu_het$common <- factor(simu_het$common, levels = c("Western", "Nigeria-Cameroon", "Central", "Eastern"))


# Import empirical Ho

## From Prado-Martinez et al. 2013
emp_data_PM <- read.table("~/H_Prado-Martinez.csv", header=TRUE)
emp_data_PM$subsp2 <- as.factor(str_split_i(emp_data_PM$Subsp, coll("_"), 3))
emp_data_PM <- emp_data_PM[, c(2,3,4)]
colnames(emp_data_PM)[3] <- "subsp"
emp_data_PM$model <- "Empirical [1]"
emp_data_PM$modeldim <- "Empi. [1]"

emp_data_PM$common <- NA
emp_data_PM[emp_data_PM$subsp=="verus", ]$common <- "Western"
emp_data_PM[emp_data_PM$subsp=="ellioti", ]$common <- "Nigeria-Cameroon"
emp_data_PM[emp_data_PM$subsp=="troglodytes", ]$common <- "Central"
emp_data_PM[emp_data_PM$subsp=="schweinfurthii", ]$common <- "Eastern"
emp_data_PM$common <- factor(emp_data_PM$common, levels = c("Western", "Nigeria-Cameroon", "Central", "Eastern"))

levels(simu_Pi$model)

## From de Manuel al. 2016

emp_data_DM <- read.table("~/H_deManuel2016.csv", header=TRUE, sep = ",")
emp_data_DM$subsp2 <- as.factor(str_split_i(emp_data_DM$Subsp, coll("_"), 3))
emp_data_DM <- emp_data_DM[, c(2,3,4)]
colnames(emp_data_DM)[3] <- "subsp"
emp_data_DM$model <- "Empirical [2]"
emp_data_DM$modeldim <- "Empi. [2]"

emp_data_DM$common <- NA
emp_data_DM[emp_data_DM$subsp=="verus", ]$common <- "Western"
emp_data_DM[emp_data_DM$subsp=="ellioti", ]$common <- "Nigeria-Cameroon"
emp_data_DM[emp_data_DM$subsp=="troglodytes", ]$common <- "Central"
emp_data_DM[emp_data_DM$subsp=="schweinfurthii", ]$common <- "Eastern"
emp_data_DM$common <- factor(emp_data_DM$common, levels = c("Western", "Nigeria-Cameroon", "Central", "Eastern"))

levels(simu_Pi$model)

lab <- c(expression("T"["CEWN"]*" = 700, T"["WN"]*" = 500 and T"["CE"]*" = 400 kya"),
         expression("T"["CEWN"]*" = 700, T"["WN"]*" = 600 and T"["CE"]*" = 500 kya"),
         expression("T"["CEWN"]*" = 800, T"["WN"]*" = 600 and T"["CE"]*" = 500 kya"),
         expression("T"["CEWN"]*" = 900, T"["WN"]*" = 500 and T"["CE"]*" = 400 kya"),
         expression("T"["CEWN"]*" = 900, T"["WN"]*" = 600 and T"["CE"]*" = 400 kya"),
         expression("T"["CEWN"]*" = 900, T"["WN"]*" = 600 and T"["CE"]*" = 500 kya"),
         expression("T"["CEWN"]*" = 900, T"["WN"]*" = 800 and T"["CE"]*" = 400 kya"),
         expression("T"["CEWN"]*" = 900, T"["WN"]*" = 800 and T"["CE"]*" = 500 kya"),
         expression("T"["CEWN"]*" = 900, T"["WN"]*" = 800 and T"["CE"]*" = 600 kya"))



# Plots ------------------------------------------------------------------------

### Ho (main text, Figure 8)
ggplot() +
  geom_point(data = simu_het[simu_het$model != "WNCE900000_WN800000_CE300000" & simu_het$model != "WNCE600000_WN500000_CE400000" & simu_het$subsp=="verus" | simu_het$subsp=="ellioti",], aes(y = simu.het*100, x=emp, color=model), position = position_dodge(width = 0.9), size=2) +
  geom_point(data = emp_data_PM[emp_data_PM$subsp2=="verus" | emp_data_PM$subsp2=="ellioti", ], aes(y = Het/10, x=model), fill = "dodgerblue3", color = "dodgerblue3", size=2) +
  geom_point(data = emp_data_DM[emp_data_DM$subsp=="verus" | emp_data_DM$subsp=="ellioti", ], aes(y = Het/10, x=model), fill = "dodgerblue4", color = "dodgerblue4", size=2) +
  facet_grid(. ~ common) +
  labs(y="Genetic diversity", color = "Simulated scenario") +
  scale_y_continuous(position = "right", limits = c(0.0, 0.2), breaks = c(0, 0.05, 0.10, 0.15, 0.20)) +
  scale_color_brewer(palette = "Set1", labels = lab) +
  theme_classic() +
  guides(color=guide_legend(nrow=3, byrow=FALSE, title.position="top"))+
  theme(strip.background    = element_blank(),
        axis.text.x.bottom = element_text(size = 14, color='black'),
        axis.text.x.top = element_text(size = 16, color='black'),
        axis.line.x         = element_blank(),
        axis.ticks.x        = element_line(colour = "black"),
        axis.ticks.length.x = unit(5, "points"),
        strip.text = element_text(size = 16),
        strip.text.y.left   = element_text(angle = 0, size = 10, hjust = 1, margin = margin(t = 2, r = 10, b = 2, l = 2)),
        panel.grid.major.y  = element_line(colour = "gray95"),
        #panel.grid.minor.y  = element_line(colour = "gray80"),
        #panel.grid.major.x  = element_line(colour = "gray90"),
        axis.line.x.bottom = element_line(colour="black"),
        legend.position = "none",
        legend.justification = "left",
        axis.title.x = element_blank(),
        axis.title.y.right = element_text(size = 14, margin = margin(t = 2, r = 2, b = 2, l = 15)),
        axis.text.y.right = element_text(size = 12, colour='black'),
        panel.spacing = unit(2, "lines"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

#ggsave(filename = "PAPER_ho_chimps_wesnig.png", width = 31, height = 11.67, dpi = 300, units = "cm", device='png') 

ggplot() +
  geom_point(data = simu_het[simu_het$model != "WNCE900000_WN800000_CE300000" & simu_het$model != "WNCE600000_WN500000_CE400000" & simu_het$subsp=="troglodytes" | simu_het$subsp=="schweinfurthii",], aes(y = simu.het*100, x=emp, color=model), position = position_dodge(width = 0.9), size=2) +
  geom_point(data = emp_data_PM[emp_data_PM$subsp2=="troglodytes" | emp_data_PM$subsp2=="schweinfurthii", ], aes(y = Het/10, x=model), fill = "dodgerblue3", color = "dodgerblue3", size=2) +
  geom_point(data = emp_data_DM[emp_data_DM$subsp=="troglodytes" | emp_data_DM$subsp=="schweinfurthii", ], aes(y = Het/10, x=model), fill = "dodgerblue4", color = "dodgerblue4", size=2) +
  facet_grid(. ~ common) +
  labs(y="Genetic diversity", color = "Simulated scenario") +
  scale_y_continuous(position = "right", limits = c(0.0, 0.2), breaks = c(0, 0.05, 0.10, 0.15, 0.20)) +
  scale_color_brewer(palette = "Set1", labels = lab) +
  theme_classic() +
  guides(color=guide_legend(nrow=3, byrow=FALSE, title.position="top"))+
  theme(strip.background    = element_blank(),
        axis.text.x.bottom = element_text(size = 14, color='black'),
        axis.text.x.top = element_text(size = 16, color='black'),
        axis.line.x         = element_blank(),
        axis.ticks.x        = element_line(colour = "black"),
        axis.ticks.length.x = unit(5, "points"),
        strip.text = element_text(size = 16),
        strip.text.y.left   = element_text(angle = 0, size = 10, hjust = 1, margin = margin(t = 2, r = 10, b = 2, l = 2)),
        panel.grid.major.y  = element_line(colour = "gray95"),
        #panel.grid.minor.y  = element_line(colour = "gray80"),
        #panel.grid.major.x  = element_line(colour = "gray90"),
        axis.line.x.bottom = element_line(colour="black"),
        legend.position = "bottom",
        legend.justification = "left",
        axis.title.x = element_blank(),
        axis.title.y.right = element_text(size = 14, margin = margin(t = 2, r = 2, b = 2, l = 15)),
        axis.text.y.right = element_text(size = 12, colour='black'),
        panel.spacing = unit(2, "lines"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

#ggsave(filename = "PAPER_ho_chimps_ceneas.png", width = 31, height = 15, dpi = 300, units = "cm", device='png') 

