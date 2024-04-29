#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(stringr)
library(scales)
library(forcats)
library(reshape)

setwd("~/stats")

# Import and format Fst from Lester et al. 2021 --------------------------------

my_data <- read.table("Fst_Lester_Camille_figures.csv", header = TRUE,  fill=TRUE, na.strings="")
rownames(my_data) <- colnames(my_data)[2:length(colnames(my_data))]
my_data <- my_data[,2:ncol(my_data)]

stats_lester <- melt(as.matrix(my_data), varnames=c("subsp1","subsp2"))
stats_lester <- na.omit(stats_lester)
stats_lester$pair <- as.factor(paste(str_split_i(stats_lester$subsp2, coll( "."), 1), str_split_i(stats_lester$subsp1, coll( "."), 1), sep = " - "))

colnames(stats_lester)[3] <- "Fst"
stats_lester$pair <- factor(stats_lester$pair, levels = c("verus - verus", "verus - ellioti", "verus - troglodytes", "verus - schweinfurthii",
                                              "ellioti - ellioti", "ellioti - troglodytes", "ellioti - schweinfurthii", 
                                              "troglodytes - troglodytes", "troglodytes - schweinfurthii", "schweinfurthii - schweinfurthii"))
stats_lester$cat <- NA
stats_lester[str_split_i(stats_lester$subsp1, coll( "."), 1)==str_split_i(stats_lester$subsp2, coll( "."), 1), ]$cat <- "within"
stats_lester[str_split_i(stats_lester$subsp1, coll( "."), 1)!=str_split_i(stats_lester$subsp2, coll( "."), 1), ]$cat <- "between"

stats_lester$commonshort <- NA
stats_lester[stats_lester$pair=="verus - verus", ]$commonshort <- "Western"
stats_lester[stats_lester$pair=="ellioti - ellioti", ]$commonshort <- "Nigeria-Cam."
stats_lester[stats_lester$pair=="troglodytes - troglodytes", ]$commonshort <- "Central"
stats_lester[stats_lester$pair=="schweinfurthii - schweinfurthii", ]$commonshort <- "Eastern"
stats_lester[stats_lester$pair=="verus - ellioti", ]$commonshort <- "W - NC"
stats_lester[stats_lester$pair=="verus - troglodytes", ]$commonshort <- "W - C"
stats_lester[stats_lester$pair=="verus - schweinfurthii", ]$commonshort <- "W - E"
stats_lester[stats_lester$pair=="ellioti - troglodytes", ]$commonshort <- "NC - C"
stats_lester[stats_lester$pair=="ellioti - schweinfurthii", ]$commonshort <- "NC - E"
stats_lester[stats_lester$pair=="troglodytes - schweinfurthii", ]$commonshort <- "C - E"
stats_lester$commonshort <- factor(stats_lester$commonshort, levels = c("Western", "Nigeria-Cam.", "Central", "Eastern", "W - NC", 
                                                                "W - C", "W - E", 
                                                                "NC - C", "NC - E", "C - E"))

stats_lester$model <- "Lester et al. 2021"


# Import and format Fst from Fischer et al. 2006 -------------------------------

stats_fisher <- data.frame(matrix(nrow = 10, ncol = 0)) 
stats_fisher$commonshort <- c("Western", "Nigeria-Cam.", "Central", "Eastern", "W - NC", "W - C", "W - E", "NC - C", 
                       "NC - E", "C - E")
stats_fisher$Fst <- c(NA, NA, NA, NA, NA, 0.29, 0.31, NA, NA, 0.09)

stats1f <- stats_fisher
for (i in 1:20) {
  stats_fisher <- rbind(stats_fisher, stats1f)
}

stats_fisher$commonshort <- factor(stats_fisher$commonshort, levels = c("Western", "Nigeria-Cam.", "Central", "Eastern", "W - NC", 
                                                                "W - C", "W - E", 
                                                                "NC - C", "NC - E", "C - E"))

stats_fisher$model <- "Fischer et al. 2006"


# Import and format simulated Fst ----------------------------------------------

simu <- read.table("~/Stats_P_troglodytes_rec07_20-12-23_genmodel_rightellioti_Ranalysis.csv", 
                      sep = ",")

simu_Fst <- simu[str_detect(simu$V1, "Fst"),] 
simu_Fst$subsp1 <- str_split_i(simu_Fst$V1, "_", 2)
simu_Fst$subsp2 <- str_split_i(simu_Fst$V1, "_", 3)
simu_Fst$model <- as.factor(paste(str_split_i(simu_Fst$V3, coll("_"), 4), str_split_i(simu_Fst$V3, coll("_"), 5), str_split_i(simu_Fst$V3, coll("_"), 6), sep="_"))
simu_Fst <- simu_Fst[, c(2,4,5,6)]

simu_Fst$subsp1 <- gsub('[[:digit:]]+', '', simu_Fst$subsp1)
simu_Fst$subsp2 <- gsub('[[:digit:]]+', '', simu_Fst$subsp2)
simu_Fst$subsp1 <- gsub("troglo", 'troglodytes', simu_Fst$subsp1)
simu_Fst$subsp1 <- gsub("schwein", 'schweinfurthii', simu_Fst$subsp1)
simu_Fst$subsp2 <- gsub("troglo", 'troglodytes', simu_Fst$subsp2)
simu_Fst$subsp2 <- gsub("schwein", 'schweinfurthii', simu_Fst$subsp2)
simu_Fst$pair <- as.factor(paste(simu_Fst$subsp1, simu_Fst$subsp2, sep=" - "))

simu_Fst$cat <- NA
simu_Fst[simu_Fst$subsp1==simu_Fst$subsp2, ]$cat <- "within"
simu_Fst[simu_Fst$subsp1!=simu_Fst$subsp2, ]$cat <- "between"

simu_Fst$common <- NA
simu_Fst[simu_Fst$pair=="verus - verus", ]$common <- "Western"
simu_Fst[simu_Fst$pair=="ellioti - ellioti", ]$common <- "Nigeria-Cameroon"
simu_Fst[simu_Fst$pair=="troglodytes - troglodytes", ]$common <- "Central"
simu_Fst[simu_Fst$pair=="schweinfurthii - schweinfurthii", ]$common <- "Eastern"
simu_Fst[simu_Fst$pair=="verus - ellioti", ]$common <- "Western - Nigeria-Cameroon"
simu_Fst[simu_Fst$pair=="verus - troglodytes", ]$common <- "Western - Central"
simu_Fst[simu_Fst$pair=="verus - schweinfurthii", ]$common <- "Western - Eastern"
simu_Fst[simu_Fst$pair=="ellioti - troglodytes", ]$common <- "Nigeria-Cameroon - Central"
simu_Fst[simu_Fst$pair=="ellioti - schweinfurthii", ]$common <- "Nigeria-Cameroon - Eastern"
simu_Fst[simu_Fst$pair=="troglodytes - schweinfurthii", ]$common <- "Central - Eastern"
simu_Fst$common <- factor(simu_Fst$common, levels = c("Western", "Nigeria-Cameroon", "Central", "Eastern", "Western - Nigeria-Cameroon", 
                                                      "Western - Central", "Western - Eastern", 
                                                              "Nigeria-Cameroon - Central", "Nigeria-Cameroon - Eastern", "Central - Eastern"))


simu_Fst$commonshort <- NA
simu_Fst[simu_Fst$pair=="verus - verus", ]$commonshort <- "Western"
simu_Fst[simu_Fst$pair=="ellioti - ellioti", ]$commonshort <- "Nigeria-Cam."
simu_Fst[simu_Fst$pair=="troglodytes - troglodytes", ]$commonshort <- "Central"
simu_Fst[simu_Fst$pair=="schweinfurthii - schweinfurthii", ]$commonshort <- "Eastern"
simu_Fst[simu_Fst$pair=="verus - ellioti", ]$commonshort <- "W - NC"
simu_Fst[simu_Fst$pair=="verus - troglodytes", ]$commonshort <- "W - C"
simu_Fst[simu_Fst$pair=="verus - schweinfurthii", ]$commonshort <- "W - E"
simu_Fst[simu_Fst$pair=="ellioti - troglodytes", ]$commonshort <- "NC - C"
simu_Fst[simu_Fst$pair=="ellioti - schweinfurthii", ]$commonshort <- "NC - E"
simu_Fst[simu_Fst$pair=="troglodytes - schweinfurthii", ]$commonshort <- "C - E"
simu_Fst$commonshort <- factor(simu_Fst$commonshort, levels = c("Western", "Nigeria-Cam.", "Central", "Eastern", "W - NC", 
                                                      "W - C", "W - E", 
                                                      "NC - C", "NC - E", "C - E"))


colnames(simu_Fst)[1] <- "Fst"

stats1f <- simu_Fst
for (i in 1:20) {
  simu_Fst <- rbind(simu_Fst, stats1f)
}


# Plots (Figure ) -----------------------------------------------------------------

### Fst (main text, Figure 9)
binw=0.01

list_mod <- c("WNCE900000_WN800000_CE600000", "WNCE800000_WN600000_CE500000", "WNCE700000_WN500000_CE400000")
leg <- c("WNCE900000_WN800000_CE600000"="orange1", "WNCE800000_WN600000_CE500000"= "green3", "WNCE700000_WN500000_CE400000"="violetred2",
         "Lester et al. 2021"=alpha("dodgerblue2", 0.8), "Fischer et al. 2006"="dodgerblue4")
lab <- c(expression("T"["CEWN"]*" = 700, T"["CEWN"]*" = 500, T"["CEWN"]*" = 400 kya"), expression("T"["CEWN"]*" = 800, T"["CEWN"]*" = 600, T"["CEWN"]*" = 500 kya"), 
                    expression("T"["CEWN"]*" = 900, T"["CEWN"]*" = 800, T"["CEWN"]*" = 600 kya"))

output <-"PAPER_Fst_MAIN.png"

ggplot() +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count.., fill=model), geom = "col", binwidth=binw) +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count..), geom = "step", color='black', position=position_nudge(x=-0.5*binw), binwidth=binw) +
  geom_bar(data = stats_fisher, aes(x = Fst, fill=model), width = 0.006, color='white', alpha=0.7) +
  geom_bar(data = simu_Fst[simu_Fst$model %in% list_mod, ], aes(x = Fst, fill = model), width = 0.004, color=alpha('white', 0)) +
  facet_grid(commonshort ~ ., switch = "y") +
  labs(y="Number of comparisons") +
  scale_y_continuous(position = "right") +
  scale_fill_manual(name="Simulated scenario                                     Empirical data", values=leg, 
                    breaks=c("WNCE700000_WN500000_CE400000", "Lester et al. 2021",  "WNCE800000_WN600000_CE500000", "Fischer et al. 2006", 
                             "WNCE900000_WN800000_CE600000"),
                    labels = c(expression("T"["CENW"]*" = 700, T"["NW"]*" = 500, T"["CE"]*" = 400 kya    "),
                               "Lester et al. 2021", expression("T"["CENW"]*" = 800, T"["NW"]*" = 600, T"["CE"]*" = 500 kya    "),
                               "Fischer et al. 2006", expression("T"["CENW"]*" = 900, T"["NW"]*" = 800, T"["CE"]*" = 600 kya    "))) +
  scale_color_manual(name="Simulated scenario", values=leg) +
  scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), limits = c(-0.1, 0.6)) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE, title.position="top"))+
  theme_classic() +
  theme(strip.background    = element_blank(),
        axis.text.x         = element_text(size = 10, color='black'),
        axis.line.x         = element_blank(),
        axis.ticks.x        = element_line(colour = "gray90"),
        axis.ticks.length.x = unit(5, "points"),
        strip.text.y.left   = element_text(angle = 0, size = 12, hjust = 1, margin = margin(t = 2, r = 0, b = 2, l = 2)),
        panel.grid.major.x  = element_line(colour = "gray90"),
        legend.position = c(0.26, -0.18),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 2, b = 2, l = 2)),
        axis.title.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        plot.margin = margin(5, 5, 100, 5),
        legend.text = element_text(size = 11))

#ggsave(filename =  output, width = 20, height = 22, dpi = 300, units = "cm", device='png')




# Plots (supplementary) --------------------------------------------------------

### Several TWNCE times (supplementary, Figure S15)

binw=0.01

list_mod <- c("WNCE700000_WN600000_CE500000", "WNCE800000_WN600000_CE500000", "WNCE900000_WN600000_CE500000")
leg <- c("WNCE700000_WN600000_CE500000"="orange1", "WNCE800000_WN600000_CE500000"= "green3", "WNCE900000_WN600000_CE500000"="violetred2",
         "Lester et al. 2021"=alpha("dodgerblue2", 0.8), "Fischer et al. 2006"="dodgerblue4")
lab <- c(expression("T"["CEWN"]*" = 700 kya"), expression("T"["CEWN"]*" = 800 kya"), expression("T"["CEWN"]*" = 900 kya    (T"["WN"]*" = 600 kya and T"["CE"]*" = 500 kya)"))

output <-"PAPER_Fst_sevCEWN.png"

ggplot() +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count.., fill=model), geom = "col", binwidth=binw) +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count..), geom = "step", color='black', position=position_nudge(x=-0.5*binw), binwidth=binw) +
  geom_bar(data = stats_fisher, aes(x = Fst, fill=model), width = 0.006, color='white', alpha=0.7) +
  geom_bar(data = simu_Fst[simu_Fst$model %in% list_mod, ], aes(x = Fst, fill = model), width = 0.004, color=alpha('white', 0)) +
  facet_grid(commonshort ~ ., switch = "y") +
  labs(y="Number of comparisons") +
  scale_y_continuous(position = "right") +
  scale_fill_manual(name="Simulated scenario                                     Empirical data", values=leg, 
                    breaks=c("WNCE700000_WN600000_CE500000", "Lester et al. 2021",  "WNCE800000_WN600000_CE500000", "Fischer et al. 2006", 
                             "WNCE900000_WN600000_CE500000"),
                    labels = c(expression("T"["CENW"]*" = 700, T"["NW"]*" = 600, T"["CE"]*" = 500 kya    "),
                               "Lester et al. 2021", expression("T"["CENW"]*" = 800, T"["NW"]*" = 600, T"["CE"]*" = 500 kya    "),
                               "Fischer et al. 2006", expression("T"["CENW"]*" = 900, T"["NW"]*" = 600, T"["CE"]*" = 500 kya    "))) +
  scale_color_manual(name="Simulated scenario", values=leg) +
  scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), limits = c(-0.1, 0.6)) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE, title.position="top"))+
  theme_classic() +
  theme(strip.background    = element_blank(),
        axis.text.x         = element_text(size = 10, color='black'),
        axis.line.x         = element_blank(),
        axis.ticks.x        = element_line(colour = "gray90"),
        axis.ticks.length.x = unit(5, "points"),
        strip.text.y.left   = element_text(angle = 0, size = 12, hjust = 1, margin = margin(t = 2, r = 0, b = 2, l = 2)),
        panel.grid.major.x  = element_line(colour = "gray90"),
        legend.position = c(0.26, -0.18),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 2, b = 2, l = 2)),
        axis.title.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        plot.margin = margin(5, 5, 100, 5),
        legend.text = element_text(size = 11))

#ggsave(filename =  output, width = 20, height = 22, dpi = 300, units = "cm", device='png')


### Several TNW times (supplementary, Figure S16)

binw=0.01

list_mod <- c("WNCE900000_WN500000_CE400000", "WNCE900000_WN600000_CE400000", "WNCE900000_WN800000_CE400000")
leg <- c("WNCE900000_WN500000_CE400000"="orange1", "WNCE900000_WN600000_CE400000"= "green3", "WNCE900000_WN800000_CE400000"="violetred2",
         "Lester et al. 2021"=alpha("dodgerblue2", 0.8), "Fischer et al. 2006"="dodgerblue4")
lab <- c(expression("T"["WN"]*" = 500 kya"), expression("T"["WN"]*" = 600 kya"), expression("T"["WN"]*" = 800 kya    (T"["CEWN"]*" = 900 kya and T"["CE"]*" = 400 kya)"))

output <-"PAPER_Fst_sevWN.png"

ggplot() +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count.., fill=model), geom = "col", binwidth=binw) +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count..), geom = "step", color='black', position=position_nudge(x=-0.5*binw), binwidth=binw) +
  geom_bar(data = stats_fisher, aes(x = Fst, fill=model), width = 0.006, color='white', alpha=0.7) +
  geom_bar(data = simu_Fst[simu_Fst$model %in% list_mod, ], aes(x = Fst, fill = model), width = 0.004, color=alpha('white', 0)) +
  facet_grid(commonshort ~ ., switch = "y") +
  labs(y="Number of comparisons") +
  scale_y_continuous(position = "right") +
  scale_fill_manual(name="Simulated scenario                                     Empirical data", values=leg, 
                    breaks=c("WNCE900000_WN500000_CE400000", "Lester et al. 2021",  "WNCE900000_WN600000_CE400000", "Fischer et al. 2006", 
                             "WNCE900000_WN800000_CE400000"),
                    labels = c(expression("T"["CENW"]*" = 900, T"["NW"]*" = 500, T"["CE"]*" = 400 kya    "),
                               "Lester et al. 2021", expression("T"["CENW"]*" = 900, T"["NW"]*" = 600, T"["CE"]*" = 400 kya    "),
                               "Fischer et al. 2006", expression("T"["CENW"]*" = 900, T"["NW"]*" = 800, T"["CE"]*" = 400 kya    "))) +
  scale_color_manual(name="Simulated scenario", values=leg) +
  scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), limits = c(-0.1, 0.6)) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE, title.position="top"))+
  theme_classic() +
  theme(strip.background    = element_blank(),
        axis.text.x         = element_text(size = 10, color='black'),
        axis.line.x         = element_blank(),
        axis.ticks.x        = element_line(colour = "gray90"),
        axis.ticks.length.x = unit(5, "points"),
        strip.text.y.left   = element_text(angle = 0, size = 12, hjust = 1, margin = margin(t = 2, r = 0, b = 2, l = 2)),
        panel.grid.major.x  = element_line(colour = "gray90"),
        legend.position = c(0.26, -0.18),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 2, b = 2, l = 2)),
        axis.title.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        plot.margin = margin(5, 5, 100, 5),
        legend.text = element_text(size = 11))

#ggsave(filename =  output, width = 20, height = 22, dpi = 300, units = "cm", device='png')


### Several TCE times (supplementary, Figure S17)

binw=0.01

list_mod <- c("WNCE900000_WN800000_CE400000", "WNCE900000_WN800000_CE500000", "WNCE900000_WN800000_CE600000")
leg <- c("WNCE900000_WN800000_CE400000"="orange1", "WNCE900000_WN800000_CE500000"= "green3", "WNCE900000_WN800000_CE600000"="violetred2",
         "Lester et al. 2021"=alpha("dodgerblue2", 0.8), "Fischer et al. 2006"="dodgerblue4")
lab <- c(expression("T"["CE"]*" = 400 kya"), expression("T"["CE"]*" = 500 kya"), expression("T"["CE"]*" = 600 kya    (T"["CEWN"]*" = 900 kya and T"["WN"]*" = 800 kya)"))
output <-"PAPER_Fst_sevCE.png"

ggplot() +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count.., fill=model), geom = "col", binwidth=binw) +
  stat_bin(data = stats_lester, aes(x = Fst, y = ..count..), geom = "step", color='black', position=position_nudge(x=-0.5*binw), binwidth=binw) +
  geom_bar(data = stats_fisher, aes(x = Fst, fill=model), width = 0.006, color='white', alpha=0.7) +
  geom_bar(data = simu_Fst[simu_Fst$model %in% list_mod, ], aes(x = Fst, fill = model), width = 0.004, color=alpha('white', 0)) +
  facet_grid(commonshort ~ ., switch = "y") +
  labs(y="Number of comparisons") +
  scale_y_continuous(position = "right") +
  scale_fill_manual(name="Simulated scenario                                     Empirical data", values=leg, 
                    breaks=c("WNCE900000_WN800000_CE400000", "Lester et al. 2021",  "WNCE900000_WN800000_CE500000", "Fischer et al. 2006", 
                             "WNCE900000_WN800000_CE600000"),
                    labels = c(expression("T"["CENW"]*" = 900, T"["NW"]*" = 800, T"["CE"]*" = 400 kya    "),
                               "Lester et al. 2021", expression("T"["CENW"]*" = 900, T"["NW"]*" = 800, T"["CE"]*" = 500 kya    "),
                               "Fischer et al. 2006", expression("T"["CENW"]*" = 900, T"["NW"]*" = 800, T"["CE"]*" = 600 kya    "))) +
  scale_color_manual(name="Simulated scenario", values=leg) +
  scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), limits = c(-0.1, 0.6)) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE, title.position="top"))+
  theme_classic() +
  theme(strip.background    = element_blank(),
        axis.text.x         = element_text(size = 10, color='black'),
        axis.line.x         = element_blank(),
        axis.ticks.x        = element_line(colour = "gray90"),
        axis.ticks.length.x = unit(5, "points"),
        strip.text.y.left   = element_text(angle = 0, size = 12, hjust = 1, margin = margin(t = 2, r = 0, b = 2, l = 2)),
        panel.grid.major.x  = element_line(colour = "gray90"),
        legend.position = c(0.26, -0.18),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 2, b = 2, l = 2)),
        axis.title.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        plot.margin = margin(5, 5, 100, 5),
        legend.text = element_text(size = 11))

#ggsave(filename =  output, width = 20, height = 22, dpi = 300, units = "cm", device='png')
