#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(stringr)
library(scales)
library(forcats)
library(tidyr)

setwd("~/INFERENCES_figures")

## Param  ---------------------------------------------------------------------

gentime <- 25
legend <- c("ellioti" = "lightblue4", "verus"="mediumpurple", "schweinfurthii"="coral2", "troglodytes"="goldenrod1")
legend_alpha <- c("ellioti" = alpha("lightblue4", 0.3), "verus"=alpha("mediumpurple", 0.3), 
                  "schweinfurthii"=alpha("coral2", 0.3), "troglodytes"=alpha("goldenrod1", 0.3))

## Import data  ---------------------------------------------------------------

subsp1 <- read.csv("verus_all_c7ebcut.csv", header = TRUE, stringsAsFactors = TRUE, sep="\t")
subsp1$sp <- rep("verus", nrow(subsp1))
subsp2 <- read.csv("eliotti_inf_12-23.csv", header = TRUE, stringsAsFactors = TRUE, sep="\t")
subsp2$sp <- rep("ellioti", nrow(subsp2))
subsp3 <- read.csv("troglo_allind_c8eb_w050.csv", header = TRUE, stringsAsFactors = TRUE, sep="\t")
subsp3$sp <- rep("troglodytes", nrow(subsp3))
subsp4 <- read.csv("schwein_all_c7eb2_w050.csv", header = TRUE, stringsAsFactors = TRUE, sep="\t")
subsp4$sp <- rep("schweinfurthii", nrow(subsp4))


## Format data ----------------------------------------------------------------

### Western
subsp1_curves <- as.data.frame(matrix(nrow = 0, ncol = 6))
for (k in c(1:nrow(subsp1))) {
  time <- c(as.numeric(unlist(select(subsp1[k,], starts_with("inf..t")))), 2e7/gentime)
  migr_rate <- as.numeric(unlist(select(subsp1[k,], starts_with("inf..M"))))
  migr_rate <- c(migr_rate, migr_rate[length(migr_rate)])
  df <- as.data.frame(matrix(data=c(rep(as.character(subsp1[k,"ind"]), length(time)), rep(subsp1[k,"id"], length(time)),
                                    time, migr_rate, rep(as.character(subsp1[k,"sp"]), length(time))), ncol = 5))
  subsp1_curves <- rbind.data.frame(subsp1_curves, setNames(df, c("ind", "rep", "time", "migr_rate", "sp")))
}
subsp1_curves <- setNames(subsp1_curves, c("ind", "rep", "time", "migr_rate", "sp"))
subsp1_curves$rep <- as.factor(subsp1_curves$rep)
subsp1_curves$ind <- as.factor(subsp1_curves$ind)
subsp1_curves$sp <- as.factor(subsp1_curves$sp)
subsp1_curves$migr_rate <- as.numeric(subsp1_curves$migr_rate)
subsp1_curves$time <- as.numeric(subsp1_curves$time)

### Nigeria-Cameroon
subsp2_curves <- as.data.frame(matrix(nrow = 0, ncol = 6))
for (k in c(1:nrow(subsp2))) {
  time <- c(as.numeric(unlist(select(subsp2[k,], starts_with("inf..t")))), 2e7/gentime)
  migr_rate <- as.numeric(unlist(select(subsp2[k,], starts_with("inf..M"))))
  migr_rate <- c(migr_rate, migr_rate[length(migr_rate)])
  df <- as.data.frame(matrix(data=c(rep(as.character(subsp2[k,"ind"]), length(time)), rep(subsp2[k,"id"], length(time)),
                                    time, migr_rate, rep(as.character(subsp2[k,"sp"]), length(time))), ncol = 5))
  subsp2_curves <- rbind.data.frame(subsp2_curves, setNames(df, c("ind", "rep", "time", "migr_rate", "sp")))
}
subsp2_curves <- setNames(subsp2_curves, c("ind", "rep", "time", "migr_rate", "sp"))
subsp2_curves$rep <- as.factor(subsp2_curves$rep)
subsp2_curves$ind <- as.factor(subsp2_curves$ind)
subsp2_curves$migr_rate <- as.numeric(subsp2_curves$migr_rate)
subsp2_curves$time <- as.numeric(subsp2_curves$time)

### Central
subsp3_curves <- as.data.frame(matrix(nrow = 0, ncol = 6))
for (k in c(1:nrow(subsp3))) {
  time <- c(as.numeric(unlist(select(subsp3[k,], starts_with("inf..t")))), 2e7/gentime)
  migr_rate <- as.numeric(unlist(select(subsp3[k,], starts_with("inf..M"))))
  migr_rate <- c(migr_rate, migr_rate[length(migr_rate)])
  df <- as.data.frame(matrix(data=c(rep(as.character(subsp3[k,"ind"]), length(time)), rep(subsp3[k,"id"], length(time)),
                                    time, migr_rate, rep(as.character(subsp3[k,"sp"]), length(time))), ncol = 5))
  subsp3_curves <- rbind.data.frame(subsp3_curves, setNames(df, c("ind", "rep", "time", "migr_rate", "sp")))
}
subsp3_curves <- setNames(subsp3_curves, c("ind", "rep", "time", "migr_rate", "sp"))
subsp3_curves$rep <- as.factor(subsp3_curves$rep)
subsp3_curves$ind <- as.factor(subsp3_curves$ind)
subsp3_curves$sp <- as.factor(subsp3_curves$sp)
subsp3_curves$migr_rate <- as.numeric(subsp3_curves$migr_rate)
subsp3_curves$time <- as.numeric(subsp3_curves$time)

### Eastern
subsp4_curves <- as.data.frame(matrix(nrow = 0, ncol = 6))
for (k in c(1:nrow(subsp4))) {
  time <- c(as.numeric(unlist(select(subsp4[k,], starts_with("inf..t")))), 2e7/gentime)
  migr_rate <- as.numeric(unlist(select(subsp4[k,], starts_with("inf..M"))))
  migr_rate <- c(migr_rate, migr_rate[length(migr_rate)])
  df <- as.data.frame(matrix(data=c(rep(as.character(subsp4[k,"ind"]), length(time)), rep(subsp4[k,"id"], length(time)),
                                    time, migr_rate, rep(as.character(subsp4[k,"sp"]), length(time))), ncol = 5))
  subsp4_curves <- rbind.data.frame(subsp4_curves, setNames(df, c("ind", "rep", "time", "migr_rate", "sp")))
}
subsp4_curves <- setNames(subsp4_curves, c("ind", "rep", "time", "migr_rate", "sp"))
subsp4_curves$rep <- as.factor(subsp4_curves$rep)
subsp4_curves$ind <- as.factor(subsp4_curves$ind)
subsp4_curves$migr_rate <- as.numeric(subsp4_curves$migr_rate)
subsp4_curves$time <- as.numeric(subsp4_curves$time)

### Bind subspecies
bindsubsp <- rbind(subsp1[,c('sp', 'ind', 'inf..n', 'inf..N_ref')], subsp2[,c('sp', 'ind', 'inf..n', 'inf..N_ref')],
                   subsp3[,c('sp', 'ind', 'inf..n', 'inf..N_ref')], subsp4[,c('sp', 'ind', 'inf..n', 'inf..N_ref')])
bindsubsp_curves <- rbind(subsp1_curves, subsp2_curves, subsp3_curves, subsp4_curves)

sub_order <- factor(bindsubsp$sp, level=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii'))
bindsubsp$sp <- factor(bindsubsp$sp, level=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii'))
bindsubsp_curves$sp <- factor(bindsubsp_curves$sp, level=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii'))

ind <- c()
for (i in str_split(bindsubsp$ind, "_")) {
  ind <- c(ind, i[length(i)])
}
ind[ind=='ind'] <- 'Kidongo'
bindsubsp$ind <- ind
bindsubsp[bindsubsp$sp=="troglodytes" & bindsubsp$ind=="Julie",]$ind <- rep("Julie2", 10)
bindsubsp$ind <- factor(bindsubsp$ind, levels = c("Bosco", "Donald", "Jimmie", "Clint", "Koby", "Akwaya-Jean", "Damian", "Julie", "Koto", "Taweh", 
                                                  "Vaillant", "Doris", "Julie2", "Clara", "Bwambale", "Kidongo", "Nakuu"))

ind_curves <- c()
for (i in str_split(bindsubsp_curves$ind, "_")) {
  ind_curves <- c(ind_curves, i[length(i)])
}
ind_curves[ind_curves=='ind'] <- 'Kidongo'
bindsubsp_curves$ind <- ind_curves
bindsubsp_curves[bindsubsp_curves$sp=="troglodytes" & bindsubsp_curves$ind=="Julie",]$ind <- rep("Julie2", 90)
bindsubsp_curves$ind <- factor(bindsubsp_curves$ind, levels = c("Bosco", "Donald", "Jimmie", "Clint", "Koby", "Akwaya-Jean", "Damian", "Julie", "Koto", "Taweh", 
                                                  "Vaillant", "Doris", "Julie2", "Clara", "Bwambale", "Kidongo", "Nakuu"))

bindsubsp$inf..n <- as.integer(bindsubsp$inf..n)
bindsubsp$inf..N_ref <- as.numeric(bindsubsp$inf..N_ref)


# Computing statistics on inferred parameters ---------------------------------

### Number of islands n
quant_n <- bindsubsp %>% group_by(sp) %>% 
  summarise(quantile(inf..n),
            .groups = 'drop') %>%
  as.data.frame()
quant_n$stat <- as.factor(rep(c("min", "1q", "med", "3q", "max"), 4))
colnames(quant_n) <- c("sp", "val", "stat")
quant_n <- quant_n %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  as.data.frame()

df <- bindsubsp %>% group_by(sp) %>% 
  summarise(mean=mean(inf..n),
            .groups = 'drop') %>%
  as.data.frame

quant_n <-cbind(quant_n, df$mean)
colnames(quant_n)[7] <- c("mean")

### Deme size N
quant_N <- bindsubsp %>% group_by(sp) %>% 
  summarise(quantile(inf..N_ref),
            .groups = 'drop') %>%
  as.data.frame()

quant_N$stat <- as.factor(rep(c("min", "1q", "med", "3q", "max"), 4))
colnames(quant_N) <- c("sp", "val", "stat")
quant_N <- quant_N %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  as.data.frame()

df <- bindsubsp %>% group_by(sp) %>% 
  summarise(mean=mean(inf..N_ref),
            .groups = 'drop') %>%
  as.data.frame

quant_N <-cbind(quant_N, df$mean)
colnames(quant_N)[7] <- c("mean")

quant_n
quant_N

#write.csv(quant_n[, c(1, 2, 3, 4, 7, 5, 6)], file="quantiles_n.csv", row.names = FALSE, quote=FALSE)
#write.csv(quant_N[, c(1, 2, 3, 4, 7, 5, 6)], file="quantiles_N.csv", row.names = FALSE, quote=FALSE)


# Plots -----------------------------------------------------------------------

### n (main text, Figure 2A)
ggplot(data = bindsubsp, aes(x=sp, y = inf..n, fill=sp, color=sp)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  geom_point(show.legend = FALSE, size=0.8) +
  scale_fill_manual(values=legend_alpha, name = "Supspecies", labels=c('W', 'NC', 'C', 'E')) +
  scale_colour_manual(values=legend, name = "Supspecies", labels=c('W', 'NC', 'C', 'E')) +
  labs(title="A", x = "Subspecies", y = "Number of islands") +
  scale_y_continuous(n.breaks = 7, lim = c(0, 60)) +
  scale_x_discrete(labels=c('W', 'NC', 'C', 'E')) +
  theme_classic() +
  theme(legend.position = "None",
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 10), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=11, color="black"), 
        legend.text = element_text(size = 12), legend.title = element_text(size=30),
        plot.title = element_blank(),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredn.png", width = 10, height = 7, dpi = 300, units = "cm", device='png') 

### N (main text, Figure 2B)
ggplot(data = bindsubsp, aes(x=sp, y = inf..N_ref, fill=sp, color=sp)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  geom_point(show.legend = FALSE, size=0.8) +
  scale_fill_manual(values=legend_alpha, name = "Supspecies", labels=c('W', 'NC', 'C', 'E')) +
  scale_colour_manual(values=legend, name = "Supspecies", labels=c('W', 'NC', 'C', 'E')) +
  labs(title="B", x = "Subspecies", y = expression('        Deme size\n (diploid individuals)')) +
  scale_y_continuous(n.breaks = 7, lim = c(0, 2200)) +
  scale_x_discrete(labels=c('W', 'NC', 'C', 'E')) +
  theme_classic() +
  theme(legend.position = "None",
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 10), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=11, color="black"), 
        legend.text = element_text(size = 12), legend.title = element_text(size=30),
        plot.title = element_blank(),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))


#ggsave(filename =  "PAPER_inferredN.png", width = 10, height = 7, dpi = 300, units = "cm", device='png') 


### CG Western and Nigeria-Cam (main text, Figure 3A)
ggplot(data= bindsubsp_curves[bindsubsp_curves$sp %in% c("verus", "ellioti"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = sp)) +
  geom_step(alpha = 0.30) +
  geom_vline(xintercept = 365000, color = alpha("#6699FF", 0.2), size=9) +
  geom_vline(xintercept = 490000, color = alpha("#FF9900", 0.2), size=7.5) +
  geom_vline(xintercept = 2700000, color = alpha("#339933", 0.8), size=1) +
  geom_vline(xintercept = 1000000, color = alpha("#66CC33", 0.8), size=1) +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "Subspecies") +
  scale_color_manual(labels = c("Western", "Nigeria-Cam."), values=legend, guide = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  scale_x_continuous(limits = c(5e3, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'right',
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=10, color="black"), 
        axis.text.x=element_text(size=9, color="black"), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = -2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

ggsave(filename =  "PAPER_inferredCG_VE2.png", width = 16, height = 6, dpi = 300, units = "cm", device='png') 

### CG Central and Eastern (main text, Figure 3B)
ggplot() +
  geom_step(data= bindsubsp_curves[bindsubsp_curves$sp %in% c("troglodytes", "schweinfurthii"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = sp), alpha = 0.20) +
  geom_vline(xintercept = 155000, color = alpha("#FF99CC", 0.2), size=17) +
  geom_vline(xintercept = 490000, color = alpha("#FF9900", 0.2), size=7) +
  geom_vline(xintercept = 2700000, color = alpha("#339933", 0.8), size=1) +
  geom_vline(xintercept = 1000000, color = alpha("#66CC33", 0.8), size=1) +
  labs(title="D", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "Subspecies") +
  scale_color_manual(labels = c("Central        ", "Eastern        "), values=legend, guide = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  scale_x_continuous(limits = c(5e3, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'right',
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=10, color="black"), 
        axis.text.x=element_text(size=9, color="black"), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = -2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredCG_CE.png", width = 16, height = 6, dpi = 300, units = "cm", device='png') 


# Plots colored by individuals ------------------------------------------------

### n (supplementary, Figure S7A)
ggplot(data = bindsubsp, aes(x=ind, y = inf..n, fill=sp, color=sp)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE, size = 0.25) +
  geom_point(show.legend = FALSE, size=0.3) +
  scale_fill_manual(values=legend_alpha, name = "Supspecies", labels=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii')) +
  scale_colour_manual(values=legend, name = "Supspecies", labels=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii')) +
  labs(title="A", x = "Subspecies", y = "Number of islands") +
  scale_y_continuous(n.breaks = 7, lim = c(0, 60)) +
  theme_classic() +
  theme(legend.position = "None",
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 10), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=11, color="black", angle = 50, vjust = 1, hjust=1), 
        legend.text = element_text(size = 12), legend.title = element_text(size=30),
        plot.title = element_blank(),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredn_perind.png", width = 13, height = 8, dpi = 300, units = "cm", device='png') 


### N (supplementary, Figure S7B)
ggplot(data = bindsubsp, aes(x=ind, y = inf..N_ref, fill=sp, color=sp)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE,  size = 0.25) +
  geom_point(show.legend = FALSE, size=0.3) +
  scale_fill_manual(values=legend_alpha, name = "Supspecies", labels=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii')) +
  scale_colour_manual(values=legend, name = "Supspecies", labels=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii')) +
  labs(title="B", x = "Subspecies", y = "Deme size") +
  scale_y_continuous(n.breaks = 7, lim = c(0, 1200)) +
  theme_classic() +
  theme(legend.position = "None",
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 10), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=10, color="black", angle = 50, vjust = 1, hjust=1), 
        legend.text = element_text(size = 12), legend.title = element_text(size=30),
        plot.title = element_blank(),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredN_perind.png", width = 13, height = 8, dpi = 300, units = "cm", device='png') 


### CG (supplementary, Figure S8A-D)

cbPalette <- c("#0072B2", "#999999", "#E69F10", "#FF689F", "#009E50")

ggplot() +
  geom_step(data= bindsubsp_curves[bindsubsp_curves$sp %in% c("verus"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = ind), alpha = 0.450, size = 0.25) +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "") +
  scale_color_manual(values=cbPalette) +
  scale_x_continuous(limits = c(1e4, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=10, color="black"), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredCG_V_perind.png", width = 13, height = 7, dpi = 300, units = "cm", device='png') 


ggplot() +
  geom_step(data= bindsubsp_curves[bindsubsp_curves$sp %in% c("ellioti"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = ind), alpha = 0.450, size = 0.25) +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "") +
  scale_color_manual(values=cbPalette) +
  scale_x_continuous(limits = c(1e4, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=10, color="black"), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 4, b = 2, l = 2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredCG_NC_perind2.png", width = 13, height = 7, dpi = 300, units = "cm", device='png') 

ggplot() +
  geom_step(data= bindsubsp_curves[bindsubsp_curves$sp %in% c("troglodytes"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = ind), alpha = 0.450, size = 0.25) +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "") +
  scale_color_manual(values=cbPalette) +
  scale_x_continuous(limits = c(1e4, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=10, color="black"), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredCG_C_perind.png", width = 13, height = 7, dpi = 300, units = "cm", device='png') 

ggplot() +
  geom_step(data= bindsubsp_curves[bindsubsp_curves$sp %in% c("schweinfurthii"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = ind), alpha = 0.450, size = 0.25) +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "") +
  scale_color_manual(values=cbPalette) +
  scale_x_continuous(limits = c(1e4, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=10, color="black"), 
        legend.text = element_text(size = 10), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "PAPER_inferredCG_E_perind.png", width = 13, height = 7, dpi = 300, units = "cm", device='png') 


