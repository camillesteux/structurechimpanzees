#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(stringr)

setwd("~/VALIDATION_figures")


# Param -----------------------------------------------------------------------

gentime <- 25
sub <- c("verus", "ellioti", "troglodytes", "schweinfurthii")
legend <- c("ellioti" = "lightblue4", "verus"="mediumpurple", "schweinfurthii"="coral2", "troglodytes"="goldenrod1")
legend_alpha <- c("ellioti" = alpha("lightblue4", 0.3), "verus"=alpha("mediumpurple", 0.3), 
                  "schweinfurthii"=alpha("coral2", 0.3), "troglodytes"=alpha("goldenrod1", 0.3))


# Import data ----------------------------------------------------------------

## Import target scenarios

scen_val <- list()
scen_val[["verus"]] <- list("n"=25, "N"=351.727, "ti"=c(0, 1560.89, 6333.65, 8476.03, 16302.6, 36261.9, 112773, 2e7/gentime),
                            "mi"=c(43.9363, 0.128708, 0.659483, 47.1583, 0.346383, 1.36904, 0.162036, 0.162036))
scen_val[["ellioti"]] <- list("n"=12, "N"=1161.93, "ti"=c(0, 548.442, 8889.64, 18111.9, 39813.9, 96375.2, 133365.0, 2e7/gentime),
                              "mi"=c(7.1041, 1.11224, 8.04395, 0.919609, 43.224, 1.29142, 0.232743, 0.232743))
scen_val[["troglodytes"]] <- list("n"=17, "N"=818.803, "ti"=c(0, 2225.83, 10806.7, 25338.3, 36991.1, 108549, 149296, 354840, 2e7/gentime),
                                  "mi"=c(43.7071, 0.205672, 0.269649, 0.731321, 40.4769, 0.678256, 0.230836, 1.43501, 1.43501))
scen_val[["schweinfurthii"]] <- list("n"=13, "N"=723.009, "ti"=c(0, 4660.82, 29970.6, 43309.8, 93549.9, 130056, 399407, 2e7/gentime),
                                     "mi"=c(49.6438, 0.207393, 0.64556, 2.84765, 0.734227, 0.152299, 0.395073, 0.395073))
for (s in sub) {
  scen_val[[s]][["df"]] <- setNames(as.data.frame(cbind(scen_val[[s]][["ti"]], scen_val[[s]][["mi"]])), c("T", "M"))
}


## Import inferred scenarios

scen_files <- list("verus"="./verus/valpsmc_verus_Donald_c7ebw050.csv", "ellioti"="./ellioti/valpsmc_ellioti_Koto_c7ebw050.csv", 
                   "troglodytes"="./troglodytes/valpsmc_troglo_Doris_c8ebw050.csv", "schweinfurthii"="./schweinfurthii/valpsmc_schwein_Kidongo_c7ebw050.csv")

## Format data

df_scen <- list()
df_scen2 <- list()

for (s in sub) {
  print(s)
  if (s=="schweinfurthii") {
    scen_snif <- read.csv(scen_files[[s]], header = TRUE, stringsAsFactors = TRUE, sep=",")
    scen_snif2 <- as.data.frame(matrix(nrow = 0, ncol = 5))
    for (k in c(1:nrow(scen_snif))) {
      time <- c(as.numeric(unlist(select(scen_snif[k,], starts_with("inf..t")))), 2e7/gentime)
      migr_rate <- as.numeric(unlist(select(scen_snif[k,], starts_with("inf..M"))))
      migr_rate <- c(migr_rate, migr_rate[length(migr_rate)])
      df <- as.data.frame(matrix(data=c(rep(s, length(time)), rep(as.character(scen_snif[k,"ind"]), length(time)), rep(scen_snif[k,"id"], length(time)), 
                                        time, migr_rate), ncol = 5))
      scen_snif2 <- rbind.data.frame(scen_snif2, setNames(df, c("sub", "ind", "rep", "time", "migr_rate")))
    }
    scen_snif2 <- setNames(scen_snif2, c("sub", "ind", "rep", "time", "migr_rate"))
    scen_snif2$rep <- as.factor(scen_snif2$rep)
    scen_snif2$ind <- as.factor(scen_snif2$ind)
    scen_snif2$migr_rate <- as.numeric(scen_snif2$migr_rate)
    scen_snif2$time <- as.numeric(scen_snif2$time)
    
    df_scen[[s]] <- scen_snif
    df_scen2[[s]] <- scen_snif2
  } else {
    scen_snif <- read.csv(scen_files[[s]], header = TRUE, stringsAsFactors = TRUE, sep="\t")
    scen_snif2 <- as.data.frame(matrix(nrow = 0, ncol = 5))
    for (k in c(1:nrow(scen_snif))) {
      time <- c(as.numeric(unlist(select(scen_snif[k,], starts_with("inf..t")))), 2e7/gentime)
      migr_rate <- as.numeric(unlist(select(scen_snif[k,], starts_with("inf..M"))))
      migr_rate <- c(migr_rate, migr_rate[length(migr_rate)])
      df <- as.data.frame(matrix(data=c(rep(s, length(time)), rep(as.character(scen_snif[k,"ind"]), length(time)), rep(scen_snif[k,"id"], length(time)), 
                                        time, migr_rate), ncol = 5))
      scen_snif2 <- rbind.data.frame(scen_snif2, setNames(df, c("sub", "ind", "rep", "time", "migr_rate")))
    }
    scen_snif2 <- setNames(scen_snif2, c("sub", "ind", "rep", "time", "migr_rate"))
    scen_snif2$rep <- as.factor(scen_snif2$rep)
    scen_snif2$ind <- as.factor(scen_snif2$ind)
    scen_snif2$migr_rate <- as.numeric(scen_snif2$migr_rate)
    scen_snif2$time <- as.numeric(scen_snif2$time)
  
    df_scen[[s]] <- scen_snif
    df_scen2[[s]] <- scen_snif2
  }
}

dataframe_n <- setNames(data.frame(matrix(nrow = 0, ncol = 3)), c("sub", "n", "val"))
dataframe_N <- setNames(data.frame(matrix(nrow = 0, ncol = 3)), c("sub", "N", "val"))
dataframe_cg <- setNames(data.frame(matrix(nrow = 0, ncol = 5)), c("sub", "ind", "rep", "time", "migr_rate"))

for (s in sub) {
  dataframe_n <- rbind.data.frame(dataframe_n, setNames(as.data.frame(cbind(rep(s, length(df_scen[[s]]$inf..n)), 
                            df_scen[[s]]$inf..n, rep(scen_val[[s]]$n, length(df_scen[[s]]$inf..n)))), c("sub", "n", "val")))
  dataframe_N <- rbind.data.frame(dataframe_N, setNames(as.data.frame(cbind(rep(s, length(df_scen[[s]]$inf..N_ref)), 
                            df_scen[[s]]$inf..N_ref, rep(scen_val[[s]]$N, length(df_scen[[s]]$inf..N_ref)))), c("sub", "N", "val")))
  dataframe_cg <- rbind.data.frame(dataframe_cg, df_scen2[[s]])
}

dataframe_n$sub <- as.factor(dataframe_n$sub)
dataframe_N$sub <- as.factor(dataframe_N$sub)
dataframe_n$n <- as.numeric(dataframe_n$n)
dataframe_N$N <- as.numeric(dataframe_N$N)
dataframe_n$val <- as.numeric(dataframe_n$val)
dataframe_N$val <- as.numeric(dataframe_N$val)
dataframe_cg$sub <- as.factor(dataframe_cg$sub)

dataframe_n$sub <- factor(dataframe_n$sub, level=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii'))
dataframe_N$sub <- factor(dataframe_N$sub, level=c('verus', 'ellioti', 'troglodytes', 'schweinfurthii'))



# Plots -----------------------------------------------------------------------

### n (supplementary, Figure S10A)
ggplot(data = dataframe_n) +
  geom_violin(aes(x=sub, y = n, fill=sub, color=sub), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  geom_point(aes(x=sub, y = n, fill=sub, color=sub), show.legend = FALSE, size=0.8) +
  geom_point(aes(x = sub, y=val), color = 'black', size = 1) +
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

#ggsave(filename =  "val_psmc_n.png", width = 10, height = 7, dpi = 300, units = "cm", device='png') 


### N (supplementary, Figure S10B)
ggplot(data = dataframe_N) +
  geom_violin(aes(x=sub, y = N, fill=sub, color=sub), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  geom_point(aes(x=sub, y = N, fill=sub, color=sub), show.legend = FALSE, size=0.8) +
  geom_point(aes(x = sub, y=val), color = 'black', size = 1) +
  scale_fill_manual(values=legend_alpha, name = "Supspecies", labels=c('W', 'NC', 'C', 'E')) +
  scale_colour_manual(values=legend, name = "Supspecies", labels=c('W', 'NC', 'C', 'E')) +
  labs(title="B", x = "Subspecies", y = expression('        Deme size\n (diploid individuals)')) +
  scale_y_continuous(n.breaks = 7, lim = c(0, 2000)) +
  scale_x_discrete(labels=c('W', 'NC', 'C', 'E')) +
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

#ggsave(filename =  "val_psmc_N.png", width = 10, height = 7, dpi = 300, units = "cm", device='png') 

### CG (supplementary, Figure S10C-F)
ggplot() +
  geom_step(data = dataframe_cg[dataframe_cg$sub %in% c("verus"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = sub), alpha = 0.20) +
  geom_step(data = scen_val[["verus"]]$df, aes(x=T*gentime, y=M), color = 'black') +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "Subspecies") +
  scale_color_manual(values=legend, guide = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  scale_x_continuous(limits = c(5e3, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=11, color="black"), 
        legend.text = element_text(size = 12), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = -2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = -2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "val_psmc_CG_W.png", width = 10, height = 6, dpi = 300, units = "cm", device='png') 

ggplot() +
  geom_step(data = dataframe_cg[dataframe_cg$sub %in% c("ellioti"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = sub), alpha = 0.20) +
  geom_step(data = scen_val[["ellioti"]]$df, aes(x=T*gentime, y=M), color = 'black') +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "Subspecies") +
  scale_color_manual(values=legend, guide = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  scale_x_continuous(limits = c(5e3, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_blank(), 
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=11, color="black"), 
        legend.text = element_text(size = 12), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = -2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = -2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "val_psmc_CG_NC.png", width = 10, height = 6, dpi = 300, units = "cm", device='png') 

ggplot() +
  geom_step(data = dataframe_cg[dataframe_cg$sub %in% c("troglodytes"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = sub), alpha = 0.20) +
  geom_step(data = scen_val[["troglodytes"]]$df, aes(x=T*gentime, y=M), color = 'black') +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "Subspecies") +
  scale_color_manual(values=legend, guide = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  scale_x_continuous(limits = c(5e3, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_text(size = 12, margin = margin(t = 7, r = 0, b = 0, l = 0)),
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=11, color="black"), 
        legend.text = element_text(size = 12), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = -2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = -2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "val_psmc_CG_C.png", width = 10, height = 7, dpi = 300, units = "cm", device='png') 

ggplot() +
  geom_step(data = dataframe_cg[dataframe_cg$sub %in% c("schweinfurthii"), ], aes(x = time*gentime, y = migr_rate, group = interaction(rep, ind), color = sub), alpha = 0.20) +
  geom_step(data = scen_val[["schweinfurthii"]]$df, aes(x=T*gentime, y=M), color = 'black') +
  labs(title="C", x = paste("Time in years (", gentime, " y/gen)", sep=""), y = "Migration rate", color = "Subspecies") +
  scale_color_manual(values=legend, guide = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  scale_x_continuous(limits = c(5e3, 2e7), trans = 'log10') +
  scale_y_continuous(trans = 'log10', n.breaks=6, limits = c(0.05, 50), breaks = c(0.1, 0.5, 1, 5, 10, 50)) +
  annotation_logticks(sides="bl") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 2), color="black"), 
        axis.title.x = element_text(size = 12, margin = margin(t = 7, r = 0, b = 0, l = 0)),
        axis.text.y=element_text(size=12, color="black"), 
        axis.text.x=element_text(size=11, color="black"), 
        legend.text = element_text(size = 12), legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = 2, r = -2, b = 2, l = 2),
        legend.margin = margin(t = 2, r = 2, b = 2, l = -2),
        panel.grid.major.y = element_line(colour = "white", size = 0.5))

#ggsave(filename =  "val_psmc_CG_E.png", width = 10, height = 7, dpi = 300, units = "cm", device='png') 

