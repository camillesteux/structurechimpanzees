# rtournebize
rm(list = ls(all = TRUE))
require(ggplot2)
require(reshape2)
require(ggdist)
require(viridis)
require(cowplot)

DIR <- "sims"
MODEL <- "chimp"

PopMap <- data.frame(
  from = c("Tw", "Tnc", "Te", "Tc", "P"),
  to = c("Western", "NC", "Eastern", "Central", "Bonobos"),
  stringsAsFactors = FALSE
)

lf <- list.files(DIR, full.names = FALSE)
lf <- lf[startsWith(lf, paste0(MODEL, ".ceu_"))]
lf <- lf[endsWith(lf, ".stats")]
lf <- paste0(DIR, "/", lf)

D <- PI <- FST <- c()
for (i in seq_along(lf)) {
  lab <- rev(strsplit(lf[i], "/")[[1]])[1]
  lab <- rev(strsplit(lab, ".", fixed = TRUE)[[1]])
  ceu <- gsub("ceu_", "", lab[3])
  yri <- gsub("yri_", "", lab[2])
  #
  x <- read.table(lf[i],
                  sep = " ",
                  fill = TRUE,
                  stringsAsFactors = FALSE)
  #
  y <- subset(x, V1 == "D_diploid")
  D <- rbind(
    D,
    data.frame(
      FILE = lf[i],
      CEU = ceu,
      YRI = yri,
      D = unlist(y[3]),
      SE = unlist(y[4]),
      stringsAsFactors = FALSE
    )
  )
  #
  y <- subset(x, V1 == "pi_CEU_per_kb")
  PI <- rbind(PI,
              data.frame(
                FILE = lf[i],
                CEU = ceu,
                PI = unlist(y[3]),
                stringsAsFactors = FALSE
              ))
  #
  y <- subset(x, V1 == "pi_YRI_per_kb")
  PI <- rbind(PI,
              data.frame(
                FILE = lf[i],
                CEU = yri,
                PI = unlist(y[3]),
                stringsAsFactors = FALSE
              ))
  #
  y <- subset(x, V1 == "Fst")
  FST <- rbind(FST,
               data.frame(
                 FILE = lf[i],
                 CEU = ceu,
                 YRI = yri,
                 FST = unlist(y[3]),
                 stringsAsFactors = FALSE
               ))
}

PI <- PI[!duplicated(PI$CEU), , drop = F]
PI$CEU <- plyr::mapvalues(PI$CEU, PopMap$from, PopMap$to, warn_missing =
                            FALSE)

FST$CEU <- plyr::mapvalues(FST$CEU, PopMap$from, PopMap$to, warn_missing =
                             FALSE)
FST$YRI <- plyr::mapvalues(FST$YRI, PopMap$from, PopMap$to, warn_missing =
                             FALSE)

# add the complements
REV.D <- D
REV.D[, c("CEU", "YRI")] <- REV.D[, c("YRI", "CEU")]
REV.D$D <- -1 * REV.D$D
D <- rbind(D, REV.D)
D$CEU <- plyr::mapvalues(D$CEU, PopMap$from, PopMap$to, warn_missing = FALSE)
D$YRI <- plyr::mapvalues(D$YRI, PopMap$from, PopMap$to, warn_missing = FALSE)
D$YRI.CEU <- apply(D[, c("YRI", "CEU")], 1, paste0, collapse = ".")
D$YRI.CEU <- factor(D$YRI.CEU, D$YRI.CEU[order(D$D)])

#===> OBS

# Human reference
OBS.FST <- rbind(
  c("Tc", "Te", 0.09),
  c("Tc", "Tw", 0.29),
  c("Te", "Tw", 0.32),
  c("Tc", "Tc", NA),
  c("Te", "Te", NA),
  c("Tw", "Tw", NA),
  c("Tnc", "Tnc", NA),
  c("Te", "Tnc", NA),
  c("Tnc", "Tw", NA),
  c("Tc", "Tnc", NA)
)
OBS.FST <- as.data.frame(OBS.FST, stringsAsFactors = FALSE)
names(OBS.FST) <- c("CEU", "YRI", "FST")
OBS.FST$FST <- as.numeric(OBS.FST$FST)
OBS.FST$CEU <- plyr::mapvalues(OBS.FST$CEU, PopMap$from, PopMap$to, warn_missing =
                                 FALSE)
OBS.FST$YRI <- plyr::mapvalues(OBS.FST$YRI, PopMap$from, PopMap$to, warn_missing =
                                 FALSE)

OBS.PI <- rbind(c("Te", 1.6), c("Tc", 2.0), c("Tw", 0.8))
OBS.PI <- as.data.frame(OBS.PI, stringsAsFactors = FALSE)
names(OBS.PI) <- c("CEU", "PI")
OBS.PI$PI <- as.numeric(OBS.PI$PI)
OBS.PI$CEU <- plyr::mapvalues(OBS.PI$CEU, PopMap$from, PopMap$to, warn_missing =
                                FALSE)

OBS.D <- rbind(
  c("Western.Central", 0.037, 0.007),
  c("Eastern.Central", 0.008, 0.004),
  c("NC.Eastern", 0.020, 0.006),
  c("NC.Central", 0.020, 0.005),
  c("Western.NC", 0.013, 0.008),
  c("Western.Eastern", 0.030, 0.007)
)
OBS.D <- as.data.frame(OBS.D, stringsAsFactors = FALSE)
names(OBS.D) <- c("YRI.CEU", "D", "SE")
OBS.D$D <- as.numeric(OBS.D$D)
OBS.D$SE <- as.numeric(OBS.D$SE)

#===> D

# YRI.CEU
FOCAL <- c(
  "Western.Central",
  "Western.Eastern",
  "NC.Central",
  "NC.Eastern",
  "Western.NC",
  "Eastern.Central"
)

DD <- subset(D, YRI.CEU %in% FOCAL)
DD$YRI.CEU <- factor(DD$YRI.CEU, levels = rev(FOCAL))
OBS.D$YRI.CEU <- factor(OBS.D$YRI.CEU, levels = rev(FOCAL))

g_D <- ggplot(DD, aes(y = YRI.CEU)) +
  geom_point(aes(x = D), size = 4.) +
  geom_vline(xintercept = 0., linetype = "dashed") +
  geom_text(aes(x = -Inf, label = YRI), hjust = 0, size = 5.) +
  geom_text(aes(x = Inf, label = CEU), hjust = 1, size = 5.) +
  xlim(c(-.06, 0.1)) +
  geom_errorbarh(aes(xmin = D - 1.96 * SE, xmax = D + 1.96 * SE),
                 height = 0,
                 size = 1.4) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 16.)
  ) +
  labs(x = "D", y = "") +
  geom_point(
    data = OBS.D,
    aes(D, YRI.CEU),
    shape = 10,
    col = "red",
    size = 6.
  ) +
  geom_errorbarh(
    data = OBS.D,
    aes(
      xmin = D - 1.96 * SE,
      xmax = D + 1.96 * SE,
      y = YRI.CEU
    ),
    col = "red",
    height = 0
  ) +
  ggtitle("(D) D-statistic")

#===> FST

diag <- data.frame(
  CEU = c("Western", "NC", "Eastern", "Central"),
  value = 0.,
  stringsAsFactors = FALSE
)
fst <- melt(FST,
            id.vars = c("CEU", "YRI"),
            measure.vars = "FST")
g_FST <- ggplot(fst) +
  geom_tile(aes(x = CEU, y = YRI, fill = value), color = "white") +
  geom_text(aes(
    x = CEU,
    y = YRI,
    label = round(value, 2)
  ), size = 6.) +
  scale_fill_viridis(begin = .7, end = .9) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 16.)
  ) +
  labs(x = "", y = "", fill = "Fst") +
  coord_equal() +
  geom_tile(data = OBS.FST,
            aes(CEU, YRI, fill = FST),
            color = "white") +
  geom_text(
    data = OBS.FST,
    aes(CEU, YRI, label = FST),
    col = "red",
    size = 6.
  ) +
  theme(legend.position = "none") +
  ggtitle(expression(paste("(C) ", F[ST]))) +
  geom_tile(
    data = diag,
    aes(x = CEU, y = CEU),
    fill = "white",
    color = "white"
  )

#===> PI

g_PI <- ggplot(PI) +
  geom_bar(
    aes(x = CEU, y = PI),
    stat = "identity",
    color = "#303030",
    fill = "#e0e0e0"
  ) +
  labs(x = "", y = expression(paste(pi, " (per kb)"))) +
  theme_bw() +
  theme(text = element_text(size = 16.), panel.grid = element_blank()) +
  geom_point(
    data = OBS.PI,
    aes(CEU, PI),
    shape = 10,
    col = "red",
    size = 6.
  ) +
  ggtitle(expression(paste("(B) ", pi)))

#===> ALL

model_plot <- paste0(DIR, "/", MODEL, ".logT.png")
m <- ggdraw() +
  draw_image(model_plot, scale = 1.2) +
  draw_figure_label("(A) Model", size = 20., position = "top.left")

G <- cowplot::plot_grid(m, g_PI, g_FST, g_D)
print(G)
ggsave(
  G,
  filename = paste0(DIR, "___", MODEL, ".pdf"),
  width = 10,
  height = 9
)

#___