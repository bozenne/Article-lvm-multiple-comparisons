## * Packages
library(data.table)
library(ggplot2)
library(ggthemes)

## * Import
## list.files("Results")
dt.figure4 <- readRDS(file = file.path("Results","dt-modelsearch-figure4.rds"))

nameStatistic.lvm <- c(Wald = "No adjustment",
                       Wald.samp = "No adjustment",
                       Wald.wild = "No adjustment",
                       Wald.bonf = "Adjustment: Bonferroni",
                       Wald.sampbonf = "Adjustment: Bonferroni",
                       Wald.wildbonf = "Adjustment: Bonferroni",
                       Wald.nummax = "Adjustment: max-test",
                       Wald.wildmax = "Adjustment: max-test",
                       Wald.sampmax = "Adjustment: max-test")
maxdist.lvm <- c(Wald = "approximation",
                 Wald.samp = "resampling",
                 Wald.wild = "wild-bootstrap",
                 Wald.bonf = "approximation",
                 Wald.sampbonf = "resampling",
                 Wald.wildbonf = "wild-bootstrap",
                 Wald.nummax = "approximation",
                 Wald.sampmax = "resampling",
                 Wald.wildmax = "wild-bootstrap")

nStatistic.lvm <- length(nameStatistic.lvm)

dt.figure4[, approach := factor(method,
                                levels = names(nameStatistic.lvm),
                                labels = as.character(nameStatistic.lvm))]
dt.figure4[, maxdist := factor(method,
                               levels = names(maxdist.lvm),
                               labels = as.character(maxdist.lvm))]


## * Figure 4
dt.figure4 <- dt.figure4[maxdist %in% c("approximation","resampling")]
vec.n <- sort(unique(dt.figure4$n))

figure4 <- ggplot(dt.figure4,
                      aes(x = medianCor, y = type1, group = n, color = n))
figure4 <- figure4 + geom_abline(intercept = 0.05, slope = 0, color = "red", size = 1.25)
figure4 <- figure4 + geom_line(size = 1.25) + geom_point(size = 3)
figure4 <- figure4 + facet_grid(maxdist~approach)
figure4 <- figure4 + theme(legend.position="bottom")

figure4 <- figure4 + theme(text = element_text(size = 25),
                                   legend.key.height = unit(0.05, "npc"),
                                   legend.key.width = unit(0.08, "npc"),
                                   axis.text.x = element_text(angle = 90, hjust = 1))
figure4 <- figure4 + labs(x = "correlation between the test statistics",
                                  y = "type 1 error (log scale)")
figure4 <- figure4 + scale_x_continuous(limits = c(-0.001,1))
figure4 <- figure4 + scale_color_manual("sample size",
                                                breaks = vec.n,
                                                label =  vec.n,
                                                values = colorblind_pal()(n = length(vec.n)))
figure4 <- figure4 + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1))
figure4 <- figure4 + coord_trans(y = 'log10') + scale_y_continuous(breaks = c(0.01,0.025,0.05,0.1,0.2,0.35,0.5)) + theme(panel.grid.minor = element_blank())

figure4

## * Export
ggsave(figure4, filename = file.path("Figures","figure4-modelsearch-type1error.pdf"),
       height = 12, width = 15)
