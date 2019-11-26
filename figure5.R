## * Packages
library(data.table)
library(ggplot2)
library(ggthemes)

## * Import
## list.files("Results")
dt.figure5 <- readRDS(file = file.path("Results","dt-modelsearch-figure5.rds"))

nameStatistic.lvm <- c(Wald = "\"No adjustment\"",
                       Wald.bonf = "\"Adjustment: Bonferroni\"",
                       Wald.max = "\"Adjustment: max-test\"",
                       diff = "\"Difference\"")

nStatistic.lvm <- length(nameStatistic.lvm)
dt.figure5[, approach := factor(method,
                                levels = names(nameStatistic.lvm),
                                labels = as.character(nameStatistic.lvm))]

## * Figure 5
vec.n <- sort(unique(dt.figure5$n))

gg.power <- ggplot(dt.figure5[approach != "\"No adjustment\""],
                   aes(x = medianCor, y = rate, group = n, color = n))
gg.power <- gg.power + geom_line(size = 1.25) + geom_point(size = 3)
gg.power <- gg.power + facet_grid(type~approach, labeller = label_parsed)
gg.power <- gg.power + theme(legend.position="bottom")

gg.power <- gg.power + theme(text = element_text(size = 25),
                             legend.key.height = unit(0.05, "npc"),
                             legend.key.width = unit(0.08, "npc"),
                             axis.text.x = element_text(angle = 90, hjust = 1))
gg.power <- gg.power + labs(x = "correlation between the test statistics",
                            y = "")
gg.power <- gg.power + scale_color_manual("sample size",
                                          breaks = vec.n,
                                          label =  vec.n,
                                          values = colorblind_pal()(n = length(vec.n)))
gg.power <- gg.power + scale_x_continuous(limits = c(-0.001,1))
gg.power <- gg.power + scale_y_continuous(labels = scales::percent, limits = c(-0.001,1.0001))
gg.power
## gg.power + coord_cartesian(ylim = c(0,0.05))

## * Export
ggsave(gg.power, filename = file.path("Results","figure5-modelsearch-power.pdf"),
       height = 12, width = 15)

