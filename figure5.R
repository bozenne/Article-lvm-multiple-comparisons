## * Packages
library(data.table)
library(ggplot2)
library(ggthemes)

## * Import
## list.files("Results")
dtAll.figure5 <- readRDS(file = file.path("Results","simulation-figure5.rds"))
## sum(is.na(dtAll.figure5))
## unique(dtAll.figure5$method)

dt.diffBonf <- dtAll.figure5[method == "Wald.bonf"]
dt.diffBonf$method <- "Difference"
dt.diffBonf$power.selected <- dtAll.figure5[method == "Wald.max",power.selected] - dtAll.figure5[method == "Wald.bonf",power.selected]
dt.diffBonf$power <- dtAll.figure5[method == "Wald.max",power] - dtAll.figure5[method == "Wald.bonf",power]
## range(dt.diffBonf$power)

dt.wide.figure5 <- rbind(dtAll.figure5[method %in% c("Wald.bonf","Wald.max")],
                         dt.diffBonf
                         )



dt.figure5 <- melt(dt.wide.figure5, id.vars = c("n","a","method","n.rep","medianCor","selection"),
                   measure.vars = c("power","power.selected"), variable.name = "type", value.name = "rate")

dt.figure5[, a := as.factor(a)]
dt.figure5[, n := as.factor(n)]
dt.figure5[, type := factor(type, levels = c("power.selected","power"), labels = c("p[eta %~% Z]<=0.05","power"))]
       
nameStatistic.lvm <- c(Wald = "\"No adjustment\"",
                       Wald.bonf = "\"Bonferroni procedure\"",
                       Wald.max = "\"Max-test procedure\"",
                       Difference = "\"Difference\"")

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
ggsave(gg.power, filename = file.path("Figures","figure5-modelsearch-power.pdf"),
       height = 12, width = 15)



## dtTempo <- dtS.power[,.(method = "diff",
##              power.selected = .SD[method == "Wald.max",power.selected]-.SD[method == "Wald.bonf",power.selected],
##              power = .SD[method == "Wald.max",power]-.SD[method == "Wald.bonf",power]),
##           by = c("n","a","n.rep","medianCor","selection")]
## setcolorder(dtTempo, neworder = names(dtS.power))

## dtSS.power <- melt(rbind(dtS.power, dtTempo), id.vars = c("n","a","method","n.rep","medianCor","selection"),
##                    measure.vars = c("power","power.selected"), variable.name = "type", value.name = "rate")

## dtSS.power[, a := as.factor(a)]
## dtSS.power[, n := as.factor(n)]
## dtSS.power[, type := factor(type, levels = c("power.selected","power"), labels = c("p[eta %~% Z]<=0.05","power"))
