### timing.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  1 2022 (11:15) 
## Version: 
## Last-Updated: feb  1 2022 (12:12) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Packages
library(data.table)
library(ggplot2)
library(ggthemes)

## * Import
## list.files("Results")
dt.timing <- readRDS(file = file.path("Results","simulation-timing.rds"))
vec.n <- sort(unique(dt.timing$n))

dt.timing[n %in% c(50,500) & a==0.0]

## * Figure
figureT <- ggplot(dt.timing,
                  aes(x = medianCor, y = time.median, group = as.factor(n), color = as.factor(n)))
figureT <- figureT + geom_line(size = 1.25) + geom_point(size = 3)
figureT <- figureT + facet_grid(~method)
figureT <- figureT + theme(legend.position="bottom")

figureT <- figureT + theme(text = element_text(size = 25),
                           legend.key.height = unit(0.05, "npc"),
                           legend.key.width = unit(0.08, "npc"),
                           axis.text.x = element_text(angle = 90, hjust = 1))
figureT <- figureT + labs(x = "correlation between the test statistics",
                          y = "computation time (s)")
figureT <- figureT + scale_x_continuous(limits = c(-0.001,1))
figureT <- figureT + scale_color_manual("sample size",
                                        breaks = vec.n,
                                        label =  vec.n,
                                        values = colorblind_pal()(n = length(vec.n)))
figureT <- figureT + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1))
figureT

## * Export
ggsave(figureT, filename = file.path("Figures","figure-timing.pdf"),
       height = 12, width = 15)

##----------------------------------------------------------------------
### timing.R ends here
