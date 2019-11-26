## * Packages
library(data.table)
library(ggplot2)
library(ggthemes)

## * Import
dt.figure2 <- readRDS(file = file.path("Results","dt-coverageMultComp-article.rds"))

## * Figure 2
seqN <- levels(dt.figure2$n)
nN <- length(seqN)
    
figure2 <- ggplot(dt.figure2,
                  aes(x = medianCor, y = type1, group = n, color = n))
figure2 <- figure2 + geom_abline(intercept = 0.05, slope = 0, color = "red", size = 1.25)
figure2 <- figure2 + geom_line(size = 1.25) + geom_point(size = 3)
figure2 <- figure2 + facet_grid(ssc~approach)
figure2 <- figure2 + theme(legend.position="bottom")

figure2 <- figure2 + theme(text = element_text(size = 25),
                                   legend.key.height = unit(0.05, "npc"),
                                   legend.key.width = unit(0.08, "npc"),
                                   axis.text.x = element_text(angle = 90, hjust = 1))
figure2 <- figure2 + labs(x = "correlation between test statistics",
                                  y = "type 1 error (log scale)")
figure2 <- figure2 + scale_color_manual("sample size",
                                                breaks = seqN,
                                                label = seqN,
                                                values = colorblind_pal()(nN))
figure2 <- figure2 + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1))
figure2 <- figure2 + coord_trans(y = 'log10') + scale_y_continuous(breaks = c(0.01,0.025,0.05,0.1,0.2,0.3,0.4)) + theme(panel.grid.minor = element_blank())
figure2

## * Export
ggsave(figure2, filename = file.path("Figures","figure2-multcomp-type1error.pdf"), height = 12, width = 15)
