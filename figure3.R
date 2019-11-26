## * Packages
library(data.table)
library(ggplot2)
library(ggthemes)

## * Import
## list.files("Results")
dtAll.figure3 <- readRDS(file = file.path("Results","dt-powerMultComp-article.rds"))

dt.diffBonf <- dtAll.figure3[approach == "Adjustment: Bonferroni"]
dt.diffBonf$approach <- "Difference"
dt.diffBonf$power <- dtAll.figure3[approach == "Adjustment: max-test",power] - dtAll.figure3[approach == "Adjustment: Bonferroni",power]

dt.figure3 <- rbind(dtAll.figure3[approach %in% c("Adjustment: Bonferroni","Adjustment: max-test")],
                         dt.diffBonf
                         )

## * Figure 3
seqN <- levels(dt.figure3$n)
nN <- length(seqN)
df.line <- expand.grid(n = seqN,
                        approach = unique(dt.figure3$approach)
                        )
df.line$intercept <- as.numeric(NA)
df.line[df.line$approach == "difference", "intercept"] <- 0
    
figure3 <- ggplot(dt.figure3,
                           aes(x = medianCor, y = power, group = n, color = n))
figure3 <- figure3 + geom_hline(data = df.line, aes(yintercept = intercept), color = "red", size = 2)
figure3 <- figure3 + geom_line(size = 1.25) + geom_point(size = 3)
figure3 <- figure3 + facet_grid(~approach)
figure3 <- figure3 + theme(legend.position="bottom")

figure3 <- figure3 + theme(text = element_text(size = 25),
                                             legend.key.height = unit(0.05, "npc"),
                                             legend.key.width = unit(0.08, "npc"),
                                             axis.text.x = element_text(angle = 90, hjust = 1))
figure3 <- figure3 + labs(x = "correlation between test statistics",
                                            y = "power")

figure3 <- figure3 + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1))

figure3 <- figure3 + scale_color_manual("sample size",
                                        breaks = seqN,
                                        label = seqN,
                                        values = colorblind_pal()(length(seqN)))
figure3

## * Export
ggsave(figure3, filename = file.path("Figures","figure3-multcomp-power.pdf"), height = 8, width = 15)
