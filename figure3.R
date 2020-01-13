## * Packages
library(data.table)
library(ggplot2)
library(ggthemes)

## * Import
## list.files("Results")
dtAll.figure3 <- readRDS(file = file.path("Results","simulation-figure3.rds"))

dt.diffBonf <- dtAll.figure3[approach == "Bonferroni adjustment"]
dt.diffBonf$approach <- "Difference"
dt.diffBonf$power <- dtAll.figure3[approach == "Max-test adjustment",power] - dtAll.figure3[approach == "Bonferroni adjustment",power]

dt.figure3 <- rbind(dtAll.figure3[approach %in% c("Bonferroni adjustment","Max-test adjustment")],
                         dt.diffBonf
                         )

## * Figure 3
seqN <- levels(dt.figure3$n)
nN <- length(seqN)
df.line <- expand.grid(n = seqN,
                       approach = unique(dt.figure3$approach)
                       )
df.line$intercept <- as.numeric(NA)
df.line[df.line$approach == "Difference", "intercept"] <- 0
    
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

## * Figure for reviewers
dt.figureR1 <- readRDS(file = file.path("Results","simulation-figureR1.rds"))

## figure R1
dtminp.figure3[, labelCor := factor(a,
                                    levels = dtminp.figure3[,round(mean(medianCor),2),by = "a"][["a"]],
                                    labels = paste0("correlation=",dtminp.figure3[,round(mean(medianCor),2),by = "a"][[2]]))]
dtminp.figure3[, labelMethod := factor(method,
                                       levels = c("none2","Fisher2","Hochberg2","Hommel2","Bonferroni2","Max2"),
                                       labels = c("No adjustment","F-test","Hochberg adjustment","Hommel adjustment","Bonferroni adjustment","Max-test adjustment"))]

figureR1 <- ggplot(dtminp.figure3[method %in% c("Bonferroni2","Hochberg2","Hommel2","Fisher2","Max2") & n == 50], aes(y = min.p.value, x = labelMethod))
figureR1 <- figureR1 + geom_violin()
figureR1 <- figureR1 + facet_wrap(~labelCor)
figureR1 <- figureR1 + ylab("smallest p-value relative to the 9 group effects") + xlab("")
figureR1 <- figureR1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
figureR1
ggsave(figureR1, filename = file.path("Figures","figureR1-multcomp-power.pdf"))

## figure R2
figureR2 <- ggplot(dtAll.figure3[method %in% c("Bonferroni2","Hochberg2","Hommel2","Fisher2","Max2")],
                     aes(x = medianCor, y = power, group = n, color = n))
figureR2 <- figureR2 + geom_line(size = 1.25) + geom_point(size = 3)
figureR2 <- figureR2 + facet_grid(~approach)
figureR2 <- figureR2 + theme(legend.position="bottom")
figureR2 <- figureR2 + theme(text = element_text(size = 25),
                                             legend.key.height = unit(0.05, "npc"),
                                             legend.key.width = unit(0.08, "npc"),
                                             axis.text.x = element_text(angle = 90, hjust = 1))
figureR2 <- figureR2 + labs(x = "correlation between test statistics",
                                            y = "power")

figureR2 <- figureR2 + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1))

figureR2 <- figureR2 + scale_color_manual("sample size",
                                        breaks = seqN,
                                        label = seqN,
                                        values = colorblind_pal()(length(seqN)))
figureR2
ggsave(figureR2, filename = file.path("Figures","figureR2-multcomp-power.pdf"))




tableR1 <- dcast(dtAll.figure3[method %in% c("Fisher2","Bonferroni2","Hochberg2","Hommel2","Max2"),.(n,a,method,power)],
                 value.var = "power", formula = n+a~method)
range(tableR1$Hommel2 - tableR1$Bonferroni2)
range(tableR1$Hochberg2 - tableR1$Bonferroni2)
