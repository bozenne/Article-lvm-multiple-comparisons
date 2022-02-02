## * Packages
library(data.table)
source("FCT.R")

path.fig2 <- "Results/simulation-multComp-type1"
path.fig3 <- "Results/simulation-multComp-power"
path.fig4 <- "Results/simulation-modelsearch-type1"
path.fig5 <- "Results/simulation-modelsearch-power"
path.psi <- "Results/simulation-selection-type1"

## * Wald, type 1 error (figure 2)
## ** load file
## list.files(path.type1LVM)
dt.all.fig2 <- sinkDirectory(path.fig2, string.keep = "(tempo)")
## saveRDS(dt.type1, file = file.path(path,"Results","dt-coverageMultComp.rds"))
## dt.type1 <- readRDS(dt.type1, file = file.path(path,"Results","dt-coverageMultComp.rds"))

## ** compute type 1 error
dt.minp.fig2 <- dt.all.fig2[,.(n.link=.N,
                               medianCor=corTest2[1],
                               min.p.value = min(p.value, na.rm = TRUE)
                               ),
                            by = c("n","a","method","rep","iFile")]

dt.type1.fig2 <- dt.minp.fig2[,.(n.rep=.N,
                                 medianCor=median(medianCor),
                                 type1=mean( min.p.value <= 0.05, na.rm=TRUE)
                                 ),
                              by = c("n","a","method")]
dt.type1.fig2[, a := as.factor(a)]
dt.type1.fig2[, n := as.factor(n)]
dt.type1.fig2[, ssc := factor(method %in% c("none2","Bonferroni2","Hochberg2","Hommel2","Max2"),
                              levels = c(FALSE,TRUE),
                              labels = c("no small sample correction","small sample correction"))]

dt.type1.fig2[, approach := gsub("2","",method)]
dt.type1.fig2[, approach := factor(approach,
                                   levels = c("none","Bonferroni","Hochberg","Hommel","Max"),
                                   labels = c("No adjustment","Bonferroni adjustment","Hochberg adjustment","Hommel adjustment","Max-test adjustment"))]

## ** export
saveRDS(dt.type1.fig2, file = "Results/simulation-figure2.rds")

## * Wald, power (figure 3)
## ** load file
## list.files(path.type1LVM)
dt.all.fig3 <- sinkDirectory(path.fig3, string.keep = "(tempo)")

## ** compute power
dt.minp.fig3 <- dt.all.fig3[,.(n.link=.N,
                               medianCor=corTest2[1],
                               min.p.value = min(p.value, na.rm = TRUE)
                               ),
                            by = c("n","a","method","rep","iFile")]

dt.power.fig3 <- dt.minp.fig3[,.(n.rep=.N,
                                 medianCor=median(medianCor),
                                 power=mean( min.p.value <= 0.05, na.rm=TRUE)
                                 ),
                              by = c("n","a","method")]
dt.power.fig3[, a := as.factor(a)]
dt.power.fig3[, n := as.factor(n)]

dt.power.fig3[, approach := gsub("2","",method)]
dt.power.fig3[, approach := factor(approach,
                                   levels = c("none","Fisher","Hochberg","Hommel","Bonferroni","Max"),
                                   labels = c("No adjustment","F-test","Hochberg adjustment","Hommel adjustment","Bonferroni adjustment","Max-test adjustment"))]

## ** export
saveRDS(dt.power.fig3, file = "Results/simulation-figure3.rds")


## * Score, type 1 error (figure 4)
## ** load file
dt.all.fig4 <- sinkDirectory(path.fig4, string.keep = "(tempo)")
n.link <- length(dt.all.fig4$statistic[[1]])
dt.all.fig4[,Wald.sampbonf := pmin(1,Wald.samp * n.link)]
dt.all.fig4[,Wald.wildbonf := pmin(1,Wald.wild * n.link)]
## names(dt.lvm)

## ** compute type 1 error
dt.long.fig4 <- melt(dt.all.fig4, id.vars = c("n","a","rep","seed","iFile","medianCor"),
                     measure.vars = c("Wald","Wald.samp","Wald.wild",
                                      "Wald.bonf","Wald.sampbonf","Wald.wildbonf",
                                      "Wald.nummax","Wald.sampmax","Wald.wildmax"),
                     variable.name = "method",
                     value.name = "p.value")

dt.type1.fig4 <- dt.long.fig4[,.(n.rep=.N,
                                 medianCor=median(medianCor),
                                 type1=mean(p.value <= 0.05, na.rm=TRUE)
                                 ),
                              by = c("n","a","method")]


dt.type1.fig4[, a := as.factor(a)]
dt.type1.fig4[, n := as.factor(n)]

dt.long.timing <- melt(dt.all.fig4, id.vars = c("n","a","rep","seed","iFile","medianCor"),
                       measure.vars = c("tps","tps.samp","tps.wild"),
                       variable.name = "method",
                       value.name = "time")
dt.long.timing[, method := factor(method, c("tps", "tps.samp", "tps.wild"), c("approximation","sample","bootstrap"))]
dtS.long.timing <- dt.long.timing[,.(n.rep=.N,
                                     medianCor=median(medianCor),
                                     time.median=median(time,na.rm =TRUE),
                                     time.q05=quantile(time,0.05,na.rm =TRUE),
                                     time.q95=quantile(time,0.95,na.rm =TRUE)
                                     ),
                                  by = c("n","a","method")]

## ** export
saveRDS(dt.type1.fig4, file = "Results/simulation-figure4.rds")
saveRDS(dtS.long.timing, file = "Results/simulation-timing.rds")

## * Score, power (figure 5)
## ** load file
dt.all.fig5 <- sinkDirectory(path.fig5, string.keep = "(tempo)")
## dt.all.fig5 <- readRDS(file = "c:/Users/hpl802/Documents/Projects/LVM/modelsearch/article/results/dt-power-max.rds")
## sum(is.na(dt.all.fig5)) 

## ** compute power
dt.long.fig5 <- melt(dt.all.fig5,
                     id.vars = c("n","a","rep","seed","iFile","medianCor","selected"),
                     variable.name = "method",
                     value.name = "p.value")
## table(dtL.lvm$method)

dt.power.fig5 <- dt.long.fig5[!is.na(p.value),.(n.rep=.N,
                                                medianCor = median(medianCor),
                                                selection = mean(selected),
                                                power = mean(p.value <= 0.05, na.rm=FALSE),
                                                power.selected = mean( (p.value <= 0.05)*(selected), na.rm=FALSE)
                                                ),
                              by = c("n","a","method")]

## ** export
saveRDS(dt.power.fig5, file = "Results/simulation-figure5.rds")


## library(ggplot2)
## ggplot(dt.power.fig5, aes(x = medianCor, y = power, color = n, group = n)) + geom_line() + geom_point() + facet_wrap(~method)
## ggplot(dt.power.fig5, aes(x = medianCor, y = power.selected, color = n, group = n)) + geom_line() + geom_point() + facet_wrap(~method)

## * Post-selection inference (section 5.3)
list.files(path.psi)
dt.all.stat <- sinkDirectory(path.psi, string.keep = "statistics-S[[:digit:]]*\\.rds")
dt.all.type1 <- sinkDirectory(path.psi, string.keep = "type1error-S[[:digit:]]*\\.rds")
## list.files(path.psi)

## ** process
dt.minP <- dt.all.type1[method=="Max",.(n.link = .N,
                                         min.p.value = min(p.value)),
                                by = c("validModel","method","n","rep","iFile")]

## ** export
saveRDS(dt.minP[,
                .(n.rep = .N, valid.pc = mean(validModel), type1error = mean(min.p.value<=0.05)),
                by = c("method","n")], file = "Results/simulation-postSelection-pvalue.rds")


## additional
if(FALSE){
    dt.all.stat[validModel == TRUE, .(nrep = .N,
                                      cor1 = cor(.SD[["max"]],.SD[["log.thalamus~~log.pallidostriatum"]]),
                                      cor2 = cor(.SD[["max"]],.SD[["log.hippocampus~~log.corpusCallosum"]])
                                      ),
                by = "n"]
    ## n  nrep          cor1          cor2
    ## 1:  50 10000 -0.0083573877  0.0105799769
    ## 2: 100 10000 -0.0007908167 -0.0082296256
    ## 3: 500 10000  0.0129545210 -0.0008700901

    dtL.plotStat <- melt(dt.all.stat, id.vars = c("max","n","validModel"),
                         measure.vars = c("log.thalamus~~log.pallidostriatum","log.hippocampus~~log.corpusCallosum"),
                         value.name = "score",
                         variable.name = "link"
                         )

    gg.stat <- ggplot(dtL.plotStat[n==50], aes(x = max, y = score))
    gg.stat <- gg.stat + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
    gg.stat <- gg.stat + facet_wrap(~link)
    gg.stat <- gg.stat + xlab("test statistic for the group effect (after FSS)") + ylab("test statistic at the first step of the FSS") + labs(fill = "density")
    gg.stat <- gg.stat +  scale_fill_gradientn(colours = terrain.colors(15)[1:13])
    gg.stat
}

