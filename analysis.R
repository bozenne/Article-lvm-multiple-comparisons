## * Setting
path.code <- file.path(path,"Code")
path.data <- file.path(path,"data")
path.source <- file.path(path,"source")

library(data.table)
library(Publish)
library(ggplot2)
library(nlme)
library(gdata)
library(butils)
library(lavaSearch2)
library(multcomp)
## library(xlsx)
library(xtable)

## * Data management 
df <- read.xls(file.path(path.source,"Input_MANOVA_primaryanalysis.xls"))
## df[,1:2]
dt <- as.data.table(df)
## dim(dt)

    ## ** Data management
    ## *** remove empty columns
    dt[,X:=NULL]
    dt[,X.1:=NULL]

    ## *** subset columns
    ## names(dt)
    old2new <- c("Subject" = "id",
                 "Genotype..1.HAB.0.MAB." = "genotype",
                 "Gender..1.M.0.F." = "gender",
                 "Type..1.concussion.0.healthy." = "group",
                 "Thalamus" = "thalamus",
                 "Pallidostriatum" = "pallidostriatum",
                 "Impact.Region.Neocortex" = "neocortex",
                 "Midbrain" = "midbrain",
                 "Pons" = "pons",
                 "Cingulate.Gyrus" = "cingulateGyrus",             
                 "Hippocampus...Parahippo" = "hippocampus",
                 "Supramarginal.Gyrus" = "supramarginalGyrus",
                 "Corpus.callosum" = "corpusCallosum"
                 )

    dt <- dt[,.SD,.SDcols = names(old2new)]

    ## *** rename columns
    setnames(dt, old = names(old2new), new = unname(old2new) )
    dt[,genotype := factor(genotype, levels = 0:1, labels = c("MAB","HAB"))]
    dt[,gender := factor(gender, levels = 0:1, labels = c("Female","Male"))]
    dt[,group := factor(group, levels = 0:1, labels = c("healthy","concussion"))]

    ## *** log transform
    regionScan <- c("thalamus",
                    "pallidostriatum",
                    "neocortex",
                    "midbrain",
                    "pons",
                    "cingulateGyrus",
                    "hippocampus",
                    "supramarginalGyrus",
                    "corpusCallosum")
    
    for(iR in c(regionScan)){
        dt[,c(paste0("log.",iR)) := log(.SD[[1]]), .SDcols = iR]
    }

    ## ** save
    if(save){
        saveRDS(dt, file = file.path(path.data, "dt.rds"))
    }

    ## ** export
    keep.col <- c("id","genotype","gender","group",paste0("log.",regionScan))
    return(dt[,.SD,.SDcols = keep.col])

## ls()
names(dt)
table(dt$genotype,dt$group)

## * initial model

## ** fit
m0.lvm <- lvm(c(log.thalamus, 
                log.pallidostriatum, 
                log.neocortex, 
                log.midbrain, 
                log.pons, 
                log.cingulateGyrus, 
                log.hippocampus, 
                log.supramarginalGyrus, 
                log.corpusCallosum)~genotype+group+eta)
latent(m0.lvm) <- ~eta

e0.lvm <- estimate(m0.lvm, data = dt)

## butils:::object2script(coef(e0.lvm), digit = 5)

## m00.lvm <- m0.lvm
## covariance(m00.lvm) <- log.supramarginalGyrus ~ log.neocortex
## covariance(m00.lvm) <- log.midbrain ~ log.pallidostriatum
## e00.lvm <- estimate(m00.lvm, data = dt, control = list(coef = coef(e0.lvm)))



## ** diagnostics
gof(e0.lvm)

## qqplot2(e0.lvm, mar = c(4,2,2,1.5))

## ## ** inference
## sCorrect(e0.lvm) <- TRUE
## ls.C <- createContrast(e0.lvm, var = "group")
## compare2(e0.lvm, contrast = ls.C$contrast, rhs = ls.C$null)

## * run modelsearch

## ** modelsearch
allNewLinks <- findNewLink(e0.lvm$model, type = "covariance")$link

system.time(
    e.search  <- modelsearch2(e0.lvm, link = allNewLinks, alpha = 0.15)
)
summary(e.search)

system.time(
    xx <- modelsearch(e0.lvm, type = "covariance")
)
## m1.lvm <- m0.lvm
## covariance(m1.lvm) <- log.supramarginalGyrus ~ log.neocortex
## e1.lvm <- estimate(m1.lvm, data = dt)
## summary2(e1.lvm)
## summary2(e1.lvm, robust = TRUE)

## ** table
tableS <- summary(e.search, print = FALSE)$table

table.search <- data.frame(matrix(NA, nrow = 3, ncol = 6))
names(table.search) <- c("covariance parameter","number of tests","max statistic","no adjustment","Bonferroni","Dunnett")

options(scipen=999)
table.search[["covariance parameter"]] <- paste0("(",gsub("\\~\\~",", ",gsub("log.","",tableS[["link"]])),")")
table.search[["number of tests"]] <- tableS[["nTests"]]
table.search[["max statistic"]] <- as.character(round(tableS[["statistic"]], 2))
table.search[["no adjustment"]] <- as.character(round(tableS[["p.value"]], digits = 6))
table.search[["Bonferroni"]] <- as.character(round(pmin(tableS[["p.value"]]*tableS[["nTests"]],1), digits = 6))
table.search[["Dunnett"]] <- as.character(round(tableS[["adjusted.p.value"]], digits = 6))

addtorow <- list(pos = list(0,0),
                 command = c("                    &                    &              & \\multicolumn{3}{c}{p-value} \\\\ \\cmidrule(r){4-6} \n",
                             paste0(paste0(names(table.search), collapse = " & "), "\\\\\n")
                             )
                 )
print(xtable(table.search, type = "latex", 
             label = "tab:Ex-modelsearchMax",
             caption = "Result of the first three steps of FSS for local dependencies. Each row corresponds to a step."),
 add.to.row = addtorow, include.rownames=FALSE, include.colnames = FALSE)

## ** comments
VarCor.stat <- e.search$sequenceSigma[[1]]
quantile(abs(VarCor.stat[upper.tri(VarCor.stat)]), probs = c(0,0.025,0.5,0.975,1))

fields::image.plot(VarCor.stat )

modelsearch(e0.lvm)

tableS[["adjusted.p.value"]]/tableS[["p.value"]]
    
## * model suggested by modelsearch
## ** fit
m1.lvm <- m0.lvm
covariance(m1.lvm) <- log.neocortex ~ log.supramarginalGyrus

m2.lvm <- m1.lvm
covariance(m2.lvm) <- log.pallidostriatum ~ log.midbrain

e1.lvm <- estimate(m1.lvm, data = dt)
e2.lvm <- estimate(m2.lvm, data = dt)

## ** gof
gof(e0.lvm)
gof(e1.lvm)
gof(e2.lvm)

qqplot2(e1.lvm, mar = c(4,2,2,1.5))
modelsearch(e1.lvm)

groupCoef <- grep("group",names(coef(e1.lvm)),value = TRUE)
iid.group <- cbind(id = dt$id, iid2(e1.lvm)[,groupCoef])
dt.inf <- melt(as.data.table(iid.group),
               id.vars = "id",
               variable.name = "region", value.name = "influence")
gg.iid <- ggplot(dt.inf, aes(x = influence))
gg.iid <- gg.iid + geom_histogram()
gg.iid <- gg.iid + facet_wrap(~region)

gg.iid



## * inference (with small sample correction)

## ** Test for no concussion effect across all regions

## chunk 20
outC <- createContrast(e1.lvm, var.test = "group")

## chunk 21
outConcussion <- compare2(e1.lvm, contrast =  outC$contrast, null = outC$null)
outConcussion[["estimate"]] <- outConcussion[["estimate"]][,1,drop=FALSE] ## hide unadjusted p-values
outConcussion

## ** Test for the concussion effect in each region separately

## chunk 22
outC <- createContrast(e1.lvm, var.test = "group")

## chunk 23
set.seed(10)
eScan.lvm.glht <- glht2(e1.lvm, linfct = outC$contrast, rhs = outC$null)
ls.pvalue <- list(none = summary(eScan.lvm.glht, test = adjusted("none"))$test,
                  Bonferroni = summary(eScan.lvm.glht, test = adjusted("bonferroni"))$test,
                  Dunnett = summary(eScan.lvm.glht, test = adjusted("single-step"))$test,
                  DunnettDown = summary(eScan.lvm.glht, test = adjusted("free"))$test
                  )

## ** table
table.multcomp <- data.frame(matrix(NA, nrow = 9, ncol = 7))
names(table.multcomp) <- c("region","concussion effect (\\%)","statistic","no adjustment","Bonferroni","Dunnett","step-down Dunnett")

table.multcomp[["region"]] <- gsub("log.|~groupconcussion","",names(ls.pvalue[[1]][["coefficients"]]))
table.multcomp[["concussion effect (\\%)"]] <- round(100*(exp(ls.pvalue[[1]][["coefficients"]]) - 1), 2)
table.multcomp[["statistic"]] <- round(ls.pvalue[[1]][["tstat"]],2)
table.multcomp[["no adjustment"]] <- round(ls.pvalue[["none"]][["pvalues"]],3)
table.multcomp[["Bonferroni"]] <- round(ls.pvalue[["Bonferroni"]][["pvalues"]],3)
table.multcomp[["Dunnett"]] <- round(ls.pvalue[["Dunnett"]][["pvalues"]],3)
table.multcomp[["step-down Dunnett"]] <- round(ls.pvalue[["DunnettDown"]][["pvalues"]],3)


addtorow <- list(pos = list(0,0),
                 command = c("      &                   &           &       \\multicolumn{4}{c}{p-value} \\\\ \\cmidrule(r){4-7} \n",
                             paste0(paste0(names(table.multcomp), collapse = " & "), "\\\\\n")
                             )
                 )
print(xtable(table.multcomp, type = "latex",
             label = "tab:Ex-multComp",
             caption = "Test of the null hypothesis of equal distribution volume in the healthy group vs. the mTBI group.
                        The column concussion indicates the percentage increase in distribution volume due to concussion.
"),
 add.to.row = addtorow, include.rownames=FALSE, include.colnames = FALSE)

## ** comments
Cov.multcomp <- eScan.lvm.glht$vcov[rownames(eScan.lvm.glht$linfct),rownames(eScan.lvm.glht$linfct)]
Cor.multcomp <- cov2cor(Cov.multcomp)
quantile(Cor.multcomp[upper.tri(Cor.multcomp)], probs = c(0,0.5,1))

fields::image.plot(cov2cor(VarCor.multcomp))

butils::ggHeatmap(cov2cor(VarCor.multcomp))

range(ls.pvalue[["Dunnett"]][["pvalues"]]/ls.pvalue[["none"]][["pvalues"]])
