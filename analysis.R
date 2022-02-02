## * Packages
library(data.table)
library(Publish)
library(ggplot2)
library(nlme)
library(readxl)
library(lavaSearch2)
library(multcomp)
library(xtable)
library(microbenchmark)

## * Data management
## list.files("Source")
if(dir.exists("Source")){
    df <- read.csv(file.path("Source","illustration-real-pet-data.csv"), header = TRUE, sep = ";", dec = ",", na.string = "")
    dt <- as.data.table(df[rowSums(!is.na(df))>0,])
    dt[,X:=NULL]
    dt[,X.1:=NULL]
    
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
    setnames(dt, old = names(old2new), new = unname(old2new) )

}else{
    
    warning("The Danish rules on data protection do not allow to freely share person-sensitive health care related data containing any kind of person-identifying information. \n",
            "Simulated data will therefore be used instead of the real data \n",
            "Numbers/figures/tables in the article may not be exactly reproduced \n")
    dt <- readRDS(file.path("Results","illustration-simulated-pet-data.rds"))   
    
}

## ** rename columns
dt[,genotype := factor(genotype, levels = 0:1, labels = c("MAB","HAB"))]
dt[,gender := factor(gender, levels = 0:1, labels = c("Female","Male"))]
dt[,group := factor(group, levels = 0:1, labels = c("healthy","concussion"))]

## ** log transform
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

## names(dt)
## table(dt$genotype,dt$group)

## * statistical modeling

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

## ** diagnostics
gof(e0.lvm)
## qqplot2(e0.lvm, mar = c(4,2,2,1.5))

allNewLinks <- findNewLink(e0.lvm$model, type = "covariance")$link

set.seed(10)
system.time(
    e.search2  <- modelsearch2(e0.lvm, link = allNewLinks, alpha = 0.15)
)

## ** final fit
## methods(class = class(e.search2))
e2.lvm <- getNewModel(e.search2, step = 1)

## ** diagnostic final fit
gof(e2.lvm)
## qqplot2(e2.lvm, mar = c(4,2,2,1.5))

groupCoef <- grep("group",names(coef(e2.lvm)),value = TRUE)
iid.group <- cbind(id = dt$id, iid2(e2.lvm)[,groupCoef])
dt.inf <- melt(as.data.table(iid.group),
               id.vars = "id",
               variable.name = "region", value.name = "influence")
gg.iid <- ggplot(dt.inf, aes(x = influence))
gg.iid <- gg.iid + geom_histogram()
gg.iid <- gg.iid + facet_wrap(~region)

## ** small sample correction
if(packageVersion("lavaSearch2")>="2.0.0"){
    e2.lvm <- estimate2(e2.lvm)
}else{
    sCorrect(e2.lvm) <- TRUE
}
## * Results for section 6
## ** number of brain regions
length(endogenous(e2.lvm))

## ** chi2 test before modelsearch
gof(e0.lvm)

## ** number of possible links
length(allNewLinks)

## ** correlation between Wald statistics (modelsearch)
VarCor.stat <- e.search2$sequenceSigma[[1]]
quantile(abs(VarCor.stat[upper.tri(VarCor.stat)]), probs = c(0,0.025,0.5,0.975,1))

## ** "benefit" of the max test
M.tests <- do.call(rbind,lapply(e.search2$sequenceTest, function(iT){iT[NROW(iT),]}))
M.tests[["adjusted.p.value"]]/M.tests[["p.value"]]

## ** chi2 test after modelsearch
gof(e2.lvm)
## qqplot2(e2.lvm)

## ** multivariate test of the concussion effect
if(packageVersion("lavaSearch2")>="2.0.0"){
    outConcussion <- compare2(e2.lvm, linfct =  "groupconcussion")
    outConcussion2 <- compare2(e2.lvm, linfct =  "groupconcussion", as.lava = FALSE)
}else{
    outC <- createContrast(e2.lvm, var.test = "group")
    outConcussion <- compare2(e2.lvm, contrast =  outC$contrast, null = outC$null)
}
outConcussion[["estimate"]] <- outConcussion[["estimate"]][,1,drop=FALSE] ## hide unadjusted confidence intervals
outConcussion
## summary(outConcussion2) ## adjusted CI


## ** covariance between the Wald statistics of the concussion effect
eScan.lvm.glht <- glht2(e2.lvm, linfct = outC$contrast, rhs = outC$null)
Cov.multcomp <- eScan.lvm.glht$vcov[rownames(eScan.lvm.glht$linfct),rownames(eScan.lvm.glht$linfct)]
Cor.multcomp <- cov2cor(Cov.multcomp)
quantile(Cor.multcomp[upper.tri(Cor.multcomp)], probs = c(0,0.5,1))

## ** "benefit" of the max-test
set.seed(10)
ls.pvalue <- list(none = summary(eScan.lvm.glht, test = adjusted("none"))$test,
                  Bonferroni = summary(eScan.lvm.glht, test = adjusted("bonferroni"))$test,
                  Dunnett = summary(eScan.lvm.glht, test = adjusted("single-step"))$test,
                  DunnettDown = summary(eScan.lvm.glht, test = adjusted("free"))$test
                  )
## summary(eScan.lvm.glht, test = adjusted("bonferroni"))
## summary(eScan.lvm.glht, test = adjusted("Shaffer"))


range(ls.pvalue[["Dunnett"]][["pvalues"]]/ls.pvalue[["none"]][["pvalues"]])
min(ls.pvalue$Dunnett$pvalues)


## * Data splitting
if(FALSE){
    dtBin <- copy(dt)
    dtBin$group <- as.numeric(dtBin$group=="concussion")
    dtBin$genotype <- as.numeric(dtBin$genotype=="HAB")

    mBin.lvm <- m0.lvm
    covariance(mBin.lvm) <- log.neocortex ~ log.supramarginalGyrus

    vec.coef <- paste0("log.",c("thalamus", "pallidostriatum", "neocortex", "midbrain", "pons", "cingulateGyrus", "hippocampus", "supramarginalGyrus", "corpusCallosum"), "~group")

    ## ** usage A
    runAnalysisSplit <- function(p, trace = FALSE){ ## p <- 0.5
    
        iIndexTrain <- sample.int(n = NROW(dtBin), size = round(p*NROW(dt)))
        iIndexTest <- setdiff(1:NROW(dtBin), iIndexTrain)
        iDataTrain <- dtBin[iIndexTrain]
        iDataTest <- dtBin[iIndexTest]
    
        iE.lvm <- estimate(mBin.lvm, data = iDataTrain)
        indexTest <- which.min(summary(iE.lvm)$coef[vec.coef,"P-value"])
        
        iT.lvm <- estimate(mBin.lvm, data = iDataTest)
        return(summary(iE.lvm)$coef[vec.coef[indexTest],,drop=FALSE])

    }

    set.seed(10)
    ls.split05 <- pbapply::pblapply(1:100,function(i){runAnalysisSplit(0.5)})
    ls.split04 <- pbapply::pblapply(1:100,function(i){runAnalysisSplit(0.4)})
    ls.split03 <- pbapply::pblapply(1:100,function(i){runAnalysisSplit(0.3)})

    ## df.split <- data.table(coef = rownames(do.call(rbind,ls.split05)), p = 0.5, do.call(rbind,ls.split05))
    df.split <- rbind(data.table(coef = rownames(do.call(rbind,ls.split05)), p = 0.5, do.call(rbind,ls.split05)),
                      data.table(coef = rownames(do.call(rbind,ls.split04)), p = 0.4, do.call(rbind,ls.split04)),
                      data.table(coef = rownames(do.call(rbind,ls.split03)), p = 0.3, do.call(rbind,ls.split03)))

    setnames(df.split, old = "P-value", new = "p.value")
    ggSplit <- ggplot(df.split, aes(x = coef, y = p.value)) + facet_wrap(~p)
    ggSplit <- ggSplit + geom_boxplot()
    ggSplit

    ## ** usage B
    estimate.lvm <- lava:::estimate.lvm
    runAnalysisSplit <- function(p, trace = FALSE){
        iIndexTrain <- sample.int(n = NROW(dtBin), size = round(p*NROW(dt)))
        iIndexTest <- setdiff(1:NROW(dtBin), iIndexTrain)
        iDataTrain <- dtBin[iIndexTrain]
        iDataTest <- dtBin[iIndexTest]
    
        iE.lvm <- estimate(m0.lvm, data = iDataTrain)
        iE.search2  <- modelsearch2(iE.lvm, link = allNewLinks, alpha = 0.05, trace = trace)
        iSelected <- which(summary(iE.search2, print = FALSE)$table$selected)
        if(length(iSelected)>0){
            iE2.lvm <- update(getNewModel(iE.search2), data = dt[iIndexTest])
        }else{
            iE2.lvm <- iE.lvm
        }
        iC <- createContrast(iE2.lvm, var.test = "group")
        return(glht2(iE2.lvm, linfct = iC$contrast, rhs = iC$null))
    }

    set.seed(10)
    ls.split05 <- pbapply::pblapply(1:10,function(i){runAnalysisSplit(0.5)})
    ls.split04 <- pbapply::pblapply(1:10,function(i){runAnalysisSplit(0.4)})
    ls.split03 <- pbapply::pblapply(1:10,function(i){runAnalysisSplit(0.3)})


    M.pvalues05 <- do.call(rbind,lapply(ls.split05, function(iSplit){summary(iSplit)$test$pvalues}))
    M.pvalues05
    M.pvalues04 <- do.call(rbind,lapply(ls.split04, function(iSplit){summary(iSplit)$test$pvalues}))
    M.pvalues04
    M.pvalues03 <- do.call(rbind,lapply(ls.split03, function(iSplit){summary(iSplit)$test$pvalues}))
    M.pvalues03
}
