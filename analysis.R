## * Packages
library(data.table)
library(Publish)
library(ggplot2)
library(nlme)
library(gdata)
library(butils)
library(lavaSearch2)
library(multcomp)
library(xtable)

## * Data management
## list.files("Source")
if(dir.exists("Source")){
    ## file.exists(file.path("Source","Input_MANOVA_primaryanalysis.xls"))
    df <- read.xls(file.path("Source","illustration-real-pet-data.xls"))
    dt <- as.data.table(df)
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
    dt <- readRDS(file.path("Results","illustration-simulated-pet-data.xls"))   
    
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
## summary(e.search2)


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
sCorrect(e2.lvm) <- TRUE

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
outC <- createContrast(e2.lvm, var.test = "group")
outConcussion <- compare2(e2.lvm, contrast =  outC$contrast, null = outC$null)
outConcussion[["estimate"]] <- outConcussion[["estimate"]][,1,drop=FALSE] ## hide unadjusted p-values
outConcussion

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


## all(eScan.lvm.glht$vcov[rownames(eScan.lvm.glht$linfct),rownames(eScan.lvm.glht$linfct)]>0)
