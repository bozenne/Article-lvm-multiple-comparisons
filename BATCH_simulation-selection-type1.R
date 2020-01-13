## cd Cluster/LVMproject/article-multipleComparisons/
## path <- "P:/Cluster/LVMproject/article-multipleComparisons"
## setwd(path)
## source("BATCH_simulation-selection-type1.R")

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 2}
if(is.na(n.iter_sim)){n.iter_sim <- 20}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n", sep = "")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","simulation-selection-type1")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-selection-type1")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
library(multcomp)
library(lavaSearch2)

## * settings 
seqN <- c(50,100,500)
n.rep <- 200
n.repMax <- 2500

## * model
## ** generative model
m.sim <- lvm(c(log.thalamus, 
               log.pallidostriatum, 
               log.neocortex, 
               log.midbrain, 
               log.pons, 
               log.cingulateGyrus, 
               log.hippocampus, 
               log.supramarginalGyrus, 
               log.corpusCallosum)~ eta + genotypeHAB, groupconcussion ~ 1)
latent(m.sim) <- ~eta
covariance(m.sim) <- log.thalamus~log.pallidostriatum
## covariance(m.sim) <- log.thalamus~log.neocortex
## covariance(m.sim) <- log.neocortex~log.pons
## covariance(m.sim) <- log.midbrain~log.cingulateGyrus
covariance(m.sim) <- log.hippocampus~log.corpusCallosum
## covariance(m.sim) <- log.hippocampus~log.corpusCallosum
name.covcoef <- names(coefType(m.sim)[coefType(m.sim)=="covariance"])
sim.coef <- setNames(rep(0.005, length(name.covcoef)), name.covcoef)

## ** investigator model
m.test <- lvm(c(log.thalamus, 
                log.pallidostriatum, 
                log.neocortex, 
                log.midbrain, 
                log.pons, 
                log.cingulateGyrus, 
                log.hippocampus, 
                log.supramarginalGyrus, 
                log.corpusCallosum)~ eta + genotypeHAB + groupconcussion)
latent(m.test) <- ~eta

## * prepare
n.N <- length(seqN)

## ** value for the simulation
sim.coef <- c(sim.coef,
              "log.pallidostriatum" = 0.07255,
              "log.neocortex" = -0.08699,
              "log.midbrain" = 0.25676,
              "log.pons" = 0.28991,
              "log.cingulateGyrus" = 0.09924,
              "log.hippocampus" = 0.09823,
              "log.supramarginalGyrus" = -0.1254,
              "log.corpusCallosum" = -0.00549,
              "eta" = 1.43044,
              "log.thalamus~genotypeHAB" = 0.60113,
              "log.thalamus~groupconcussion" = 0.11999,
              "log.pallidostriatum~eta" = 0.83452,
              "log.pallidostriatum~genotypeHAB" = 0.57197,
              "log.pallidostriatum~groupconcussion" = 0.11358,
              "log.neocortex~eta" = 0.85006,
              "log.neocortex~genotypeHAB" = 0.57948,
              "log.neocortex~groupconcussion" = 0.04283,
              "log.midbrain~eta" = 0.84739,
              "log.midbrain~genotypeHAB" = 0.61591,
              "log.midbrain~groupconcussion" = 0.09895,
              "log.pons~eta" = 0.86516,
              "log.pons~genotypeHAB" = 0.53958,
              "log.pons~groupconcussion" = 0.01545,
              "log.cingulateGyrus~eta" = 0.76682,
              "log.cingulateGyrus~genotypeHAB" = 0.65551,
              "log.cingulateGyrus~groupconcussion" = 0.15936,
              "log.hippocampus~eta" = 0.76147,
              "log.hippocampus~genotypeHAB" = 0.57525,
              "log.hippocampus~groupconcussion" = 0.11901,
              "log.supramarginalGyrus~eta" = 0.87999,
              "log.supramarginalGyrus~genotypeHAB" = 0.57436,
              "log.supramarginalGyrus~groupconcussion" = 0.05089,
              "log.corpusCallosum~eta" = 0.67779,
              "log.corpusCallosum~genotypeHAB" = 0.57192,
              "log.corpusCallosum~groupconcussion" = 0.17416,
              "log.thalamus~~log.thalamus" = 0.01308,
              "log.pallidostriatum~~log.pallidostriatum" = 0.00987,
              "log.neocortex~~log.neocortex" = 0.00603,
              "log.midbrain~~log.midbrain" = 0.00402,
              "log.pons~~log.pons" = 0.01053,
              "log.cingulateGyrus~~log.cingulateGyrus" = 0.00451,
              "log.hippocampus~~log.hippocampus" = 0.01247,
              "log.supramarginalGyrus~~log.supramarginalGyrus" = 0.00612,
              "log.corpusCallosum~~log.corpusCallosum" = 0.01602,
              "eta~~eta" = 0.06319)

## ** give appropriate name
dfType <- coefType(lava::estimate(m.sim,lava::sim(m.sim,1e2)),
                   as.lava = FALSE)[,c("name","param","lava")]
## name2lava <- setNames(dfType[!is.na(dfType$lava),"lava"],dfType[!is.na(dfType$lava),"name"])
sim.coefLava <- sim.coef[setdiff(names(sim.coef), dfType[is.na(dfType$lava),"name"])]
## dfType[is.na(dfType$lava),"name"]

## ** null hypotheses
name.test <- paste0(c("log.thalamus", 
                      "log.pallidostriatum", 
                      "log.neocortex", 
                      "log.midbrain", 
                      "log.pons", 
                      "log.cingulateGyrus", 
                      "log.hippocampus", 
                      "log.supramarginalGyrus",
                      "log.corpusCallosum"),"~groupconcussion")
n.test <- length(name.test)

newLinks <- findNewLink(m.test, type = "covariance")$links


## dt.data <- lava::sim(m.sim, n = 1e4, p = sim.coefLava, latent = FALSE)
## m.test2 <- m.test
## covariance(m.test2) <- log.thalamus ~ log.pallidostriatum
## covariance(m.test2) <- log.hippocampus ~ log.supramarginalGyrus
## e.lvm <- estimate(m.test,  data = dt.data)
## modelsearch2(e.lvm, link = newLinks, method.p.adjust = "bonferroni")

## * loop
dt.stat <- NULL
dt.res <- NULL

for(iN in 1:n.N){ # iN <- 1
    cat("sample size=",seqN[iN],": ", sep = "")
    n.tempo <- seqN[iN]
        
    for(iRep in 1:n.repMax){ # iRep <- 1
        cat(iRep," ")
        ls.max <- list()

        ## ** Simulate data
        dt.data <- lava::sim(m.sim, n = n.tempo, p = sim.coefLava, latent = FALSE)
        
        ## ** fit lvm
        e.lvm <- estimate(m.test,  data = dt.data)
        n.coef0 <- length(coef(e.lvm))
        
        if (e.lvm$opt$convergence == 1) {
            next
        }
        if (any(eigen(getVarCov2(e.lvm))$values <= 0)) {
            next
        }

        ## ** search for extentions
        eS.lvm <- try(modelsearch2(e.lvm, link = newLinks, method.p.adjust = "bonferroni", trace = 0, nStep = 2),
                      silent = TRUE)
        if(inherits(eS.lvm,"try-error")){
            next
        }
        ## m.test2 <- m.test
        ## covariance(m.test2) <- log.thalamus~log.pallidostriatum
        
        ## name.covcoef 
        e.lvm <- getNewModel(eS.lvm)

        if (e.lvm$opt$convergence == 1) {
            next
        }
        if (any(eigen(getVarCov2(e.lvm))$values <= 0)) {
            next
        }       
        
        name.coef <- names(coef(e.lvm))
        n.coef <- length(name.coef)

        test.validModel <- all(name.covcoef %in% name.coef)
        vec.statSearch <- c("log.thalamus~~log.pallidostriatum" = eS.lvm$sequenceTest[[1]][eS.lvm$sequenceTest[[1]]$link=="log.thalamus~~log.pallidostriatum","statistic"],
                            "log.hippocampus~~log.corpusCallosum" = eS.lvm$sequenceTest[[1]][eS.lvm$sequenceTest[[1]]$link=="log.hippocampus~~log.corpusCallosum","statistic"])
                         
        ## ** create contrast matrix
        Ccontrast <- matrix(0, ncol = n.coef, nrow = n.test, 
                            dimnames = list(name.test,name.coef))
        diag(Ccontrast[name.test,name.test]) <- 1
            
        ## ** adjustment for multiple comparison
        e.glht2 <- try(glht2(e.lvm, linfct = Ccontrast, rhs = rep(0, n.test)),
                      silent = TRUE)
        if("try-error" %in% class(e.glht2) || is.na(e.glht2$df) || (e.glht2$df<0)){
            p.value.none2 <- NA
            p.value.bonf2 <- NA
            p.value.max2 <- NA
            medianCor.test2 <- NA                
            p.value.max <- rep(NA,length(name.test))
        }else{            
            cor.test <- cov2cor(e.glht2$vcov[name.test,name.test])
            medianCor.test <- median(abs(cor.test[lower.tri(cor.test)]))
            
            e0.glht <- summary(e.glht2, test = univariate())
            p.value.none <- as.double(e0.glht$test$pvalues)
            name.X <- names(e0.glht$test$coef)
            vec.statTest <- e0.glht$test$tstat
        
            eB.glht <- summary(e.glht2, test = adjusted(type = "bonferroni"))
            p.value.bonf <- as.double(eB.glht$test$pvalues)

            if (any(eigen(cor.test)$values <= 0)) {
                p.value.max <- rep(NA,length(name.test))
            } else{
                eS.glht <- summary(e.glht2, test = adjusted(type = "single-step"))
                p.value.max <- as.double(eS.glht$test$pvalues)
            }

        }

        ## ** store results
        dt.tempo <- rbind(data.table(method = "none", variable = name.X, p.value = p.value.none),
                          data.table(method = "Bonferroni", variable = name.X, p.value = p.value.bonf),
                          data.table(method = "Max", variable = name.X, p.value = p.value.max)
                          )        
        dt.tempo[, n := n.tempo]
        dt.tempo[, validModel := test.validModel]
        dt.tempo[, rep := iRep]
        dt.tempo[, seed := iSeed]
        dt.tempo[, corTest := medianCor.test]

        dt.tempo2 <- data.table(rbind(c(vec.statSearch, vec.statTest)))
        dt.tempo2[, max := max(abs(vec.statTest))]
        dt.tempo2[, n := n.tempo]
        dt.tempo2[, validModel := test.validModel]
        dt.tempo2[, rep := iRep]
        dt.tempo2[, seed := iSeed]

        dt.res <- rbind(dt.res, dt.tempo)
        dt.stat <- rbind(dt.stat,dt.tempo2)
        
        if(iRep %% 100 == 0){
            cat("\n")
            filename <- paste0("type1error-S",iter_sim,"(tempo).rds")
            saveRDS(dt.res, file = file.path(path.res,filename))
            filename <- paste0("statistics-S",iter_sim,"(tempo).rds")
            saveRDS(dt.stat, file = file.path(path.res,filename))
        }
        
        if(sum(dt.res[n==n.tempo & method == "none" & variable == "log.thalamus~groupconcussion"]$validModel) == n.rep){
            break
        }
        
    }
}



## * export
filename <- paste0("type1error-S",iter_sim,".rds")
saveRDS(dt.res, file = file.path(path.res,filename))
filename <- paste0("statistics-S",iter_sim,".rds")
saveRDS(dt.stat, file = file.path(path.res,filename))

## * display
print(sessionInfo())

