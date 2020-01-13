## path <- "P:/Cluster/LVMproject/article-multipleComparisons"
## setwd(path)
## source("BATCH_simulation-multComp-type1.R")

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 97}
if(is.na(n.iter_sim)){n.iter_sim <- 100}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n", sep = "")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","simulation-multComp-type1")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-multComp-type1")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
library(multcomp)
library(lavaSearch2)

## * settings 
seqN <- c(30,50,75,100,150,200,300,500)
seqCor <-  c(0.1,0.2,0.35,0.65,1,5)##c(0,0.5,1,1.5,3,5)
n.Cor <- length(seqCor)
n.rep <- 100

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
               log.corpusCallosum)~ a * eta + genotypeHAB + b * groupconcussion)
latent(m.sim) <- ~eta

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
n.Cor <- length(seqCor)

## ** value for the simulation
sim.coef <- c("log.pallidostriatum" = 0.07255,
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

## * loop
dt <- NULL
pb <- txtProgressBar(max = n.Cor)
dt.res <- NULL

for(iN in 1:n.N){ # iN <- 5
    for(iCor in 1:n.Cor){ # iCor <- 1
        cat("sample size=",seqN[iN],", correlation=",seqCor[iCor],": ", sep = "")
        n.tempo <- seqN[iN]
        a.tempo <- seqCor[iCor]
        
        for(iRep in 1:n.rep){ # iRep <- 1
            cat(iRep," ")
            ls.max <- list()

            ## ** Simulate data
            dt.data <- lava::sim(m.sim, n = n.tempo, p = c(a = a.tempo, sim.coefLava, b = 0), latent = FALSE)

            ## ** models
            e.lvm <- estimate(m.test,  data = dt.data)
            ## summary(e.lvm)
            ##setdiff(c(endogenous(m.test),exogenous(m.test)), names(dt.data))            
            ## e.lvm <- estimate(m.test,  data = dt.data, control = list(constrain = TRUE, starterfun = "startvalues"))
            ## e.lvm <- estimate(m.test,  data = dt.data, control = list(constrain = TRUE, starterfun = "startvalues0"))
            ## e.lvm <- estimate(m.test,  data = dt.data, control = list(constrain = TRUE, starterfun = "startvalues1"))
            ## e.lvm <- estimate(m.test,  data = dt.data, control = list(starterfun = "startvalues"))
            ## e.lvm <- estimate(m.test,  data = dt.data, control = list(starterfun = "startvalues0"))
            ## e.lvm <- estimate(m.test,  data = dt.data, control = list(starterfun = "startvalues1"))
            
            if (e.lvm$opt$convergence == 1) {
                next
            }
            if (any(eigen(getVarCov2(e.lvm))$values <= 0)) {
                next
            }
            
            name.coef <- names(coef(e.lvm))
            n.coef <- length(name.coef)
            
            ## ** create contrast matrix
            Ccontrast <- matrix(0, ncol = n.coef, nrow = n.test, 
                                dimnames = list(name.test,name.coef))
            diag(Ccontrast[name.test,name.test]) <- 1
            
            ## ** adjustment for multiple comparison
            ## *** lava 
            e.glht <- glht(e.lvm, linfct = Ccontrast)
            cor.test <- cov2cor(e.glht$vcov[name.test,name.test])
            medianCor.test <- median(abs(cor.test[lower.tri(cor.test)]))

            e0.glht <- summary(e.glht, test = univariate())
            eB.glht <- summary(e.glht, test = adjusted(type = "bonferroni"))
            eHochberg.glht <- summary(e.glht, test = adjusted(type = "hochberg"))
            eHommel.glht <- summary(e.glht, test = adjusted(type = "hommel"))
            eS.glht <- summary(e.glht, test = adjusted(type = "single-step"))

            name.X <- names(e0.glht$test$coef)
            p.value.none <- as.double(e0.glht$test$pvalues)
            p.value.bonf <- as.double(eB.glht$test$pvalues)
            p.value.hoch <- as.double(eHochberg.glht$test$pvalues)
            p.value.homm <- as.double(eHommel.glht$test$pvalues)
            p.value.max <- as.double(eS.glht$test$pvalues)

            ### *** lavaSearch 2
            e.glht2 <- try(glht2(e.lvm, linfct = Ccontrast, rhs = rep(0, n.test)),
                           silent = TRUE)
            ## e.glht2 <- try(glht2(e.lvm, linfct = Ccontrast), silent = TRUE)            
            if("try-error" %in% class(e.glht2) || is.na(e.glht2$df) || (e.glht2$df<0)){
                p.value.none2 <- NA
                p.value.bonf2 <- NA
                p.value.max2 <- NA
                medianCor.test2 <- NA                
            }else{            
                cor.test2 <- cov2cor(e.glht2$vcov[name.test,name.test])
                medianCor.test2 <- median(abs(cor.test2[lower.tri(cor.test)]))

                e0.glht2 <- summary(e.glht2, test = univariate())
                eB.glht2 <- summary(e.glht2, test = adjusted(type = "bonferroni"))
                eHochberg.glht2 <- summary(e.glht2, test = adjusted(type = "hochberg"))
                eHommel.glht2 <- summary(e.glht2, test = adjusted(type = "hommel"))
                eS.glht2 <- summary(e.glht2, test = adjusted(type = "single-step"))
                p.value.none2 <- as.double(e0.glht2$test$pvalues)
                p.value.bonf2 <- as.double(eB.glht2$test$pvalues)
                p.value.hoch2 <- as.double(eHochberg.glht2$test$pvalues)
                p.value.homm2 <- as.double(eHommel.glht2$test$pvalues)
                p.value.max2 <- as.double(eS.glht2$test$pvalues)
            }
            
            
            ## ** store results
            dt.tempo <- rbind(data.table(method = "none", variable = name.X, p.value = p.value.none),
                              data.table(method = "Bonferroni", variable = name.X, p.value = p.value.bonf),
                              data.table(method = "Hochberg", variable = name.X, p.value = p.value.hoch),
                              data.table(method = "Hommel", variable = name.X, p.value = p.value.homm),
                              data.table(method = "Max", variable = name.X, p.value = p.value.max),
                              data.table(method = "none2", variable = name.X, p.value = p.value.none2),                              
                              data.table(method = "Bonferroni2", variable = name.X, p.value = p.value.bonf2),                              
                              data.table(method = "Hochberg2", variable = name.X, p.value = p.value.hoch2),                              
                              data.table(method = "Hommel2", variable = name.X, p.value = p.value.homm2),                              
                              data.table(method = "Max2", variable = name.X, p.value = p.value.max2)                              
                              )
            dt.tempo[, n := n.tempo]
            dt.tempo[, a := a.tempo]
            dt.tempo[, rep := iRep]
            dt.tempo[, seed := iSeed]
            dt.tempo[, corTest := medianCor.test]
            dt.tempo[, corTest2 := medianCor.test2]

            dt.res <- rbind(dt.res, dt.tempo)
        }
        cat("\n")
    }
    filename <- paste0("type1error-S",iter_sim,"(tempo).rds")
    saveRDS(dt.res, file = file.path(path.res,filename))

}

## * export
filename <- paste0("type1error-S",iter_sim,".rds")
saveRDS(dt.res, file = file.path(path.res,filename))

## * display
print(sessionInfo())
