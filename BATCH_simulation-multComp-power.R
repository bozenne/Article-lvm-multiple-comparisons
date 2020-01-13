## path <- "P:/Cluster/LVMproject/article-multipleComparisons"
## setwd(path)
## source("BATCH_simulation-multComp-power.R")

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
path.res <- file.path(path,"Results","simulation-multComp-power")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-multComp-power")
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
seqCor <-  c(0.25,1,2,3,4,5)
n.Cor <- length(seqCor)
n.rep <- 50

## * model
## ** generative model
m.sim <- lvm(log.thalamus[0:sigma] ~ a * eta + genotypeHAB + b * groupconcussion,
             log.pallidostriatum[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion, 
             log.neocortex[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion,
             log.midbrain[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion,
             log.pons[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion,
             log.cingulateGyrus[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion,
             log.hippocampus[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion,
             log.supramarginalGyrus[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion,
             log.corpusCallosum[0:sigma] ~ a * eta + genotypeHAB + 0 * groupconcussion,
             eta[0:psi] ~ 1)
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

## ** give appropriate name
dfType <- coefType(lava::estimate(m.sim,lava::sim(m.sim,1e2)),
                   as.lava = FALSE)[,c("name","param","lava")]

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

for(iN in 1:n.N){ # iN <- 8 
    for(iCor in 1:n.Cor){ # iCor <- 6
        cat("sample size=",seqN[iN],", correlation=",seqCor[iCor],": ", sep = "")
        n.tempo <- seqN[iN]
        a.tempo <- seqCor[iCor]
        
        for(iRep in 1:n.rep){ # iRep <- 1
            cat(iRep," ")
            ls.max <- list()

            ## ** Simulate data
            ## n.tempo <- 500; a.tempo <- 5
            ## m.sim$par
            dt.data <- lava::sim(m.sim, n = n.tempo, p = c(a = sqrt(a.tempo), "sigma" = 5.25-a.tempo, "eta" = 1, b = 0.4), latent = FALSE)
            ## print(apply(dt.data[,grep("log.",names(dt.data))],2,sd))
            
            ## ** models
            e.lvm <- estimate(m.test,  data = dt.data)
            ## cat("pons :",coef(e.lvm)["eta~~eta"] + coef(e.lvm)["log.pons~eta"] * coef(e.lvm)["log.pons~~log.pons"],
                ## " pallidostriatum :",coef(e.lvm)["eta~~eta"] + coef(e.lvm)["log.pallidostriatum~eta"] * coef(e.lvm)["log.pallidostriatum~~log.pallidostriatum"],"\n")
            ## summary(e.lvm)$coef
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
            
            ## *** lavaSearch 2
            eS.lvm <- e.lvm
            mess <- try({
                ## summary(eS.lvm)$coef[name.test,]
                ## summary2(eS.lvm)$coef[name.test,]
                sCorrect(eS.lvm) <- TRUE
                e.glht2 <- try(glht2(eS.lvm, linfct = Ccontrast, rhs = rep(0, n.test)),
                               silent = TRUE)},
                silent = TRUE)
            
            ## e.glht2 <- try(glht2(e.lvm, linfct = Ccontrast), silent = TRUE)            
            if("try-error" %in% class(mess) || is.na(e.glht2$df) || (e.glht2$df<0)){
                p.value.none2 <- NA
                p.value.bonf2 <- NA
                p.value.max2 <- NA
                medianCor.test2 <- NA
            }else{            
                cor.test2 <- cov2cor(e.glht2$vcov[name.test,name.test])
                medianCor.test2 <- median(abs(cor.test2[lower.tri(cor.test2)]))

                e0.glht2 <- summary(e.glht2, test = univariate())
                eB.glht2 <- summary(e.glht2, test = adjusted(type = "bonferroni"))
                eHochberg.glht2 <- summary(e.glht2, test = adjusted(type = "hochberg"))
                eHommel.glht2 <- summary(e.glht2, test = adjusted(type = "hommel"))
                eS.glht2 <- suppressWarnings(summary(e.glht2, test = adjusted(type = "single-step")))
                p.value.none2 <- as.double(e0.glht2$test$pvalues)
                p.value.bonf2 <- as.double(eB.glht2$test$pvalues)
                p.value.hoch2 <- as.double(eHochberg.glht2$test$pvalues)
                p.value.homm2 <- as.double(eHommel.glht2$test$pvalues)
                p.value.max2 <- as.double(eS.glht2$test$pvalues)
                p.value.F2 <- compare2(eS.lvm, par = rownames(Ccontrast))$p.value                
            }
            
            
            ## ** store results
            dt.tempo <- rbind(data.table(method = "none2", variable = name.test, p.value = p.value.none2),                              
                              data.table(method = "Bonferroni2", variable = name.test, p.value = p.value.bonf2),                              
                              data.table(method = "Hochberg2", variable = name.test, p.value = p.value.hoch2),                              
                              data.table(method = "Hommel2", variable = name.test, p.value = p.value.homm2),                              
                              data.table(method = "Max2", variable = name.test, p.value = p.value.max2),                              
                              data.table(method = "Fisher2", variable = "all", p.value = p.value.F2)                              
                              )
            dt.tempo[, n := n.tempo]
            dt.tempo[, a := a.tempo]
            dt.tempo[, rep := iRep]
            dt.tempo[, seed := iSeed]
            dt.tempo[, corTest2 := medianCor.test2]

            dt.res <- rbind(dt.res, dt.tempo)
        }
        cat("\n")
    }
    filename <- paste0("type1error-S",iter_sim,"(tempo).rds")
    saveRDS(dt.res, file = file.path(path.res,filename))

}

## * export
## dt.res[method=="Max2", .(min.p.value = min(p.value)), by = c("a","rep")][,mean(min.p.value<=0.05),by="a"]
filename <- paste0("type1error-S",iter_sim,".rds")
saveRDS(dt.res, file = file.path(path.res,filename))

## * display
print(sessionInfo())
