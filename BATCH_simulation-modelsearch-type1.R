## path <- "P:/Cluster/LVMproject/article-multipleComparisons/"
## setwd(path)
## source("BATCH_simulation-modelsearch-type1.R")

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 10}
if(is.na(n.iter_sim)){n.iter_sim <- 100}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n", sep = "")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","simulation-modelsearch-type1")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-modelsearch-type1")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
library(lavaSearch2)
lava.options(search.calc.quantile.int = FALSE)

## * settings
seqN <- c(30,50,75,100,150,200,300,500) ## c(100,300,500)
seqCor <- c(0,0.6,1,1.5,2.5,5)
n.rep <- 100
p <- 15

## * prepare
n.N <- length(seqN)
n.Cor <- length(seqCor)
X.all <- paste0("X",1:p)

m <- lvm(c(Y1,Y2,Y3,Y4,Y5)~eta, eta~Treatment)
m <- addvar(m, var = c(X.all))
latent(m) <- ~eta
link.all <- paste0("eta","~",X.all)

covariance(m) <- Y1 ~ Y2
covariance(m) <- Y3 ~ Y4
covariance(m) <- Y3 ~ Y5

m.sim <- m
for(iConfounder in X.all){
  regression(m.sim) <- as.formula(paste0(iConfounder," ~ a * confounder"))
}
latent(m.sim) <- ~confounder

## * loop
dt.res <- NULL

for(iN in 1:n.N){ # iN <- 5
    for(iCor in 1:n.Cor){ # iCor <- 3
        cat("sample size=",seqN[iN],", correlation=",seqCor[iCor],": ", sep = "")
        n.tempo <- seqN[iN]
        a.tempo <- seqCor[iCor]
      
         for(iRep in 1:n.rep){ # iRep <- 1
            cat(iRep," ")

            ls.max <- list()
      
            ## ** Simulate data
            dt.data <- lava::sim(m.sim, n = n.tempo, p = c(a = a.tempo), latent = FALSE)
            ## head(dt.data) ; dt.data$eta <- NULL
            
            ## ** lvm
            e.lvm <- estimate(m,  data = dt.data)
      
            ## ** run score test
            tps <- system.time(
                res.Wald <- try(modelsearch2(e.lvm,
                                             data = dt.data,
                                             link = link.all,
                                             nStep = 1,
                                             method.p.adjust = "fastmax",
                                             method.maxdist = "approximate",
                                             trace = 0), silent = TRUE)
            )
            if(inherits(res.Wald,"try-error")){next}
            
            iTable <- getStep(res.Wald, slot = "sequenceTest")
            ## iTable[match(paste0("eta~X",1:15),iTable$link),]
            ls.max$Wald <- iTable[NROW(iTable),"p.value"]
            ls.max$Wald.bonf <- p.adjust(iTable[,"p.value"], method = "bonferroni")[NROW(iTable)]
            ls.max$Wald.hoch <- p.adjust(iTable[,"p.value"], method = "hochberg")[NROW(iTable)]
            ls.max$Wald.hommel <- p.adjust(iTable[,"p.value"], method = "hommel")[NROW(iTable)]
            ls.max$Wald.nummax <- iTable[NROW(iTable),"adjusted.p.value"]
            ls.max$tps <- tps["elapsed"]


            cor.test <- getStep(res.Wald, slot = "sequenceSigma")
            medianCor.test <- median(abs(cor.test[lower.tri(cor.test, diag = FALSE)]))

            ls.max$Wald.samp <- as.numeric(NA)
            ls.max$Wald.sampmax <- as.numeric(NA)
            ls.max$tps.samp <- as.numeric(NA)
            ls.max$Wald.wild <- as.numeric(NA)
            ls.max$Wald.wildmax <- as.numeric(NA)
            ls.max$tps.wild <- as.numeric(NA)
            
            tps2 <- system.time(
                res.Wald2 <- try(modelsearch2(e.lvm, data = dt.data, link = link.all,
                                              nStep = 1, method.p.adjust = "fastmax",
                                              method.maxdist = "resampling",
                                              trace = 0, n.sample = 1e4), silent = TRUE)
            )
            if(!inherits(res.Wald2,"try-error")){
                iTable2 <- getStep(res.Wald2, slot = "sequenceTest")
                ls.max$Wald.samp <- iTable2[NROW(iTable2),"p.value"]
                ls.max$Wald.sampmax <- iTable2[NROW(iTable2),"adjusted.p.value"]
                ls.max$tps.samp <- tps2["elapsed"]
            }else{
                print(res.Wald2)
            }

            tps3 <- system.time(
                res.Wald3 <- try(modelsearch2(e.lvm, data = dt.data, link = link.all,
                                              nStep = 1, method.p.adjust = "fastmax",
                                              method.maxdist = "bootstrap",
                                              trace = 0, n.sample = 1e4), silent = TRUE)
            )
            if(!inherits(res.Wald3,"try-error")){
                iTable3 <- getStep(res.Wald3, slot = "sequenceTest")
                ls.max$Wald.wild <- iTable3[NROW(iTable3),"p.value"]
                ls.max$Wald.wildmax <- iTable3[NROW(iTable3),"adjusted.p.value"]
                ls.max$tps.wild <- tps3["elapsed"]
            }else{
                print(res.Wald3)
            }

            ### ** merge
            dt.mean <- cbind(n = n.tempo,
                             a = a.tempo,
                             rep = iRep,
                             seed = iSeed,                       
                             medianCor = medianCor.test,
                             as.data.table(ls.max))

            dt.mean[,statistic := .(.(iTable[match(link.all,iTable$link),"statistic"]))]
            dt.mean[,Sigma := .(.(cor.test))]
            dt.res <- rbind(dt.res, dt.mean)
        }
        cat("\n")
    }
    ### ** export (tempo)
    filename <- paste0("type1error-S",iter_sim,"(tempo).rds")
    saveRDS(dt.res, file = file.path(path.res,filename))
}

## * export
filename <- paste0("type1error-S",iter_sim,".rds")
saveRDS(dt.res, file = file.path(path.res,filename))

## * display
print(sessionInfo())
