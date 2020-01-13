## path <- "P:/Cluster/LVMproject/article-multipleComparisons/"
## setwd(path)
## source("BATCH_simulation-modelsearch-power.R")

## * seed
iter_sim <- as.numeric(Sys.getenv("SGE_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SGE_TASK_LAST"))
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 1}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n", sep = "")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","simulation-modelsearch-power")
if(dir.exists(path.res)==FALSE){
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulation-modelsearch-power")
if(dir.exists(path.output)==FALSE){
    dir.create(path.output)
}

## * libraries
library(lava)
library(data.table)
library(lavaSearch2)
lava.options(search.calc.quantile.int = FALSE)

## * settings
seqN <- c(30,50,75,100,150,200,300,500)
seqCor <- c(0,0.75,1,1.5,2.5,5)
n.rep <- 50
p <- 15

## * prepare
n.N <- length(seqN)
n.Cor <- length(seqCor)
X.all <- paste0("X",1:p)
link.all <- c("eta~Z",paste0("eta~",X.all))

m <- lvm(c(Y1,Y2,Y3,Y4,Y5)~eta, eta~Treatment)
covariance(m) <- Y1 ~ Y2
covariance(m) <- Y3 ~ Y4
covariance(m) <- Y3 ~ Y5
latent(m) <- ~eta

m.sim <- m
for(iConfounder in X.all){
    regression(m.sim) <- as.formula(paste0(iConfounder," ~ a * confounder"))
}
regression(m.sim) <- eta ~ 0.25 * Z
latent(m.sim) <- ~confounder
## plot(m.sim)

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
            dt.data <- lava::sim(m.sim, n = n.tempo, p = c("a" = a.tempo))
            
            ## ** lvm
            e.lvm <- estimate(m,  data = dt.data)
      
            ## ** run score test
            res.Wald <- try(modelsearch2(e.lvm,
                                         data = dt.data,
                                         nStep = 1,
                                         link = link.all,
                                         method.p.adjust = "max",
                                         method.maxdist = "approximate",                                         
                                         trace = 0),
                            silent = TRUE)
            if(inherits(res.Wald,"try-error")){
                print(res.Wald)
                next
            }            ## res.Wald$sequenceModel

            iTable <- getStep(res.Wald, slot = "sequenceTest")
            iTable$bonf.p.value <- p.adjust(iTable$p.value, method = "bonferroni")
            iTable$hoch.p.value <- p.adjust(iTable$p.value, method = "hochberg")
            iTable$hommel.p.value <- p.adjust(iTable$p.value, method = "hommel")

            ls.max$Wald <- iTable[iTable$link=="eta~Z","p.value"]
            ls.max$Wald.bonf <- iTable[iTable$link=="eta~Z","bonf.p.value"]
            ls.max$Wald.max <- iTable[iTable$link=="eta~Z","adjusted.p.value"]

            ls.max$selected <- which(iTable$link=="eta~Z") %in% which.max(abs(iTable$statistic))
            
            cor.test <- getStep(res.Wald, slot = "sequenceSigma")
            medianCor.test <- median(abs(cor.test[lower.tri(cor.test, diag = FALSE)]))

            ## ** merge
            dt.tempo <- cbind(n = n.tempo,
                              a = a.tempo,
                              rep = iRep,
                              seed = iSeed,
                              medianCor = medianCor.test,
                              as.data.table(ls.max))
            dt.res <- rbind(dt.res, dt.tempo)
        }
        cat("\n")

        ## ** export (tempo)
        filename <- paste0("type1error-S",iter_sim,"(tempo).rds")
        saveRDS(dt.res, file = file.path(path.res,filename))
    }
}

## * export
filename <- paste0("type1error-S",iter_sim,".rds")
saveRDS(dt.res, file = file.path(path.res,filename))

## * display
print(sessionInfo())

