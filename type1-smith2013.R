### type1-smith2013.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 19 2019 (16:36) 
## Version: 
## Last-Updated: dec 19 2019 (17:12) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(ggplot2)

## * settings
cpus <- 25
n.sim <- 1e4
type <- 2

vec.n <- c(100,200,500,1000)

## * generative model
if(type==1){
    Sigma <- diag(3.25,6,6)
    Sigma[4:6,1:3] <- 1
    Sigma[1:3,4:6] <- 1
    ## eigen(Sigma)
}else if(type==2){
    Sigma <- diag(1.1,6,6)
    Sigma[4,1] <- Sigma[1,4] <- 1
    Sigma[5,2] <- Sigma[2,5] <- 1
    Sigma[6,3] <- Sigma[3,6] <- 1
    ## eigen(Sigma)
}

## * warper
warper <- function(n){ ## n <- 1000

    M <- mvtnorm::rmvnorm(n = n, mean = rep(0,6), sigma = Sigma)
    d <- as.data.frame(M)
    names(d) <- c("X1","X2","X3","Z1","Z2","Z3")
    d$Y <- rnorm(n)
    ## ggHeatmap(cov2cor(cor(d)))
 
    e <- lm(Y~0+X1+X2+X3+Z1+Z2+Z3, data = d)
    ## ggHeatmap(cov2cor(vcov(e)))

    vec.p <- summary(e)$coef[paste0("X",1:3),"Pr(>|t|)"]
    iM  <- matrix(NA, nrow = 3, ncol = 9,
                  dimnames = list(names(vec.p),
                                  c("none","bonferroni","hochberg","AB1","AB2","single-step","free","rho","calpha")))

    iM[,"none"] <- vec.p
    iM[,"bonferroni"] <- p.adjust2(vec.p, method = "bonferroni")
    iM[,"hochberg"] <- p.adjust2(vec.p, method = "hochberg")

    AB1 <- p.adjust2(vec.p, method = "AB1", vcov.param = vcov(e))
    iM[,"AB1"] <- AB1
    iM[,"AB2"] <- p.adjust2(vec.p, method = "AB2", vcov.param = vcov(e))
    iM[,"rho"] <- attr(AB1,"r")

 
    C <- lavaSearch2::createContrast(e, c("X1","X2","X3"), add.variance = FALSE)
    e.glht <- multcomp::glht(e, linfct = C$contrast)
 
    iM[,"single-step"] <- summary(e.glht, test = multcomp::adjusted("single-step"))$test$pvalues
    iM[,"free"] <- summary(e.glht, test = multcomp::adjusted("free"))$test$pvalues
    iM[,"calpha"] <- attr(confint(e.glht)$confint,"calpha")
    return(iM)
}
## warper(100)

## * simulation

## ** sequential
if(cpus==1){
    set.seed(10)
    ls.sim <- pblapply(1:n.sim, function(i){
        iOut <- NULL
        for(iN in vec.n){
            iOut <- rbind(iOut,
                          data.table(n = iN, i = i, warper(iN))
                          ) 
        }
        return(iOut)
    })
}

## ** parallel computation
if(cpus>1){
    cl <- snow::makeSOCKcluster(cpus)
    doSNOW::registerDoSNOW(cl)

    pb <- txtProgressBar(max = n.sim, style=3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))

    ls.sim <- foreach::`%dopar%`(
                           foreach::foreach(i=1:n.sim, .options.snow=opts, .packages = "data.table"), {
                               iOut <- NULL
                               for(iN in vec.n){
                                   iOut <- rbind(iOut,
                                                 data.table(n = iN, i = i, warper(iN))
                                                 ) 
                               }
                               return(iOut)
                           })
}

## * process results

dt.sim <- do.call(rbind,ls.sim)
dtL.sim <- melt(dt.sim, id.vars = c("n","i"), value.name = "p", variable.name = "method")
dtLS.sim <- dtL.sim[method %in% c("calpha","rho") == FALSE, .(rep = .N, type1 = any(p<=0.05)), by = c("method","i","n")]
dtLSS.sim <- dtLS.sim[, .(rep = .N, type1 = mean(type1)), by = c("method","n")]

## butils::object2script(dtLSS.sim)
dtLSS.sim <- data.table("method" = c("none", "none", "none", "none", "bonferroni", "bonferroni", "bonferroni", "bonferroni", "hochberg", "hochberg", "hochberg", "hochberg", "AB1", "AB1", "AB1", "AB1", "AB2", "AB2", "AB2", "AB2", "single-step", "single-step", "single-step", "single-step", "free", "free", "free", "free"),
                        "n" = c( 100,  200,  500, 1000,  100,  200,  500, 1000,  100,  200,  500, 1000,  100,  200,  500, 1000,  100,  200,  500, 1000,  100,  200,  500, 1000,  100,  200,  500, 1000),
                        "rep" = c(10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000),
                        "type1" = c(0.1357, 0.1466, 0.1456, 0.1446, 0.0492, 0.0501, 0.0531, 0.0501, 0.0498, 0.0503, 0.0535, 0.0504, 0.0674, 0.0749, 0.0746, 0.0713, 0.0773, 0.0860, 0.0856, 0.0805, 0.0504, 0.0507, 0.0540, 0.0511, 0.0504, 0.0507, 0.0540, 0.0511))


## * display results
 
gg.curve <- ggplot(dtLSS.sim, aes(x = n, y = type1, group = method, color = method))
gg.curve <- gg.curve + geom_point() + geom_line()
gg.curve <- gg.curve + geom_hline(yintercept = 0.05)

## gg.uninf <- ggplot(dtL.sim[method %in% c("calpha","rho") == FALSE], aes(p, group = method, color = method))
## gg.uninf <- gg.uninf + geom_density(size = 2)
## gg.uninf <- gg.uninf + facet_wrap(~n, labeller = label_both)

Sigma <- diag(1,6,6)


mvtnorm::qmvnorm(p = 0.975, mu = rep(0,6), sigma = cov2cor(vcov(e)))
mvtnorm::qmvnorm(p = 0.975, mu = rep(0,6), sigma = diag(1,6,6))


mvtnorm::qmvnorm(p = 0.975, mu = rep(0,3), sigma = diag(1,3,3))$quantile
r <- 0.2
GS <- qnorm(1-0.025/3)
AB1 <- qnorm(1-0.025/(3-(3-1)*sqrt(abs(r))))
AB2 <- qnorm(1-0.025/(3^(1-sqrt(abs(r)))))


1-mvtnorm::pmvnorm(lower = rep(-AB1,3), upper = rep(AB1,3), mean = rep(0,3), sigma = diag(1,3,3))
1-mvtnorm::pmvnorm(lower = rep(-AB2,3), upper = rep(AB2,3), mean = rep(0,3), sigma = diag(1,3,3))
1-mvtnorm::pmvnorm(lower = rep(-GS,3), upper = rep(GS,3), mean = rep(0,3), sigma = diag(1,3,3))

######################################################################
### type1-smith2013.R ends here
