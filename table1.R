## * Run analysis
source("analysis.R")

## * table1
tableS <- summary(e.search2, print = FALSE)$table

table1 <- data.frame(matrix(NA, nrow = 3, ncol = 6))
names(table1) <- c("covariance parameter","number of tests","max statistic","no adjustment","Bonferroni","max-test")

options(scipen=999)
table1[["covariance parameter"]] <- paste0("(",gsub("\\~\\~",", ",gsub("log.","",tableS[["link"]])),")")
table1[["number of tests"]] <- tableS[["nTests"]]
table1[["max statistic"]] <- as.character(round(tableS[["statistic"]], 2))
table1[["no adjustment"]] <- as.character(round(tableS[["p.value"]], digits = 6))
table1[["Bonferroni"]] <- as.character(round(pmin(tableS[["p.value"]]*tableS[["nTests"]],1), digits = 6))
table1[["Dunnett"]] <- as.character(round(tableS[["adjusted.p.value"]], digits = 6))

addtorow <- list(pos = list(0,0),
                 command = c("                    &                    &              & \\multicolumn{3}{c}{p-value} \\\\ \\cmidrule(r){4-6} \n",
                             paste0(paste0(names(table1), collapse = " & "), "\\\\\n")
                             )
                 )
print(xtable(table1, type = "latex", 
             label = "tab:Ex-modelsearchMax",
             caption = "Result of the first three steps of FSS for local dependencies. Each row corresponds to a step."),
      add.to.row = addtorow, include.rownames=FALSE, include.colnames = FALSE)
