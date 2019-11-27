## * Run analysis
source("analysis.R")

## * table2
table2 <- data.frame(matrix(NA, nrow = 9, ncol = 7))
names(table2) <- c("region","concussion effect (\\%)","statistic","no adjustment","Bonferroni","max-test","step-down max-test")

table2[["region"]] <- gsub("log.|~groupconcussion","",names(ls.pvalue[[1]][["coefficients"]]))
table2[["concussion effect (\\%)"]] <- round(100*(exp(ls.pvalue[[1]][["coefficients"]]) - 1), 2)
table2[["statistic"]] <- round(ls.pvalue[[1]][["tstat"]],2)
table2[["no adjustment"]] <- round(ls.pvalue[["none"]][["pvalues"]],3)
table2[["Bonferroni"]] <- round(ls.pvalue[["Bonferroni"]][["pvalues"]],3)
table2[["Dunnett"]] <- round(ls.pvalue[["Dunnett"]][["pvalues"]],3)
table2[["step-down Dunnett"]] <- round(ls.pvalue[["DunnettDown"]][["pvalues"]],3)


addtorow <- list(pos = list(0,0),
                 command = c("      &                   &           &       \\multicolumn{4}{c}{p-value} \\\\ \\cmidrule(r){4-7} \n",
                             paste0(paste0(names(table2), collapse = " & "), "\\\\\n")
                             )
                 )
print(xtable(table2, type = "latex",
             label = "tab:Ex-multComp",
             caption = "Test of the null hypothesis of equal distribution volume in the healthy group vs. the mTBI group.
                        The column concussion indicates the percentage increase in distribution volume due to concussion.
"),
 add.to.row = addtorow, include.rownames=FALSE, include.colnames = FALSE)



