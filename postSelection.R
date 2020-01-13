## * Packages
library(data.table)

## * Import
## list.files("Results")
dt.postSelection <- readRDS(file = file.path("Results","simulation-postSelection-pvalue.rds"))

print(dt.postSelection)
##    method   n  n.rep   valid.pc type1error
## 1:    Max  50 106536 0.09386498 0.05595292
## 2:    Max 100  18806 0.53174519 0.05370626
## 3:    Max 500  10000 1.00000000 0.04930000
