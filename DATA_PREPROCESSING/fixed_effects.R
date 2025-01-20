library(data.table)

ped <- fread("pedigree.csv")
setorder(ped, V4)   # sort by date
X <- ped[,c(4,5)]  # when the ped file is sorted, just extract the 2 last columns as the 2 fixed effects
write.csv(X, "TRANSFORMACJE_DANYCH/FixedEffects_X.csv",
          row.names = FALSE)