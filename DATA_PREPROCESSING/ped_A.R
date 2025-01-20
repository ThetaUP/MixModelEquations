library(data.table)
library(AGHmatrix)

path_ped <- "pedigree.csv"
ped_orig <- fread(path_ped)
setorder(ped_orig, V4)
View(ped_orig)
unique_sire <- unique(ped_orig$V2)
unique_dam <- unique(ped_orig$V3)
ped_new <- matrix(rep(0,(length(unique_sire)+length(unique_dam))*3), nrow = length(unique_dam)+length(unique_sire), ncol = 3)
unique_parents <- c(unique_sire, unique_dam)
ped_new[,1] <- unique_parents
ped_new[,c(2,3)] <- 0
cows <- ped_orig[,c(1,2,3)]
ped_final <- rbind(ped_new, cows)

A <- Amatrix(ped_final)
write.csv(A,"TRANSFORMACJE_DANYCH/A_matrix.csv",
          row.names = TRUE)
