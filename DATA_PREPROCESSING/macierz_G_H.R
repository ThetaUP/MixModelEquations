library(data.table)
library(AGHmatrix)
library(MASS)

# przygotowywanie macierzy G____________________________________________________________________________________________________-
genotypes <- fread("Genotypes.csv")

# sorting genotypes according to the order of the ped file
ped <- fread("pedigree.csv")
genotypes <- merge(ped, genotypes, by = "V1")
setorder(genotypes, "V4.x")
genotypes <- genotypes[,5:ncol(genotypes)]
genotypes <- as.matrix(genotypes)
unique(as.vector(genotypes))  # there are values which are > 2, we need to recode them to -9
genotypes[genotypes > 2] <- -9
G <- Gmatrix(SNPmatrix = as.matrix(genotypes),method = "VanRaden",missingValue = -9,maf = 0.01)

# write the G matrix to a file just in case
write.csv(G,"TRANSFORMACJE_DANYCH/G_matrix.csv",
          row.names = FALSE)


# tworzenie macierzy H_____________________________________________________________________________________________________
A <- fread("TRANSFORMACJE_DANYCH/A_matrix.csv")
A$V1 <- NULL
N_total <- dim(A)[1]
N_PhenoGeno <- dim(ped)[1]
delta_N <- N_total - N_PhenoGeno
A_inv <- solve(A)
G_inv <- solve(as.matrix(G))
A22 <- A[1:N_PhenoGeno, 1:N_PhenoGeno]
A22_inv <- solve(as.matrix(A22))

TopLeft <- matrix(0, nrow = delta_N, ncol = delta_N)
BottomLeft <- matrix(0, nrow = N_PhenoGeno, ncol = delta_N)
TopRight <- matrix(0, nrow = delta_N, ncol = N_PhenoGeno)
BottomRight <- G_inv - A22_inv

Top <- cbind(TopLeft, TopRight)
Bottom <- cbind(BottomLeft, BottomRight)
H_inv <- rbind(Top, Bottom)

# write to file
write.csv(H_inv,"TRANSFORMACJE_DANYCH/H_inv_matrix.csv",
          row.names = FALSE)
