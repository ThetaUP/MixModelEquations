library(data.table)

# load one of the two response variables (y)
y1 <- fread("TRANSFORMACJE_DANYCH/Y1_FirstLastInsem.txt")

# load the A matrix
A_mat <- fread("TRANSFORMACJE_DANYCH/A_matrix.csv")
A_mat[,1] <- NULL  # remove 1st column, because it indicates the row names

N_phenotyped_cows <- dim(y1)[1]  # number of cows for which we have phenotypes/genotypes
N_total_cows <- dim(A_mat)[1]    # total number of cows
N_unphenotyped_cows <- N_total_cows - N_phenotyped_cows  # number of cows without phenotype/genotype


# create the Z matrix
diag_part <- diag(N_phenotyped_cows)
Z <- matrix(rep(0, N_phenotyped_cows*N_total_cows),
            nrow = N_phenotyped_cows, ncol = N_total_cows)  # dim(Z) = N_phenotyped x N_total

Z[ , (N_unphenotyped_cows+1):ncol(Z)] <- diag_part

write.csv(Z, "TRANSFORMACJE_DANYCH/Z_matrix.csv",
          row.names = FALSE)
