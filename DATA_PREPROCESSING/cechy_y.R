library(data.table)

lakt <- fread("laktacje.csv")
lakt$calvingDate <- as.Date(lakt$calvingDate, format = "%d.%m.%Y")
lakt$inseminationDate <- as.Date(lakt$inseminationDate, format = "%d.%m.%Y")
unique_cows <- unique(lakt$cowId)

# first - last insemination_______________________________________________________________________________________
FirstLastInsem <- 1:length(unique_cows)
CowID <- 1:length(unique_cows)

for(i in 1:length(FirstLastInsem)){
  current_cow <- lakt[lakt$cowId == unique_cows[i],]
  first_insemination <- as.Date(current_cow$inseminationDate[1], format = "%d.%m.%Y")
  last_insemination <- as.Date(current_cow$inseminationDate[nrow(current_cow)], format = "%d.%m.%Y")
  days_diff <- as.numeric(last_insemination - first_insemination)
  FirstLastInsem[i] <- days_diff
  CowID[i] <- current_cow$cowId[1]
}

FirstLastInsem <- data.frame(cbind(CowID, FirstLastInsem))

# set the order to the same as in the pedigree file
ped <- fread("pedigree.csv")
setorder(ped, V4)  # sort ped by year
colnames(ped)[1] <- "CowID"
FirstLastInsem <- merge(FirstLastInsem, ped, by = "CowID")
setorder(FirstLastInsem, V4)
FirstLastInsem <- as.vector(FirstLastInsem$FirstLastInsem)
write(FirstLastInsem, "TRANSFORMACJE_DANYCH/Y1_FirstLastInsem.txt")



# calving date - last insemination_______________________________________________________________________________________
LastInsemCalving <- 1:length(unique_cows)
CowID <- 1:length(unique_cows)

for(i in 1:length(LastInsemCalving)){
  current_cow <- lakt[lakt$cowId == unique_cows[i],]
  last_insemination <- as.Date(current_cow$inseminationDate[nrow(current_cow)], format = "%d.%m.%Y")
  calving <- as.Date(current_cow[current_cow$success == 1,]$calvingDate, format = "%d.%m.%Y")
  LastInsemCalving[i] <- calving - last_insemination
  CowID[i] <- current_cow$cowId[1]
}


LastInsemCalving <- data.frame(cbind(CowID, LastInsemCalving))

# set the order to the same as in the pedigree file
setorder(ped, V4)
LastInsemCalving <- merge(LastInsemCalving, ped, by = "CowID")
setorder(LastInsemCalving, V4)
LastInsemCalving <- as.vector(LastInsemCalving$LastInsemCalving)
write(LastInsemCalving, "TRANSFORMACJE_DANYCH/Y2_LastInsemCalving.txt")

