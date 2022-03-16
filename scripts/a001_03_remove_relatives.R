#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
options(digits=6)
library(data.table)


#----
# Load
#----

iids <- fread("st001_03_iids.txt",header=FALSE,col.names="iid")
rel <- fread("/exports/igmm/eddie/UK-BioBank-proj19655/processing/ukb19655_rel_s488363.dat")
rel1 <- rel[ID1 %in% iids$iid & ID2 %in% iids$iid,]


#----
# Get most related individuals
#----

rel2 <- copy(rel1) 
rel2[,c("ID1","ID2"):=.(ID2,ID1)] # Reverse order of ID1 and ID2
rel3 <- rbind(rel1, rel2) 
ID1n <- rel3[,.N,by="ID1"][order(-N)] # Count ID1 pairs, and sort so highest is on top


#----
# Remove sequentially
#----

# Create array for results
remove <- array(NA, length(unique(ID1n$ID1)))


i <- 1
repeat {
	# Remove ID1 that has most pairs
	remove[i] <- ID1n[1,ID1] 
	rel3 <- rel3[ID1 != remove[i] & ID2 != remove[i], ]
	cat(sprintf("\r%i",i))

	if(nrow(rel3) == 0) break # Stop if no more related pairs

	i <- i+1
	ID1n <- rel3[,.N,by="ID1"][order(-N)] # Recount ID1 pairs, and sort so new highest is on top 
}

# Format results
remove <- remove[!is.na(remove)] 



#----
# Output
#----

# Create new list of iids
iids1 <- iids[!iid %in% remove,]

# Update QC statistics
stats <- fread("st001_03_stats.txt", header=FALSE, col.names=c("stat","n"))
stats <- rbind(stats, data.table(stat="Relatedness", n=nrow(iids1)))

# Write to file
fwrite(stats, "st001_03_stats.txt", col.names=FALSE, sep="\t", quote=FALSE, na="NA")
fwrite(iids1, "st001_03_iids.txt", col.names=FALSE, sep="\t", quote=FALSE, na="NA")
