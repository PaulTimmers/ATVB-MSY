#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
options(digits=6)
library(data.table)
library(ggplot2)

#----
# Load data
#----

haplo <- fread("st002_02_haplogroups.tsv")
ethnicity <- fread("st002_03_ethnicity.tsv")
df <- ethnicity[haplo,,on="iid"]

#----
# Most common haplogroup by ethnicity
#----

# Subset to well-defined ethnic groupings
df1 <- df[!is.na(ethnicity) & !grepl("and|or|not|other|mixed",ethnicity,ignore.case=T),]

# Keep only most abundant haplogroup by ethnicity
df2 <- df1[,.(N_haplo=.N),by=c("ethnicity","haplogroup_short")][order(-N_haplo,ethnicity)]
df3 <- df2[!duplicated(ethnicity),]

# Calculate percentage of people carrying the most common haplogroup
total <- df2[,.(N_total=sum(N_haplo)),by="ethnicity"]
df4 <- df3[total,,on="ethnicity"]
df4[,`%`:=sprintf("%2.2f", 100*N_haplo/N_total)]

# Order table nicely
df4[,haplogroup_short:=factor(haplogroup_short,levels=unique(haplogroup_short))]
df5 <- df4[order(haplogroup_short, -N_haplo)]

# Write to file
fwrite(df5, "st002_04_most_common_by_ethnicity.tsv", sep="\t", quote=FALSE, na="NA")



# Summarise
#----------

ph <- copy(df1)

n <- ph[,.(N_total=.N),by=c("ethnicity")]
hg <- ph[,.(N_haplo=.N),by=c("ethnicity","haplogroup_short")]
hgn <- n[hg,,on="ethnicity"]
hgn[,pct:=N_haplo/N_total]
hgn[,pct_se:=sqrt(pct * (1 - pct) / N_total)]


pdf("st002_04_most_common_by_ethnicity.pdf", width=8, height=6)
ggplot(hgn[pct > 0.01, ], aes(x=pct, y=haplogroup_short, fill=haplogroup_short)) + 
geom_errorbarh(aes(xmin=pct - qnorm(0.975)*pct_se, xmax=pct + qnorm(0.975)*pct_se)) +
geom_bar(stat="identity") + guides(fill="none") +
facet_wrap(.~ethnicity, scales="free_y") + 
labs(x="Percentage with Y-Haplogroup",y="") +
theme_minimal()
dev.off()