library(rentrez)
library(DECIPHER)
library(tidyverse)
library(stringi)
library(ape)
library(RSQLite)
library(dplyr)
library(muscle)
library(ggplot2)
library(ggtree)
library(Biostrings)
library(phytools)
library(phangorn)
library(msa)
library(viridisLite)
library(scico)
library(igraph)

SALMONIDAE_search <- entrez_search(db = "nuccore", term = ("Salmonidae[ORGN] AND CYTB[GENE] AND 1000:1200[SLEN]"), retmax = 1498, use_history = TRUE)

SALMONIDAE_search

Salm_cytb <- entrez_fetch(db = "nuccore", web_history = SALMONIDAE_search$web_history, rettype = "fasta")

class(Salm_cytb)

write(Salm_cytb, "Salm_cytb.fasta", sep = "\n", )

Salm_string <- readDNAStringSet("Salm_cytb.fasta")

Salm_string

class(Salm_string)

names(Salm_string)

df <- data.frame(Salm_Title = names(Salm_string), Salm_Sequence = paste(Salm_string))

df$Species <- word(df$Salm_Title, 2L, 3L)

df <- df[, c("Species", "Salm_Sequence")]

dim(df)

names(df)

summary(nchar(df$Salm_Sequence))

filtered_df <- df %>%
 filter(!grepl("N+", df$Salm_Sequence))

sum(str_count(filtered_df$Salm_Sequence, "-"))
sum(str_count(filtered_df$Salm_Sequence, "N"))

hist(nchar(filtered_df$Salm_Sequence))

sequence_data <- filtered_df %>%
  mutate(sequence_length = nchar(Salm_Sequence))

ggplot(sequence_data, aes(x = sequence_length, fill = Species)) + geom_histogram(binwidth = 25, position = "dodge")+ labs( x = "Sequence Length", y = "Frequency, title = Sequence Length Distribution by Species") + scale_fill_discrete(name = "Species")
       
less_than_1100 <- filtered_df[nchar(filtered_df$Salm_Sequence) < 1100, ]

more_than_1100 <- filtered_df[nchar(filtered_df$Salm_Sequence) > 1100, ]

common_names <- intersect(less_than_1100$Species, more_than_1100$Species)

common_names 

filtered_Salm <- filtered_df[nchar(df$Salm_Sequence) >= 1140, ]
filtered_Salm <- filtered_Salm[nchar(filtered_Salm$Salm_Sequence) <= 1170, ]



filtered_species <- filtered_Salm[!duplicated(filtered_Salm$Species), ]

cleaned_species <- na.omit(filtered_species)
cleaned_Salm <- na.omit(filtered_Salm)

sum(is.na(cleaned_species))
sum(is.na(cleaned_Salm))

##alignment and phylogeny creation using single species

summary(cleaned_Salm)
summary(cleaned_species)

class(cleaned_Salm)
class(cleaned_species)

Salmstringset <- DNAStringSet(cleaned_species$Salm_Sequence)

names(Salmstringset) <- cleaned_species$Species
class(Salmstringset)

Salm_alignment <- DNAStringSet(muscle::muscle(Salmstringset, quiet = TRUE, maxiters = 5, gapopen = -1000), use.names = TRUE)

Salm_alignment

BrowseSeqs(Salm_alignment)


Salm_Bin <- as.DNAbin(Salm_alignment)
Salm_Dist <- as.dist(Salm_Matrix)
Salm_hclust <- hclust(dist(Salm_Dist))

Salm_hclust

plot(Salm_hclust)


treeUPGMA = upgma(Salm_Dist)
treeNJ = NJ(Salm_Dist)
layout(matrix(c(1,2)))
par(mar = c(.1,.1,.1,.1), nrow = 1)      
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "fan", main="NJ")

treeUPGMA <- as.phylo(treeUPGMA)

Salm_phydat <- as.phyDat(Salm_Bin)

model_list <- c("JC", "F81", "K80", "HKY", "GTR")

likelihoods <- sapply(model_list, function(model) { fit <- pml(Salm_Bin, tree, k = 4, model = model)
-fit$logLik})

class(tree)
class(treeNJ)
class(treeUPGMA)



###CODING PT 2

Salm_BOLD <- read_tsv("http://boldsystems.org/index.php/API_Public/combined?taxon=Salmonidae&format=tsv")

