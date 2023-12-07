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
library(msa)
library(phytools)

SALMONIDAE_search <- entrez_search(db = "nuccore", term = ("Salmonidae[ORGN] AND CYTB[GENE] AND 1000:1200[SLEN]"), retmax = 1498, use_history = TRUE)

SALMONIDAE_search

Salm_cytb <- entrez_fetch(db = "nuccore", web_history = SALMONIDAE_search$web_history, rettype = "fasta")

class(Salm_cytb)

head(Salm_cytb)

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

hist(nchar(df$Salm_Sequence))

filtered_df <- df %>%
 filter(!grepl("N+", df$Salm_Sequence))


hist(nchar(filtered_df$Salm_Sequence))

less_than_1100 <- filtered_df[nchar(filtered_df$Salm_Sequence) < 1100, ]

more_than_1100 <- filtered_df[nchar(filtered_df$Salm_Sequence) > 1100, ]

common_names <- intersect(less_than_1100$Species, more_than_1100$Species)


hist(nchar(more_than_1100$Salm_Sequence))

filtered_Salm <- filtered_df[nchar(df$Salm_Sequence) >= 1140, ]
filtered_Salm <- filtered_Salm[nchar(filtered_Salm$Salm_Sequence) <= 1145, ]

hist(nchar(filtered_Salm$Salm_Sequence))

filtered_species <- filtered_Salm[!duplicated(filtered_Salm$Species), ]

cleaned_species <- na.omit(filtered_species)
cleaned_Salm <- na.omit(filtered_Salm)

sum(is.na(cleaned_species))
sum(is.na(cleaned_Salm))
##alignment and phylogeny creation using single species

summary(cleaned_Salm)
summary(cleaned_species)

Salmstringset <- DNAStringSet(cleaned_species$Salm_Sequence)

names(Salmstringset) <- cleaned_species$Species

Salmstringset

Salm_alignment <- AlignSeqs(Salmstringset)
                                

BrowseSeqs(Salm_alignment)

AdjustAlignment(Salm_alignment, perfectMatch = 5, misMatch = 0, gapLetter = -3, gapOpening = -0.1, gapExtension = 0, substitutionMatrix = NULL, shiftPenalty = -0.2, threshold = 0.1, weight = 1, processors = 1)

BrowseSeqs(Salm_alignment)

##Phylogeny creation 

