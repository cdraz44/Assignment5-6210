#####################CODING PART 1 #######################
######Data Acquisition#########
##Load all required packages for use in this script
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

##Search organism and fetch using web history from NCBI
SALMONIDAE_search <- entrez_search(db = "nuccore", term = ("Salmonidae[ORGN] AND CYTB[GENE] AND 1000:1200[SLEN]"), retmax = 1498, use_history = TRUE)

Salm_cytb <- entrez_fetch(db = "nuccore", web_history = SALMONIDAE_search$web_history, rettype = "fasta")
class(Salm_cytb)

##Write to fasta file in order to create stringset and data frame
write(Salm_cytb, "Salm_cytb.fasta", sep = "\n", )


Salm_string <- readDNAStringSet("Salm_cytb.fasta")

class(Salm_string)
names(Salm_string)

df <- data.frame(Salm_Title = names(Salm_string), Salm_Sequence = paste(Salm_string))

##Cleaning up names and data frame to only contain information I want, species and sequence
df$Species <- word(df$Salm_Title, 2L, 3L)
df <- df[, c("Species", "Salm_Sequence")]

dim(df)
names(df)

######Cleaning up data ##########

##Checking for extremes in my data
summary(nchar(df$Salm_Sequence))

##Removing any sequence with incomplete bases.  I want to have a very clean dataset for this alignment and phylogenetic analysis
filtered_df <- df %>%
  filter(!grepl("N+", df$Salm_Sequence))

##Ensuring I have clean data
sum(str_count(filtered_df$Salm_Sequence, "-"))
sum(str_count(filtered_df$Salm_Sequence, "N"))

##Viewing distribution of sequence lengths
hist(nchar(filtered_df$Salm_Sequence))

##Mutating to be used with ggplot
sequence_data <- filtered_df %>%
  mutate(sequence_length = nchar(Salm_Sequence))

##I want to visualize the distribution of sequence length with species, so I can gauge where to cut outliers out.  This histogram has colored bars by species showing distribution.  
ggplot(sequence_data, aes(x = sequence_length, fill = Species)) + geom_histogram(binwidth = 25, position = "dodge")+ labs( x = "Sequence Length", y = "Frequency, title = Sequence Length Distribution by Species") + scale_fill_discrete(name = "Species")

##More looking into my data, seeing which species are below 1100 bases and which are above/
less_than_1100 <- filtered_df[nchar(filtered_df$Salm_Sequence) < 1100, ]
more_than_1100 <- filtered_df[nchar(filtered_df$Salm_Sequence) > 1100, ]

common_names <- intersect(less_than_1100$Species, more_than_1100$Species)

common_names 

##Filtering out sequences containing less than 1140 and more than 1170 nucleotides.
filtered_Salm <- filtered_df[nchar(df$Salm_Sequence) >= 1140, ]
filtered_Salm <- filtered_Salm[nchar(filtered_Salm$Salm_Sequence) <= 1170, ]


##Creating a dataframe containing only one sequence per species for phylogenetic analysis
filtered_species <- filtered_Salm[!duplicated(filtered_Salm$Species), ]

##More clean up, removing any NAs from both dataframes
cleaned_species <- na.omit(filtered_species)
cleaned_Salm <- na.omit(filtered_Salm)
##Making sure this worked
sum(is.na(cleaned_species))
sum(is.na(cleaned_Salm))

######Alignment and Phylogenetic Analysis#######

##Checking summaries and class of dataframes before moving forward
summary(cleaned_Salm)
summary(cleaned_species)

class(cleaned_Salm)
class(cleaned_species)

##Creating my single sequence species data frame into a DNA string set
Salmstringset <- DNAStringSet(cleaned_species$Salm_Sequence)
class(Salmstringset)

##Setting names to carry through 
names(Salmstringset) <- cleaned_species$Species

##MSA using muscle, increasing iterations and changing gap penalties
Salm_alignment <- DNAStringSet(muscle::muscle(Salmstringset, quiet = TRUE, maxiters = 5, gapopen = -500), use.names = TRUE)

##Viewing my alignment
Salm_alignment
BrowseSeqs(Salm_alignment)

##Creating distance matrixes and DNA bins for further clustering analysis
Salm_Bin <- as.DNAbin(Salm_alignment)
Salm_Dist <- as.dist(Salm_Matrix)

##Creating a phylogenetic tree using the distance matrix
Salm_hclust <- hclust(dist(Salm_Dist))
plot(Salm_hclust)

##Creating two different phylogenetic trees and plotting them together
treeUPGMA = upgma(Salm_Dist)
treeNJ = NJ(Salm_Dist)
layout(matrix(c(1,2), nrow = 1))
par(mar = c(.1,.1,.1,.1))      
plot(treeUPGMA, main="UPGMA", cex = 0.5)
plot(treeNJ, "fan", main="NJ", cex = 0.7)


##Converting the tree to phylo format for phangorn analysis.  (This tree was already in phylo format, but for some reason phangorn is not recognizing any of my trees as class phylo even though they all are)
treeUPGMA <- as.phylo(treeUPGMA)


##Creating a phydat from my DNA bin
Salm_phydat <- as.phyDat(Salm_Bin)


##Setting models to run phangorn likelihood analysis
model_list <- c("JC", "F81", "K80", "HKY", "GTR")

##Running likelihood analysis, it did not work because it won't recognize any of my trees. I tried for hours and hours and eventually had to give up but I wasted so much time on making this my "novel component" I have left it in.
likelihoods <- sapply(model_list, function(model) { fit <- pml(Salm_Bin, tree, k = 4, model = model)
-fit$logLik})

##This is showing that all my trees are class phylo and I don't know why I kept getting the error "tree must be of class phylo"
class(tree)
class(treeNJ)
class(treeUPGMA)


##########CODING PART 2#############

