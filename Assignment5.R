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

stats <- function(x) {list(class = class(x), names = if (is.data.frame(x)) names(x) else NULL, summary = if (is.numeric(x)) summary(x) else NULL)}


df <- data.frame(Salm_Title = names(Salm_string), Salm_Sequence = paste(Salm_string))

##Cleaning up names and data frame to only contain information I want, species and sequence

df$Species <- word(df$Salm_Title, 2L, 3L)
df <- df[, c("Species", "Salm_Sequence")]

stats(df)

######Cleaning up data ##########

##Checking for extremes in my data

stats(nchar(df$Salm_Sequence))

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

##More looking into my data, seeing which species are below 1100 bases and which are above

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

##Making sure this worked and all NAs removed

sum(is.na(cleaned_species))
sum(is.na(cleaned_Salm))

######Alignment and Phylogenetic Analysis#######

##Checking summaries and class of dataframes before moving forward

stats(cleaned_Salm)
stats(cleaned_species)
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

#THIS WAS THE CODE::::::
#likelihoods <- sapply(model_list, function(model) { fit <- pml(Salm_Bin, treeUPGMA, k = 4, model = model)-fit$logLik})

##This is showing that all my trees are class phylo and I don't know why I kept getting the error "tree must be of class phylo"

stats(treeNJ)
stats(treeUPGMA)

##########CODING PART 2#############

##Loading packages needed for this section

library(geojsonio)
library(sf)
library(maps)
library(plotly)

##Importing data from BOLD for Salmonidae

Salm_BOLD <- read_tsv("http://boldsystems.org/index.php/API_Public/combined?taxon=Salmonidae&format=tsv")

##Specifying data wanted in data frame from BOLD data

BOLDdf <- Salm_BOLD[, c("species_name", "lat", "lon")]

##Removing NAs from data frame 

BOLDdf <- na.omit(BOLDdf)

##Ensuring that NAs are removed

sum(is.na(BOLDdf))

##Setting column names of new data frame

colnames(BOLDdf) <- c("Species", "Latitude", "Longitude")

##Joining data frames, data frame containing sequences of all species from part 1, and BOLD data frame

dfRightJoin <- cleaned_Salm %>%
  right_join(BOLDdf, by = "Species", relationship = "many-to-many")

stats(dfRightJoin)

##Fetching a world map to plot over

world <- map_data("world")

##Plotting occurences of species onto a world map to visualize distribution

ggplot() + geom_polygon(data = world, aes(x = long, y= lat, group = group), fill = "white", color = "black") +
  geom_point(data = dfRightJoin, aes(x = Longitude, y = Latitude, color = Species), size = 2) + scale_color_discrete(name = "Species") + labs(title = "Salmonidae Distribution Map of the World", xlab = "Longitude", ylab = "Latitude") + theme(legend.position = "bottom")

##Fetching a geojson file which is a map of Canada and making it an sf object

canada <- geojson_read("https://raw.githubusercontent.com/codeforgermany/click_that_hood/main/public/data/canada.geojson", what = "sp")

## Setting these as sf objects to use for plotting choropleth 

canada <- st_as_sf(canada, crs = 4326)
species <- st_as_sf(dfRightJoin, coords = c("Longitude", "Latitude"), crs = 4326)

##Checking class to make sure I have the right type of object

stats(canada)
stats(species)

## Joining both sf data frames 

combined_data <- st_join(canada, species, join = st_intersects)

##Grouping by province name so I can plot species occurrences 
combined_data <- combined_data %>%
  group_by(name) %>%
  summarize(Species = n())

## There was an issue creating my interactive plot so I had to circle back and change the geometry of this object

combined_data <- sf::st_cast(combined_data, "GEOMETRYCOLLECTION") %>%
  st_collection_extract("POLYGON")

## Plotting choropleth plot

choro_map <- combined_data %>%
  ggplot(aes(fill = Species, text = str_c(name, ": ", Species))) + 
  geom_sf() + scale_fill_scico("Salmonidae Occurences", breaks = c(0, 100, 200, 300, 400,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)) + theme_void() + ggtitle("Species Occurences of Salmonidae spp. Across Canada") + theme(legend.key.size = unit(3, "cm"))

plot(choro_map)

## Plotting interactive map of same data, just building on the previous plot 

interactive_map <- ggplotly(choro_map, tooltip = "text")

interactive_map


