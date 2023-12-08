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
par(mar = c(.5,.5,.5,.5), nrow = 1)      
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "fan", main="NJ")

treeUPGMA <- as.phylo(treeUPGMA)

Salm_phydat <- as.phyDat(Salm_Bin)

model_list <- c("JC", "F81", "K80", "HKY", "GTR")

##likelihoods <- sapply(model_list, function(model) { fit <- pml(Salm_Bin, tree, k = 4, model = model)-fit$logLik})

class(tree)
class(treeNJ)
class(treeUPGMA)



###CODING PT 2

library(geojsonio)
library(sf)
library(maps)
library(plotly)


Salm_BOLD <- read_tsv("http://boldsystems.org/index.php/API_Public/combined?taxon=Salmonidae&format=tsv")

BOLDdf <- Salm_BOLD[, c("species_name", "lat", "lon")]

BOLDdf <- na.omit(BOLDdf)

sum(is.na(BOLDdf))

colnames(BOLDdf) <- c("Species", "Latitude", "Longitude")

dfRightJoin <- cleaned_Salm %>%
  right_join(BOLDdf, by = "Species", relationship = "many-to-many")


world <- map_data("world")

ggplot() + geom_polygon(data = world, aes(x = long, y= lat, group = group), fill = "white", color = "black") +
  geom_point(data = dfRightJoin, aes(x = Longitude, y = Latitude, color = Species), size = 2) + scale_color_discrete(name = "Species") + labs(title = "Salmonidae Distribution Map of the World", xlab = "Longitude", ylab = "Latitude") + theme(legend.position = "bottom")
                        
canada <- geojson_read("https://raw.githubusercontent.com/codeforgermany/click_that_hood/main/public/data/canada.geojson", what = "sp")

canada <- st_as_sf(canada, crs = 4326)
species <- st_as_sf(dfRightJoin, coords = c("Longitude", "Latitude"), crs = 4326)

species <- sf::st_cast(species, "MULTIPOLYGONAL")

class(species)
names(species)

combined_data <- st_join(canada, species, join = st_intersects)

combined_data <- combined_data %>%
  group_by(name) %>%
  summarize(Species = n())

combined_data <- sf::st_cast(combined_data, "GEOMETRYCOLLECTION") %>%
  st_collection_extract("POLYGON")

choro_map <- combined_data %>%
  ggplot(aes(fill = Species, text = str_c(name, ": ", Species))) + 
  geom_sf() + scale_fill_scico("Salmonidae Occurences", breaks = c(0, 100, 200, 300, 400,500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)) + theme_void() + ggtitle("Species Occurences of Salmonidae spp. Across Canada") + theme(legend.key.size = unit(3, "cm"))

plot(choro_map)


interactive_map <- ggplotly(choro_map, tooltip = "text")

interactive_map

