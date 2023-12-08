# This R script is intended to show examples for handling spatial data in R. We will use various kinds of data and geometries like polygons, point occurrences and raster files. The objective is to present some of the most commonly used packages and geospatial operations. Some examples are from the Wiki page of the epm R package that can be found in (https://github.com/ptitle/epm/wiki)

#Author:Jessica Castellanos-Labarcena (jcaste01@uoguelph.ca)
#Most updated version: November 19, 2023.

#First we are going to start by creating a grid cell. 
#The R package epm (short for EcoPhyloMapper) allows one to combine individual species' geographic ranges into a self-contained R object that includes a full accounting of which species occur in each grid cell. The epmGrid object also has the ability to contain morphological and phylogenetic data.

#Title, P. O., Swiderski, D. L., & Zelditch, M. L. (2022). EcoPhyloMapper: An r package for integrating geographical ranges, phylogeny and morphology. Methods in Ecology and Evolution, 13, 1912–1922. https://doi.org/10.1111/2041-210X.13914

#installing required packages and loading the libraries. 

# sf library is the workhorse library for all spatial vector data (points, polygons, lines)
install.packages("sf")
library(sf)
install.packages("epm")
library(epm)
install.packages("geodata")
library(geodata)
install.packages("terra")
library(terra)
install.packages("raster")
library(raster)

#The first step will be to download the data. We will download species range polygons from the IUCN Spatial Data repository. The IUCN Red List of Threatened Species contains global assessments for more than 150300 species. More than 82% of these (>123600 species) have spatial data. The spatial data provided are mostly for comprehensively assessed taxonomic groups and selected freshwater groups.

#Go to the following site and download the TERRESTRIAL MAMMALS polygon https://www.iucnredlist.org/resources/spatial-data-download. It is necessary to unzip the downloaded folder and leave the extracted files in the folder of the original name. Although we are reading the .shp file, all the other files in the folder are needed to read the data correctly. Set your working directory to the folder where the data is located, and then read the data. I am storing the path in a variable for then reading. This path contains the location of the .shp file. 
# filename and location of IUCN mammals shapefile.
IUCNfile <- "/Users/jessica/Downloads/data_Goespatial /MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp"

# Load the shapefile as a simple features object. Here we are using the function st_read. The function read simple features from file or database, or retrieve layer names and their geometry type(s). 

mammals <- st_read(IUCNfile, stringsAsFactors = FALSE)
mammals

# 1. Geographic data. 
#Here, we will extract the squirrels from the IUCN mammals dataset and create a list, where each list element is the spatial data for a species. Importantly, the names of the list are the names of the species.

# We can see that there is an associated data table, and a geometry column.
# Each row in this table is a geometric feature with associated attributes (the species, whether this is the breeding range, other metadata...)

# We can actually separate the data table from the geographic features

## Notice that the geometry column has disappeared
head(st_drop_geometry(mammals))

# Conversely, we can pull out the geographic features, separating them from the attribute data
st_geometry(mammals)

#As we are focusing on North American squirrels, we will identify all the relevant species using a bit of regular expression. Here, any taxon names that contains any of the listed terms will be returned (the | implies "or").

allsp <- unique(mammals$sci_name)

# Because there is a data table, we can select parts of it as we would in R with dataframes.

squirrelSp <- grep("Spermophilus\\s|Citellus\\s|Tamias\\s|Sciurus\\s|Glaucomys\\s|Marmota\\s|Cynomys\\s", allsp, value = TRUE, ignore.case = TRUE)

head(squirrelSp)

#The mammal ranges are currently all combined into a multipolygon object. We will now extract the squirrel species and store them as a list of separate species range polygons.

spList <- vector('list', length(squirrelSp))
names(spList) <- squirrelSp

for (i in 1:length(squirrelSp)) {
  	ind <- which(mammals$sci_name == squirrelSp[i])
    spList[[i]] <- mammals[ind,]
}

#These range polygons are unprojected, in longitude/latitude. We would like to work with equal area grid cells so that we will transform these range polygons to an equal area projection. The North America Albers EqualArea projection works well for North America, so we will use the following coordinate system definition:
#You can find more information about Coordinate Reference Systems here (https://docs.qgis.org/3.28/en/docs/gentle_gis_introduction/coordinate_reference_systems.html#:~:text=Geographic%20Coordinate%20Systems,popular%20is%20called%20WGS%2084.)

EAproj <- '+proj=aea +lat_1x=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

# In this line we are transforming the CRS for all the species in the list. We are iterating over all the elements in the list and applying the st_transform function. 
spListEA <- lapply(spList, function(x) st_transform(x, st_crs(EAproj)))

## Here, we will subset all rows where the genus is Neotamias
neotamias <- mammals[mammals$genus == 'Neotamias', ]

# Let's pull out a specific species
neotamias_dorsalis <- mammals[mammals$sci_name == 'Neotamias dorsalis',]

# Let's plot just the polygon of the species
plot(st_geometry(neotamias_dorsalis), col = 'orange')

# We need some context. Here, we will load a world landmass dataset. The folder with the data should be unzipped in your working directory. 
land <- st_read('ne_50m_land/ne_50m_land.shp')

plot(st_geometry(neotamias_dorsalis), col = 'white', border = NA)
plot(st_geometry(land), add = TRUE, col = 'pale green', lwd = 0.5)
plot(neotamias_dorsalis, col = 'orange', lwd = 0.5, add = TRUE)

#Here we just checked that the polygon is plotted correctly in the land, continental area. 

#The function createEPMgrid() will take range polygons and create an epmGrid object. We will set the resolution to 50km, we will opt to retain small-ranged species that would otherwise be lost, we will use hexagonal cells. This function intersects the range polygons with the grid cell and counts each species recorded in each cell. 

extentPoly <- "POLYGON ((-2906193 4690015, -3110223 4376122, -3377032 4172092, -3769399 3873894, -4083291 3340276, -4067597 2728184, -2702163 -1509370, 782048.5 -4036208, 1425529 -3973429, 1802200 -3785093, 2163177 -2639384, 3622779 -2027293, 2900826 5396274, 578018.1 5207938, -2906193 4690015))"

#Creates an epmGrid object from a range of species-specific inputs.epmGrid objects can be generated either from species geographic range polygons or from point occurrence records. An epmGrid object is a list that contains the following components:
# 1.a grid made up of either hexagonal or square cells
# 2. a list of the unique species communities found in these grid cells
# 3. an indexing vector that links each unique species community to the appropriate grid cells
# 4. a vector of the unique species found throughout the cell communities
# 5. the counts of cells that make up each species' geographic range (for calculating range-weighted metrics)
# 6. ecological / morphological data (if included)
# 7. a phylogeny (if included)

squirrelEPM <- createEPMgrid(spListEA, resolution = 50000, retainSmallRanges = TRUE, extent = extentPoly, method = 'centroid', cellType = 'hex')

#Let's plot it 

plot(squirrelEPM)

#Let's also create an epmGrid object with square grid cells for comparison. You may see the Warning message Failed to compute min/max, no valid pixels found in sampling. (GDAL error 1) . This just means that a species did not register in any grid cells. If you specified retainSmallRanges = TRUE, then those species will be included in a subsequent step. Therefore, this message can be ignored

squareEPM <- createEPMgrid(spListEA, resolution = 50000, retainSmallRanges = TRUE, extent = extentPoly, method = 'centroid', cellType = 'square')

plot(squareEPM)

#We can also see how a particular species is encoded in the epmGrid object. This can be helpful for confirming that this worked as expected.

map1 <- plotSpRange(squirrelEPM, taxon = "Neotamias dorsalis", lwd = 0.2)
map2 <- plotSpRange(squareEPM, taxon = "Neotamias dorsalis")

tmap::tmap_arrange(map1, map2)

#Now that we've created our epmGrid object, we can add additional data, such as trait data and/or a phylogeny. By using the grid cell geometry, we can retrieve more data and add it to the grid object. 

#Climatic data. In this section we are going to see how to retrieve climatic data from the WorldClim repository by using grographic coordinates. 

#Downloading  BioClimatic data from WorldClim website (https://www.worldclim.org/data/bioclim.html)
#in the path argument write the path to the folder where the data will be saved. This should be the working directory. 

bioCl <- worldclim_global(var = "bio", res = 5, path = "/Users/jessica/Downloads/data_Goespatial ")

#Now, we need the geographic coordinates to retrieve climatic data. 
#We will use the grid object part of the epm object and retrieve the coordinates for each cell. 

squirrel_grid <- squirrelEPM [[1]]

#get the centroid values for each cell. We will use the centroids to get the climatic values. What is happening here is that, using the coordinates, e are extracting the values of the raster files with the climatic variables. 

templateCentroids <- st_centroid(sf::st_geometry(squirrel_grid$gridTemplate))

#transform to sf object by retrieving coordinates in matrix form
coord_grid <- sf::st_coordinates(templateCentroids)

#create the spatial points of the grid cells centroids using the same projection of the epm object. This is important because the points must be in the correct projection to extract the values. 

points<- SpatialPoints(coord_grid, proj4string = CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

#spTransform for map projection: this is creating a spatial object with the points to extract the values.   

points <- spTransform(points, CRS("+proj=longlat +datum=WGS84 +no_defs"))

#Finally, transform the object to a SpatVector to use with the terra package. 

points <- vect(points)

#extract values from the raster object.
values_bioCl <- terra::extract(bioCl, points)

#Here, we got the climatic values for each cell for the 19 bioclimatic variables. You might notice that for some of the cells no values are retrieved. That is most likely related to cells close to the coastline and that we are using the centroid value for the cell. 


