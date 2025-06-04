
# Meerkats ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs)

theme_set(theme_cowplot())

FocalSystem <- "KalahariMeerkats"

# Importing and Cleaning ####

Root <- paste0("Datasets/",FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

Burrows <- read.csv(paste0(Root, "/Data/tblBurrows.csv"))

# Burrows %>% ggplot(aes(Longitude, Latitude)) + geom_point() + coord_sf()# +
# facet_wrap(~GroupRef)

# Burrows %>% 
#   filter(GroupRef == 3) %>% 
#   ggplot(aes(Longitude, Latitude)) + geom_point(aes(colour = )) + coord_sf()# +
# # facet_wrap(~GroupRef)

# Burrows %>% ggplot(aes(Easting, Northing)) + geom_point()

Burrows %<>% mutate_at("FirstSeenDate", 
                       ~ifelse(.x == "", "1900-01-01", .x) %>% 
                         lubridate::ymd())

Burrows %<>% arrange(GroupRef, FirstSeenDate)

Censuses <- read.csv(paste0(Root, "/Data/tblIndividualSightings.csv"))

Censuses %<>% filter(!is.na(GroupRef))

Censuses %<>% mutate_at("SeenDate", ~lubridate::ymd(.x))

Censuses[,c("X", "Y")] <- NA

Censuses %>% 
  separate(SeenDate, "-", into = c("Year", "Month", "Day")) %>% 
  dplyr::select(c("Day", "Month", "Year")) %>% 
  mutate_all(as.numeric) ->
  Censuses[,c("Day", "Month", "Year")]

Censuses %<>% mutate_at("SeenDate", as.numeric)
Burrows %<>% mutate_at("FirstSeenDate", as.numeric)

# Burrows %<>% dpplyr::select(FirstSeenDate, GroupRef, Latitude, Longitude) %>% na.omit

# Non Loop ####

DateDF <- 
  Censuses %>% 
  group_by(GroupRef) %>% 
  summarise(SeenDate = list(unique(sort(SeenDate)))) %>% 
  unnest(SeenDate) %>% 
  as.data.frame

DateDF[,c("X", "Y")] <- NA

# Loop ####

i <- 1

for(i in i:nrow(DateDF)){
  
  print(i)
  
  # Replacement <-
  #   Burrows %>% filter(GroupRef == DateDF[i, c("GroupRef")],
  #                      FirstSeenDate < DateDF[i, c("SeenDate")]) %>%
  #   slice(n()) %>%
  #   dplyr::select("Longitude", "Latitude")
  
  Replacement <-
    Burrows[last(which(Burrows$GroupRef == DateDF[i, c("GroupRef")] &
                         Burrows$FirstSeenDate < DateDF[i, c("SeenDate")])),
            c("Longitude", "Latitude")]
  
  if(nrow(Replacement) > 0){
    
    DateDF[i, c("X", "Y")] <-
      Replacement[nrow(Replacement),]
    
  }
  
}

Censuses %<>% 
  dplyr::select(-c(X, Y)) %>% 
  left_join(DateDF, by = c("GroupRef", "SeenDate"))

# Censuses %>% RandomSlice(10000) %>% ggplot(aes(X, Y)) + geom_point()

# Old Loop ####

# i <- 1
# 
# for(i in i:nrow(Censuses)){
#   
#   print(i)
#   
#   # Replacement <- 
#   #   Burrows %>% filter(GroupRef == Censuses[i, c("GroupRef")],
#   #                      FirstSeenDate < Censuses[i, c("SeenDate")]) %>% 
#   #   slice(n()) %>% 
#   #   dplyr::select("Longitude", "Latitude")
#   
#   Replacement <- 
#     Burrows[Burrows$GroupRef == Censuses[i, c("GroupRef")] &
#               Burrows$FirstSeenDate < Censuses[i, c("SeenDate")],
#             c("Longitude", "Latitude")]
#   
#   if(nrow(Replacement) > 0){
#     
#     Censuses[i, c("X", "Y")] <- 
#       Replacement[nrow(Replacement),]
#     
#   }
#   
# }

Censuses %>% 
  filter(Year > 1996, Year < 2019) %>% 
  arrange(Year) %>% 
  pull(Year) %>% unique %>% 
  sort %>% as.character -> FocalYears

FocalYears %<>% #c(min(as.numeric(FocalYears))-1, .) %>% 
  as.numeric

# FocalYears <- FocalYears[2:length(FocalYears)]

Records = 5

Censuses %<>% mutate(GroupDate = paste0(GroupRef, ",", SeenDate))

FocalYears %>% map(~Censuses %>% filter(Year == .x)) -> GroupList

ObservationDF <- 
  Censuses %>% 
  # dplyr::select(-Year) %>% 
  dplyr::select(Id = IndividID, GroupDate, X, Y, Year) %>% 
  na.omit

ObservationDF %<>% mutate_at(c("X", "Y"), ~.x*100) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x))

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

# Social Phase ####

AMList <- list()

i <- 1

for(i in i:length(GroupList)){
  
  print(i)
  
  GroupList[[i]] -> Groups
  
  Groups %>% 
    dplyr::select(Id = IndividID, GroupDate) %>% unique %>% droplevels %>% 
    # unite("Assoc", Group, Date, sep = "_") %>%
    # mutate(Assoc = GroupDate)
    table() -> M1
  
  SocGraph <- graph_from_incidence_matrix(M1)
  
  Proj <- bipartite.projection(SocGraph)$proj1
  
  AM <- Proj %>% get.adjacency(attr = "weight") %>% as.matrix
  
  Observations = rowSums(M1)
  
  Matrix <- AM # %>% ProportionalMatrix(Observations)
  
  diag(Matrix) <- Observations
  
  AMList[[i]] <- Matrix
  
}

AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

# Spatial Phase ####

# Making Centroids ####

ObservationDF %>% group_by(Year, Id) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  AnnualCentroids

AnnualCentroids %>% group_by(Id) %>%
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  LifetimeCentroids

ObservationDF %>% group_by(Year, Id) %>% count(name = "Obs") %>% 
  left_join(AnnualCentroids, .) -> AnnualCentroids

ObservationDF %>% group_by(Id) %>% count(name = "Obs") %>%
  left_join(LifetimeCentroids, .) -> LifetimeCentroids

AnnualCentroids %<>% mutate_at(c("X", "Y"), ~.x / 10)

LifetimeCentroids %<>% mutate_at(c("X", "Y"), ~.x / 10)

# Making Densities 

SPDF <- SpatialPointsDataFrame(data = LifetimeCentroids[,c("X", "Y")], 
                               coords = LifetimeCentroids[,c("X", "Y")])

MetreDims <- 
  LifetimeCentroids[,c("X", "Y")] %>% 
  map_dbl(c(range, diff)) %>% 
  multiply_by(1000) %>% ceiling

LifetimeCentroids[,c("X", "Y")] %>% 
  rename_all(tolower) %>% raster::extent() %>% 
  raster::raster(ncols = MetreDims[["X"]], 
                 nrows = MetreDims[["Y"]]) -> BlankRaster

# LifetimeKUDL <- kernelUD(SPDF, 
#                          extent = 0,
#                          same4all = TRUE, 
#                          grid = 500)
# 
# KUDLRaster <- LifetimeKUDL %>% raster::raster()
# 
# KUDLRaster <- raster::resample(KUDLRaster, BlankRaster)
# 
# raster::values(KUDLRaster) <- 
#   raster::values(KUDLRaster)*nrow(LifetimeCentroids)/
#   sum(raster::values(KUDLRaster), na.rm = T)

AnnualCentroids %>% 
  filter(Year > 1996, Year < 2019) %>% 
  arrange(Year) %>% pull(Year) %>% unique %>% as.character -> FocalYears

# FocalYears %<>% c(min(as.numeric(FocalYears))-1, .)

SPDF <- SpatialPointsDataFrame(data = AnnualCentroids[,c("X", "Y", "Year")], 
                               coords = AnnualCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Year"]

KUDL <- kernelUD(SPDF, 
                 extent = 0,
                 same4all = TRUE, 
                 grid = 500)

1:length(FocalYears) %>% lapply(function(a){
  
  print(FocalYears[a])
  
  DF <- AnnualCentroids %>% filter(Year == FocalYears[a])
  
  KUDL2 <- KUDL[[FocalYears[a]]]
  
  KUDLRaster <- KUDL2 %>% raster::raster() 
  
  KUDLRaster <- raster::resample(KUDLRaster, BlankRaster)
  
  raster::values(KUDLRaster) <-
    raster::values(KUDLRaster)*nrow(DF)/
    sum(raster::values(KUDLRaster), na.rm = T)
  
  KUDLRaster %>% 
    raster::extract(DF[,c("X", "Y")]) ->
    DF$Density.Annual
  
  return(DF)
  
}) -> DensityList

DensityList %>% bind_rows -> AnnualCentroids

NodeDF <- 
  AnnualCentroids

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Adding Social nodes to Spatial nodes ####

names(AMList) <- as.character(FocalYears)

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% mutate_at("Year", as.numeric) -> SocialNodeTraits

NodeDF %<>% left_join(SocialNodeTraits %>% mutate_at("Id", as.numeric))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

# LifetimeDistance <-
#   LifetimeCentroids[,c("X", "Y")] %>% 
# GristMatrix(LifetimeCentroids$Id)

FocalYears %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Year == a) 
  
  SubMat <- SubDF %>% 
    dplyr::select(X, Y) %>% 
    GristMatrix(SubDF$Id)
  
  SubMat
  
}) -> AnnualDistances

# Making Space Use Matrices ####

# Annual

HROList <- AreaList <- MCPListList  <- list()

for(i in FocalYears){
  
  print(i)
  
  MCPIndividuals <- AnnualCentroids %>% 
    filter(Year == i) %>% 
    filter(Obs > 5) %>% 
    pull(Id)
  
  MCPIndividuals %>% 
    
    lapply(function(a){
      
      # print(a)
      ObservationDF %>% 
        filter(Id == a, Year == i) %>% 
        dplyr::select(X, Y) %>% na.omit %>% 
        SpatialPoints %>% mcp
      
    }) -> MCPList
  
  names(MCPList) <- MCPIndividuals
  
  MCPListList[[i]] <- MCPList
  
  AnnualAreas <- 
    MCPList %>% map_dbl(~.x$area) %>% as.data.frame() %>% 
    rownames_to_column("Id") %>% rename(Area = 2)
  
  AreaList[[i]] <- AnnualAreas
  
  ZeroAreas <- AnnualAreas %>% filter(Area == 0) %>% pull(Id)
  
  Dyads <- expand.grid(Id1 = MCPIndividuals %>% setdiff(ZeroAreas), 
                       Id2 = MCPIndividuals %>% setdiff(ZeroAreas))
  
  if(nrow(Dyads) > 0){
    
    Dyads %>% #t %>% data.frame %>% 
      # dplyr::select(1:2)
      # slice(1:1000) %>% 
      t %>% data.frame %>% 
      map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas) %>% as.character], .x %>% unlist %>% as.character, 
                      Symmetrical = T)) %>% 
      unlist %>% 
      matrix(nrow = length(MCPIndividuals %>% setdiff(ZeroAreas))) ->
      AnnualHRO
    
    dimnames(AnnualHRO) <- list(MCPIndividuals %>% setdiff(ZeroAreas), 
                                MCPIndividuals %>% setdiff(ZeroAreas))
    
    HROList[[i]] <- AnnualHRO
    
  }
  
}

MCPListList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/MCPList.rds"))

FullSpatialList <-
  list(AnnualDistance = AnnualDistances,
       AnnualHRO = HROList,
       AnnualAreas = AreaList)

FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

NodeDF <-
  AreaList %>% bind_rows(.id = "Year") %>% 
  mutate_at("Year", as.numeric) %>% 
  mutate_at("Id", as.numeric) %>% 
  full_join(NodeDF, ., by = c("Year", "Id"))

NodeDF %>% dplyr::select(Year, Id) %>% 
  # unique %>%
  nrow

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <- readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))
