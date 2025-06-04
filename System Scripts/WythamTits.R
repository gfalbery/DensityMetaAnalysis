
# WythamTits ####

# https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork)
library(fs)

theme_set(theme_cowplot())

FocalSystem <- "WythamTits"

Root <- paste0("Datasets/",FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Importing and Cleaning ####

"Datasets/WythamTits/Input" %>% list.files(pattern = ".zip", full.names = T) %>% 
  unzip(exdir = "Datasets/WythamTits/Input")

"Datasets/WythamTits/Input/csvs for greg" %>% list.files(pattern = "csv", full.names = T) %>% 
  map(fread) -> TitList

TitList %>% bind_rows(.id = "Year") %>% 
  rename_all(CamelConvert) -> 
  FullTitDF

FullTitDF %>% dplyr::select(Year, ID = Id, Time, Location, X, Y) -> ObservationDF

ObservationDF %<>% arrange(Time)

ObservationDF %<>% 
  mutate_at(c("X", "Y"), ~.x/1000) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x))

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

ObservationDF %>% Grelect(X, Y, Location) %>% unique -> NestLocations

load(glue("Datasets/{FocalSystem}/Input/gt_group_information_for_simple_greg.RData"))

GroupList <- yr.gd.list %>% 
  map(~.x %>% dplyr::select(-c(x, y))) %>% 
  map(~unnest(.x, ids)) %>% 
  map(~rename_all(.x, CamelConvert)) %>% 
  map(~rename(.x, ID = Ids) %>% unite(Assoc, Time:Location))

# Social Phase ####

AMList <- list()

i <- 1

for(i in i:length(GroupList)){
  
  print(i)
  
  GroupList[[i]] -> Groups
  
  Groups %>% 
    dplyr::select(ID, Assoc) %>% unique %>% droplevels %>% 
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

ObservationDF %>% group_by(Year, ID) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  AnnualCentroids

AnnualCentroids %>% group_by(ID) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  LifetimeCentroids

ObservationDF %>% group_by(Year, ID) %>% count(name = "Obs") %>% 
  left_join(AnnualCentroids, .) -> AnnualCentroids

ObservationDF %>% group_by(ID) %>% count(name = "Obs") %>% 
  left_join(LifetimeCentroids, .) -> LifetimeCentroids

# Making Densities ####

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

LifetimeKUDL <- kernelUD(SPDF, 
                         extent = 0,
                         same4all = TRUE, 
                         grid = 500)

KUDLRaster <- LifetimeKUDL %>% raster::raster()

KUDLRaster <- raster::resample(KUDLRaster, BlankRaster)

raster::values(KUDLRaster) <- 
  raster::values(KUDLRaster)*nrow(LifetimeCentroids)/
  sum(raster::values(KUDLRaster), na.rm = T)

AnnualCentroids %>% arrange(Year) %>% pull(Year) %>% unique %>% as.character -> FocalYears

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
        rownames_to_column("ID")) %>% 
  bind_rows(.id = "Year") -> SocialNodeTraits

NodeDF %<>% left_join(SocialNodeTraits)

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

LifetimeDistance <-
  LifetimeCentroids[,c("X", "Y")] %>% 
  GristMatrix(LifetimeCentroids$ID)

FocalYears %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Year == a) 
  
  SubMat <- SubDF %>% 
    dplyr::select(X, Y) %>% 
    GristMatrix(SubDF$ID)
  
  SubMat
  
}) -> AnnualDistances

# Making Space Use Matrices ####

# Annual

HROList <- AreaList <- MCPListList <- list()

for(i in FocalYears){
  
  print(i)
  
  MCPIndividuals <- AnnualCentroids %>% 
    filter(Year == i) %>% 
    filter(Obs > 5) %>% 
    pull(ID)
  
  MCPIndividuals %>% 
    
    lapply(function(a){
      
      # print(a)
      ObservationDF %>% 
        filter(ID == a, Year == i) %>% 
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
  
  Dyads <- expand.grid(ID1 = MCPIndividuals %>% setdiff(ZeroAreas), 
                       ID2 = MCPIndividuals %>% setdiff(ZeroAreas))
  
  if(nrow(Dyads) > 0){
    
    Dyads %>% #t %>% data.frame %>% 
      # dplyr::select(1:2)
      # slice(1:1000) %>% 
      t %>% data.frame %>% 
      map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas) %>% as.character], .x %>% unlist %>% as.character, Symmetrical = T)) %>% 
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
  AreaList %>% bind_rows(.id = "Year") %>% #mutate_at("Year", as.numeric) %>% 
  full_join(NodeDF, ., by = c("Year", "ID" = "Id"))

NodeDF %>% dplyr::select(Year, ID) %>% 
  # unique %>%
  nrow

NodeDF %<>% mutate_at(vars(Associations:Observations), ~replace_na(.x, 0))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))
