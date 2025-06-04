
# FieldingCows ####

rm(list = ls())

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "FieldingCows"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Importing data ####

Files <- 
  Root %>% 
  list.files(full.names = T,
             pattern = ".rds")

FileList <- 
  Files %>% map(readRDS)

Contacts <- FileList[[2]]

Locations <- FileList[[1]] %>% rename_all(~str_remove_all(.x, "min_adj_"))

# Locations %<>% mutate_at("cow_id", ~str_remove_all(.x, "[a-zA-Z]") %>% 
#                           as.numeric)
# 
# Cows %<>% mutate_at("GPS", ~str_remove_all(.x, "[a-zA-Z]") %>% 
#                       as.numeric)

Locations %<>% 
  rename_all(~str_replace_all(.x, " ", "_")) %>% 
  mutate(X = longitude*100, Y = latitude*1000) %>% 
  mutate(X = X - min(X), Y = Y - min(Y)) # Remember to put it in kilometres!!!!!!!!

FocalSites <- Locations %>% group_by(farm_id) %>% summarise(N = nunique(cow_id)) %>% 
  filter(N > 50) %>% pull(farm_id)

Locations %<>% filter(farm_id %in% FocalSites)
Contacts %<>% filter(farm_id %in% FocalSites)

# Social Phase ####

AMList <- list()

Contacts$Site <- Contacts$farm_id

GraphList <- 
  Contacts$farm_id %>% unique %>% sort %>% 
  map(function(a){
    
    Contacts %>% 
      filter(farm_id == a) %>% 
      dplyr::select(id1, id2) %>%
      # dplyr::select(cow_id1, cow_id2) %>% 
      # as.matrix %>% 
      graph_from_data_frame(directed = F) %>% 
      as_tbl_graph
    
  })

AdjList <- GraphList %>% 
  map(get.adjacency) %>% 
  map(function(a){
    
    a <- as.matrix(a)
    
    diag(a) <- 0
    
    Observations <- rowSums(a)
    
    diag(a) <- Observations
    
    a %>% return
    
  })

AdjList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

# Deriving stuff ####

Censuses <- Locations# %>% left_join(Cows %>% rename(Id = prox_id))

# Making data frame

Censuses$Site <- Censuses$farm_id

Censuses %<>% 
  semi_join(Contacts, by = c("cow_id" = "id1", "Site")) %>% 
  semi_join(Contacts, by = c("cow_id" = "id2", "Site"))

ObservationDF <- Censuses %>% 
  # dplyr::select(-Site) %>% 
  dplyr::select(Id = cow_id, X, Y, Site) %>% 
  na.omit

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

FocalSites <- 
  ObservationDF$Site %>% unique %>% sort

Censuses %>% group_by(Site) %>% count %>% filter(n>50) %>% 
  dplyr::select(-c("n")) -> FocalSubdivisions

1:nrow(FocalSubdivisions) %>% 
  map(~ObservationDF %>% semi_join(FocalSubdivisions[.x,], by = c("Site"))) -> 
  GroupList

GroupList %>% bind_rows() %>% nrow

# Spatial Phase ####

# Making Centroids ####

ObservationDF %>% group_by(Site, Id) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  AnnualCentroids

AnnualCentroids %>% group_by(Id) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  LifetimeCentroids

ObservationDF %>% group_by(Site, Id) %>% count(name = "Obs") %>% 
  left_join(AnnualCentroids, .) -> AnnualCentroids

ObservationDF %>% group_by(Id) %>% count(name = "Obs") %>% 
  left_join(LifetimeCentroids, .) -> LifetimeCentroids

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

LifetimeKUDL <- kernelUD(SPDF, 
                         extent = 0,
                         same4all = TRUE, 
                         grid = 500)

KUDLRaster <- LifetimeKUDL %>% raster::raster()

KUDLRaster <- raster::resample(KUDLRaster, BlankRaster)

raster::values(KUDLRaster) <- 
  raster::values(KUDLRaster)*nrow(LifetimeCentroids)/
  sum(raster::values(KUDLRaster), na.rm = T)

AnnualCentroids %>% 
  arrange(Site) %>% pull(Site) %>% 
  unique %>% as.character -> FocalSites

AnnualCentroids %>% group_by(Site) %>% count %>% 
  # filter(n>50) %>% 
  pull(Site) %>% intersect(FocalSites) ->
  FocalSites

SubCentroids <- AnnualCentroids %>% filter(Site %in% FocalSites)

SPDF <- SpatialPointsDataFrame(data = SubCentroids[,c("X", "Y", "Site")], 
                               coords = SubCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Site"]

KUDL <- kernelUD(SPDF, 
                 extent = 0,
                 same4all = TRUE, 
                 grid = 500)

1:length(FocalSites) %>% lapply(function(a){
  
  print(FocalSites[a])
  
  DF <- AnnualCentroids %>% filter(Site == FocalSites[a])
  
  KUDL2 <- KUDL[[FocalSites[a]]]
  
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

AMList <- AdjList

names(AMList) <- FocalSubdivisions$Site

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Site") %>% 
  mutate_at("Id", as.numeric) %>% 
  mutate_at("Site", as.numeric) -> SocialNodeTraits

NodeDF %<>% left_join(SocialNodeTraits, by = c("Id", "Site"))

NodeDF$Associations

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

LifetimeDistance <-
  LifetimeCentroids[,c("X", "Y")] %>% 
  GristMatrix(LifetimeCentroids$Id)

FocalSites %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Site == a) 
  
  SubMat <- SubDF %>% 
    dplyr::select(X, Y) %>% 
    GristMatrix(SubDF$Id)
  
  SubMat
  
}) -> AnnualDistances

# Making Space Use Matrices ####

# Annual

HROList <- AreaList <- MCPListList <- list()

for(i in FocalSites){
  
  print(i)
  
  MCPIndividuals <- ObservationDF %>% 
    filter(Site == i) %>% dplyr::select(X, Y, Id) %>% 
    unique %>% 
    group_by(Id) %>% 
    summarise(Obs = n(), 
              X = mean(X, na.rm = T), 
              Y = mean(Y, na.rm = T)) %>% 
    filter(Obs > 5) %>% 
    pull(Id)
  
  if(length(MCPIndividuals)>0){
    
    MCPIndividuals %>% 
      
      lapply(function(a){
        
        print(a)
        ObservationDF %>% 
          filter(Id == a, Site == i) %>% 
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
    
    Dyads %>% 
      t %>% data.frame %>% 
      map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas) %>% as.character], .x %>% unlist %>% as.character, Symmetrical = T)) %>% 
      unlist %>% 
      matrix(nrow = length(MCPIndividuals %>% setdiff(ZeroAreas))) ->
      AnnualHRO
    
    dimnames(AnnualHRO) <- list(MCPIndividuals %>% setdiff(ZeroAreas), 
                                MCPIndividuals %>% setdiff(ZeroAreas))
    
    HROList[[i]] <- AnnualHRO
    
  }else{
    
    HROList[[i]] <- list()
    AreaList[[i]] <- list()
    
  }
}

names(MCPList) <-
  names(AreaList) <- 
  names(HROList) <- FocalSites

saveRDS(MCPListList, glue("Datasets/{FocalSystem}/Intermediate/MCPList.rds"))

FullSpatialList <-
  list(AnnualDistance = AnnualDistances,
       AnnualHRO = HROList,
       AnnualAreas = AreaList)

FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

AreaList <- AreaList[which((AreaList %>% map(length))>0)]

if(length(AreaList)>0){
  
  NodeDF <-
    AreaList %>% bind_rows(.id = "Site") %>% mutate_at(c("Site", "Id"), as.numeric) %>% 
    full_join(NodeDF %>% mutate_at(c("Site", "Id"), as.numeric), ., by = c("Site", "Id"))
  
}

NodeDF %>% 
  dplyr::select(Site, Id)
# unique %>%
nrow

NodeDF %<>% 
  mutate(Year = "Year", Site_Year = paste0(Site, "_Year"))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <- readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))
