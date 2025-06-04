
# GermanBoar ####

rm(list = ls())

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "EuropeanBoar"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Importing data ####

Files <- 
  Root %>% 
  list.files(full.names = T)

# Root %>% dir_ls(regex = ".zip") %>% unzip(exdir = Root)

FileList <- Root %>% dir_ls(regex = ".csv$") %>% map(read.csv)

Contacts <- FileList[[1]]

Locations <- FileList[[2]] %>% 
  dplyr::select(-X) %>% 
  rename_all(CamelConvert) %>% 
  rename(X = Longitude, Y = Latitude)

Locations %<>% 
  rename_all(~str_replace_all(.x, " ", "_"))
Locations %<>% mutate(Year = str_split(Acquisition_time, "-") %>% map_chr(1))

Locations %<>% 
  group_by(Year, Study_area) %>% summarise_at("Animals_id", nunique) %>% 
  filter(Animals_id > 30) %>% 
  semi_join(Locations, ., by = c("Year", "Study_area"))

Contacts %<>% 
  rename_all(CamelConvert)

Contacts %<>% 
  inner_join(Locations %>% dplyr::select(ID1 = Animals_id, Study_area) %>% unique,
             by = c("ID1"))

Contacts %<>% 
  inner_join(Locations %>% dplyr::select(Study_area, Timegroup, Year) %>% unique,
             by = c("Study_area", "Timegroup"))

Contacts <- 
  Locations %>% 
  group_by(Year, Study_area) %>% summarise_at("Animals_id", nunique) %>% 
  filter(Animals_id > 30) %>% 
  semi_join(Contacts, ., by = c("Year", "Study_area"))

Locations %<>% 
  mutate(Site = Study_area) %>% 
  mutate(X = X*100, Y = Y*100) %>% 
  group_by(Site) %>% 
  mutate(X = X - min(X), Y = Y - min(Y))  # Remember to put it in kilometres!!!!!!!!

# Locations %<>% filter(Study_area == "HainichNP")

# Social Phase ####

AMList <- list()

Contacts$Site <- Contacts$Study_area

Contacts %<>% mutate(Site_Year = paste0(Site, "_", Year))

GraphList <- 
  Contacts$Site_Year %>% unique %>% sort %>% 
  map(function(a){
    
    Contacts %>% 
      filter(Site_Year == a) %>% 
      dplyr::select(ID1, ID2) %>%
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

Censuses$Site <- Censuses$Study_area

Censuses %<>% mutate(Site_Year = paste0(Site, "_", Year))

ObservationDF <- Censuses %>% 
  # dplyr::select(-Site) %>% 
  dplyr::select(Id = Animals_id, X, Y, Site, Year, Site_Year) %>% 
  na.omit

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

FocalSites <- 
  ObservationDF$Site_Year %>% unique %>% sort

Censuses %>% group_by(Site_Year) %>% count %>% filter(n>30) %>% 
  dplyr::select(-c("n")) -> FocalSubdivisions

1:nrow(FocalSubdivisions) %>% 
  map(~ObservationDF %>% semi_join(FocalSubdivisions[.x,], by = c("Site_Year"))) -> 
  GroupList

GroupList %>% bind_rows() %>% nrow

# Spatial Phase ####

# Making Centroids ####

ObservationDF %>% group_by(Site_Year, Id) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  AnnualCentroids

AnnualCentroids %>% group_by(Id) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
  LifetimeCentroids

ObservationDF %>% group_by(Site_Year, Id) %>% count(name = "Obs") %>% 
  left_join(AnnualCentroids, .) -> AnnualCentroids

ObservationDF %>% group_by(Id) %>% count(name = "Obs") %>% 
  left_join(LifetimeCentroids, .) -> LifetimeCentroids

# Making Densities 

SPDF <- SpatialPointsDataFrame(data = LifetimeCentroids[,c("X", "Y")], 
                               coords = LifetimeCentroids[,c("X", "Y")])

MetreDims <- 
  LifetimeCentroids[,c("X", "Y")] %>% 
  map_dbl(c(range, diff)) %>% 
  multiply_by(10) %>% ceiling

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
  arrange(Site_Year) %>% pull(Site_Year) %>% 
  unique %>% as.character -> FocalSites

AnnualCentroids %>% group_by(Site_Year) %>% count %>% 
  filter(n>30) %>%
  pull(Site_Year) %>% intersect(FocalSites) ->
  FocalSites

SubCentroids <- AnnualCentroids %>% filter(Site_Year %in% FocalSites)

SPDF <- SpatialPointsDataFrame(data = SubCentroids[,c("X", "Y", "Site_Year")], 
                               coords = SubCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Site_Year"]

KUDL <- kernelUD(SPDF, 
                 extent = 0,
                 same4all = TRUE, 
                 grid = 500)

1:length(FocalSites) %>% lapply(function(a){
  
  print(FocalSites[a])
  
  DF <- AnnualCentroids %>% filter(Site_Year == FocalSites[a])
  
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
  AnnualCentroids# %>% 
# left_join(LifetimeCentroids, by = c("Id"), 
#           suffix = c(".Annual", ".Lifetime"))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Adding Social nodes to Spatial nodes ####

AMList <- AdjList

names(AMList) <- FocalSubdivisions$Site_Year

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Site_Year") %>% 
  mutate_at("Id", as.numeric) -> SocialNodeTraits

NodeDF %<>% left_join(SocialNodeTraits, by = c("Id", "Site_Year"))

NodeDF$Associations %<>% replace_na(0)

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

LifetimeDistance <-
  LifetimeCentroids[,c("X", "Y")] %>% 
  GristMatrix(LifetimeCentroids$Id)

FocalSites %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Site_Year == a) 
  
  SubMat <- SubDF %>% 
    dplyr::select(X, Y) %>% 
    GristMatrix(SubDF$Id)
  
  SubMat
  
}) -> AnnualDistances

# Making Space Use Matrices ####

HROList <- AreaList <- MCPListList <- list()

for(i in FocalSites){
  
  print(i)
  
  MCPIndividuals <- ObservationDF %>% ungroup %>% 
    filter(Site_Year == i) %>% dplyr::select(X, Y, Id) %>% 
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
        ObservationDF %>% ungroup %>% 
          filter(Id == a, Site_Year == i) %>% 
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
      map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas) %>% as.character], 
                      .x %>% unlist %>% as.character, Symmetrical = T)) %>% 
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

names(MCPListList) <- 
  
  names(AreaList) <- 
  
  names(HROList) <- FocalSites

saveRDS(MCPListList, glue("Datasets/{FocalSystem}/Intermediate/MCPList.rds"))

FullSpatialList <-
  list(#LifetimeDistance = LifetimeDistance,
    AnnualDistance = AnnualDistances,
    # LifetimeHRO = LifetimeHRO,
    AnnualHRO = HROList,
    # LifetimeAreas = LifetimeAreas,
    AnnualAreas = AreaList)

FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

AreaList <- AreaList[which((AreaList %>% map(length))>0)]

if(length(AreaList)>0){
  
  NodeDF <-
    AreaList %>% bind_rows(.id = "Site_Year") %>% mutate_at("Id", as.numeric) %>% 
    # left_join(LifetimeAreas, suffix = c(".Annual", ".Lifetime"), by = c("Id")) %>% 
    full_join(NodeDF, ., by = c("Site_Year", "Id"))
  
}

NodeDF %<>% mutate_at(c("Associations", "Strength", "Degree"), ~replace_na(.x, 0))

NodeDF %>% dplyr::select(Site_Year, Id) %>% 
  unique %>%
  nrow

# NodeDF %<>% mutate(Observations = Obs)

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <- readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% ggplot(aes(Density.Annual, Associations)) + geom_point() + geom_smooth(method = lm)
