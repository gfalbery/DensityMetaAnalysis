
# BeakedWhales ####

rm(list = ls())

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "BeakedWhales"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Importing data ####

Files <- 
  Root %>% 
  dir_ls(regex = ".xlsx")

SheetList <- 
  Root %>% 
  dir_ls(regex = ".xlsx") %>% 
  excel_sheets() %>% 
  purrr::set_names() %>% 
  map(~read_excel(Files, sheet = .x))

Whales <- 
  SheetList[[1]] %>%
  # mutate(Year = str_split(Date, "-") %>% map_chr(1)) %>% arrange(Date)
  rename_all(CamelConvert)

Whales %<>% 
  # filter(START_LAT < 57.35) %>% 
  filter(Latitude > 43.75)

Whales %<>% 
  # rename_all(tolower) %>% rename_all(CamelConvert) %>% 
  mutate(X = Longitude, Y = Latitude) %>% 
  mutate(X = X*100, Y = Y*100) %>% 
  mutate(X = X - min(X), Y = Y - min(Y)) # Remember to put it in kilometres!!!!!!!!

Whales %<>% as.data.frame

Whales %<>% #mutate(GroupDate = as.numeric(factor(paste(Date, X, Y))))
  rename(GroupDate = GroupID)

# Deriving stuff ####

Censuses <- 
  Whales

# Making data frame

ObservationDF <- Censuses %>% 
  # dplyr::select(-Year) %>% 
  dplyr::select(Id = Individual, GroupDate, X, Y, 
                Year) %>% 
  na.omit

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

FocalYears <- 
  ObservationDF$Year %>% unique %>% sort

Censuses %>% group_by(Year) %>% count %>% filter(n>50) %>% 
  dplyr::select(-c("n")) -> FocalSubdivisions

1:nrow(FocalSubdivisions) %>% 
  map(~ObservationDF %>% semi_join(FocalSubdivisions[.x,], by = c("Year"))) -> 
  GroupList

GroupList %>% bind_rows() %>% nrow

# Social Phase ####

AMList <- list()

i <- 1

for(i in i:length(GroupList)){
  
  print(i)
  
  GroupList[[i]] -> Groups
  
  Groups %>% 
    dplyr::select(Id, GroupDate) %>% #unique %>% 
    droplevels %>% 
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

# Making Densities 

SPDF <- SpatialPointsDataFrame(data = LifetimeCentroids[,c("X", "Y")], 
                               coords = LifetimeCentroids[,c("X", "Y")])

MetreDims <- 
  LifetimeCentroids[,c("X", "Y")] %>% 
  map_dbl(c(range, diff)) %>% 
  multiply_by(100) %>%
  ceiling

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
  arrange(Year) %>% pull(Year) %>% 
  unique %>% as.character -> FocalYears

# AnnualCentroids %>% group_by(Year) %>% count %>% 
#   filter(n>50) %>% 
#   pull(Year) %>% intersect(FocalYears) ->
#   FocalYears

SubCentroids <- AnnualCentroids %>% filter(Year %in% FocalYears)

SPDF <- SpatialPointsDataFrame(data = SubCentroids[,c("X", "Y", "Year")], 
                               coords = SubCentroids[,c("X", "Y")])

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

names(AMList) <- FocalSubdivisions$Year

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% mutate_at("Year", as.numeric) %>% 
  mutate_at("Id", as.numeric) -> SocialNodeTraits

NodeDF %<>% left_join(SocialNodeTraits, by = c("Id", "Year"))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

FocalYears %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Year == a) 
  
  SubMat <- SubDF %>% 
    dplyr::select(X, Y) %>% 
    GristMatrix(SubDF$Id)
  
  SubMat
  
}) -> AnnualDistances

# Making Space Use Matrices ####

# Annual

HROList <- AreaList <- MCPListList <- list()

for(i in FocalYears){
  
  print(i)
  
  MCPIndividuals <- ObservationDF %>% 
    filter(Year == i) %>% dplyr::select(X, Y, Id) %>% 
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
      
      Dyads %>% 
        t %>% data.frame %>% 
        map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas) %>% as.character], .x %>% unlist %>% as.character, Symmetrical = T)) %>% 
        unlist %>% 
        matrix(nrow = length(MCPIndividuals %>% setdiff(ZeroAreas))) ->
        AnnualHRO
      
      dimnames(AnnualHRO) <- list(MCPIndividuals %>% setdiff(ZeroAreas), 
                                  MCPIndividuals %>% setdiff(ZeroAreas))
      
      HROList[[i]] <- AnnualHRO
      
    }
    
  }else{
    
    HROList[[i]] <- list()
    AreaList[[i]] <- list()
    
  }
}

MCPListList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/MCPList.rds"))

FullSpatialList <-
  list(AnnualDistance = AnnualDistances,
       AnnualHRO = HROList,
       AnnualAreas = AreaList)

FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

AreaList <- AreaList[which((AreaList %>% map(length))>0)]

if(length(AreaList)<0){
  
  NodeDF <-
    AreaList %>% bind_rows(.id = "Year") %>% mutate_at(c("Year"), as.numeric) %>% 
    full_join(NodeDF, ., by = c("Year", "Id"))
  
}

NodeDF %>% dplyr::select(Year, Id) %>% 
  # unique %>%
  nrow

# NodeDF %<>% filter(Density.Annual < 2e-4)

NodeDF %<>% filter(!is.na(Associations))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))
