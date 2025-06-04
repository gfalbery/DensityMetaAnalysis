
# AcornWoodpeckers ####

# rm(list = ls())

# STILL NEEDS RERUNNING #####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "AcornWoodpeckers"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Importing data ####

# Root %>% 
#   list.files(full.names = T,
#              pattern = ".zip") %>% 
#   unzip(exdir = Root)

FileList <- 
  Root %>% 
  dir_ls(regex = ".csv") %>% 
  map(read.csv) %>% 
  map(~rename_all(.x, CamelConvert))

Birds <- FileList[[1]]
Locations <- FileList[[2]]

Birds %<>% left_join(Locations, by = c("GROUP"))

Birds %<>% rename_all(tolower) %>%  rename_all(CamelConvert)

Birds %<>% mutate(Obs.id = paste(Group, Year, Month, Day, Bin, sep = "_"))

Contacts <- Birds

# Deriving stuff ####

FocalYears <- 
  Contacts$Year %>% unique %>% sort

Contacts %>% group_by(Year) %>% count %>% filter(n>50) %>% 
  dplyr::select(-c("n")) -> FocalSubdivisions

1:nrow(FocalSubdivisions) %>% 
  map(~Contacts %>% semi_join(FocalSubdivisions[.x,], by = c("Year"))) -> 
  GroupList

GroupList %>% bind_rows() %>% nrow

GroupList %>% map(nrow)

# Social Phase ####

AMList <- list()

i <- 1

# for(i in i:length(GroupList)){
# 
#   # print(i)
# 
#   GroupList[[i]] -> Groups
# 
#   FocalBirds <- Groups %>% count(Bird) %>% pull(Bird)
# 
#   FocalGroups <- Groups %>% count(Obs.id) %>% filter(n > 1) %>% pull(Obs.id)
# 
#   Groups %>%
#     dplyr::select(GroupDate = Obs.id, Id = Bird) %>% #unique %>%
#     # filter(Bird %in% FocalBirds) %>%
#     # filter(GroupDate %in% FocalGroups) %>%
#     droplevels %>%
#     # unite("Assoc", Group, Date, sep = "_") %>%
#     # mutate(Assoc = GroupDate)
#     table() -> M1
# 
#   Matrix <- matrix(0, length(FocalBirds), length(FocalBirds))
# 
#   dimnames(Matrix) <- list(FocalBirds, FocalBirds)
# 
#   LongDF <- Matrix %>% reshape2::melt()
# 
#   LongDF %>% nrow %>% print
# 
#   x <- 1
# 
#   for(x in 1:nrow(LongDF)){
# 
#     print(x)
# 
#     LongDF[x, "value"] <-
#       M1[,LongDF[x, c("Var1", "Var2")] %>% unlist] %>% rowSums %>% as.data.frame() %>%
#       rename(Value = 1) %>% filter(Value > 1) %>% nrow
# 
#   }
# 
#   # Observations = Groups %>% count(Bird) %>% pull(n)
#   #
#   # Matrix <- matrix(0, length(FocalBirds), length(FocalBirds))
#   #
#   # dimnames(Matrix) <- list(FocalBirds, FocalBirds)
#   #
#   # Matrix[colnames(AM), colnames(AM)] <- AM
#   #
#   # diag(Matrix) <- Observations
# 
#   Matrix[] <- LongDF$value
# 
#   AMList[[i]] <- Matrix
# 
# }


for(i in i:length(GroupList)){
  
  # print(i)
  
  GroupList[[i]] -> Groups
  
  FocalBirds <- Groups %>% count(Bird) %>% pull(Bird)
  
  FocalGroups <- Groups %>% count(Obs.id) %>% filter(n > 1) %>% pull(Obs.id)
  
  Matrix <- matrix(0, length(FocalBirds), length(FocalBirds))
  
  dimnames(Matrix) <- list(FocalBirds, FocalBirds)
  
  LongDF2 <- Matrix %>% reshape2::melt()
  
  LongDF2 %>% nrow %>% print
  
  x <- 1
  
  for(x in x:nrow(LongDF2)){
    
    print(x)
    
    LongDF2[x, "value"] <-
      Groups %>% filter(Bird %in% (LongDF2[x, c("Var1", "Var2")] %>% unlist)) %>%
      count(Obs.id) %>%
      filter(n > 1) %>%
      nrow
    
  }
  
  Matrix[] <- LongDF2$value
  
  Observations = Groups %>% count(Bird) %>% pull(n)
  
  diag(Matrix) <- Observations
  
  AMList[[i]] <- Matrix
  
}

AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

# Spatial Phase ####

# Making Centroids ####

ObservationDF <- 
  Birds %>% rename(Id = Bird) %>% 
  rename(X = Easting, Y = Northing) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x, na.rm = T)) %>% 
  mutate_at(c("X", "Y"), ~.x/1000) %>%
  # dplyr::select(-contains("euclidean")) %>% 
  filter(!is.na(X), !is.na(Y))

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

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

FocalYears <- AnnualCentroids$Year %>% unique %>% sort

# Making Densities

MetreDims <-
  AnnualCentroids[,c("X", "Y")] %>%
  map_dbl(c(range, diff)) %>%
  multiply_by(1000) %>% ceiling

AnnualCentroids[,c("X", "Y")] %>%
  rename_all(tolower) %>% raster::extent() %>%
  raster::raster(ncols = MetreDims[["X"]],
                 nrows = MetreDims[["Y"]]) -> BlankRaster

AnnualCentroids %>%
  arrange(Year) %>% pull(Year) %>%
  unique %>% as.character ->
  FocalYears

AnnualCentroids %>% group_by(Year) %>% count %>%
  # filter(n>50) %>%
  pull(Year) %>% intersect(FocalYears) ->
  FocalYears

SubCentroids <- AnnualCentroids %>% filter(Year %in% FocalYears)

SPDF <- SpatialPointsDataFrame(data = SubCentroids[,c("X", "Y", "Year")],
                               coords = SubCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Year"]

KUDL <- kernelUD(SPDF,
                 # extent = 0,
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
  bind_rows(.id = "Year") %>% mutate_at("Year", as.numeric) -> SocialNodeTraits

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
    
    Dyads %>% 
      t %>% data.frame %>% 
      map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas)], .x %>% unlist %>% as.character, Symmetrical = T)) %>% 
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

saveRDS(MCPListList, glue("Datasets/{FocalSystem}/Intermediate/MCPList.rds"))

FullSpatialList <-
  list(AnnualDistance = AnnualDistances, 
       AnnualHRO = HROList,
       AnnualAreas = AreaList)

FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

AreaList <- AreaList[which((AreaList %>% map(length))>0)]

if(length(AreaList)<0){
  
  NodeDF <-
    AreaList %>% bind_rows(.id = "Year") %>% 
    full_join(NodeDF, ., by = c("Year", "Id"))
  
}

NodeDF %>% 
  dplyr::select(Year, Id) %>% 
  # unique %>%
  nrow

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <- readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))
