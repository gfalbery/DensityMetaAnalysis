### MultimammateMice ### 

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph); library(igraph); library(ggraph)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "MultimammateMice"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Importing ####

Spouce <- read.csv(paste0(Root, "/data_greg_multimammate.csv"))

Spouce <- read.csv(paste0(Root, "/DATA_GREG_DENSITY.csv"))

Spouce %>% rename_all(CamelConvert) %>% 
  mutate(Date = as.Date(DATE, format="%d/%m/%Y"), 
         Month = month(Date) %>% as.factor(), 
         Year = year(Date) %>% as.factor()) ->
  Spouce

Spouce %<>% mutate(Y = sapply(Y, function(a) which(LETTERS == a)))

Spouce %<>% mutate_at(c("X", "Y"), ~.x/1000)

# Running pipeline ####

# Deriving stuff ####

Censuses <- Spouce %>% 
  # mutate(GroupDate = paste(Site, Burrow_number, sep = ",")) %>% 
  mutate(GroupDate = paste(Date, X, Y, sep = ","))

# Making data frame

ObservationDF <- Censuses %>% 
  # dplyr::select(-Year) %>% 
  dplyr::select(Id = TOE, Date, GroupDate, X, Y, Year) %>% 
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
  
  Dates <- Groups$Date %>% unique %>% sort
  
  Dates %>% map(function(j){
    
    SubMat <- Groups %>% filter(Date == j) %>% 
      dplyr::select(X, Y) %>% 
      dist %>% 
      as.matrix
    
    dimnames(SubMat) <- Groups %>% filter(Date == j) %>% 
      pull(Id) %>% list(., .)
    
    SubMat[upper.tri(SubMat, diag = T)] <- NA
    
    ifelse(SubMat < (((10^2 + 10^2)^(1/2))/1000), 1, 0) %>% 
      reshape2::melt() %>% filter(value == 1) %>% 
      dplyr::select(-value) %>% 
      mutate(Date = j)
    
  }) %>% bind_rows() %>% filter(!Var1 == Var2) -> EdgeList
  
  SocGraph <- 
    graph_from_data_frame(EdgeList, directed = F)
  
  AM <- SocGraph %>% get.adjacency(sparse = F)
  
  Observations = rowSums(AM)
  
  Matrix <- AM # %>% ProportionalMatrix(Observations)
  
  diag(Matrix) <- Observations
  
  AMList[[i]] <- Matrix
  
}

AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

# Spatial Phase ####

# Making Centroids ####

ObservationDF %>% 
  filter(Year %in% c(0:2018)) %>% 
  group_by(Year, Id) %>% 
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
  multiply_by(1000) %>%
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

AnnualCentroids %>% 
  mutate_at("Year", as.character) %>% 
  filter(Year < 2019) %>% 
  group_by(Year) %>% 
  summarise(n = n(), 
            NSites = nunique(paste0(X, Y))) %>% 
  filter(n>50, NSites > 50) %>% pull(Year) %>% intersect(FocalYears) ->
  FocalYears

SubCentroids <- AnnualCentroids %>% filter(Year %in% FocalYears) %>% droplevels()

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
  bind_rows(.id = "Year") %>% 
  mutate_at("Year", as.factor) -> SocialNodeTraits

NodeDF %<>% left_join(SocialNodeTraits, by = c("Id", "Year"))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

LifetimeDistance <-
  LifetimeCentroids[,c("X", "Y")] %>% 
  GristMatrix(LifetimeCentroids$Id)

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
    AreaList %>% bind_rows(.id = "Year") %>% #mutate_at("Year", as.numeric) %>% 
    full_join(NodeDF, ., by = c("Year", "Id"))
  
}

NodeDF %>% dplyr::select(Year, Id) %>% 
  # unique %>%
  nrow

NodeDF %<>% mutate_at(vars(Associations:Observations), ~replace_na(.x, 0))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))
