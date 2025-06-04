# 0_Liverpool Data Import ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph); library(igraph); library(ggraph)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "LiverpoolMice"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Final tidied data 

Liverpool <- read.csv(paste0(Root, "/allyears_CLEANED_211118.csv"))

Liverpool <- Liverpool[!is.na(Liverpool$id),]

Liverpool <-  Liverpool %>% subset(species=="WM")

Liverpool %>% rename_all(CamelConvert) %>%
  droplevels %>%
  mutate(Sex = Sex_within_week,
         Age = Age_std,
         Month = as.factor(Month_num),
         Year = as.factor(Year)) ->
  
  Liverpool

# Liverpool %>% ggplot(aes(Long, Lat)) + geom_point() + facet_wrap(~Grid + Site)

BadGrids <- c("HA2")

Liverpool %>% filter(!Grid %in% BadGrids) %>% 
  filter(!Age == "J") %>% droplevels %>%
  mutate(X = Long/1000, Y = Lat/1000) ->
  
  Liverpool

# Deriving stuff ####

MouseCaptures <- Liverpool %>% rename(Date = Capture_date)

Censuses <- 
  MouseCaptures %>% 
  mutate(GroupDate = paste(Date, X, Y, Site, sep = ","))

Censuses %<>% mutate_at("Id", ~paste0("X", .x))

# Making data frame

ObservationDF <- Censuses %>% 
  dplyr::select(Id, Date, GroupDate, X, Y, Year, Site) %>% 
  mutate(Site_Year = paste0(Site, "_", Year))  %>% 
  na.omit

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

(FocalYears <- 
    ObservationDF$Site_Year %>% unique %>% sort)

ObservationDF %>% group_by(Site_Year) %>% count %>% filter(n>30) %>% 
  dplyr::select(-c("n")) -> FocalSubdivisions

1:nrow(FocalSubdivisions) %>% 
  map(~ObservationDF %>% semi_join(FocalSubdivisions[.x,], by = c("Site_Year"))) -> 
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
    
    EdgeDF <- 
      ifelse(SubMat < (((10^2 + 10^2)^(1/2))/1000), 1, 0) %>% 
      reshape2::melt() %>% filter(value == 1) %>%
      mutate_at(c("Var1", "Var2"), as.character) %>% 
      dplyr::select(-value) %>% 
      mutate(Date = j)
    
    if(nrow(EdgeDF) > 0) return(EdgeDF) else NULL
    
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
  group_by(Site_Year, Id) %>% 
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

# SPDF <- SpatialPointsDataFrame(data = LifetimeCentroids[,c("X", "Y")], 
#                                coords = LifetimeCentroids[,c("X", "Y")])

MetreDims <-
  LifetimeCentroids[,c("X", "Y")] %>%
  map_dbl(c(range, diff)) %>%
  multiply_by(1000) %>%
  ceiling

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

# raster::values(KUDLRaster) <- 
#   raster::values(KUDLRaster)*nrow(LifetimeCentroids)/
#   sum(raster::values(KUDLRaster), na.rm = T)

AnnualCentroids %>% 
  arrange(Site_Year) %>% pull(Site_Year) %>% 
  unique %>% as.character -> FocalYears

AnnualCentroids %>% 
  mutate_at("Site_Year", as.character) %>% 
  group_by(Site_Year) %>% 
  summarise(n = n(), 
            NSites = nunique(paste0(X, Y))) %>% 
  filter(n>30) %>% pull(Site_Year) %>% intersect(FocalYears) ->
  FocalYears

SubCentroids <- AnnualCentroids %>% filter(Site_Year %in% FocalYears) %>% droplevels()

SPDF <- SpatialPointsDataFrame(data = SubCentroids[,c("X", "Y", "Site_Year")], 
                               coords = SubCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Site_Year"]

KUDL <- kernelUD(SPDF, 
                 extent = 0,
                 same4all = TRUE, 
                 grid = 500)

1:length(FocalYears) %>% lapply(function(a){
  
  print(FocalYears[a])
  
  DF <- AnnualCentroids %>% filter(Site_Year == FocalYears[a])
  
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

names(AMList) <- FocalSubdivisions$Site_Year

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Site_Year") %>% 
  mutate_at("Site_Year", as.factor) -> SocialNodeTraits

NodeDF %<>% 
  mutate_at("Id", as.character) %>% 
  left_join(SocialNodeTraits, by = c("Id", "Site_Year"))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

# LifetimeDistance <-
#   LifetimeCentroids[,c("X", "Y")] %>% 
#   GristMatrix(LifetimeCentroids$Id)

FocalYears %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Site_Year == a) 
  
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
  
  MCPIndividuals <- 
    ObservationDF %>% 
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
        
        # print(a)
        ObservationDF %>% 
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
    
    if(nrow(Dyads) > 0){
      
      Dyads %>% 
        t %>% data.frame %>% 
        map(function(a){
          
          print(a)
          
          MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas) %>% as.character], 
                     a %>% unlist %>% as.character, Symmetrical = T)
          
        }) %>% 
        
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

NodeDF %>% dplyr::select(Site_Year, Id) %>% 
  # unique %>%
  nrow

NodeDF %<>% filter(!Id == "*")

NodeDF %>% group_by(Site_Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))

NodeDF %<>% mutate_at(vars(Associations:Observations), ~replace_na(.x, 0))

NodeDF %>% 
  # rename(Site_Year = Site_Year) %>% 
  saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <- readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums()
