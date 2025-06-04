
# PrairieDogs ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "PrairieDogs"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

Root %>% list.files(full.names = T, pattern = "zip") %>% 
  unzip(exdir = Root)

# Importing data ####

Files <- 
  Root %>% 
  list.files(full.names = T,
             pattern = ".xlsx$")

FileList <- Files %>% map(read_xlsx)

Contacts <- FileList[[2]]

Locations <- FileList[[3]]

Locations %<>% 
  rename_all(~str_replace_all(.x, " ", "_")) %>% 
  # mutate(X = x*100, Y = y*100) %>% 
  mutate(X = UTM_Easting, Y = UTM_Northing) %>% 
  mutate(X = X - min(X), Y = Y - min(Y)) # Remember to put it in kilometres!!!!!!!!

Locations$Year <- 2017# Locations$season

Locations %<>% mutate_at("Animal", ~str_replace_all(.x, " ", "_"))

Contacts %<>% mutate_at(c("Actor", "Receiver"), ~str_replace_all(.x, " ", "_"))

Contacts$Actor %>% setdiff(Contacts$Receiver, .)
Contacts$Actor %>% setdiff(Locations$Animal, .)
Contacts$Receiver %>% setdiff(Locations$Animal)

FocalIndividuals <- 
  Locations$Animal %>% 
  unique %>% sort

Contacts %<>% filter(Actor %in% FocalIndividuals, 
                     Receiver %in% FocalIndividuals)

# for(FocalSite in Sites){
#   
#   print(FocalSite)
#   
#   # Deriving stuff ####

Censuses <- Locations #%>% filter(Site == FocalSite)

# Making data frame

ObservationDF <- Censuses %>% 
  # dplyr::select(-Year) %>% 
  dplyr::select(Id = Animal, X, Y, Date, Year) %>% 
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

# Doing the same with contacts ####

Contacts$Year <- 2017

1:nrow(FocalSubdivisions) %>% 
  map(~Contacts %>% 
        filter(Contacts$Actor %in% unique(ObservationDF$Id),
               Contacts$Actor %in% unique(ObservationDF$Id)) %>% 
        semi_join(FocalSubdivisions[.x,], by = c("Year"))) -> 
  GroupList

# Social Phase ####

AMList <- list()

i <- 1

# for(i in i:length(GroupList)){

print(i)

GroupList[[i]] -> Groups

EdgeList <- Groups %>% 
  arrange(Actor, Receiver) %>% 
  dplyr::select(Actor, Receiver, Social_Behavior) %>% #unique %>% 
  droplevels# %>% as.matrix

Behaviours <- Groups %>% 
  count(Social_Behavior) %>% 
  filter(n > 40) %>% 
  pull(Social_Behavior)

library(tidygraph)

FocalBehaviour <- Behaviours[1]

Behaviours <- "Proximal foraging"

for(i in seq_along(Behaviours)){
  
  FocalBehaviour <- Behaviours[i]
  
  print(FocalBehaviour)
  
  SocGraph <- graph_from_data_frame(EdgeList, directed = F) %>% as_tbl_graph
  
  SocGraph %<>% activate(edges) %>% filter(Social_Behavior == FocalBehaviour)
  
  Layout <- ObservationDF %>% 
    filter(Id %in% V(SocGraph)$name) %>% 
    group_by(Id) %>% 
    summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) %>% 
    ungroup %>% 
    dplyr::select(X, Y) %>% as.matrix
  
  library(ggraph)
  
  SocGraph %>%
    activate(nodes) %>%
    ggraph(Layout) +
    # ggraph() +
    geom_edge_link(alpha = 0.2) +
    geom_node_point() +
    coord_fixed()
  
  AM <- SocGraph %>% get.adjacency() %>% as.matrix
  
  Observations <- Groups %>% dplyr::select(Actor, Receiver) %>% unlist %>% table
  
  Observations <- Observations[colnames(AM)]
  
  Matrix <- AM # %>% ProportionalMatrix(Observations)
  
  diag(Matrix) <- Observations
  
  AMList[[i]] <- Matrix
  
  AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalBehaviour}AMList.rds"))
  
}

# Spatial Phase ####

# Making daily Centroids ####

DailyCentroids <- 
  ObservationDF %>% group_by(Year, Date, Id) %>% 
  summarise_at(c("X", "Y"), ~mean(.x, na.rm = T))

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

SPDF <- SpatialPointsDataFrame(data = DailyCentroids[,c("X", "Y")],
                               coords = DailyCentroids[,c("X", "Y")])

MetreDims <-
  DailyCentroids[,c("X", "Y")] %>%
  map_dbl(c(range, diff)) %>%
  multiply_by(1) %>% ceiling

DailyCentroids[,c("X", "Y")] %>%
  rename_all(tolower) %>% raster::extent() %>%
  raster::raster(ncols = MetreDims[["X"]],
                 nrows = MetreDims[["Y"]]) -> BlankRaster

DailyKUDL <- kernelUD(SPDF,
                      extent = 0,
                      same4all = TRUE,
                      grid = 500)

KUDLRaster <- DailyKUDL %>% raster::raster()

KUDLRaster <- raster::resample(KUDLRaster, BlankRaster)

raster::values(KUDLRaster) <-
  raster::values(KUDLRaster)*nrow(DailyCentroids)/
  sum(raster::values(KUDLRaster), na.rm = T)

KUDLRaster %>%
  raster::extract(AnnualCentroids[,c("X", "Y")]) ->
  AnnualCentroids$Density.Daily

# Making Densities 

SPDF <- SpatialPointsDataFrame(data = LifetimeCentroids[,c("X", "Y")],
                               coords = LifetimeCentroids[,c("X", "Y")])

MetreDims <-
  LifetimeCentroids[,c("X", "Y")] %>%
  map_dbl(c(range, diff)) %>%
  multiply_by(1) %>% ceiling

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

AnnualCentroids %>% group_by(Year) %>% count %>%
  filter(n>10) %>% pull(Year) %>% intersect(FocalYears) ->
  FocalYears

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

names(AMList) <- Behaviours

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  reduce(~full_join(.x, .y, by = "Id")) -> SocialNodeTraits

names(SocialNodeTraits)[2:ncol(SocialNodeTraits)] <- 
  paste0(rep(c("Associations", "Strength", "Degree", "Observations")))#, "_", rep(Behaviours, each = 4))

NodeDF %<>% 
  left_join(SocialNodeTraits, by = c("Id"))

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

HROList <- AreaList <- list()

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

# NodeDF$Associations <- NodeDF %>% ungroup %>% 
#   dplyr::select(contains("Associations_")) %>% rowSums
# 
# NodeDF$Observations <- 
#   NodeDF %>% ungroup %>% 
#   dplyr::select(contains("Observations")) %>% apply(1, max)
# 
# NodeDF$Degree <- 
#   NodeDF %>% ungroup %>% 
#   dplyr::select(contains("Degree")) %>% apply(1, max)
# 
# NodeDF$Degree <- 
#   NodeDF %>% ungroup %>% 
#   dplyr::select(contains("Strength")) %>% 
#   rowSums %>% divide_by(length(Behaviours))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))
