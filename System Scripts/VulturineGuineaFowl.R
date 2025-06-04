
# VulturineGuineaFowl ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph)

theme_set(theme_cowplot())

FocalSystem <- "VulturineGuineafowl"

# Importing and Cleaning ####

Root <- paste0("Datasets/",FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Data

VGF <- Root %>% paste0("/Data") %>% dir_ls(regex = "csv$") %>% fread

VGF %<>% as.data.frame

VGF %>% 
  separate(Date, sep = "-", into = c("Year", "Month", "Day")) %>% 
  dplyr::select(c("Year", "Month", "Day")) ->
  VGF[,c("Year", "Month", "Day")]

VGF %>% filter(Date == Date[1104]) %>% 
  ggplot(aes(East_UTM, North_UTM)) + 
  geom_point(aes(colour = ID)) + facet_wrap(~timestamp)

VGF %<>% 
  mutate_at(c("East_UTM", "North_UTM"), ~(.x - min(.x))/1000)

ThinVGF <- VGF %>% dplyr::select(c("East_UTM", "North_UTM", 
                                   "Group",
                                   "timestamp"))

TimeStamps <- VGF$timestamp %>% unique

DistList <- 
  TimeStamps %>% 
  map(function(a){
    
    print(a)
    
    SubVGF <- ThinVGF[ThinVGF$timestamp == a,]
    
    Dist <- SubVGF[, 1:2] %>% dist %>% as.matrix
    
    dimnames(Dist) <- list(SubVGF$Group, SubVGF$Group)
    
    Dist %>% return
    
  })

DistList %>% map_dbl(function(a) sum(a[a>0]<0.02)) %>% divide_by(2) %>% table

Threshold <- 0.2

ContactList <- 
  DistList %>% 
  map(function(a){
    
    a %>% 
      reshape2::melt() %>% 
      slice(which(lower.tri(a))) %>% 
      filter(value < Threshold) 
    
  })

names(ContactList) <- TimeStamps

EdgeList <- 
  ContactList %>% 
  bind_rows(.id = "TimeStamp")

EdgeList %>% 
  separate(TimeStamp, sep = "-", into = c("Year", "Month", "Day")) %>% 
  separate(Day, sep = " ", into = c("Day", "Time")) %>% 
  dplyr::select(c("Year", "Month", "Day")) ->
  EdgeList[,c("Year", "Month", "Day")]

EdgeList %<>% dplyr::select(-value)

GroupSizes <- 
  VGF %>% 
  dplyr::select(Group, Group_size, Year, Month, Day) %>% 
  unique

EdgeList %<>% 
  left_join(GroupSizes, by = c("Var1" = "Group", "Year", "Month", "Day"))

EdgeList %<>% 
  left_join(GroupSizes, by = c("Var2" = "Group", "Year", "Month", "Day"), 
            suffix = c("1", "2"))

# EdgeList %<>% 
#   mutate(Weight = Group_size1 + Group_size2)

EdgeList %<>% 
  mutate(Weight = Group_size1 * Group_size2)

EdgeList %>% saveRDS(paste0(Root, "/Intermediate/Edgelist.rds"))
EdgeList <- readRDS(paste0(Root, "/Intermediate/Edgelist.rds"))

EdgeList %<>% dplyr::select(Var1, Var2, Year:Group_size2) %>% unique

Contacts <- EdgeList

# Next step ####

GraphList <- 
  Contacts$Year %>% unique %>% sort %>% 
  map(function(a){
    
    Contacts %>% 
      filter(Year == a) %>% mutate(Weight = Group_size1 * Group_size2) %>% 
      dplyr::select(Var1, Var2, Weight) %>%
      graph_from_data_frame(directed = F) %>% 
      as_tbl_graph
    
  })

AdjList <- GraphList %>% 
  map(~get.adjacency(.x, attr = "Weight")) %>% 
  map(function(a){
    
    a <- as.matrix(a)
    
    diag(a) <- 0
    
    Observations <- rowSums(a)
    
    diag(a) <- Observations
    
    a %>% return
    
  })

AdjList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

# Deriving stuff ####

Censuses <- VGF %>% mutate(X = location_lat * 100, 
                           Y = location_long * 100) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x))

# Making data frame

ObservationDF <- Censuses %>% 
  # dplyr::select(-Year) %>% 
  dplyr::select(Id = Group, X, Y, Year, Month, Day) %>% 
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

# Spatial Phase ####

# Making Centroids ####

ObservationDF %<>% 
  left_join(GroupSizes, by = c("Id" = "Group", "Year", "Month", "Day"))

ObservationDF %>% 
  group_by(Year, Id) %>%
  summarise_at(c("X", "Y", "Group_size"), ~mean(.x, na.rm = T)) %>%
  mutate_at(c("Group_size"), ~ceiling(.x)) ->
  AnnualCentroids

AnnualCentroids <- AnnualCentroids[rep(1:length(AnnualCentroids$Group_size), AnnualCentroids$Group_size), ]

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
  arrange(Year) %>% pull(Year) %>% 
  unique %>% as.character -> FocalYears

AnnualCentroids %>% group_by(Year) %>% count %>% 
  # filter(n>50) %>% 
  pull(Year) %>% intersect(FocalYears) ->
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

AMList <- AdjList

names(AMList) <- FocalSubdivisions$Year

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% 
  mutate_at("Id", as.numeric) -> SocialNodeTraits

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
        
        print(a)
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

names(AreaList) <- 
  names(HROList) <- 
  names(MCPListList) <- FocalYears

MCPListList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/MCPList.rds"))

FullSpatialList <-
  list(AnnualDistance = AnnualDistances,
       AnnualHRO = HROList,
       AnnualAreas = AreaList)

FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

AreaList <- AreaList[which((AreaList %>% map(length))>0)]

if(length(AreaList)>0){
  
  NodeDF <-
    AreaList %>% bind_rows(.id = "Year") %>% mutate_at(c("Year", "Id"), as.numeric) %>% 
    full_join(NodeDF %>% mutate_at(c("Year", "Id"), as.numeric), ., by = c("Year", "Id"))
  
}

NodeDF %>% 
  dplyr::select(Year, Id) %>% unique() %>% data.frame

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))


