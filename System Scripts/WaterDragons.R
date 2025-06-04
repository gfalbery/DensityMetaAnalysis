
# WaterDragons ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph); library(igraph); library(ggraph)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "WaterDragons"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Their script ####

library(adehabitatHR)
library(raster)
library(rgdal)
library(maptools)
library(geosphere)

# Read the data, previously saved as csv

data <- read.csv(paste0(Root, "/Latest_database_2012-2022.csv"))

data[6] <- NULL # remove Location column
data <- data[, c(3, 1, 4, 5, 2)]
data$Name <- as.factor(data$Name)

# Read the date and time properly in R

data$Time <- gsub("7am-10am", "07:00", data$Time)
data$Time <- gsub("10am-1pm", "07:00", data$Time)
data$Time <- gsub("1pm-3pm", "13:00", data$Time)
data$datetime <- paste(data$Date, data$Time, sep=" ")
str(data)
data$datetime <- strptime(data$datetime, format=c("%d/%m/%Y %H:%M"))
data$datetime <- as.POSIXct(data$datetime)
head(data$datetime)

# Convert latitude/longitude coordinates into X/Y coordinates

data$Lat <-  ((((data$Lat/1000)+27)/60)+27)*-1
data$Long <- as.numeric(data$Long)
data$Long <- (((data$Long/1000)+1)/60)+153

names(data)[3:4] <- c("Y","X")
data <- data[c(1,2,5,6,4,3)]
head(data)
# plot(data$X, data$Y, col = data$Name) #Is there any obvious outliers?

##need to map this, and clean it if necessary - make a spatial points data frame
#####################

# Split per date

data <- na.omit(data)
surveys <- split(data, data$datetime)
length(surveys) #is this the number of surveys we have? check!
surveys2 <- list()
surveys1 <- list()
for(i in 1:length(surveys)){
  if(dim(surveys[[i]])[1] == 1){
    surveys1[[i]] <- surveys[[i]]
  }else{
    surveys2[[i]] <- surveys[[i]]
  }}
surveys1 <- surveys1[!sapply(surveys1, is.null)] 
surveys2 <- surveys2[!sapply(surveys2, is.null)] 


# Calculate geographic distances per date

distances <- list()

l <- 1

for(l in 1:length(surveys2)){
  
  track <- tr <- surveys2[[l]]
  date.time <- c(track[1, "datetime"])
  coordinates(track) <- c("X", "Y")
  dm <- distm(track)
  
  tr <- surveys2[[l]]
  rownames(dm) <- #paste(tr[,1])
    colnames(dm) <- 
    tr[,1]
  
  distances[[l]] <- dm %>% 
    reshape2::melt() %>% # Turns matrix into data frame
    slice(which(lower.tri(dm))) %>% # Only takes the lower triangle
    rename(names1 = 1, names2 = 2, dm.vec.tri = 3) %>% 
    mutate(datetime = date.time) %>% 
    dplyr::select(datetime, names1, names2, dm.vec.tri)
  
}

names(distances) <- surveys2 %>% map_chr(~.x$datetime %>% unique)

#### The list called distances will give you geographic proximity between pairs on each day #####

# Set min and max distance to include in grouping
min = 0 
max = 1.85 # dragons within 1.85m are considered to be associating

EdgeList <- 
  distances %>% 
  map(function(a){
    
    a %>% filter(between(dm.vec.tri, min, max))
    
  }) %>% bind_rows

EdgeList %<>% mutate(Year = str_split(datetime, "-") %>% map_chr(1))

GraphList <- 
  EdgeList$Year %>% unique %>% sort %>% 
  map(function(a){
    
    EdgeList %>% 
      filter(Year == a) %>% 
      dplyr::select(names1, names2, datetime) %>%
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

names(AdjList) <- EdgeList$Year %>% unique %>% sort

AdjList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

ObservationDF <- data %>% rename(Id = Name) %>% 
  mutate(Year = str_split(datetime, "-") %>% map_chr(1))

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

# Making Densities 

SPDF <- SpatialPointsDataFrame(data = LifetimeCentroids[,c("X", "Y")], 
                               coords = LifetimeCentroids[,c("X", "Y")])

MetreDims <- 
  LifetimeCentroids[,c("X", "Y")] %>% 
  map_dbl(c(range, diff)) %>% 
  multiply_by(100000) %>%
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
  unique %>% as.character -> FocalSites

AnnualCentroids %>% group_by(Year) %>% count %>% 
  # filter(n>50) %>% 
  pull(Year) %>% 
  intersect(FocalSites) ->
  FocalSites

SubCentroids <- AnnualCentroids %>% filter(Year %in% FocalSites)

SPDF <- SpatialPointsDataFrame(data = SubCentroids[,c("X", "Y", "Year")], 
                               coords = SubCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Year"]

KUDL <- kernelUD(SPDF, 
                 extent = 0,
                 same4all = TRUE, 
                 grid = 500)

1:length(FocalSites) %>% lapply(function(a){
  
  print(FocalSites[a])
  
  DF <- AnnualCentroids %>% filter(Year == FocalSites[a])
  
  KUDL2 <- KUDL[[FocalSites[a]]]
  
  KUDLRaster <- KUDL2 %>% raster::raster() 
  
  # KUDLRaster <- raster::resample(KUDLRaster, BlankRaster)
  
  # raster::values(KUDLRaster) <-
  #   raster::values(KUDLRaster)*nrow(DF)/
  #   sum(raster::values(KUDLRaster), na.rm = T)
  
  KUDLRaster %>% 
    raster::extract(DF[,c("X", "Y")]) ->
    DF$Density.Annual
  
  return(DF)
  
}) -> DensityList

DensityList %>% bind_rows -> AnnualCentroids

NodeDF <- 
  AnnualCentroids

# NodeDF %<>% mutate_at(vars(contains("Density")), ~.x/(10^6))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Adding Social nodes to Spatial nodes ####

AMList <- AdjList

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% mutate_at("Year", as.character) -> SocialNodeTraits

SocialNodeTraits %<>%
  full_join(NodeDF, by = c("Id", "Year"))

NodeDF <- SocialNodeTraits

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))
NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

# Making Distance Matrices ####

LifetimeDistance <-
  LifetimeCentroids[,c("X", "Y")] %>% 
  GristMatrix(LifetimeCentroids$Id)

FocalSites %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Year == a) 
  
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

names(AreaList) <- 
  names(MCPListList) <-
  names(HROList) <- FocalSites

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
    AreaList %>% bind_rows(.id = "Year") %>% mutate_at("Year", as.numeric) %>% 
    mutate_at("Id", as.numeric) %>% 
    # left_join(LifetimeAreas, suffix = c(".Annual", ".Lifetime"), by = c("Id")) %>% 
    full_join(NodeDF, ., by = c("Year", "Id"))
  
}

NodeDF %>% dplyr::select(Year, Id) %>% 
  # unique %>%
  nrow

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

# Repeating with fight data ####

# Read the data, previously saved as csv

data <- read.csv(paste0(Root, "/Latest_database_2012-2022(2).csv"))

# data[6] <- NULL # remove Location columns
data <- data[, c("Name", "Date", "Lat", "Long", "Time", "Fighting")]
data$Name <- as.factor(data$Name)

# Read the date and time properly in R

data$Time <- gsub("7am-10am", "07:00", data$Time)
data$Time <- gsub("10am-1pm", "07:00", data$Time)
data$Time <- gsub("1pm-3pm", "13:00", data$Time)
data$datetime <- paste(data$Date, data$Time, sep=" ")
str(data)
data$datetime <- strptime(data$datetime, format=c("%d/%m/%Y %H:%M"))
data$datetime <- as.POSIXct(data$datetime)
head(data$datetime)

# Convert latitude/longitude coordinates into X/Y coordinates

data$Lat <-  ((((data$Lat/1000)+27)/60)+27)*-1
data$Long <- as.numeric(data$Long)
data$Long <- (((data$Long/1000)+1)/60)+153

names(data)[3:4] <- c("Y","X")
head(data)
# plot(data$X, data$Y, col = data$Name) #Is there any obvious outliers?

##need to map this, and clean it if necessary - make a spatial points data frame
#####################

# Split per date

data <- na.omit(data)
surveys <- split(data %>% filter(Fighting == "Yes"), 
                 data %>% filter(Fighting == "Yes") %>% pull(datetime))
length(surveys) #is this the number of surveys we have? check!
surveys2 <- list()
surveys1 <- list()
for(i in 1:length(surveys)){
  if(dim(surveys[[i]])[1] == 1){
    surveys1[[i]] <- surveys[[i]]
  }else{
    surveys2[[i]] <- surveys[[i]]
  }}
surveys1 <- surveys1[!sapply(surveys1, is.null)] 
surveys2 <- surveys2[!sapply(surveys2, is.null)] 


# Calculate geographic distances per date

distances <- list()

l <- 1

for(l in 1:length(surveys2)){
  
  track <- tr <- surveys2[[l]]
  date.time <- c(track[1, "datetime"])
  coordinates(track) <- c("X", "Y")
  dm <- distm(track)
  
  tr <- surveys2[[l]]
  rownames(dm) <- #paste(tr[,1])
    colnames(dm) <- 
    tr[,1]
  
  distances[[l]] <- dm %>% 
    reshape2::melt() %>% # Turns matrix into data frame
    slice(which(lower.tri(dm))) %>% # Only takes the lower triangle
    rename(names1 = 1, names2 = 2, dm.vec.tri = 3) %>% 
    mutate(datetime = date.time) %>% 
    dplyr::select(datetime, names1, names2, dm.vec.tri)
  
}

names(distances) <- surveys2 %>% map_chr(~.x$datetime %>% unique)

#### The list called distances will give you geographic proximity between pairs on each day #####

# Set min and max distance to include in grouping
min = 0 
max = 1.85 # dragons within 1.85m are considered to be associating

EdgeList <- 
  distances %>% 
  map(function(a){
    
    a %>% filter(between(dm.vec.tri, min, max))
    
  }) %>% bind_rows

EdgeList %<>% mutate(Year = str_split(datetime, "-") %>% map_chr(1))

GraphList <- 
  EdgeList$Year %>% unique %>% sort %>% 
  map(function(a){
    
    EdgeList %>% 
      filter(Year == a) %>% 
      dplyr::select(names1, names2, datetime) %>%
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

names(AdjList) <- EdgeList$Year %>% unique %>% sort

AdjList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

ObservationDF <- data %>% rename(Id = Name) %>% 
  mutate(Year = str_split(datetime, "-") %>% map_chr(1))

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
  multiply_by(100000) %>%
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
  unique %>% as.character -> FocalSites

AnnualCentroids %>% group_by(Year) %>% count %>% 
  # filter(n>50) %>% 
  pull(Year) %>% 
  intersect(FocalSites) ->
  FocalSites

SubCentroids <- AnnualCentroids %>% filter(Year %in% FocalSites)

SPDF <- SpatialPointsDataFrame(data = SubCentroids[,c("X", "Y", "Year")], 
                               coords = SubCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Year"]

KUDL <- kernelUD(SPDF, 
                 extent = 0,
                 same4all = TRUE, 
                 grid = 500)

1:length(FocalSites) %>% lapply(function(a){
  
  print(FocalSites[a])
  
  DF <- AnnualCentroids %>% filter(Year == FocalSites[a])
  
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

# NodeDF %<>% mutate_at(vars(contains("Density")), ~.x/(10^6))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Adding Social nodes to Spatial nodes ####

AMList <- AdjList

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% mutate_at("Year", as.character) -> SocialNodeTraits

SocialNodeTraits %<>%
  full_join(NodeDF, by = c("Id", "Year"))

NodeDF <- SocialNodeTraits

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))
NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

# Making Distance Matrices ####

LifetimeDistance <-
  LifetimeCentroids[,c("X", "Y")] %>% 
  GristMatrix(LifetimeCentroids$Id)

FocalSites %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% filter(Year == a) 
  
  SubMat <- SubDF %>% 
    dplyr::select(X, Y) %>% 
    GristMatrix(SubDF$Id)
  
  SubMat
  
}) -> AnnualDistances

# Making Space Use Matrices ####

# Annual

HROList <- AreaList <- list()

for(i in FocalSites){
  
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
        map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas) %>% as.character], 
                        .x %>% unlist %>% as.character, Symmetrical = T)) %>% 
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

names(AreaList) <- FocalSites

names(HROList) <- FocalSites

FullSpatialList <-
  list(AnnualDistance = AnnualDistances,
       AnnualHRO = HROList,
       AnnualAreas = AreaList)

FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

AreaList <- AreaList[which((AreaList %>% map(length))>0)]

if(length(AreaList)>0){
  
  NodeDF <-
    AreaList %>% bind_rows(.id = "Year") %>% mutate_at("Year", as.character) %>% 
    # mutate_at("Id", as.numeric) %>% 
    full_join(NodeDF, ., by = c("Year", "Id"))
  
}

NodeDF %>% dplyr::select(Year, Id) %>% 
  # unique %>%
  nrow

NodeDF[,c("X", "Y")] <- NodeDF[,c("X", "Y")]*100

NodeDF[,c("X", "Y")] <- NodeDF[,c("X", "Y")] %>% mutate_all(~.x - min(.x))

NodeDF %<>% mutate_at(vars(Associations:Observations), ~replace_na(.x, 0))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))
