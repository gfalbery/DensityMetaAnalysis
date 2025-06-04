
# SoaySheep ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)

theme_set(theme_cowplot())

FocalSystem <- "SoaySheep"

# Importing and Cleaning ####

Root <- paste0("Datasets/", FocalSystem)

Sheep <- read.csv(paste0(Root, "/SheepIndividuals.csv"), header = T)
SheepCensuses <- read.csv(paste0(Root, "/SheepCensuses.csv"), header = T)
SheepParasites <- read.csv(paste0(Root, "/SheepParasites.csv"), header = T)
SheepBirths <- read.csv(paste0(Root, "/SheepBirths.csv"), header = T)
HirtaVillPop <- read.csv(paste0(Root, "/HirtaVillPop.csv"), header = T)

Dates <- read.csv(paste0(Root, "/Dates.csv"), header = T)

# Sorting out censuses ####

# SheepCensuses$ID <- as.factor(SheepCensuses$Id)

# SheepCensuses <- merge(SheepCensuses, Dates[,c("DateUK","Ndate")], all.x = T, by.x = "Date", by.y = "DateUK")
# 
# YearStarts <- Dates[substr(Dates$DateUK,1,5) == "01/01","Ndate"]
# SpringStarts <- Dates[substr(Dates$DateUK,1,5) == "01/05","Ndate"][2:52]
# 
# Timing <- data.frame(Year = 1969:2019, YearStarts, SpringStarts)
# 
# SheepCensuses$Year <- as.numeric(as.character(cut(SheepCensuses$Ndate,breaks=YearStarts,labels=1969:2018)))
# SheepCensuses$SheepYear <- as.numeric(as.character(cut(SheepCensuses$Ndate,breaks=SpringStarts,labels=1969:2018)))

SheepCensuses %<>% mutate_at("Date", ~lubridate::dmy(.x))

SheepCensuses %>% 
  separate(Date, sep = "-", into = c("Year", "Month", "Day")) %>% 
  dplyr::select(c("Year", "Month", "Day")) %>% mutate_all(as.numeric) ->
  SheepCensuses[,c("Year", "Month", "Day")]

SheepCensuses %<>% mutate(SheepYear = ifelse(Month < 5, Year - 1, Year))

# YearLocations <- with(SheepCensuses, cbind(reshape2::melt(tapply(Easting, list(Id, SheepYear), function(m) mean(m, na.rm=T))),
#                                            reshape2::melt(tapply(Northing, list(Id, SheepYear), function(m) mean(m, na.rm=T)))$value))
# 
# names(YearLocations) <- c("ID", "SheepYear", "Easting", "Northing")
# 
# YearLocations <- na.omit(YearLocations[YearLocations$Northing>975,])

SheepCensuses %<>% dplyr::rename(Id = ID)

SheepCensuses %<>% mutate_at("Northing", ~ifelse(.x < 90, .x + 1000, .x))

YearLocations <- 
  SheepCensuses %>% 
  group_by(Id, SheepYear) %>% 
  summarise_at(c("Easting", "Northing"), ~mean(.x, na.rm = T)) %>% 
  na.omit %>% 
  filter(Northing >= 975)

FocalYears <- SheepCensuses$SheepYear %>% unique %>% sort

SheepCensuses %<>% unite(Assoc, GroupNo, Date, sep = "_")

ObservationDF <- 
  SheepCensuses %>% filter(Northing >= 975) %>% 
  dplyr::select(-Year) %>% 
  dplyr::select(Id, X = Easting, Y = Northing, Year = SheepYear, Assoc) %>% 
  na.omit

ObservationDF %<>% mutate_at(c("X", "Y"), ~.x/10) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x))

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

FocalYears %>% map(~ObservationDF %>% filter(Year == .x)) -> GroupList

# Social Phase ####

AMList <- list()

i <- 1

for(i in i:length(GroupList)){
  
  print(i)
  
  GroupList[[i]] -> Groups
  
  Groups %>% 
    dplyr::select(Id, Assoc) %>% unique %>% droplevels %>% 
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
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% 
  mutate_at("Year", as.numeric) ->
  SocialNodeTraits

NodeDF %<>% 
  mutate_at("Id", as.character) %>% 
  left_join(SocialNodeTraits)

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
  
  MCPIndividuals <- AnnualCentroids %>% 
    filter(Year == i) %>% 
    filter(Obs > 5) %>% 
    pull(Id)
  
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
  
  AreaList[[as.character(i)]] <- AnnualAreas
  
  ZeroAreas <- AnnualAreas %>% filter(Area == 0) %>% pull(Id)
  
  Dyads <- expand.grid(ID1 = MCPIndividuals %>% setdiff(ZeroAreas), 
                       ID2 = MCPIndividuals %>% setdiff(ZeroAreas))
  
  if(nrow(Dyads) > 0){
    
    Lower <- lower.tri(matrix(NA, 
                              nrow = length(MCPIndividuals %>% setdiff(ZeroAreas)), 
                              ncol = length(MCPIndividuals %>% setdiff(ZeroAreas))
    ), diag = F) %>% which
    
    Dyads <- Dyads[Lower,]
    
    AnnualHROList <- list()
    
    AnnualHROList[1:nrow(Dyads)] <- list()
    
    j <- 1
    
    NID <- length(MCPIndividuals %>% setdiff(ZeroAreas))
    
    for(j in j:nrow(Dyads)){
      
      # print(round(j/nrow(Dyads), 3))
      
      MCPOverlap(MCPList, 
                 Dyads[j,] %>% unlist %>% as.character, 
                 Symmetrical = T) -> # %>% 
        # unlist %>% 
        # matrix(nrow = NID) ->
        # LifetimeHRO
        
        AnnualHROList[[j]] # <- LifetimeHRO
      
    }
    
    AnnualHRO <- matrix(0, nrow = NID, ncol = NID)
    
    dimnames(AnnualHRO) <- list(MCPIndividuals %>% setdiff(ZeroAreas), 
                                MCPIndividuals %>% setdiff(ZeroAreas))
    
    AnnualHRO[lower.tri(AnnualHRO)] <- AnnualHROList %>% unlist
    
    AnnualHRO %>% 
      EnforceSymmetry() ->
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

FullSpatialList %>% 
  saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

# Joining Spatial and Social Node Traits #####

AreaList <- AreaList[AreaList %>% map_lgl(~!is.null(.x))]
names(AreaList) <- FocalYears

NodeDF <-
  AreaList %>% 
  bind_rows(.id = "Year") %>% mutate_at("Year", as.numeric) %>% 
  full_join(NodeDF, ., by = c("Year", "Id"))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))
