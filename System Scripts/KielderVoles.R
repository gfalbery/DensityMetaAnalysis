
### KielderVoles ### 

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph); library(igraph); library(ggraph)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "KielderVoles"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

# Importing ####

# load(paste0(Root, "/Kielder_SCP_2009_2010_SCP_tags&locations"))

# Voles <- obs_cent

VoleCaptures <- read.csv(paste0(Root, "/VoleCaptures.csv"))

VoleCaptures %>% rename_all(CamelConvert) ->
  VoleCaptures

VoleLocations <- read.csv(paste0(Root, "/VoleLocations.csv"))

VoleLocations %>% rename_all(CamelConvert) ->
  VoleLocations

VoleLocations %>% mutate(Long = substr(Trap, 1, 1),
                         Lat = substr(Trap, 2, 3)) %>%
  mutate(X = factor(Long, levels = LETTERS[1:10]) %>% as.factor %>% as.numeric %>% multiply_by(5),
         Y = factor(Lat, levels = 1:10) %>% as.factor %>% as.numeric %>% multiply_by(5)) ->
  
  VoleLocations

VoleLocations %>% 
  arrange(Date) %>%
  group_by(Tagid, Date) %>% 
  mutate(N = 1:n()) %>% filter(N == 1) ->
  FirstLocations

VoleCaptures %>% 
  inner_join(FirstLocations, by = c("Tagid", "Date", "Site")) ->
  VoleCaptures

VoleCaptures %>% separate(Date, "-", into = c("Day", "Month", "Year")) %>% 
  mutate(Date = paste(Day, Month, Year, sep = "-") %>% lubridate::dmy()) ->
  VoleCaptures

VoleCaptures %<>% filter(Year != "07")

VoleCaptures %<>% mutate_at(c("X", "Y"), ~.x/1000)

# Running pipeline ####

# Klara data too ####

KlaraVoles <- read.csv(paste0(Root, "/", "kielder_site_check_date_loc_ID.csv"))

KlaraVoles %>% 
  separate(Date, sep = "/", into = c("Day", "Month", "Year")) %>% 
  dplyr::select(Day:Year) ->
  KlaraVoles[,c("Day", "Month", "Year")]

KlaraVoles %<>% 
  mutate_at("Location", str_trim) %>% 
  mutate(X = substr(Location, 1, 1)) %>% 
  mutate(Y = substr(Location, 2, 4))

KlaraVoles %<>% 
  mutate_at("Y", ~as.numeric(as.character(.x))) %>% 
  mutate_at("X", ~as.numeric(factor(.x, levels = LETTERS))) %>% 
  mutate_at("X", ~ifelse(.x > 20, .x - 26, .x)) %>% 
  group_by(Site) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x, na.rm = T)) %>% 
  mutate_at(c("X", "Y"), ~.x*5) %>% 
  mutate_at(c("X", "Y"), ~.x/1000)

KlaraVoles %>% 
  ggplot(aes(X, Y)) + geom_point() + facet_wrap(~Site) + coord_fixed()

KlaraVoles %<>% rename(Tagid = Tag_ID) %>% 
  mutate_at("Tagid", ~.x %>% str_trim %>% str_replace_all(" ", "_")) %>% 
  mutate_at("Date", ~lubridate::dmy(.x))

KlaraVoles %<>% filter(!Year == "2010")

VoleCaptures %<>% 
  mutate_at("Tagid", as.character) %>% 
  bind_rows(KlaraVoles)

# AND ANOTHER OOOOOONE ####

Kielder3 <- read.csv(paste0(Root, "/kielder_2015_2017.csv"))

Kielder3 %<>% rename_all(CamelConvert) %>% 
  rename(Tagid = Chip, Date = DateCapt)

Kielder3 %<>% mutate(Long = substr(Trap, 1, 1),
                     Lat = substr(Trap, 2, 3)) %>%
  mutate(X = factor(Long, levels = LETTERS) %>% as.factor %>% as.numeric %>% multiply_by(5),
         Y = factor(Lat, levels = 1:100) %>% as.factor %>% as.numeric %>% multiply_by(5)) %>% 
  group_by(Site) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x, na.rm = T)) %>% 
  # mutate_at(c("X", "Y"), ~.x*5) %>% 
  mutate_at(c("X", "Y"), ~.x/1000)

Kielder3 %>% ungroup %>% 
  separate(Date, sep = "/", into = c("Day", "Month", "Year")) %>% 
  dplyr::select(Day:Year) ->
  Kielder3[,c("Day", "Month", "Year")]

Kielder3 %<>% 
  mutate_at("Date", ~lubridate::dmy(.x))

VoleCaptures %<>% 
  bind_rows(Kielder3 %>% mutate_at("Tagid", as.character))

# Deriving stuff ####

Censuses <- 
  VoleCaptures %>% 
  mutate(GroupDate = paste(Date, X, Y, Site, sep = ","))

# Making data frame

ObservationDF <- Censuses %>% 
  dplyr::select(Id = Tagid, Date, GroupDate, X, Y, Year, Site) %>% 
  mutate(Site_Year = paste0(Site, "_", Year))  %>% 
  na.omit

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

(FocalYears <- 
    ObservationDF$Site_Year %>% unique %>% sort)

ObservationDF %>% group_by(Site_Year) %>% count %>% filter(n>50) %>% 
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
  filter(n>50) %>% pull(Site_Year) %>% intersect(FocalYears) ->
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

# NodeDF %<>% 
#   # mutate()
#   filter(is.na(Density.Annual))

NodeDF %<>% mutate_at(vars(Associations:Degree), ~replace_na(.x, 0))

NodeDF %>% 
  saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <- 
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% 
  ggplot(aes(Density.Annual, Associations)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~Site, scales = "free")
