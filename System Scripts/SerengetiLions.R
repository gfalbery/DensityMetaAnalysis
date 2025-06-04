
# 0_Lion Data Import ####

rm(list = ls())

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(TDPanalysis)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "SerengetiLions"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

load(paste0("Datasets/", FocalSystem, "/Sightings_and_demo_for_Nets.Rdata"))

Lions <- dat.sights.u %>% rename_all(CamelConvert)

Lions %<>% filter(between(North, 9675000, 9742000))

Lions %<>% filter(between(East, 670000, 760000))

Lions %<>% mutate_at(vars(contains("date")), ~lubridate::ymd(.x))

Lions %<>% mutate(Age = as.numeric(Full_date - Birth_date)/365)

Lions %<>% mutate(Date = Full_date) %>% arrange(Date)

Lions %<>% mutate_at("Season", CamelConvert)

Lions %<>% 
  mutate(LionYear = ifelse(Month < 6, Year - 1, Year)) %>% 
  mutate(Year = ifelse(Month < 6, Year - 1, Year)) %>% 
  mutate(SeasonNo = paste0(LionYear, ", ", Season))

Lions %<>% filter(Sex %in% c("F", "M"))

# Lions$SeasonNo <- Lions$LionYear

# Lions[,c("Year", "Month", "Season", "LionYear", "SeasonNo")] %>% unique

Lions %<>% 
  separate(Full_date, "-", into = c("Year", "Month", "Day")) %>% 
  mutate(DOY = date.to.DOY(paste(Day, Month, Year, sep = "/")))

Lions %<>% 
  mutate(Full_date = Birth_date) %>% 
  separate(Full_date, "-", into = c("Year", "Month", "Day")) %>% 
  mutate(BirthDOY = date.to.DOY(paste(Day, Month, Year, sep = "/")))

LionTraits <- 
  Lions %>% 
  dplyr::select(Id, #Birth_date,
                Collared,
                Sex) %>% 
  unique

FocalYears <- 
  Lions$LionYear %>% 
  unique %>% sort

FY <- FocalYears[1]

AgeMatList <- AMList <- list()

GraphList <- 
  FocalYears %>% 
  map(function(FY){
    
    print(FY)
    
    FocalLions <- 
      Lions %>% filter(LionYear == FY) %>% 
      pull(Id) %>% unique %>% sort 
    
    IM <- 
      Lions %>% filter(LionYear == FY) %>% 
      dplyr::select(Groupn, Id) %>% arrange(Id) %>% 
      table
    
    AM <- 
      IM %>% 
      graph_from_incidence_matrix %>% 
      bipartite.projection %>%
      extract2(2) %>% 
      get.adjacency(attr = "weight", sparse = F)
    
    Matrix <- AM[FocalLions, FocalLions]
    
    Observations = 
      Lions %>% filter(LionYear == FY) %>% 
      count(Id) %>% 
      slice(order(Id)) %>% 
      pull(n)
    
    diag(Matrix) <- Observations
    
    AMList[[which(FocalYears == FY)]] <<-
      Matrix
    
  })

# names(GraphList) <- 
  names(AMList) <- FocalYears

# AMList <- GraphList

AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList.rds"))

# Spatial Phase ####

# Making Centroids ####

ObservationDF <- 
  Lions %>% #rename(Id = Bird) %>% 
  rename(X = East, Y = North) %>% 
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
  # multiply_by(1000) %>% 
  ceiling

AnnualCentroids[,c("X", "Y")] %>%
  rename_all(tolower) %>% raster::extent() %>%
  raster::raster(ncols = MetreDims[["X"]],
                 nrows = MetreDims[["Y"]]) -> BlankRaster

AnnualCentroids %>%
  arrange(Year) %>% pull(Year) %>%
  unique %>% as.character ->
  FocalYears

AnnualCentroids %>% group_by(Year) %>% count %>%
  filter(n>20) %>%
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

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") -> SocialNodeTraits

NodeDF %<>% left_join(SocialNodeTraits, by = c("Id", "Year"))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/NodeData.rds"))

# Making Distance Matrices ####

FocalYears %>% map(function(a){
  
  SubDF <- AnnualCentroids %>% 
    ungroup %>% 
    filter(Year == a) 
  
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
