
# DesertTortoises ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs)

theme_set(theme_cowplot())

FocalSystem <- "DesertTortoises"

# Importing and Cleaning ####

Root <- paste0("Datasets/",FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

DateParse <- function(Vector, Sep){
  
  String <- str_split(Vector, Sep)
  
  DateDF <- String %>% bind_cols %>% t %>% data.frame %>% 
    mutate_all(~.x %>% as.character %>% as.numeric)
  
  Ranges <- DateDF %>% apply(2, range)
  
  Year <- DateDF %>% apply(2, function(a) max(a)>31) %>% which()
  
  if(length(Year) > 0){
    
    Day <- DateDF %>% 
      apply(2, function(a) max(a)>12) %>% 
      which() %>% setdiff(Year)
    
  }
  
  if(length(Year) > 0 & length(Day) > 0){
    
    Month <- 1:3 %>% setdiff(Day) %>% setdiff(Year)
    
  }
  
  colnames(DateDF)[c(Day, Month, Year)] <- 
    c("Day", "Month", "Year")
  
  DateDF[DateDF$Year < 100, "Year"] <- 
    DateDF[DateDF$Year < 100, "Year"] + 1900
  
  DateDF[DateDF$Year < 1921, "Year"] <- 
    DateDF[DateDF$Year < 1921, "Year"] + 100
  
  DateDF %>% dplyr::select(Day, Month, Year) %>% 
    return
  
}

if(!file.exists(paste0(Root, "/Data/", ".TortoiseDF.csv"))){
  
  list.files(paste0(Root, "/Data"), full.names = T, pattern = ".csv") %>% 
    map(fread) -> DFList
  
  names(DFList) <- 
    list.files(paste0(Root, "/Data"), full.names = F, pattern = ".csv") %>% 
    str_split("_") %>% map_chr(1)
  
  for(i in 1:length(DFList)){
    
    print(i)
    
    DFList[[i]]$Date %>% DateParse(Sep = "/|-") ->
      DFList[[i]][,c("Day", "Month", "Year")]
    
  }
  
  DFList %>% 
    map(~.x %>% mutate_all(as.character)) %>%
    bind_rows(.id = "Site") -> Tortoises
  
  Tortoises %>% write.csv(paste0(Root, "/Data/", ".TortoiseDF.csv"), row.names = F)
  
}

Tortoises <- read.csv(paste0(Root, "/Data/", ".TortoiseDF.csv"))

Tortoises %<>% 
  rename_all(~str_remove(.x, "^UTM_") %>% CamelConvert) %>% 
  rename(Code = Tortoise_number)

Tortoises %<>% mutate_at(c("Easting", "Northing"), ~.x/100000) %>% 
  mutate_at(c("Easting", "Northing"), ~.x - min(.x, na.rm = T))

Tortoises %<>% filter(Easting < 100) %>% 
  filter(Easting > 4)

# Making data frame

Censuses <- Tortoises

Censuses %<>% 
  mutate_at("Site", 
            
            ~case_when(.x == "FI" & Year %in% 2005:2008 ~ "FI1",
                       .x == "FI" & Year %in% 2009:2011 ~ "FI2",
                       .x == "FI" & Year %in% 2013:2015 ~ "FI3",
                       TRUE ~ .x)
            
  )

Censuses %<>% 
  # mutate(GroupDate = paste(Site, Burrow_number, sep = ",")) %>% 
  mutate(GroupDate = paste(Site, Date, Burrow_number, sep = ",")) %>% 
  mutate(Site_Year = paste(Site, Year, sep = "_"))

Sites <- Censuses$Site %>% unique %>% sort %>% 
  setdiff(c("HW", "MC", "PV", "SL", "FI"))

FocalSite <- Sites[1]

for(FocalSite in Sites){
  
  print(FocalSite)
  
  ObservationDF <- Censuses %>% 
    filter(Site == FocalSite) %>% 
    # dplyr::select(-Year) %>% 
    dplyr::select(Id = Code, X = Easting, Y = Northing, Year, Site, Site_Year) %>% 
    na.omit
  
  if(FocalSite == "CS"){
    
    ObservationDF %<>% filter(Y < 1.95)
    
  }
  
  ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}ObservationDF.rds"))
  
  FocalYears <- 
    ObservationDF$Year %>% unique %>% sort
  
  Censuses %>% 
    filter(Site == FocalSite) %>% 
    group_by(Site, Year) %>% count %>% filter(n>50) %>% 
    mutate(Site_Year = paste(Site, Year, sep = "_")) %>%
    dplyr::select(-c("n")) -> FocalSubdivisions
  
  1:nrow(FocalSubdivisions) %>% 
    map(~Censuses %>% semi_join(FocalSubdivisions[.x,], by = c("Site", "Year"))) -> 
    GroupList
  
  GroupList %>% bind_rows() %>% nrow
  
  # Social Phase ####
  
  AMList <- list()
  
  i <- 1
  
  for(i in i:length(GroupList)){
    
    print(i)
    
    GroupList[[i]] -> Groups
    
    Groups %>% 
      dplyr::select(Code, GroupDate) %>% #unique %>% 
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
  
  AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}AMList.rds"))
  
  # Spatial Phase ####
  
  # Making Centroids ####
  
  ObservationDF %>% group_by(Site_Year, Id) %>% 
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
    arrange(Site_Year) %>% pull(Site_Year) %>% 
    unique %>% as.character -> FocalYears
  
  AnnualCentroids %>% group_by(Site_Year) %>% count %>% 
    filter(n>50) %>% pull(Site_Year) %>% intersect(FocalYears) ->
    FocalYears
  
  SubCentroids <- AnnualCentroids %>% filter(Site_Year %in% FocalYears)
  
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
  
  NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}NodeData.rds"))
  
  # Adding Social nodes to Spatial nodes ####
  
  names(AMList) <- FocalSubdivisions$Site_Year
  
  AMList %>% 
    map(~data.frame(Associations = rowSums(.x), 
                    Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                    Degree = rowSums(DegreeGet(.x)),
                    Observations = diag(.x)) %>% 
          rownames_to_column("Id")) %>% 
    bind_rows(.id = "Site_Year") -> SocialNodeTraits
  
  NodeDF %<>% left_join(SocialNodeTraits, by = c("Id", "Site_Year"))
  
  NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}NodeData.rds"))
  
  # Making Distance Matrices ####
  
  LifetimeDistance <-
    LifetimeCentroids[,c("X", "Y")] %>% 
    GristMatrix(LifetimeCentroids$Id)
  
  FocalYears %>% map(function(a){
    
    SubDF <- AnnualCentroids %>% filter(Site_Year == a) 
    
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
    
    MCPIndividuals <- AnnualCentroids %>% 
      filter(Site_Year == i) %>% 
      filter(Obs > 5) %>% 
      pull(Id)
    
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
    
    Dyads %>% #t %>% data.frame %>% 
      # dplyr::select(1:2)
      # slice(1:1000) %>% 
      t %>% data.frame %>% 
      map(~MCPOverlap(MCPList[MCPIndividuals %>% setdiff(ZeroAreas)], .x %>% unlist %>% as.character, Symmetrical = T)) %>% 
      unlist %>% 
      matrix(nrow = length(MCPIndividuals %>% setdiff(ZeroAreas))) ->
      AnnualHRO
    
    dimnames(AnnualHRO) <- list(MCPIndividuals %>% setdiff(ZeroAreas), 
                                MCPIndividuals %>% setdiff(ZeroAreas))
    
    HROList[[i]] <- AnnualHRO
    
  }
  
  saveRDS(MCPListList, glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}MCPList.rds"))
  
  FullSpatialList <-
    list(#LifetimeDistance = LifetimeDistance,
         AnnualDistance = AnnualDistances,
         AnnualHRO = HROList,
         AnnualAreas = AreaList)
  
  FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}FullSpatialList.rds"))
  
}

NodeDataList <- Sites %>% 
  map(~readRDS(glue("Datasets/{FocalSystem}/Intermediate/{.x}NodeData.rds")))

names(NodeDataList) <- Sites

NodeDataList %>% bind_rows(.id = "Site") %>% 
  # mutate(Site_Year = paste0(Site, "_", Year)) %>% 
  pull(Site_Year) %>% table

NodeDataList %<>% bind_rows(.id = "Site") %>% 
  # mutate(Site_Year = paste0(Site, "_", Year)) %>% 
  mutate(Year = str_remove(Site_Year, paste0(Site, "_")))

NodeDataList %>% 
  saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

SpatialDataList <- Sites[-3] %>% map(~readRDS(glue("Datasets/{FocalSystem}/Intermediate/{.x}FullSpatialList.rds")))

# names(SpatialDataList) <- Sites[-3]

SpatialDataList %>% map("AnnualHRO") %>% unlist(recursive = F) %>% list(AnnualHRO = .) %>% 
  saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

NodeDF <- readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% ggplot(aes(Density.Annual, Associations)) + geom_point() + geom_smooth() +
  facet_wrap(~Site, scales = "free")

NodeDF %>% is.na %>% colSums
