
# AfricanDogs ####

rm(list = ls())

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "AfricanDogs"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

Root %>% list.files(full.names = T, pattern = "zip") %>% 
  unzip(exdir = Root)

# Importing data ####

Files <- 
  Root %>% 
  list.files(full.names = T,
             pattern = "csv$")

FileList <- Files %>% map(read.csv)

Contacts <- FileList[[1]]

Locations <- FileList[[2]]

Locations %<>% 
  rename_all(~str_replace_all(.x, " ", "_")) %>% 
  mutate(X = x*100, Y = y*100) %>% 
  mutate(X = X - min(X), Y = Y - min(Y)) # Remember to put it in kilometres!!!!!!!!

Locations$Year <- Locations$season

Locations %<>% mutate(Date = dateTime %>% str_split(" ") %>% map_chr(1))

# Locations %<>% filter(str_detect(villageID, "MEDEGUE"))

FocalIndividuals <- 
  Contacts %>% dplyr::select(ID1, ID2) %>% unlist %>% 
  unique %>% sort %>% 
  list %>% 
  append(Locations["dogID"]) %>% 
  reduce(intersect) %>% sort

Contacts %<>% filter(ID1 %in% FocalIndividuals, 
                     ID2 %in% FocalIndividuals)

Locations %<>% filter(dogID %in% FocalIndividuals)

Locations %<>% rename(Site = villageID)

Sites <- Locations$Site %>% unique %>% sort %>% setdiff("TARANGARA")

FocalSite <- Sites[1]

for(FocalSite in Sites){
  
  print(FocalSite)
  
  # Deriving stuff ####
  
  Censuses <- Locations %>% filter(Site == FocalSite)
  
  # Making data frame
  
  ObservationDF <- Censuses %>% 
    # dplyr::select(-Year) %>% 
    dplyr::select(Id = dogID, X, Y, Date, Year) %>% 
    na.omit
  
  ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}ObservationDF.rds"))
  
  FocalYears <- 
    ObservationDF$Year %>% unique %>% sort
  
  Censuses %>% group_by(Year, Site) %>% count %>% filter(n>50) %>% 
    dplyr::select(-c("n")) -> FocalSubdivisions
  
  1:nrow(FocalSubdivisions) %>% 
    map(~ObservationDF %>% semi_join(FocalSubdivisions[.x,], by = c("Year"))) -> 
    GroupList
  
  GroupList %>% bind_rows() %>% nrow
  
  # Doing the same with contacts ####
  
  Contacts$Year <- Contacts$season
  
  1:nrow(FocalSubdivisions) %>% 
    map(~Contacts %>% 
          filter(Contacts$ID1 %in% unique(ObservationDF$Id),
                 Contacts$ID2 %in% unique(ObservationDF$Id)) %>% 
          semi_join(FocalSubdivisions[.x,], by = c("Year"))) -> 
    GroupList
  
  # Social Phase ####
  
  AMList <- list()
  
  i <- 1
  
  for(i in i:length(GroupList)){
    
    print(i)
    
    GroupList[[i]] -> Groups
    
    EdgeList <- Groups %>% 
      arrange(ID1, ID2) %>% 
      dplyr::select(ID1, ID2) %>% #unique %>% 
      droplevels %>% as.matrix
    
    SocGraph <- graph_from_edgelist(EdgeList, directed = F) %>% as_tbl_graph
    
    Layout <- ObservationDF %>% 
      filter(Id %in% V(SocGraph)$name) %>% 
      group_by(Id) %>% 
      summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) %>% 
      ungroup %>% 
      dplyr::select(X, Y) %>% as.matrix
    
    # SocGraph %>% 
    #   activate(nodes) %>% 
    #   left_join(Locations %>% dplyr::select(dogID, villageID) %>% unique,
    #             by = c("name" = "dogID")) %>% 
    #   activate(edges) %>% slice(sample(1:length(E(SocGraph)), 5000)) %>% 
    #   ggraph(Layout) +
    #   # ggraph() + 
    #   geom_edge_link(alpha = 0.2) + 
    #   geom_node_point(aes(colour = villageID)) + 
    #   coord_sf()
    
    AM <- SocGraph %>% get.adjacency() %>% as.matrix
    
    Observations <- Groups %>% dplyr::select(ID1, ID2) %>% unlist %>% table
    
    Observations <- Observations[colnames(AM)]
    
    Matrix <- AM # %>% ProportionalMatrix(Observations)
    
    diag(Matrix) <- Observations
    
    AMList[[i]] <- Matrix
    
  }
  
  AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}AMList.rds"))
  
  # Spatial Phase ####
  
  # Making daily Centroids ####
  
  DailyCentroids <- 
    ObservationDF %>% group_by(Year, Date, Id) %>% 
    summarise_at(c("X", "Y"), ~mean(.x, na.rm = T))
  
  # Making Centroids ####
  
  ObservationDF %>% group_by(Year, Id) %>% 
    summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
    AnnualCentroids
  
  if(FocalSite == "KIRA"){
    
    AnnualCentroids %<>% filter(X > 290)
    
  }
  
  if(FocalSite == "BEMBAYA"){
    
    AnnualCentroids %<>% filter(Y < 17.4)
    
  }
  
  AnnualCentroids %>% group_by(Id) %>% 
    summarise_at(c("X", "Y"), ~mean(.x, na.rm = T)) ->
    LifetimeCentroids
  
  ObservationDF %>% group_by(Year, Id) %>% count(name = "Obs") %>% 
    left_join(AnnualCentroids, .) -> AnnualCentroids
  
  ObservationDF %>% group_by(Id) %>% count(name = "Obs") %>% 
    left_join(LifetimeCentroids, .) -> LifetimeCentroids
  
  # Making Densities 
  
  SPDF <- SpatialPointsDataFrame(data = AnnualCentroids[,c("X", "Y")],
                                 coords = AnnualCentroids[,c("X", "Y")])
  
  MetreDims <-
    AnnualCentroids[,c("X", "Y")] %>%
    map_dbl(c(range, diff)) %>%
    multiply_by(1000) %>% ceiling
  
  AnnualCentroids[,c("X", "Y")] %>%
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
  
  NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}NodeData.rds"))
  
  # Adding Social nodes to Spatial nodes ####
  
  names(AMList) <- FocalSubdivisions$Year
  
  AMList %>% 
    map(~data.frame(Associations = rowSums(.x), 
                    Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                    Degree = rowSums(DegreeGet(.x)),
                    Observations = diag(.x)) %>% 
          rownames_to_column("Id")) %>% 
    bind_rows(.id = "Year") -> SocialNodeTraits
  
  NodeDF %<>% left_join(SocialNodeTraits, by = c("Id", "Year"))
  
  NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}NodeData.rds"))
  
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
  
  FullSpatialList <-
    list(#AnnualDistance = AnnualDistances,
      AnnualHRO = HROList,
      AnnualAreas = AreaList)
  
  saveRDS(MCPListList, glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}MCPList.rds"))
  
  FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}FullSpatialList.rds"))
  
  # Joining Spatial and Social Node Traits #####
  
  AreaList <- AreaList[which((AreaList %>% map(length))>0)]
  
  if(length(AreaList)<0){
    
    NodeDF <-
      AreaList %>% bind_rows(.id = "Year") %>% 
      full_join(NodeDF, ., by = c("Year", "Id"))
    
  }
  
  NodeDF %>% 
    dplyr::select(Year, Id) %>% 
    nrow
  
  NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/{FocalSite}CleanData.rds"))
  
}

NodeDataList <- Sites %>% map(~readRDS(glue("Datasets/{FocalSystem}/Intermediate/{.x}CleanData.rds")))

names(NodeDataList) <- Sites

NodeDataList %>% bind_rows(.id = "Site") %>% 
  mutate(Site_Year = paste0(Site, "_", Year)) %>% pull(Site_Year) %>% table

NodeDataList %>% bind_rows(.id = "Site") %>% 
  mutate(Site_Year = paste0(Site, "_", Year)) %>% 
  mutate_at(vars(Associations:Observations), ~replace_na(.x, 0)) %>% 
  saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

SpatialDataList <- Sites %>% map(~readRDS(glue("Datasets/{FocalSystem}/Intermediate/{.x}FullSpatialList.rds")))

names(SpatialDataList) <- Sites

SpatialDataList %>% map("AnnualHRO") %>% unlist(recursive = F) %>% list(AnnualHRO = .) %>% 
  saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))

NodeDF <- readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))
