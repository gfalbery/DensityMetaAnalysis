
# RockHyraxes ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(tidygraph); library(readxl)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

FocalSystem <- "RockHyraxes"

Root <- paste0("Datasets/", FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

Hyraxes <- Root %>% dir_ls(regex = "xlsx") %>% extract(1) %>% read_xlsx %>% as.data.frame

Hyraxes %>% dim

Hyraxes$Site %>% table

Hyraxes[1,]

if(0){
  
  (Locations <- 
     # Hyraxes$Location %>% table %>% as.data.frame %>% 
     # rename(Location = 1, Count = 2) %>% 
     Hyraxes %>% 
     mutate_at("Location", tolower) %>% 
     filter(!Location %in% c("unknown", "?", "<NA>")) %>% 
     group_by(Site, Location) %>% count %>% rename(Count = n) %>% arrange(desc(Count)) %>% 
     as.data.frame)
  
  (Locations <- 
      Locations %>% 
      # mutate_at("Location", ~str_remove_all(.x, "west of |east of |west of |west to |south east in |east to |south of |next to |trail |beside ")) %>% 
      # mutate_at("Location", ~str_remove_all(.x, "near |area | area$|above |below |under |behind |riverbed under |entrance to |[.]|[,]")) %>% 
      # mutate_at("Location", str_trim) %>% 
      # mutate_at("Location", ~str_replace_all(.x, c("twf" = "taily weed field", "jujube " = "jujubes ", "jujube$" = "jujubes",
      #                                              "erect rocks" = "erected rocks",
      #                                              "entrence" = "entrance", "highst" = "highest", "  " = " ", " m " = "m ",
      #                                              "eastern" = "east", "western" = "west", "northern" = "north", "southern" = "south",
      #                                              "sticky rock$" = "sticky rocks"))) %>% 
      # mutate_at("Location", ~.x %>% str_split("by |of |around ") %>% map_chr(last)) %>% 
      # mutate_at("Location", ~.x %>% str_split(" - ") %>% map_chr(first)) %>% 
      # arrange(desc(Count)) %>% filter(str_detect(Location, "sukkot moringa")) %>% 
    group_by(Site, Location) %>% 
      summarise_at("Count", sum) %>% 
      arrange(desc(Count)) %>% 
      as.data.frame())
  
  Locations %>% filter(str_detect(Location, " - "))
  
  Locations %>% filter(str_detect(Location, "duet rock"))
  Locations %>% filter(str_detect(Location, "jujube"))
  
  Locations %>% filter(Count > 50) %>% pull(Count) %>% sum
  
  library(gsheet)
  
  NameSheet <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1JgK3FIihgOgNvgvkq0--5AfgSAF_y6O8Wj0i9j9gdh8/edit#gid=0")
  
  NameSheet$Count <- 
    NameSheet$Location %>% 
    map_dbl(function(a){
      
      # sum(as.numeric(str_detect(Locations$Location, a)),
      #     na.rm = T)
      
      Locations %>% 
        filter(str_detect(Location, a)) %>% 
        pull(Count) %>% sum
      
    })
  
  NameSheet %>% filter(Count == 0)
  
}

# 

Hyraxes %<>% 
  mutate_at("Location", ~str_remove_all(.x, "west of |east of |west of |west to |south east in |east to |south of |next to |trail |beside ")) %>%
  mutate_at("Location", ~str_remove_all(.x, "near |area | area$|above |below |under |behind |riverbed under |entrance to |[.]|[,]")) %>%
  mutate_at("Location", str_trim) %>%
  mutate_at("Location", ~str_replace_all(.x, c("twf" = "taily weed field", "jujube " = "jujubes ", "jujube$" = "jujubes",
                                               "erect rocks" = "erected rocks",
                                               "entrence" = "entrance", "highst" = "highest", "  " = " ", " m " = "m ",
                                               "eastern" = "east", "western" = "west", "northern" = "north", "southern" = "south",
                                               "sticky rock$" = "sticky rocks"))) %>%
  mutate_at("Location", ~.x %>% str_split("by |of |around ") %>% map_chr(last)) %>%
  mutate_at("Location", ~.x %>% str_split(" - ") %>% map_chr(first))


# Hyraxes$Location

Locations <- read_xlsx(paste0(Root, "/HyraxLocations.xlsx")) %>% dplyr::select(Location, X, Y)

Hyraxes %<>% 
  left_join(Locations) %>% 
  filter(!is.na(X))

Hyraxes %<>% mutate_at(c("X", "Y"), as.numeric)

# Hyraxes %>% 
#   ggplot(aes(X, Y)) + 
#   geom_point(aes(colour = Site)) + coord_fixed()

Hyraxes %<>% filter(!(Site == "arugot" & X > 35.37),
                    !(Site == "david" & X < 35.37),
                    !(Site == "outside"))

Hyraxes %<>% 
  dplyr::select(Date) %>% 
  separate(Date, sep = "-", into = c("Year", "Month", "Day")) %>% 
  bind_cols(Hyraxes, .)

# Making long ####

LongHyraxes <- 
  Hyraxes %>% 
  rename(Id = `Focal /initatior (chip)`) %>% 
  filter(!is.na(Id)) %>% 
  mutate(GroupDate = 1:n())

LongHyraxes %<>% 
  mutate_at("Id", ~.x %>% 
              str_replace_all(",", " ") %>% 
              str_replace_all("  ", " ") %>% 
              str_split(" ")) %>% 
  unnest(Id)

LongHyraxes %<>% 
  filter(!str_detect(Id, "Un"))

LongHyraxes <- 
  LongHyraxes[LongHyraxes$Id %>% 
                # paste0("a", .) %>% 
                map_lgl(~str_detect(.x, "^[1-9]")&str_detect(.x, "[1-9]$")) %>% which(),]

LongHyraxes %>% count(Id) %>% arrange(desc(n)) %>% data.frame

LongHyraxes %<>% filter(str_count(Id) > 1)

LongHyraxes %<>% 
  mutate_at(c("X", "Y"), ~.x*100) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x)) 

LongHyraxes %<>% filter(Site == "arugot")

# Pipeline ####

Censuses <- LongHyraxes

Sites <- Censuses$Site %>% unique %>% sort

FocalSite <- Sites[1]

print(FocalSite)

ObservationDF <- Censuses %>% 
  filter(Site == FocalSite) %>% 
  # dplyr::select(-Year) %>% 
  dplyr::select(Id, X, Y, Year, Site) %>% 
  na.omit

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

FocalYears <- 
  ObservationDF$Year %>% unique %>% sort

Censuses %>% 
  filter(Site == FocalSite) %>% 
  group_by(Site, Year) %>% count %>% filter(n>50) %>% 
  # mutate(Year = paste(Site, Year, sep = "_")) %>%
  dplyr::select(-c("n")) -> FocalSubdivisions

1:nrow(FocalSubdivisions) %>% 
  map(~Censuses %>% semi_join(FocalSubdivisions[.x,], by = c("Year"))) -> 
  GroupList

GroupList %>% bind_rows() %>% nrow

# Social Phase ####

AMList <- list()

i <- 1

for(i in i:length(GroupList)){
  
  print(i)
  
  GroupList[[i]] -> Groups
  
  Groups %>% 
    dplyr::select(Id, GroupDate) %>% #unique %>% 
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
  filter(n>20) %>% pull(Year) %>% intersect(FocalYears) ->
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

names(AMList) <- FocalSubdivisions$Year

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

if(0){
  
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
  
  MCPListList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/MCPList.rds"))
  
  FullSpatialList <-
    list(#LifetimeDistance = LifetimeDistance,
      AnnualDistance = AnnualDistances,
      AnnualHRO = HROList,
      AnnualAreas = AreaList)
  
  FullSpatialList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/FullSpatialList.rds"))
  
}

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))
