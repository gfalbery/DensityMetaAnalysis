
# RumDeer ####

# https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs)

theme_set(theme_cowplot())

FocalSystem <- "RumDeer"

SummaryTable <- 
  gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1fTcR8cefcihKxO277VNXci_SKDPXEbV5exg6HVg42ks/edit#gid=1894078012") %>% 
  filter(!is.na(RName))

ContactTypes <- 
  SummaryTable %>% 
  filter(RName == FocalSystem) %>% 
  pull(ContactTypes) %>% 
  str_split("_") %>% 
  unlist

# Importing and Cleaning ####

Root <- paste0("Datasets/",FocalSystem)

dir_create(paste0(Root, "/Intermediate"))
dir_create(paste0(Root, "/Output"))

Censuses <- read.csv(paste0(Root, "/Behaviour/FullCensuses.csv"))

Censuses %>% 
  separate(Date, "/", into = c("Day", "Month", "Year")) %>% 
  dplyr::select(c("Day", "Month", "Year")) %>% 
  mutate_all(as.numeric) ->
  Censuses[,c("Day", "Month", "Year")]

Censuses %<>% 
  mutate(DeerYear = ifelse(Month<5, Year - 1, Year)) %>% 
  mutate(Season = ifelse(Year == DeerYear, "Rut", "Spring")) %>% 
  mutate(GroupDate = paste(Date, Group, sep=",")) %>% 
  filter(!is.na(Code), !Code == "") %>% 
  dplyr::select(c("Date","Code","Easting","Northing","GroupSize","Year","DeerYear","Season","GroupDate", GrazeType))

Individuals <- read_xlsx(paste0(Root, "/Phenotypes/tblLife.xlsx"))
Names <- read_xlsx(paste0(Root, "/Phenotypes/tblNames.xlsx"))

Individuals <- merge(Individuals, Names, by="Code",all.x=TRUE)

Individuals %<>% mutate_at(c("GivenName", "FamilyName"), as.character)

Individuals[is.na(Individuals$GivenName),"GivenName"] <- Individuals[is.na(Individuals$GivenName),"FamilyName"]
Individuals$Animal <- Individuals$GivenName
colnames(Individuals)[colnames(Individuals)=="Birth.Date"]<-"BirthYear"

Individuals$Sex<-cut(Individuals$Sex,breaks=c(0,1.5,2.5,3.5),labels=c("F","M","3"))

for(x in which(colnames(Individuals)=="BirthDay"):which(colnames(Individuals)=="DeathYear")){
  Individuals[,x]<-factor(as.factor(Individuals[,x]),levels=c(0:10000,paste(0,1:10,sep="")))
}

for(x in which(colnames(Individuals)=="BirthDay"):which(colnames(Individuals)=="DeathYear")){
  Individuals[(sapply(Individuals[,x],function(y) nchar(substr(y,1,2)))==1)&!is.na(Individuals[,x]),x]<-paste(0,Individuals[(sapply(Individuals[,x],function(y) nchar(substr(y,1,2)))==1)&!is.na(Individuals[,x]),x],sep="")
}

Individuals$Birth.Date<-as.character(factor(with(Individuals,paste(BirthDay,BirthMonth,BirthYear,sep="/"))))
Individuals$Death.Date<-as.character(factor(with(Individuals,paste(DeathDay,DeathMonth,DeathYear,sep="/"))))

Individuals[Individuals$Birth.Date=="NA/NA/NA","Birth.Date"]<-NA
Individuals[Individuals$Death.Date=="NA/NA/NA","Death.Date"]<-NA

Individuals$BirthYear<-as.numeric(as.character(Individuals$BirthYear))
Individuals$DeathYear<-as.numeric(as.character(Individuals$DeathYear))
Individuals$Name = Individuals$Code

Censuses <- merge(Censuses, Individuals[,c("Sex", "Name", "BirthYear")], 
                  by.x = c("Code"), by.y = c("Name"))

Censuses$Age <- with(Censuses, DeerYear - BirthYear)
Censuses$Hind <- with(Censuses, ifelse(Age > 2 & Sex == "F", "Y", "N"))

Censuses %<>% arrange(lubridate::dmy(Date))

Censuses %<>% filter(Hind == "Y")

FocalSeason <- c("Rut", "Spring")

Censuses %<>% filter(Season %in% FocalSeason)

Censuses %>% arrange(DeerYear) %>% 
  pull(DeerYear) %>% unique %>% 
  sort %>% as.character -> FocalYears

FocalYears %<>% #c(min(as.numeric(FocalYears))-1, .) %>% 
  as.numeric

# FocalYears <- FocalYears[2:length(FocalYears)]

Records = 5

FocalYears %>% map(~Censuses %>% filter(DeerYear == .x)) -> GroupList

ObservationDF <- Censuses %>% 
  dplyr::select(-Year) %>% 
  dplyr::select(Id = Code, X = Easting, Y = Northing, Year = DeerYear) %>% 
  na.omit

ObservationDF %<>% mutate_at(c("X", "Y"), ~.x/10) %>% 
  mutate_at(c("X", "Y"), ~.x - min(.x))

ObservationDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/ObservationDF.rds"))

Dominance <- read_xlsx(paste0(Root, "/Behaviour/frmDominance.xlsx"))

Dominance %<>% 
  rename(From = cboDomCode, To = cboSubCode) %>% 
  separate("txtDate", sep = "-", into = c("Year", "Month", "Day")) %>% 
  mutate_at(c("Year", "Month", "Day"), as.numeric) %>% 
  mutate(DeerYear = ifelse(Month<5, Year - 1, Year)) %>% 
  mutate(Season = ifelse(Year == DeerYear, "Rut", "Spring")) 

# Social Phase: GOG ####

AMList <- list()

i <- 1

for(i in i:length(GroupList)){
  
  print(i)
  
  GroupList[[i]] -> Groups
  
  Groups %>% 
    dplyr::select(Code, GroupDate) %>% unique %>% droplevels %>% 
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

names(AMList) <- 
  FocalYears

AMList %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList_GOG.rds"))

# Social Phase: Dominance ####

Dominance$DeerYear %>% 
  unique %>% sort %>% 
  map(~Dominance %>% filter(DeerYear == .x)) -> GroupList

AMList2 <- list()

i <- 1

for(i in i:length(GroupList)){
  
  print(i)
  
  GroupList[[i]] -> Groups
  
  Groups %>% 
    dplyr::select(From, To) %>% as.matrix %>% 
    # unite("Assoc", Group, Date, sep = "_") %>%
    # mutate(Assoc = GroupDate)
    graph_from_edgelist(directed = F) %>% 
    as_tbl_graph %>% activate(edges) %>% mutate(weight = 1) %>% 
    simplify(edge.attr.comb = "sum") -> 
    
    SocGraph
  
  AM <- SocGraph %>% get.adjacency(attr = "weight") %>% as.matrix
  
  Observations = rowSums(AM)
  
  Matrix <- AM # %>% ProportionalMatrix(Observations)
  
  diag(Matrix) <- Observations
  
  AMList2[[i]] <- Matrix
  
}

names(AMList2) <- 
  Dominance$DeerYear %>% 
  unique %>% sort %>% 
  as.character

AMList2 %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/AMList_Dominance.rds"))

# AMList <- AMList2

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

AMList %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% mutate_at("Year", as.numeric) -> SocialNodeTraits

AMList2 %>% 
  map(~data.frame(Associations = rowSums(.x), 
                  Strength = .x %>% ProportionalMatrix(Observations = diag(.x)) %>% rowSums,
                  Degree = rowSums(DegreeGet(.x)),
                  Observations = diag(.x)) %>% 
        rownames_to_column("Id")) %>% 
  bind_rows(.id = "Year") %>% mutate_at("Year", as.numeric) -> SocialNodeTraits2

SocialNodeTraits %<>%
  full_join(NodeDF) %>% 
  mutate(ContactType = "GOG") %>% 
  # mutate_at(paste0(c("Associations", "Strength", "Degree")),
  # ~replace_na(.x, 0)
  # ) %>% 
  bind_rows(SocialNodeTraits2 %<>%
              full_join(NodeDF) %>% 
              mutate(ContactType = "Dominance")# %>% 
            # mutate_at(paste0(c("Associations", "Strength", "Degree")),
            # ~replace_na(.x, 0))
  )

NodeDF <- SocialNodeTraits

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
  
  AreaList[[i]] <- AnnualAreas
  
  ZeroAreas <- AnnualAreas %>% filter(Area == 0) %>% pull(Id)
  
  Dyads <- expand.grid(Id1 = MCPIndividuals %>% setdiff(ZeroAreas), 
                       Id2 = MCPIndividuals %>% setdiff(ZeroAreas))
  
  if(nrow(Dyads) > 0){
    
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

NodeDF %>% dplyr::select(Year, Id) %>% 
  # unique %>%
  nrow

NodeDF %<>% 
  filter(!(ContactType == "Dominance" & Year > 2000))

NodeDF %>% saveRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF <-
  readRDS(glue("Datasets/{FocalSystem}/Intermediate/CleanData.rds"))

NodeDF %>% filter(ContactType == "GOG") %>% is.na %>% colSums

NodeDF %>% filter(ContactType == "Dominance") %>% is.na %>% colSums

NodeDF %>% group_by(Year) %>% summarise_at(c("Density.Annual", "Associations"), ~Prev(is.na(.x)))
