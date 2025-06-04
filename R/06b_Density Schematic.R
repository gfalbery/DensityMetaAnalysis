
# X_Density Schematic ####

library(tidyverse); library(zip); library(data.table); library(ggregplot); library(magrittr); library(igraph)
library(adehabitatHR); library(rgeos); library(glue); library(cowplot); library(patchwork); library(readxl)
library(fs); library(colorspace)

theme_set(theme_cowplot())

FocalSystem <- "RumDeer"

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

FocalYears2 <- FocalYears[c(1, 11, 21, 41)]

SPDF <- SpatialPointsDataFrame(data = AnnualCentroids[,c("X", "Y", "Year")], 
                               coords = AnnualCentroids[,c("X", "Y")])

SPDF <- SPDF[,"Year"]

KUDL <- kernelUD(SPDF, 
                 extent = 0,
                 same4all = TRUE, 
                 grid = 500)

a <- 1

1:length(FocalYears) %>% lapply(function(a){
  
  print(FocalYears[a])
  
  if(FocalYears[a] %in% FocalYears2){
    
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
    
  }else NULL
  
}) -> DensityList

DensityList %>% bind_rows -> AnnualCentroids

NodeDF <- 
  AnnualCentroids 

# Puting the figure together ####

FocalYears2 <- FocalYears[c(1, 11, 21, 41)]

XLim <- c(min(NodeDF$X) - 0.1, max(NodeDF$X) + 0.1)
YLim <- c(min(NodeDF$Y) + 0.5, max(NodeDF$Y) + 0.25)

PanelA <- 
  ObservationDF %>% 
  # RandomSlice(2500) %>% 
  filter(Year %in% FocalYears2) %>% 
  ggplot(aes(X, Y)) + geom_point(alpha = 0.2) +
  coord_fixed(xlim = XLim, ylim = YLim)

PanelB <- 
  AnnualCentroids %>% 
  filter(Year %in% FocalYears2) %>% 
  ggplot(aes(X, Y)) + geom_point(alpha = 0.4) +
  coord_fixed(xlim = XLim, ylim = YLim)

(PanelC <- 
    KUDL[FocalYears2] %>% 
    map(GetKUDL) %>% bind_rows(.id = "Year") %>% 
    group_by(Year) %>% 
    mutate_at("Density", ~scales::rescale(.x, to = c(0, 1))) %>% 
    ungroup %>% 
    filter(Density > 0.01) %>% 
    ggplot(aes(X, Y)) + geom_tile(aes(fill = Density)) +
    coord_fixed(xlim = XLim, ylim = YLim) +
    geom_point(data = NodeDF %>% 
                 filter(Year %in% FocalYears2),
               aes(shape = as.factor(Year)),
               alpha = 0.2) +
    guides(shape = "none") +
    facet_wrap(~Year) + labs(shape = "Year"))

PanelD <- 
  NodeDF %>% rename(Density = Density.Annual) %>% 
  filter(Year %in% FocalYears2) %>% 
  mutate_at("Density", ~scales::rescale(.x, to = c(0, 1))) %>% 
  ggplot(aes(X, Y)) + 
  geom_point(aes(colour = Density,
                 shape = as.factor(Year)),
             alpha = 0.6) +
  coord_fixed(xlim = XLim, ylim = YLim) + 
  labs(shape = "Year")

PlotList <- list(PanelA, PanelB, PanelC, PanelD)

PlotList[[1]] <- PlotList[[1]] + ggtitle("Raw observations")
PlotList[[2]] <- PlotList[[2]] + ggtitle("Annual centroids")
PlotList[[3]] <- PlotList[[3]] + ggtitle("Density distributions")
PlotList[[4]] <- PlotList[[4]] + ggtitle("Individual values")

PlotList %>% 
  map(~.x + theme_void() + 
        theme(plot.title = element_text(hjust = 0.5)) +
        # scale_fill_continuous_sequential(palette = "Terrain", limits = c(0, 1)) +
        # scale_colour_continuous_sequential(palette = "Terrain", limits = c(0, 1)) +
        scale_fill_continuous_sequential(palette = AlberPalettes[[2]], limits = c(0, 1)) +
        scale_colour_continuous_sequential(palette = AlberPalettes[[2]], limits = c(0, 1)) +
        theme(strip.text = element_blank()) +
        theme(strip.background = element_rect(fill = "white", colour = "white"))) %>% 
  ArrangeCowplot + 
  plot_layout(nrow = 1) + 
  plot_layout(guides = "collect") + 
  # & theme(legend.position = "bottom")# +
  # plot_annotation(tag_levels = "A") +
  NULL

ggsave("Figures/DensitySchematic.jpeg", units = "mm", height = 100, width = 250, dpi = 300)
