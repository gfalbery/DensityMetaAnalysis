
# 01_Summarising System Replicates ####

{
  
  rm(list = ls())
  
  library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
  library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)
  
  theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))
  
  dir_create("Intermediate")
  
  if(file.exists("SavedSummaryTable.csv")) SummaryTable <- read.csv("SavedSummaryTable.csv")
  
  SummaryTable <- 
    gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1fTcR8cefcihKxO277VNXci_SKDPXEbV5exg6HVg42ks/edit#gid=1894078012") %>% 
    filter(!is.na(RName))
  
  SummaryTable %>% write.csv("SavedSummaryTable.csv", row.names = F)
  
  FocalSystems <- 
    "Datasets" %>% list.files %>% 
    intersect(SummaryTable$RName)
  
  a <- FocalSystems[2]
  
  NodeData <-
    FocalSystems %>%
    map(function(a){
      
      print(a)
      
      RepType <- SummaryTable %>% filter(RName == a) %>% pull(ReplicateType)
      
      ContactTypes <- SummaryTable %>% filter(RName == a) %>% pull(ContactTypes) %>% 
        str_split("_") %>% unlist
      
      DF <- glue("Datasets/{a}/Intermediate/CleanData.rds") %>% readRDS %>% 
        data.frame %>% 
        mutate_at(vars(contains("Year")), as.character) 
      
      if(RepType == "Site_Year" & (!any(str_detect(colnames(DF), "Site_Year")))) DF %<>% rename(Site_Year = SiteYear)
      
      DF$Rep <- DF[,RepType]
      
      if(RepType == "Site_Year"){
        
        DF %<>% separate(Site_Year, sep = "_", into = c("Site", "Year"))
        
      }else  if(RepType == "Year"){
        
        DF$Site <- "Site"
        
      }
      
      
      DF %<>% 
        rename_at(vars(contains("ID")), ~"Id") %>% 
        mutate_at("Id", as.character)
      
      if(is.null(DF$ContactType)){
        
        DF$ContactType <- ContactTypes
        
      }
      
      DF %<>%
        mutate_at(vars(contains("Year")), as.character) %>% 
        rename_at(vars(contains("ID")), ~"Id") %>% 
        mutate(Behaviour = ContactType) %>% 
        mutate(System = paste0(a)) %>% 
        mutate(SystemReplicate = paste0(as.character(a), "_", ContactType)) %>% 
        mutate(SpatialBehaviour = case_when(
          str_detect(Behaviour, "GOG|CoTrapping|Proximity") ~ "SpatialParameter",
          str_detect(Behaviour, "HRO") ~ "Spatial",
          TRUE ~ "NonSpatial"))
      
      DF %>% 
        return
      
    })
  
  names(NodeData) <- FocalSystems
  
  NodeData %<>%
    map(~mutate_at(.x, vars(matches("Site|Id|Rep")), as.character)) %>% 
    bind_rows()
  
  NodeData %<>% 
    mutate_at("Associations", ~.x - Observations)# %>% 
    # mutate_at("Degree", ~.x - 1) %>% 
    # mutate_at("Strength", ~.x - 1)
  
  ObsCorrect <- F
  MinObs <- 50
  
  # Combining with HRO pipeline ####
  
  a <- FocalSystems[1]
  
  SpatialData <- 
    FocalSystems %>%
    map(function(a){
      
      print(a)
      
      RepType <- SummaryTable %>% filter(RName == a) %>% pull(ReplicateType)
      
      DF <- glue("Datasets/{a}/Intermediate/CleanData.rds") %>% readRDS %>% 
        data.frame %>% 
        mutate_at(vars(contains("Year")), as.character) 
      
      if(RepType == "Site_Year" & (!any(str_detect(colnames(DF), "Site_Year")))) DF %<>% rename(Site_Year = SiteYear)
      
      DF$Rep <- DF[,RepType] %>% as.character()
      
      if(RepType == "Site_Year"){
        
        DF %<>% separate(Site_Year, sep = "_", into = c("Site", "Year"))
        
      }else if(RepType == "Year"){
        
        DF$Site <- "Site"
        
      }
      
      if(all(is.na(DF$Year))) DF$Year <- "Year"
      
      DF %<>% 
        rename_at(vars(contains("ID")), ~"Id") %>% 
        mutate_at("Id", as.character)
      
      SpatialList <- glue("Datasets/{a}/Intermediate/FullSpatialList.rds") %>% readRDS
      
      HRODF <- 
        SpatialList$AnnualHRO[SpatialList$AnnualHRO %>% map_lgl(~any(class(.x) == "matrix"))]
      
      HRODF %<>% 
        map(function(a){
          
          diag(a) <- 0
          
          a
          
        })
      
      names(HRODF) %<>% str_replace_all("[.]", "_")
      
      if(ObsCorrect){
        
        for(i in names(HRODF)){
          
          IncludeNodes <- DF %>% filter(Rep == i) %>% 
            filter(Obs > MinObs) %>% pull(Id) %>% intersect(rownames(HRODF[[i]]))
          
          HRODF[[i]] <- 
            HRODF[[i]][IncludeNodes, IncludeNodes]
          
        }
        
      }
      
      HRODF <- 
        # HRODF[HRODF %>% map_lgl(~(nrow(.x) > 0)|(length(.x)>1))]
        HRODF[HRODF %>% map_lgl(~(length(.x)>1))]
      
      HRODF %<>% 
        map(~data.frame(Associations = rowSums(.x), 
                        Strength = .x %>% rowSums,
                        Degree = rowSums(DegreeGet(.x))) %>% 
              rownames_to_column("Id")) %>% 
        bind_rows(.id = RepType)
      
      if(nrow(HRODF) > 30){
        
        HRODF$Rep <- HRODF[,RepType] %>% as.character()
        
        if(RepType == "Site_Year"){
          
          HRODF %<>% separate(Site_Year, sep = "_", into = c("Site", "Year"))
          
        }else  if(RepType == "Year"){
          
          HRODF$Site <- "Site"
          
        }
        
        if(all(is.na(HRODF$Year))) HRODF$Year <- "Year"
        
        DF %<>% 
          dplyr::select(-c(matches("Associations|Degree|Strength"))) %>% 
          left_join(HRODF %>% dplyr::select(Id, Rep, 
                                            matches("Associations|Degree|Strength"),
                                            all_of(setdiff(colnames(HRODF), colnames(DF)))),
                    .,
                    by = c("Id", "Rep"))
        
        DF <- 
          DF %>%
          mutate_at(vars(contains("Year")), as.character) %>% 
          mutate(Behaviour = "HRO") %>% 
          mutate(System = paste0(a)) %>% 
          mutate(SystemReplicate = paste0(as.character(a), "_", "HRO"))
        
        DF %>% 
          return
        
      }else NULL
      
    })
  
  SpatialData <- 
    SpatialData[!SpatialData %>% map_lgl(is.null)] %>%
    map(~mutate_at(.x, vars(matches("Site|Id|Rep")), as.character)) %>% 
    bind_rows()
  
  NodeData <- NodeData %>% bind_rows(SpatialData)
  
  NodeData %<>%
    mutate_at("Density.Annual", ~.x*(10^6))
  
  NodeData <- 
    NodeData %>% 
    count(System) %>% 
    mutate_at("n", ~(1/log10(n))) %>% 
    rename(SystemSizeAlpha = n) %>% 
    left_join(NodeData, ., by = "System")
  
  NodeData %<>% mutate_at("Associations", 
                          ~ifelse(SystemReplicate %in% 
                                    c("WaterDragons_CensusProximity", 
                                      "WildCrickets_Fight",
                                      "WildCrickets_Mate") & 
                                    is.na(.x), 0, .x))
  
  NodeData %<>% filter_at(vars(Density.Annual:Observations), ~(.x != Inf))
  
  NodeData %<>% mutate(UniqueID = paste(System, "_", Id))
  
  # NodeData %>%
  #   group_by(SystemReplicate) %>% 
  #   summarise(NAssoc = Prev(is.na(Associations)),
  #             NObs = Prev(is.na(Observations)),
  #             NX = Prev(is.na(X)),
  #             NY = Prev(is.na(Y)),
  #             NDensity = Prev(is.na(Density.Annual)),
  #             MinAssocs = min(Associations, na.rm = T), 
  #             Total = n()) %>% 
  #   data.frame
  
  NodeData %<>% filter(!SystemReplicate %in% c("PrairieDogs_HRO"))
  NodeData %<>% filter(!(SystemReplicate %in% c("RumDeer_Dominance") & Year > 2000))
  
  NodeData %>% 
    group_by(SystemReplicate, Site) %>% filter(!is.na(Site)) %>% 
    summarise(N = n(),
              NID = nunique(Id), 
              NYears = nunique(Year), 
              Area = diff(range(X, na.rm = T))*diff(range(Y, na.rm = T))) %>% 
    mutate(PopDensity = NID/Area/NYears) %>% 
    ungroup %>% group_by(SystemReplicate) %>% mutate(Sites = 1) %>% 
    summarise_at(c("N", "NID", "Area", "NYears", "PopDensity", "Sites"), 
                 list(Sum = ~sum(.x, na.rm = T), 
                      Mean = ~mean(.x, na.rm = T), 
                      Max = ~max(.x, na.rm = T))) %>% 
    dplyr::select(SystemReplicate, N = N_Sum, NID = NID_Sum, 
                  Sites = Sites_Sum,
                  Area = Area_Sum, 
                  NYears = NYears_Max, 
                  SiteYears = NYears_Sum,
                  PopDensity = PopDensity_Mean) %>% 
    write.csv("DerivedSystemTraits.csv", row.names = F)
  
  NodeData %<>% 
    group_by(System, Behaviour, Site) %>%
    count %>% filter(n > 10) %>% 
    semi_join(NodeData, ., by = c("System", "Behaviour", "Site"))
  
}

NodeData %<>% filter(!is.na(Density.Annual), !is.na(Strength))

NodeData %>% nrow

NodeData$UniqueID %>% nunique

NodeData$SystemReplicate %>% unique %>% sort

NodeData %>% saveRDS("Intermediate/NodeData.rds")

Factorise <- function(Number){
  
  1:Number %>% 
    
    map(~if(Number %% .x == 0){
      
      .x
      
    }) %>% unlist
  
}
