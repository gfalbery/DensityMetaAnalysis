
# 02b_GAMMs #####

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

dir_create("Output/GAMM")

NodeData <- readRDS("Intermediate/NodeData.rds")

Responses <- c(#"Associations", "Degree", 
               "Strength")

(FocalSystems <- 
    NodeData %>% 
    dplyr::select(SystemReplicate, 
                  all_of(Responses),
                  Density.Annual,
                  X, Y) %>% 
    # na.omit %>% 
    pull(SystemReplicate) %>% unique %>% sort)

(FocalSystems <-
    "Output/GAMM" %>% list.files %>% str_remove(".rds$") %>% setdiff(FocalSystems, .))

(FocalSystem <- FocalSystems[1])

for(FocalSystem in FocalSystems){
  
  print(FocalSystem)
  
  # if(FocalSystem == "FieldingCows_HRO"){
  #   
  #   SubResponses <- Responses[c(1, 3)]
  #   
  # }else 
  
  SubResponses <- Responses
  
  SubsetDF <- 
    NodeData %>% 
    filter(SystemReplicate == FocalSystem) %>% 
    dplyr::select(
      all_of(SubResponses),
      Density.Annual,
      Year, Site,
      X, Y) %>% 
    na.omit %>%
    group_by(Site) %>% na.omit %>% filter_all(~ (!.x == Inf)) %>% 
    mutate_at(vars(1:Density.Annual), ~c(scale(.x))) %>% 
    ungroup %>% na.omit %>% 
    data.frame
  
  SubsetDF %>% nrow %>% print
  
  ModelList <- 
    SubResponses %>% 
    map(function(a){
      
      BAM1 <- BAMModelAdd(Data = SubsetDF, 
                          Explanatory = c("Density.Annual"),
                          Response = a,
                          Base = T,
                          Beep = F,
                          Gamma = 1.4,
                          Scale = F) %>% 
        list()
      
      glue::glue("s(Density.Annual, k = {3:7})") %>%
        
        map(function(E){
          
          BAMModelAdd(Data = SubsetDF,
                      Explanatory = E,
                      Response = a,
                      Base = T,
                      Beep = F,
                      Gamma = 1.4,
                      Scale = F)
          
        }) %>% append(BAM1, .)
      
    })
  
  names(ModelList) <- SubResponses
  
  ModelList %>% saveRDS(glue("Output/GAMM/{FocalSystem}.rds"))
  
}
