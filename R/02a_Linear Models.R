
# 02a_Linear Models #####

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

dir_create("Output/Linear")
dir_create("Output/LinearRank")
dir_create("Output/Spearman")

dir_create("Output/Spatial")

NodeData <- readRDS("Intermediate/NodeData.rds")

FocalSystem <- "RumDeer_Dominance"

Responses <- c("Associations", "Degree", "Strength")

(FocalSystems <- 
    NodeData %>% 
    dplyr::select(SystemReplicate, 
                  all_of(Responses),
                  Density.Annual,
                  X, Y) %>% 
    na.omit %>% 
    pull(SystemReplicate) %>% unique %>% sort)

(FocalSystems <-
  "Output/Linear" %>% list.files %>% str_remove(".rds$") %>% setdiff(FocalSystems, .))

# LM Models ####

for(FocalSystem in FocalSystems){
  
  print(FocalSystem)
  
  if(FocalSystem == "FieldingCows_HRO"){
    
    SubResponses <- Responses[c(1, 3)]
    
  }else SubResponses <- Responses
  
  SubsetDF <- 
    NodeData %>% 
    filter(SystemReplicate == FocalSystem) %>% 
    dplyr::select(all_of(SubResponses), 
                  Density.Annual,
                  Year, Site,
                  X, Y) %>% 
    na.omit %>% 
    group_by(Site) %>% 
    mutate_at(vars(1:Density.Annual), ~c(scale(.x))) %>%
    ungroup %>% 
    data.frame %>% na.omit
  
  SubsetDF %>% nrow %>% print
  
  if(nrow(SubsetDF) > 0){
    
    ModelList <- 
      
      SubResponses %>% 
      
      map(function(a){
        
        SubsetDF$Response <- SubsetDF[,a]
        
        lm(data = SubsetDF, Response ~ Density.Annual)
        
      })
    
    names(ModelList) <- SubResponses
    
    ModelList$Data <- SubsetDF
    
    ModelList %>% saveRDS(glue("Output/Linear/{FocalSystem}.rds"))
    
  }
}

# Log-Linear Models ####

dir_create("Output/LogLinear")

for(FocalSystem in FocalSystems){
  
  print(FocalSystem)
  
  if(FocalSystem == "FieldingCows_HRO"){
    
    SubResponses <- Responses[c(1, 3)]
    
  }else SubResponses <- Responses
  
  SubsetDF <- 
    NodeData %>% 
    filter(SystemReplicate == FocalSystem) %>% 
    dplyr::select(all_of(SubResponses), 
                  Density.Annual,
                  Year, Site,
                  X, Y) %>% 
    na.omit %>% 
    group_by(Site) %>% 
    # mutate_at(vars(1:Density.Annual), ~c(scale(.x))) %>% 
    mutate_at(vars(1:Density.Annual), ~c(log10(.x+1))) %>% 
    ungroup %>% 
    data.frame %>% na.omit
  
  SubsetDF %>% nrow %>% print
  
  if(nrow(SubsetDF) > 0){
    
    ModelList <- 
      
      SubResponses %>% 
      
      map(function(a){
        
        SubsetDF$Response <- SubsetDF[,a]
        
        lm(data = SubsetDF, Response ~ Density.Annual)
        
      })
    
    names(ModelList) <- SubResponses
    
    ModelList$Data <- SubsetDF
    
    ModelList %>% saveRDS(glue("Output/LogLinear/{FocalSystem}.rds"))
    
  }
}


# Rank Models ####

for(FocalSystem in FocalSystems){
  
  print(FocalSystem)
  
  SubsetDF <- 
    NodeData %>% 
    filter(SystemReplicate == FocalSystem) %>% 
    dplyr::select(all_of(SubResponses), 
                  Density.Annual,
                  Year, Site,
                  X, Y) %>% 
    na.omit %>% group_by(Site) %>% 
    mutate_at(vars(1:Density.Annual), ~c(scale(.x))) %>% 
    ungroup %>% 
    data.frame %>% na.omit
  
  SubsetDF %>% nrow %>% print
  
  if(nrow(SubsetDF) > 0){
    
    ModelList <- 
      SubResponses %>% 
      map(function(a){
        
        SubsetDF$Response <- SubsetDF[,a]
        
        SubsetDF %<>% group_by(Site) %>% 
          mutate_at("Response", ~.x %>% rank %>% scale %>% c) %>% 
          ungroup
        
        lm(data = SubsetDF, Response ~ Density.Annual)
        
      })
    
    names(ModelList) <- SubResponses
    
    ModelList$Data <- SubsetDF
    
    ModelList %>% saveRDS(glue("Output/LinearRank/{FocalSystem}.rds"))
    
  }
}

# Spearman's Rank ####

for(FocalSystem in FocalSystems){
  
  if(FocalSystem == "FieldingCows_HRO"){
    
    SubResponses <- Responses[c(1, 3)]
    
  }else SubResponses <- Responses
  
  print(FocalSystem)
  
  SubsetDF <- 
    NodeData %>% 
    filter(SystemReplicate == FocalSystem) %>% 
    dplyr::select(all_of(SubResponses), 
                  Density.Annual,
                  Year, Site,
                  X, Y) %>% 
    na.omit %>% group_by(Site) %>% 
    mutate_at(vars(1:Density.Annual), ~c(scale(.x))) %>% 
    ungroup %>% 
    data.frame %>% na.omit
  
  SubsetDF %>% nrow %>% print
  
  if(nrow(SubsetDF) > 0){
    
    ModelList <- 
      SubResponses %>% 
      map(function(a){
        
        SubsetDF$Response <- SubsetDF[,a]
        
        SubsetDF %<>% group_by(Site) %>% 
          mutate_at("Response", ~.x %>% rank %>% scale %>% c) %>% 
          ungroup
        
        lm(data = SubsetDF, Response ~ Density.Annual)
        
        cor(SubsetDF$Density.Annual, SubsetDF$Response, method = "spearman")
        
      })
    
    names(ModelList) <- SubResponses
    
    ModelList$Data <- SubsetDF
    
    ModelList %>% saveRDS(glue("Output/Spearman/{FocalSystem}.rds"))
    
  }
}
