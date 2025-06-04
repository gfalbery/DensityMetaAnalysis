
# 2c_InflectionModels ####

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

dir_create("Output/InflectionModels")
dir_create("Output/InflectionRank")

NodeData <- readRDS("Intermediate/NodeData.rds")

FocalSystem <- "RumDeer_Dominance"

Responses <- c("Associations", "Degree", "Strength")

(FocalSystems <- 
    NodeData %>% 
    dplyr::select(SystemReplicate, 
                  # Response,
                  all_of(Responses),
                  Density.Annual,
                  # Observations,
                  X, Y) %>% 
    na.omit %>% 
    pull(SystemReplicate) %>% unique %>% sort)

SplitPoints <- 
  FocalSystems %>% #setdiff("KenyanElephants_GOG") %>% 
  map(function(FocalSystem){
    
    print(FocalSystem)
    
    LMList <- 
      Responses %>% map(function(FocalResponse){
        
        print(FocalResponse)
        
        TestDF <- 
          NodeData %>% 
          filter(SystemReplicate == FocalSystem) %>% 
          # mutate_at("Observations", ~c(scale(log(.x + 1)))) %>% 
          mutate_at(c("Density.Annual", FocalResponse), ~c(scale((.x)))) %>% 
          arrange(Density.Annual)
        
        TestDF$Response <- TestDF[,FocalResponse]
        
        TestDF1 <- TestDF %>% 
          slice(1:round((n()/2)))
        
        TestDF2 <- TestDF %>% 
          slice(round((n()/4)*1):round(n()/4*3))
        
        TestDF3 <- TestDF %>% 
          slice(round((n()/2)*1):n())
        
        list(TestDF1, TestDF2, TestDF3) %>% 
          map_dbl(nrow) %>% print
        
        list(TestDF1, TestDF2, TestDF3) %>%
          # map(function(a) a %>% mutate_at(c("Associations"), ~c(scale(rank(.x))))) %>%
          # map(function(a) a %>% mutate_at(c("Associations"), ~c(scale((.x))))) %>%
          # map(function(a) a %>% mutate_at(c("Associations"), ~c(rank((.x))))) %>%
          map(~lm(data = .x, Response ~ Density.Annual))
        
      })
    
    names(LMList) <- Responses
    
    LMList %>% saveRDS(glue("Output/InflectionModels/{FocalSystem}.rds"))
    
  })

SplitPoints <- 
  FocalSystems %>% setdiff("KenyanElephants_GOG") %>% 
  map(function(FocalSystem){
    
    print(FocalSystem)
    
    TestDF <- 
      NodeData %>% 
      filter(SystemReplicate == FocalSystem) %>% 
      mutate_at(c("Density.Annual", "Associations"), ~c(scale((.x)))) %>% 
      arrange(Density.Annual)
    
    TestDF1 <- TestDF %>% 
      slice(1:round((n()/2)))
    
    TestDF2 <- TestDF %>% 
      slice(round((n()/4)*1):round(n()/4*3))
    
    TestDF3 <- TestDF %>% 
      slice(round((n()/2)*1):n())
    
    list(TestDF1, TestDF2, TestDF3) %>% 
      map_dbl(nrow) %>% print
    
    LMList <- 
      list(TestDF1, TestDF2, TestDF3) %>%
      map(function(a) a %>% mutate_at(c("Associations"), ~c(scale(rank(.x))))) %>%
      # map(function(a) a %>% mutate_at(c("Associations"), ~c(scale((.x))))) %>%
      # map(function(a) a %>% mutate_at(c("Associations"), ~c(rank((.x))))) %>%
      map(~lm(data = .x, Associations ~ Density.Annual))
    
    LMList %>% saveRDS(glue("Output/InflectionRank/{FocalSystem}.rds"))
    
  })

