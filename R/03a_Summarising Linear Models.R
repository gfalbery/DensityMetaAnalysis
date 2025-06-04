
# 03a_Summarising Linear Models #####

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

NodeData <- readRDS("Intermediate/NodeData.rds")

Responses <- c("Associations", "Degree", "Strength")

# Response <- Responses[3]

(FocalSystems <- "Output/Linear" %>% list.files %>% str_remove(".rds$"))

LMList <- "Output/Linear" %>% dir_ls %>% map(readRDS)

names(LMList) <- FocalSystems

# For LM meta-analysis

System <- LMList[[1]]

LinearEstimates <- 
  LMList %>% 
  map(function(System){
    
    SubList <- 
      Responses %>% intersect(names(System)) %>% 
      map(function(FocalResponse){
        
        preds <- System[[FocalResponse]] %>% fitted
        actual <- System$Data[,FocalResponse]
        rss <- sum((preds - actual) ^ 2)  ## residual sum of squares
        tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
        rsq <- 1 - rss/tss
        
        SubDF <-
          System[[FocalResponse]] %>% summary() %>% extract2("coefficients") %>% data.frame %>% 
          slice(2) %>% 
          dplyr::select(Estimate, SE = 2, t = 3) %>% 
          mutate(R2 = rsq)
        
        System[[FocalResponse]] %>% confint %>% data.frame %>% slice(2) %>% 
          rename(LinearLower = 1, LinearUpper = 2) %>% 
          bind_cols(SubDF, .)
        
      })
    
    names(SubList) <- Responses %>% intersect(names(System))
    
    SubList %>% bind_rows(.id = "Response")# %>% 
    # mutate_at("Response", ~Responses[as.numeric(.x)])
    
  }) %>% bind_rows(.id = "SystemReplicate")

saveRDS(LinearEstimates, "Output/LinearEstimates.rds") # Next run 05_Meta-analysis.R

# Check the models ####

if(0){
  
  LMList %>% 
    map(1) %>% 
    map(resid) %>% unlist %>% data.frame %>% rownames_to_column("SystemReplicate") %>% rename(Res = 2) %>% 
    ggplot(aes(Res)) + geom_density(aes(group = SystemReplicate)) +
    theme(legend.position = "none")
  
}
