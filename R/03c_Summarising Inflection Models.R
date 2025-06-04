
# 3c_Summarising Inflection Models.R ####

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)
library(mgcv); library(MASS); library(coda); library(MCMCglmm)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

NodeData <- readRDS("Intermediate/NodeData.rds")

# MetaDF <- read.csv("MetaDF.csv")

Responses <- c("Associations", "Degree", "Strength")
# Response <- Responses[3]

(FocalSystems <- "Output/InflectionModels" %>% list.files %>% str_remove(".rds$"))

LMList <- "Output/InflectionModels" %>% dir_ls %>% map(readRDS)

names(LMList) <- FocalSystems

FocalSystem <- FocalSystems[1]

ReturnList <- 
  FocalSystems %>% 
  map(function(FocalSystem){
    
    names(LMList[[FocalSystem]]) %>% intersect(Responses) %>% 
      
      map(function(FocalResponse){
        
        ReturnDF <- 
          
          LMList[[FocalSystem]][[FocalResponse]] %>% 
          
          map(function(Model){
            
            Model %>% summary %>% extract2("coefficients") %>% 
              data.frame %>% dplyr::select(Density.Annual = 1, SE = 2, P = 4) %>% slice(2) %>% 
              rownames_to_column %>% 
              bind_cols(
                Model %>% confint %>% data.frame %>% rename(Lower = 1, Upper = 2) %>% slice(2)
              )
            
          }) %>% reduce(~left_join(.x, .y, by = c("rowname")))
        
        ReturnDF %<>% rename_all(~str_replace_all(.x, c("[.]x" = "_First", "[.]y" = "_Middle")))
        
        ReturnDF %<>% rename_at(vars(Density.Annual:Upper), ~paste0(.x, "_Last"))
        
        ReturnDF
        
      }) %>% bind_rows(.id = "Response") %>% 
      
      mutate_at("Response", ~Responses[as.numeric(.x)])
    
  })

InflectionPointDF <- ReturnList %>% bind_rows(.id = "SystemReplicate")

InflectionPointDF <- 
  1:nrow(InflectionPointDF) %>% 
  map(function(i){
    
    Differences <- 
      rnorm(1000, InflectionPointDF[i, "Density.Annual_First"], InflectionPointDF[i, "SE_First"]) - 
      rnorm(1000, InflectionPointDF[i, "Density.Annual_Last"], InflectionPointDF[i, "SE_Last"]) %>% 
      as.mcmc
    
    HPD <- HPDinterval(Differences)
    
    data.frame(Mean_Saturation = posterior.mode(Differences), 
               Lower_Saturation = HPD[1], 
               Upper_Saturation = HPD[2])
    
  }) %>% bind_rows %>% bind_cols(InflectionPointDF, .)

InflectionPointDF %>% saveRDS("Output/InflectionPointDF.rds")

System <- ModelList[[1]]

FocalSystem <- FocalSystems[1]

SaturationEstimates <- 
  
  FocalSystems %>% 
  
  map(function(FocalSystem){
    
    print(FocalSystem)
    
    SubList <- 
      
      names(LMList[[FocalSystem]]) %>% intersect(Responses) %>% 
      
      map(function(FocalResponse){
        
        preds <- LMList[[FocalSystem]][[FocalResponse]][[1]] %>% fitted
        actual <- LMList[[FocalSystem]][[FocalResponse]][[1]]$model$Response
        rss <- sum((preds - actual) ^ 2)  ## residual sum of squares
        tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
        rsq <- 1 - rss/tss
        
        SubDF <-
          actual <- LMList[[FocalSystem]][[FocalResponse]][[1]] %>% summary() %>% extract2("coefficients") %>% data.frame %>% 
          slice(2) %>% 
          dplyr::select(Estimate, SE = 2, t = 3) %>% 
          mutate(R2 = rsq)
        
        ReturnDF1 <- 
          LMList[[FocalSystem]][[FocalResponse]][[1]] %>% confint %>% data.frame %>% slice(2) %>% 
          rename(LinearLower = 1, LinearUpper = 2) %>% 
          bind_cols(SubDF, .) %>%
          # rename_all(~paste0(.x, "_First"))
          mutate(Portion = "First")
        
        preds <- LMList[[FocalSystem]][[FocalResponse]][[3]] %>% fitted
        actual <- LMList[[FocalSystem]][[FocalResponse]][[3]]$model$Response
        rss <- sum((preds - actual) ^ 2)  ## residual sum of squares
        tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
        rsq <- 1 - rss/tss
        
        SubDF <-
          LMList[[FocalSystem]][[FocalResponse]][[3]] %>% summary() %>% extract2("coefficients") %>% data.frame %>% 
          slice(2) %>% 
          dplyr::select(Estimate, SE = 2, t = 3) %>% 
          mutate(R2 = rsq)
        
        ReturnDF2 <- 
          LMList[[FocalSystem]][[FocalResponse]][[3]] %>% confint %>% data.frame %>% slice(2) %>% 
          rename(LinearLower = 1, LinearUpper = 2) %>% 
          bind_cols(SubDF, .) %>%
          # rename_all(~paste0(.x, "_Last"))
          mutate(Portion = "Last")
        
        ReturnDF1 %>% bind_rows(ReturnDF2)
        
      })
    
    names(SubList) <- names(LMList[[FocalSystem]]) %>% intersect(Responses)
    
    SubList %<>% bind_rows(.id = "Response")
    
  })

names(SaturationEstimates) <- FocalSystems

SaturationEstimates %<>% bind_rows(.id = "SystemReplicate")

saveRDS(SaturationEstimates, "Output/SaturationEstimates.rds") # Next run 05_Meta-analysis.R

# Rank models ####


if(0){
  
  ModelList <- "Output/InflectionRank" %>% dir_ls %>% map(readRDS)
  
  names(ModelList) <- FocalSystems
  
  ReturnList <- 
    ModelList %>% 
    map(function(LMList){
      
      ReturnDF <- 
        LMList %>% 
        map(function(Model){
          
          Model %>% summary %>% extract2("coefficients") %>% data.frame %>% dplyr::select(Density.Annual = 1, SE = 2, P = 4) %>% slice(2) %>% 
            rownames_to_column %>% 
            bind_cols(
              Model %>% confint %>% data.frame %>% rename(Lower = 1, Upper = 2) %>% slice(2)
            )
          
        }) %>% reduce(~left_join(.x, .y, by = c("rowname")))
      
      ReturnDF %<>% rename_all(~str_replace_all(.x, c("[.]x" = "_First", "[.]y" = "_Middle")))
      
      ReturnDF %<>% rename_at(vars(Density.Annual:Upper), ~paste0(.x, "_Last"))
      
      ReturnDF
      
    })
  
  InflectionPointDF <- ReturnList %>% bind_rows(.id = "SystemReplicate")
  
  InflectionPointDF <- 
    1:nrow(InflectionPointDF) %>% 
    map(function(i){
      
      Differences <- 
        rnorm(1000, InflectionPointDF[i, "Density.Annual_First"], InflectionPointDF[i, "SE_First"]) - 
        rnorm(1000, InflectionPointDF[i, "Density.Annual_Last"], InflectionPointDF[i, "SE_Last"]) %>% 
        as.mcmc
      
      HPD <- HPDinterval(Differences)
      
      data.frame(Mean_Saturation = posterior.mode(Differences), 
                 Lower_Saturation = HPD[1], 
                 Upper_Saturation = HPD[2])
      
    }) %>% bind_rows %>% bind_cols(InflectionPointDF, .)
  
  InflectionPointDF %>% saveRDS("Output/InflectionPointRankDF.rds")
  
  
  
  if(0){
    
    # The rest is plotting
    
    InflectionPointDF %>% left_join(MetaDF) %>% 
      mutate(Sig = as.numeric((Lower_Saturation*Upper_Saturation) > 0 )) %>%
      arrange(Mean_Saturation) %>% 
      mutate_at("SystemReplicate", ~factor(.x, levels = unique(.x))) %>% 
      ggplot(aes(SystemReplicate, Mean_Saturation)) +
      geom_hline(yintercept = 0, lty = 2, colour = "light grey") +
      geom_errorbar(aes(colour = Behaviour,
                        ymin = Lower_Saturation, ymax = Upper_Saturation, alpha = as.factor(Sig)), 
                    width = 0.3)  +
      geom_point(colour = "white", size = 3.5, alpha = 0.8) +
      geom_point(aes(colour = Behaviour, shape = Behaviour), size = 2.5) +
      # coord_flip() + 
      theme(axis.line.y = element_blank(), 
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank()) +
      theme(legend.position = c(0, 1), legend.justification = c(0, 1)) +
      labs(y = "Density effect estimate") +
      scale_alpha_manual(values = c(0.3, 1), guide = "none") +
      # scale_shape_manual(#name = "Contact") +
      scale_colour_manual(#name = "Contact", 
        values = rev(c(alpha(ParasiteColours[[3]], 1), AlberColours[[2]], AlberColours[[1]]))) +
      NULL
    
    InflectionPointDF %>% 
      mutate(Sig = as.numeric((Lower_Saturation*Upper_Saturation) > 0 )) %>%
      arrange(Mean_Saturation) %>% 
      mutate_at("SystemReplicate", ~factor(.x, levels = unique(.x))) %>% 
      ggplot(aes(Behaviour, Mean_Saturation)) +
      geom_boxplot() +
      geom_sina()
    
    # These are all plots ####
    
    PlotDF <- 
      c("_First", "_Middle", "_Last") %>% 
      map(~InflectionPointDF %>% 
            pivot_longer(contains(.x))) %>% 
      bind_rows() %>% data.frame %>% 
      separate(name, sep = "_", into = c("Metric", "Portion")) %>% 
      dplyr::select(-contains("_"), -c(rowname:Upper, SE)) %>% 
      pivot_wider(names_from = "Metric", values_from = "value") %>% 
      data.frame
    
    PlotDF %>% 
      ggplot(aes(SystemReplicate, Density.Annual, colour = Portion)) + 
      geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(w = 0.4)) +
      geom_point(position = position_dodge(w = 0.4)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    PlotDF %>%
      filter(Density.Annual < 3) %>% 
      group_by(SystemReplicate) %>% #mutate_at("Density.Annual", ~scales::rescale(.x, to = c(0, 1))) %>% 
      mutate_at("Portion", ~factor(.x, levels = c("First", "Middle", "Last")) %>% as.numeric) %>% 
      ggplot(aes(Portion, Density.Annual)) + 
      geom_point() +
      geom_smooth(method = lm)
    
    PlotDF %>%
      # filter(Density.Annual < 3) %>% 
      group_by(SystemReplicate) %>% #mutate_at("Density.Annual", ~scales::rescale(.x, to = c(0, 1))) %>% 
      mutate_at("Portion", ~factor(.x, levels = c("First", "Middle", "Last"))) %>%
      ggplot(aes(Portion, Density.Annual)) + 
      geom_boxplot() +
      geom_sina() +
      labs(x = "Data portion", y = "Density slope")
    
    XList <- list(0:1/2, 1:2/2)
    
    InflectedLines <- 
      1:nrow(InflectionPointDF) %>% 
      map(function(i){
        
        First <- InflectionPointDF[i, "Density.Annual_First"]*XList[[1]]
        
        Second <- InflectionPointDF[i, "Density.Annual_Last"]*XList[[1]]
        
        Second <- Second - Second[1] + First[2]
        
        ReturnDF <- 
          data.frame(X = XList %>% unlist, 
                     Y = c(First, Second),
                     Portion = rep(c("First",
                                     "Last"), XList %>% map_dbl(length)))
        
        ReturnDF[-2, ] %>% 
          return
        
      })
    
    names(InflectedLines) <- FocalSystems 
    
    InflectedLines %<>% bind_rows(.id = "SystemReplicate")
    
    PolygonDF <- data.frame(X = c(0, 0.5, 1), 
                            Y = c(0, 1, 1))
    
    InflectedLines %>% 
      group_by(SystemReplicate) %>% mutate_at(c("X", "Y"), ~scales::rescale(.x, to = c(0, 1))) %>% 
      # slice(c(1, 2, 4)) %>% 
      # ggplot(aes(jitter(X, amount = 0.1),  Y)) + #jitter(Y, amount = 0.1))) + 
      ggplot(aes(X, Y)) +
      geom_polygon(data = PolygonDF, 
                   aes(X, Y),
                   fill = "white", colour = "black", linewidth = 1) +
      geom_line(alpha = 0.3, 
                aes(group = paste0(SystemReplicate), colour = Portion)) +
      # theme_void() + 
      theme(legend.position = "none") +
      labs(x = "Relative density", y = "Relative contact rates") + 
      theme(axis.text = element_blank()) +
      scale_colour_manual(values = c(AlberColours[[2]], AlberColours[[1]]))
    
    ggsave("Figures/RelativeInflections.jpeg", units = "mm", height = 80, width = 80)
    
    # Trying this with a fan ####
    
    InflectedLines %>% 
      group_by(SystemReplicate) %>% mutate_at(c("X", "Y"), ~scales::rescale(.x, to = c(0, 1)))
    
    map(FocalSystems, function(i){
      
      expand.grid(From = c("First", "Middle", "Last"), 
                  To = c("First", "Middle", "Last")) %>% 
        slice(4, 8) %>% 
        mutate(SystemReplicate = i)
      
    }) %>% bind_rows() -> EdgeDF
    
    library(igraph); library(tidygraph); library(ggraph)
    
    # EdgeDF %>% left_join(InflectedLines)
    
    EdgeDF %>% graph_from_data_frame() %>% as_tbl_graph() %>% 
      ggraph(layout = as.matrix(InflectedLines[,c("X", "Y")])) + geom_edge_fan()
    
    # Once more with the middle ####
    
    XList <- list(0:1/3, 1:2/3, 2:3/3)
    
    InflectedLines <- 
      1:nrow(InflectionPointDF) %>% 
      map(function(i){
        
        First <- InflectionPointDF[i, "Density.Annual_First"]*XList[[1]]
        
        Second <- InflectionPointDF[i, "Density.Annual_Middle"]*XList[[2]]
        
        Second <- Second - Second[1] + First[2]
        
        Third <- InflectionPointDF[i, "Density.Annual_Last"]*XList[[3]]
        
        Third <- Third - Third[1] + Second[2]
        
        ReturnDF <- 
          data.frame(X = XList %>% unlist, 
                     Y = c(First, Second, Third),
                     Portion = rep(c("First",
                                     "Middle",
                                     "Last"), XList %>% map_dbl(length)))
        
        ReturnDF[-2, ] %>% 
          return
        
      })
    
    names(InflectedLines) <- FocalSystems 
    
    InflectedLines %<>% bind_rows(.id = "SystemReplicate")
    
    InflectedLines %>% 
      group_by(SystemReplicate) %>% mutate_at(c("X", "Y"), ~scales::rescale(.x, to = c(0, 1))) %>% 
      # slice(c(1, 2, 4)) %>% 
      # ggplot(aes(jitter(X, amount = 0.1),  Y)) + #jitter(Y, amount = 0.1))) + 
      ggplot(aes(X, Y)) +
      # geom_polygon(data = PolygonDF, 
      #              aes(X, Y),
      #              fill = "white", colour = "black", linewidth = 1) +
      geom_line(alpha = 0.3, 
                aes(group = paste0(SystemReplicate), colour = Portion)) +
      # theme_void() + 
      theme(legend.position = "none") +
      labs(x = "Relative density", y = "Relative contact rates") + 
      theme(axis.text = element_blank()) +
      scale_colour_manual(values = c(AlberColours[[2]], AlberColours[[1]], AlberColours[[3]]))
    
  }
  
}