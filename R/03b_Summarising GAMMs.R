
# 03c_Summarising GAMMs with Derivatives

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)
library(mgcv); library(MASS); library(MCMCglmm)

# source("https://gist.githubusercontent.com/gavinsimpson/ca18c9c789ef5237dbc6/raw/295fc5cf7366c831ab166efaee42093a80622fa8/derivSimulCI.")

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

NodeData <- readRDS("Intermediate/NodeData.rds")

Responses <- c(#"Associations", "Degree", 
               "Strength")

(FocalSystems <- "Output/GAMM" %>% list.files %>% str_remove(".rds$"))

BAMList <- "Output/GAMM" %>% dir_ls %>% map(readRDS)

names(BAMList) <- FocalSystems

BestFitDF <- 
  BAMList %>% 
  map(function(b){
    
    Responses %>% map(function(a){
      
      Devs <- b[[a]] %>% map("FinalModel") %>% map_dbl("aic")
      
      if(length(Devs) > 0){
        
        SubDF <- data.frame(BestFit = c("Linear", glue::glue("K{3:7}"))[which(Devs == min(Devs))],
                            Delta = abs(Devs[1] - min(Devs)))
        
        SubDF %<>% mutate_at("BestFit", ~ifelse((Delta < 2 & BestFit != "Linear"), "Linear", BestFit))
        
        SubDF[,c("LinearDeviance", glue::glue("K{3:7}Deviance"))] <- Devs
        
        SubDF %>% mutate(Response = a) %>% return
        
      }
      
    }) %>% bind_rows#(.id = "Response")
    
  }) %>% bind_rows(.id = "SystemReplicate")

BestFitDF %<>% 
  mutate(NumBestFit = as.numeric(BestFit != "Linear"))

ChosenModels <- 
  1:nrow(BestFitDF) %>% 
  map(function(a){
    
    BAMList[[BestFitDF[a, "SystemReplicate"]]][[BestFitDF[a, "Response"]]][[3]] # Selecting K = 4
    
  })

names(ChosenModels) <- FocalSystems

# Looking at the splines ####

SaturationList <- list()
SaturationList[1:length(ChosenModels)] <- list(list())

GAMInflection <- list()
GAMInflection[1:length(ChosenModels)] <- list(list())

a <- b <- 1

PredictedData <-
  1:length(ChosenModels) %>%
  map(function(b){
    
    print(paste0(BestFitDF$SystemReplicate[b], "_", BestFitDF$Response[b]))
    
    # a <- ChosenModels[[b]]
    a <- b
    
    ReturnList <- 
      BAMList[[BestFitDF[a, "SystemReplicate"]]][[BestFitDF[a, "Response"]]][2:6] %>% 
      map(function(a){
        
        PredDF <- 
          MakePredictDF(a$Data %>% 
                          na.omit %>% 
                          dplyr::select(Density.Annual)) %>%
          bind_cols() %>% data.frame
        
        mod <- a$FinalModel
        
        eps = 1e-7
        
        Xp <- (predict(mod, PredDF + eps, type = "lpmatrix") -
                 predict(mod, PredDF, type = "lpmatrix")) / eps
        
        Rbeta <- 
          mvrnorm(n = 10000, 
                  coef(mod), 
                  vcov(mod)) %>% 
          t
        
        i <- 1
        
        want <- grep("Density", colnames(Xp))
        
        Simulations = Xp[, want] %*% Rbeta[want, ]
        
        # Calculating inflection from GAMMs ####
        
        GAMInflection[[b]][[length(GAMInflection[[b]]) + 1]] <<-
          Simulations[1:100,] %>%
          data.frame %>%
          bind_cols(PredDF %>% slice(1:100) %>%
                      mutate(Portion = ifelse(Density.Annual < 0, "First", "Last")), .) %>%
          # rename(Density.Annual = 1) %>%s
          pivot_longer(matches("X")) %>%
          group_by(Portion) %>% slice(sample(1:n(), 1000)) %>%
          mutate_at("value", as.mcmc) %>%
          summarise(Mean = mean(value),
                    Lower = quantile(value, prob = 0.025),
                    Upper = quantile(value, prob = 0.975)) %>%
          data.frame()
        
        # PredDF2 <- 
        #   data.frame(Density.Annual = c(runif(1000, min(a$Data$Density.Annual), 0),
        #                                 runif(1000, 0, max(a$Data$Density.Annual))),
        #              Portion = rep(c("First", "Last"), each = 1000))
        # 
        # Xp <- predict(mod, PredDF2, type = "lpmatrix")
        # 
        # Rbeta2 <- 
        #   mvrnorm(n = 1000, 
        #           coef(mod), 
        #           vcov(mod)) %>% 
        #   t
        # 
        # Simulations2 = Xp[, want] %*% Rbeta2[want,]
        # 
        # Simulations2 %>% 
        #   apply(1, mean) %>% t %>% 
        #   as.data.frame %>% 
        #   rename(DerivativeLower = 1, 
        #          DerivativeMean = 2, 
        #          DerivativeUpper = 3)
        
        # QuickGAMPredict(Data = a$Data %>% na.omit, 
        #                 Model = a$Base, 
        #                 Covariates = "Density.Annual", 
        #                 OutputCovariate = "Density.Annual") %>% 
        #   
        #   bind_cols(
        #     
        #     Simulations %>% 
        #       apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% 
        #       as.data.frame %>% 
        #       rename(DerivativeLower = 1, 
        #              DerivativeMean = 2, 
        #              DerivativeUpper = 3)
        #     
        #   )
        
        # Returning predicted ####
        
        ReturnDF <- 
          
          QuickGAMPredict(Data = a$Data %>% na.omit, 
                          Model = a$Base, 
                          Covariates = c("Density.Annual"),
                          OutputCovariate = "Density.Annual") %>% 
          
          bind_cols(
            
            Simulations %>% 
              apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% 
              as.data.frame %>% 
              rename(DerivativeLower = 1, 
                     DerivativeMean = 2, 
                     DerivativeUpper = 3)
            
          ) %>% 
          mutate(Slope = case_when(
            
            DerivativeLower > 0 & DerivativeUpper > 0 ~ "P",
            DerivativeLower < 0 & DerivativeUpper < 0 ~ "N",
            TRUE ~ "Z"
            
          ))
        
        ReturnDF$Group <- 
          ReturnDF$Slope %>% 
          paste(collapse = " ") %>% 
          scan(what = character(), text = .) %>%
          rle %>% 
          extract2("lengths") %>% 
          # return %>% 
          rep(1:length(.), .)
        
        ReturnDF %>% 
          slice(1:100) %>% 
          return
        
      }) %>% reduce(~full_join(.x, .y, by = "Density.Annual"))
    
    ReturnList %>% 
      rename_at(vars(matches(".x.x$")), ~str_replace(.x, ".x.x$", ".K5")) %>% 
      rename_at(vars(matches(".y.y$")), ~str_replace(.x, ".y.y$", ".K6")) %>% 
      rename_at(vars(matches(".x$")), ~str_replace(.x, ".x$", ".K3")) %>% 
      rename_at(vars(matches(".y$")), ~str_replace(.x, ".y$", ".K4")) %>%
      rename_at(46:56, ~paste0(.x, ".K7")) %>% 
      return
    
  })

names(PredictedData) <- names(GAMInflection) <- 
  paste0(BestFitDF$SystemReplicate, "_", BestFitDF$Response)

GAMInflectionDF <- 
  GAMInflection %>% map(~bind_rows(.x, .id = "K") %>% mutate_at("K", ~as.numeric(.x) + 2)) %>% 
  bind_rows(.id = "SystemReplicate")

GAMInflectionDF %>% filter(K == 4) %>% 
  rename(Estimate = Mean, LinearLower = Lower, LinearUpper = Upper) %>% 
  # pivot_wider(names_from = "Portion", values_from = c("Density.Annual", "Lower", "Upper")) %>% 
  saveRDS("Output/GAMInflectionDF.rds")

InterpretDF <-
  glue::glue(".K{3:7}") %>% 
  map(function(FocalK){
    
    PredictedData %>% 
      map(~.x %>% dplyr::select(matches(FocalK)) %>% rename_all(~str_remove(.x, FocalK))) %>% 
      map("Slope") %>% 
      map(~paste(.x, collapse = " ")) %>%  
      map(~.x %>% scan(what = character(), text = .) %>%
            rle %>% 
            extract2("values")) %>% 
      map(~paste(.x, collapse = " ")) %>% 
      map(function(a){
        
        print(a)
        
        # b <- str_remove_all(a, " |Z") %>% str_split("") %>% unlist %>% unique %>% 
        #   paste(collapse = " ")
        
        b <- a %>%
          str_remove_all(" |Z") %>% str_split("") %>% unlist %>%
          scan(what = character(), text = .) %>%
          rle %>% 
          extract2("values") %>% 
          paste(collapse = " ")
        
        if(length(b) == 0) b <- " "
        
        data.frame(
          
          Curves = a,
          Interpretation = case_when(a == "P" ~ "Positive", 
                                     a == "P Z" ~ "Saturating", 
                                     a == "Z P" ~ "Accelerating", 
                                     a == "Z" ~ "Flat", 
                                     a == "Z P Z" ~ "Logistic",
                                     a == "P Z N"|a == "Z P Z N Z"|b == "P N" ~ "Quadratic",
                                     b == "P" ~ "NLPositive",
                                     TRUE ~ "Nonlinear")
          
        ) %>% rename_all(~paste0(.x, FocalK))
        
      }) %>% bind_rows#(.id = "SystemReplicate")
    
  }) %>% bind_cols(SystemReplicate = paste0(BestFitDF$SystemReplicate, "_", BestFitDF$Response), .)

InterpretDF %<>% 
  mutate_at(vars(matches("Interpretation")), 
            ~factor(.x, levels = c("Accelerating", "NLPositive", "Positive", "Saturating", 
                                   "Logistic", "Quadratic", "Flat", "Nonlinear"))) %>% 
  separate(SystemReplicate, remove = F,
           sep = "_", into = c("System", "ContactType", "Response"))

PredictedData %<>% 
  bind_rows(.id = "SystemReplicate") %>% 
  left_join(InterpretDF, by = "SystemReplicate") %>% 
  mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>% 
  mutate_at("SystemReplicate", ~factor(.x, levels = unique(.x))) 

# PredictedData %<>% 
#   separate(SystemReplicate, remove = F,
#            sep = "\n", into = c("System", "ContactType", "Response"))

PredictedData %>% write.csv("Output/FullPredictedSmooths.csv", row.names = F)

InterpretDF %<>% # Cutting down to only a middling K
  dplyr::select(SystemReplicate, matches("K4$")) %>% 
  rename_all(~str_remove(.x, ".K4"))

InterpretDF %<>% 
  separate(SystemReplicate, sep = "_", into = c("System", "Behaviour", "Response"))

PredictedData$Behaviour <- PredictedData$SystemReplicate %>% str_split("\n") %>% map_chr(2)

InterpretDF <- 
  NodeData %>%
  group_by(SystemReplicate, Site) %>%
  mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) %>% 
  ungroup %>% 
  group_by(SystemReplicate, System, Behaviour) %>% 
  mutate(Strength2 = Strength + diff(range(Strength, na.rm = T))/15) %>% 
  mutate_at("Strength", ~.x + diff(range(.x, na.rm = T))/4.5) %>% 
  summarise(MinDensity.Annual = min(Density.Annual, na.rm = T),
            MaxDensity.Annual = max(Density.Annual - diff(range(Density.Annual, na.rm = T))/10, na.rm = T),
            Strength = max(Strength, na.rm = T),
            Strength2 = max(Strength2, na.rm = T),
            N = n()) %>% 
  rename(Fit = Strength) %>% 
  left_join(InterpretDF %>% filter(Response == "Strength"), ., by = c("System", "Behaviour")) %>%
  mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
  # mutate_at("SystemReplicate", ~factor(.x, levels = levels(PredictedData$SystemReplicate))) %>% 
  # left_join(SlopeCatDF) %>% 
  left_join(BestFitDF %>% filter(Response == "Strength") %>% dplyr::select(SystemReplicate, BestFit) %>%
              mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")))

InterpretDF %<>% 
  arrange(Interpretation) %>%
  mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
  # mutate_at("SystemReplicate", ~factor(.x, levels = levels(PredictedData$SystemReplicate)))
  mutate_at("SystemReplicate", ~factor(.x))

# InterpretDF %<>% 
#   mutate(System = str_split(SystemReplicate, "\n") %>% map_chr(1),
#          Behaviour = str_split(SystemReplicate, "\n") %>% map_chr(2))

InterpretDF %<>% mutate(Label = paste0(System, "\n", Interpretation, ", ", Behaviour))

InterpretDF %>% saveRDS("Output/NonLinearInterpretations.rds")

# Adding phylopics ####

if(0){
  
  source("R/00_Phylopic Setup.R")
  
  InterpretDF %<>% 
    # mutate(System = str_split(SystemReplicate, "\n") %>% map_chr(1)) %>% 
    left_join(SummaryTable[,c("RName", "ImageFile")], by = c("System" = "RName"))
  
  library(ggimage)
  
  PredictedData %>% 
    ggplot(aes(Density.Annual, Fit)) +
    geom_point(data = NodeData %>%
                 mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                 group_by(SystemReplicate, Site) %>%
                 mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) %>%
                 mutate_at("SystemReplicate", ~factor(.x, levels = levels(PredictedData$SystemReplicate))) ,
               aes(Density.Annual, Strength), inherit.aes = F, alpha = 0.1) +
    geom_line(size = 2, colour = "grey") + 
    geom_line(data = PredictedData %>% filter(Slope == "P"), 
              aes(group = Group),
              size = 2, colour = AlberColours[[2]]) +
    geom_line(data = PredictedData %>% filter(Slope == "N"), 
              aes(group = Group),
              size = 2, colour = AlberColours[[1]]) +
    # facet_wrap(~System + Behaviour, scales = "free", ncol = 5) +
    facet_wrap(~SystemReplicate, scales = "free", nrow = 10) +
    theme(axis.text = element_blank()) +
    scale_linetype_manual(values = c(2, 1), guide = "none") +
    # theme(legend.position = "top") +
    geom_text(data = InterpretDF %>% arrange(Interpretation) %>%
                mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                mutate_at("SystemReplicate", ~factor(.x, levels = levels(PredictedData$SystemReplicate))), 
              hjust = 0, size = 3,
              # colour = "white", label.colour = "black",
              aes(label = Interpretation, x = MinDensity.Annual)) +
    # geom_image(data = InterpretDF %>% arrange(Interpretation) %>%
    #              mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
    #              mutate_at("SystemReplicate", ~factor(.x, levels = levels(PredictedData$SystemReplicate))), 
    #            hjust = 1, vjust = 1,
    #            asp = 1.5,
    #            size = 0.1,
    #            # colour = "white", label.colour = "black",
    #            aes(image = ImageFile, x = MaxDensity.Annual)) + # %>% substr(1, 4))) +
    labs(x = "Density", y = Response)
  
  ggsave("Figures/PreliminaryFigure.jpeg", 
         width = 350*8.27/11.69, height = 350, 
         units = "mm", dpi = 600)
  
  # Multiple curves ####
  
  PredictedData %>% 
    # bind_rows(.id = "SystemReplicate") %>% 
    pivot_longer(contains("Fit")) %>% arrange(desc(name)) %>% 
    mutate_at("name", ~factor(.x, levels = unique(name))) %>% 
    ggplot(aes(Density.Annual, value)) +
    geom_point(data = NodeData %>%
                 # mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                 group_by(SystemReplicate, Site) %>%
                 mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) ,
               aes(Density.Annual, Strength), inherit.aes = F, alpha = 0.1) +
    geom_line(aes(colour = name),
              size = 2) +
    facet_wrap(~SystemReplicate, scales = "free", nrow = 10) +
    scale_colour_discrete_sequential(palette = AlberPalettes[[1]]) +
    theme(legend.position = "none") +
    theme(strip.text = element_blank())
  
  ggsave("Figures/MultipleFits.jpeg", 
         width = 350*8.27/11.69, height = 350, 
         units = "mm", dpi = 600)
  
  # Poster print ####
  
  PosterPlotData <- 
    PredictedData %>% 
    # bind_rows(.id = "SystemReplicate") %>% 
    pivot_longer(contains("Fit")) %>% arrange(desc(name)) %>% filter(name == "Fit.K4") %>%  
    # mutate_at("name", ~factor(.x, levels = unique(name))) %>% 
    arrange(Interpretation)
  
  SystemOrder <- PosterPlotData$SystemReplicate %>% unique
  
  PosterPlotData %<>% 
    mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder))
  
  InterpretDF %<>% 
    mutate(System = str_split(SystemReplicate, "_") %>% map_chr(1),
           Behaviour = str_split(SystemReplicate, "_") %>% map_chr(2))
  
  InterpretDF %<>% mutate(Label = paste0(System, "\n", Interpretation, ", ", Behaviour))
  
  InterpretDF %<>% 
    # mutate(System = str_split(SystemReplicate, "\n") %>% map_chr(1)) %>% 
    left_join(SummaryTable[,c("RName", "ImageFile")], by = c("System" = "RName")) %>% 
    left_join(NodeData %>%
                # mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                group_by(SystemReplicate, Site) %>%
                # mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
                mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) %>% 
                ungroup %>% group_by(SystemReplicate) %>% 
                summarise(Strength = max(Strength),
                          MaxDensity.Annual = max(Density.Annual))
    )
  
  PosterPlotData %>% 
    ggplot(aes(Density.Annual, value)) +
    geom_point(data = NodeData %>%
                 mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                 group_by(SystemReplicate, Site) %>%
                 mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
                 mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) ,
               # colour = AlberColours[[1]],
               aes(Density.Annual, Strength, alpha = SystemSizeAlpha), inherit.aes = F, alpha = 0.1) +
    geom_line(size = 2, colour = "grey") + 
    geom_line(data = PosterPlotData %>% 
                filter(Slope.K4 == "P"), 
              aes(group = Group.K4),
              size = 2, colour = AlberColours[[2]]) +
    geom_line(data = PosterPlotData %>% 
                filter(Slope.K4 == "N"), 
              aes(group = Group.K4),
              size = 2, colour = AlberColours[[1]]) +
    facet_wrap(~SystemReplicate, scales = "free", nrow = 10) +
    scale_colour_discrete_sequential(palette = AlberPalettes[[1]]) +
    theme(legend.position = "none") +
    theme(strip.text = element_blank()) +
    theme(axis.text = element_blank()) +
    labs(x = "Local population density", y = "Contact rate") +
    geom_image(data = InterpretDF %>% arrange(Interpretation) %>%
                 mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                 mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)),
               hjust = 1, #vjust = 1,
               asp = 1.5,
               size = 0.1,
               # colour = "white", label.colour = "black",
               aes(image = ImageFile, y = Strength, x = MaxDensity.Annual))  # %>% substr(1, 4)))
  
  ggsave("Figures/PosterFitsOrdered.jpeg", 
         width = 350*8.27/11.69, height = 350, 
         units = "mm", dpi = 1200)
  
  PosterPlotData %>% 
    ggplot(aes(Density.Annual, value)) +
    geom_point(data = NodeData %>%
                 mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                 group_by(SystemReplicate, Site) %>%
                 mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
                 mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) ,
               # colour = AlberColours[[1]],
               aes(Density.Annual, Strength, alpha = SystemSizeAlpha), inherit.aes = F, alpha = 0.1) +
    geom_line(size = 2, colour = "grey") + 
    geom_line(data = PosterPlotData %>% 
                filter(Slope.K4 == "P"), 
              aes(group = Group.K4),
              size = 2, colour = AlberColours[[2]]) +
    geom_line(data = PosterPlotData %>% 
                filter(Slope.K4 == "N"), 
              aes(group = Group.K4),
              size = 2, colour = AlberColours[[1]]) +
    facet_wrap(~SystemReplicate, scales = "free", nrow = 6) +
    scale_colour_discrete_sequential(palette = AlberPalettes[[1]]) +
    theme(legend.position = "none") +
    theme(strip.text = element_blank()) +
    theme(axis.text = element_blank()) +
    labs(x = "Local population density", y = "Contact rate") +
    geom_image(data = InterpretDF %>% arrange(Interpretation) %>%
                 mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
                 mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)),
               hjust = 1, #vjust = 1,
               asp = 1/1.25,
               size = 0.1,
               # colour = "white", label.colour = "black",
               aes(image = ImageFile, y = Strength, x = MaxDensity.Annual))  # %>% substr(1, 4)))
  
  ggsave("Figures/PosterFitsHorizontalasp2.jpeg", 
         height = 350*8.27/11.69, width = 350, 
         units = "mm", 
         
         dpi = 1200)
  
}
