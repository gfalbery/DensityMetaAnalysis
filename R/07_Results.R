
# 07_Results ####

{
  
  library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
  library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce)
  library("rnaturalearth"); library("rnaturalearthdata"); library(rotl); library(ape); library(ggtree)
  
  theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))
  
  source("R/03b_Summarising GAMMs.R")
  
  source("R/04_Combining Summary Data.R")
  
  source("R/00_Phylopic Setup.R")
  
  # MetaDF ####
  
  MetaDF <- read.csv("Output/MetaDFOutput.csv")
  
  # "Output/MetaAnalyses" %>% dir_ls(regex = "Strength_Meta") %>% 
  #   map(read.csv) %>% bind_rows#(.id = "Behaviour")
  
  MetaDF %<>% mutate_at("Behaviour", ~factor(.x, levels = rev(c("Indirect", "Proximity"))))
  
  MetaDF <- SummaryTable %>% dplyr::select(RName, ImageFile) %>% left_join(MetaDF, ., by = c("System" = "RName"))
  
  InterpretDF %<>% 
    # mutate(System = str_split(SystemReplicate, "\n") %>% map_chr(1)) %>% 
    left_join(SummaryTable[,c("RName", "HostGroup", "ImageFile")], by = c("System" = "RName"))
  
  LinearPredDF <- 
    "Output/MetaAnalyses" %>% dir_ls(regex = "Strength_LinearMetaPredictions") %>% map(readRDS) %>% 
    map(~.x %>% slice(n())) %>% 
    bind_rows(.id = "Rep") %>% 
    mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_LinearMetaPredictions.rds")) %>% 
    filter(Rep != "Indirect_Associations") %>% 
    separate(Rep, sep = "_", into = c("Behaviour", "Response")) %>% 
    mutate_at("Behaviour", ~factor(.x, levels = levels(MetaDF$Behaviour)))
  
  MetaAnalysisModelList <- 
    "Output/MetaAnalyses" %>% dir_ls(regex = "StrengthLinearMeta") %>% map(readRDS)
  
  StrengthMetaAnalysis <- 
    readRDS(paste0("Output/StrengthMetaAnalysis/", "StrengthLinearMetaModel.rds"))
  
  SaturationPredDF <- 
    "Output/SaturationMetaAnalyses" %>% 
    dir_ls(regex = "Strength_LinearMetaPredictions") %>% map(readRDS) %>% 
    map(~.x %>% group_by(Portion) %>% slice(n())) %>% 
    bind_rows(.id = "Rep") %>% 
    mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_LinearMetaPredictions.rds")) %>% 
    filter(Rep != "Indirect_Associations") %>% 
    mutate_at("Rep", ~factor(.x, levels = c("Proximity_Associations", paste0(rep(c("Proximity", "Indirect")), 
                                                                             "_",
                                                                             rep(c("Degree", "Strength"), each = 2))))) %>% 
    separate(Rep, sep = "_", into = c("Behaviour", "Response"), remove = F) %>%
    mutate_at("Behaviour", ~factor(.x, levels = c("Proximity", "Indirect"))) %>% 
    filter(Response != "Degree")
  
  SaturationMetaDF <- read.csv("Output/SaturationMetaDFOutput.csv")
  
  SaturationModelList <- 
    "Output/SaturationMetaAnalyses" %>% dir_ls(regex = "StrengthLinearMeta") %>% map(readRDS)
  
}


# Getting the findings ####

# Abstract ####
# we found XXX% of density-contact relationships were significantly positive, 

LinearEstimates %>% group_by(Response) %>% 
  mutate(Sig = as.numeric(LinearLower*LinearUpper > 0 & Estimate > 0)) %>% 
  summarise_at("Sig", ~Prev(.x))

# XXX% nonlinear, 

BestFitDF %>% 
  filter(Response == "Strength") %>% 
  pull(NumBestFit) %>% Prev

# and XXXX% saturating

SaturationMetaDF %>% 
  group_by(SystemReplicate, Behaviour, Response) %>% 
  summarise(Saturation = as.numeric(Estimate[1] > Estimate[2])) %>% 
  group_by(Response) %>% 
  summarise_at("Saturation", Prev)


# of our 60 replicates, XX  (XX%) were significantly positive when analysed using linear models

LinearEstimates %>% group_by(Response) %>% 
  mutate(Sig = as.numeric(LinearLower*LinearUpper > 0 & Estimate > 0)) %>% 
  summarise_at("Sig", ~sum(.x))

LinearEstimates %>% group_by(Response) %>% 
  mutate(Sig = as.numeric(LinearLower*LinearUpper > 0 & Estimate > 0)) %>% 
  summarise_at("Sig", ~Prev(.x))

# Additionally, ). Meta-analyses identified a highly significant positive mean correlation between density and connectedness,
# both for social networks (XXX)

LinearPredDF %>% 
  filter(Behaviour == "Proximity") %>% 
  dplyr::select(Fit, Lower, Upper)

# and spatial networks (XXXX)

LinearPredDF %>% 
  filter(Behaviour == "Indirect") %>% 
  dplyr::select(Fit, Lower, Upper)

# Despite strong positive density-contact relationships overall, slopes were highly variable across systems
# for both direct (XXXX)

MetaAnalysisModelList$`Output/MetaAnalyses/Proximity_StrengthLinearMetaModelList.rds`$Phylogeny

# and indirect contacts (XXX)

MetaAnalysisModelList$`Output/MetaAnalyses/Indirect_StrengthLinearMetaModelList.rds`$Phylogeny

# Density’s effect was remarkably more positive for indirect than direct contacts 
# (Figure 3B; XXXX 

LinearPredDF %>% 
  filter(Behaviour == "Proximity") %>% 
  dplyr::select(Fit, Lower, Upper) 

# versus XXX)

LinearPredDF %>% 
  filter(Behaviour == "Indirect") %>% 
  dplyr::select(Fit, Lower, Upper)

# XXX for the difference when modelled together

StrengthMetaAnalysis

# X/60 (X%) of systems had a steeper slope at low density values than at high ones. 

SaturationMetaDF %>% 
  group_by(SystemReplicate, Behaviour, Response) %>% 
  summarise(Saturation = as.numeric(Estimate[1] > Estimate[2])) %>% 
  group_by(Response) %>% 
  filter(Saturation == 1) %>% count

SaturationMetaDF %>% 
  group_by(SystemReplicate, Behaviour, Response) %>% 
  summarise(Saturation = as.numeric(Estimate[1] > Estimate[2])) %>% 
  group_by(Response) %>% 
  summarise_at("Saturation", Prev)

# This effect was highly significant for direct contact (XXXX) 

SaturationModelList$`Output/SaturationMetaAnalyses/Proximity_StrengthLinearMetaModelList.rds`$`Portion + NID + Area + NYears`

# and less so for indirect contact (XXXX)

SaturationModelList$`Output/SaturationMetaAnalyses/Indirect_StrengthLinearMetaModelList.rds`$`Portion + NID + Area + NYears`

# the latter half of the indirect contact effect was still higher than the first half of the direct contact effect, 
# and the second half of the direct effect overlapped with zero (Figure 3C XXXX), 

# SaturationPredDF %>% filter(Behaviour == "Proximity", Portion == "Last")

# Beyond these general trends, our GAMs revealed that XX/60 relationships (XX%) were significantly nonlinear (DeltaAIC<-2)

BestFitDF %>% 
  filter(Response == "Strength") %>% 
  pull(NumBestFit) %>% table

BestFitDF %>% 
  filter(Response == "Strength") %>% 
  pull(NumBestFit) %>% Prev

# 

# MetaEstimates <- read.csv(paste0("Output/StrengthMetaAnalysis/", FocalRep, "_Table SX.csv"))

MetaEstimates
MetaEstimates$AICc %>% diff

# Old findings ####

BestFitDF$BestFit %>% table
(BestFitDF$BestFit=="Non-Linear") %>% Prev

41/60

InterpretDF$Saturating20 %>% table
InterpretDF$Saturating50 %>% table

39/60

pbinom(41, 60, 0.5, lower.tail = F)

pbinom(39, 60, 0.5, lower.tail = F)

MetaDF$Area %>% median

MetaDF %>% filter(as.logical(Sig), Estimate > 0) %>% nrow

MetaDF %>% filter(as.logical(Sig), Estimate > 0) %>% nrow %>% divide_by(nrow(MetaDF))

InterpretDF %>% 
  group_by(System) %>% 
  summarise(HROYes = as.numeric(any(Behaviour == "HRO"))) %>% 
  pull(2) %>% table

InterpretDF %>% 
  group_by(System) %>% 
  summarise(HROYes = as.numeric(any(Behaviour == "HRO"))) %>% 
  pull(2) %>% Prev

InterpretDF %>% 
  group_by(Behaviour, Interpretation) %>% 
  count %>% 
  data.frame

InterpretDF %>% 
  group_by(Behaviour) %>% 
  summarise_at("Interpretation", nunique)

# Methods ####

SummaryTable$Species %>% nunique

#HRO networks were restricted to only individuals with five or more observations in a given year to allow us to create convex polygons effectively; 
# XX/36 (XX%) systems did not have sufficient sampling for this analysis. 

NodeData %>% 
  
  filter(Behaviour == "HRO") %>% 
  pull(System) %>% nunique %>% subtract(36, .) %>% 
  divide_by(36)

NodeData %>% 
  
  filter(Behaviour != "HRO") %>% 
  pull(System) %>% nunique

NodeData %>% 
  dplyr::select(Behaviour, System) %>% unique %>% 
  count(System) %>% 
  filter(n>1) %>% 
  nrow



