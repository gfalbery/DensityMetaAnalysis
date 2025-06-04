
# 04_Combining Summary Data ####

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)
library(metafor)

SummaryTable <- read.csv("SavedSummaryTable.csv") # Produced in 01_

SummaryTable %>% 
  mutate_at("ContactTypes", ~str_replace(.x, " ", ", ")) %>% 
  dplyr::select(System, Data.Type, Species, ContactTypes) %>% 
  write.csv("Output/SuppTableSpecies.csv", row.names = F)

DerivedSystemTraits <- read.csv("DerivedSystemTraits.csv") # Produced in 01_

NodeData <- readRDS("Intermediate/NodeData.rds") # Produced in 01_

LinearEstimates <- readRDS("Output/LinearEstimates.rds") # Produced in 03a_

NonLinearInterpretations <- readRDS("Output/NonLinearInterpretations.rds") # Produced in 03b_

InflectionPoints <- readRDS("Output/InflectionPointDF.rds") # Produced in 03c_

SummaryDF <- 
  NodeData %>% 
  mutate_at("SystemReplicate", ~str_replace(.x, "_", "\n")) %>% 
  dplyr::select(System, SystemReplicate) %>% unique %>% 
  left_join(SummaryTable %>% rename(SystemLabel = System), by = c("System" = "RName")) %>% 
  left_join(NonLinearInterpretations, by = "SystemReplicate") %>% dplyr::select(-N) %>% 
  left_join(DerivedSystemTraits %>% 
              mutate_at("SystemReplicate", ~str_replace(.x, "_", "\n")), 
            by = c("SystemReplicate"))

SummaryDF %<>% mutate(Behaviour = str_split(SystemReplicate, "\n") %>% map_chr(last))

SummaryDF %<>% 
  mutate_at("Behaviour", 
            ~case_when(.x %in% c("Fight", "Mate", "Dominance") ~ "Direct",
                       .x == "HRO" ~ "Indirect", 
                       TRUE ~ "Proximity")) %>% 
  mutate_at("HostGroup", 
            ~case_when(.x %in% c("Ungulate", "Elephant") ~ "LargeHerbivore",
                       .x %in% c("Rodent") ~ "SmallMammal",
                       .x %in% c("Insect", "Reptile") ~ "Ectotherm",
                       TRUE ~ .x))

MetaDF <- 
  SummaryDF %>% 
  mutate_at("SystemReplicate", ~str_replace(.x, "\n", "_")) %>% 
  left_join(LinearEstimates, by = c("SystemReplicate", "Response")) %>% 
  dplyr::select(SystemReplicate, Species, Behaviour, HostGroup, 
                N:NYears, 
                Response,
                Estimate:LinearUpper) %>% 
  mutate(Sig = as.numeric((LinearLower*LinearUpper) > 0 ))

MetaDF %>% # MetaDF is the data frame for running the meta-analyses, not their outputs
  
  write.csv("Output/MetaDF.csv", row.names = F)

MetaDF %>% 
  rename(BehaviourClass = Behaviour) %>% dplyr::select(-Sig) %>% 
  separate(SystemReplicate, sep = "_", into = c("System", "ContactType")) %>% 
  arrange(System, ContactType) %>% 
  left_join(SummaryTable %>% dplyr::select(RName, SystemLabel = System), by = c("System" = "RName")) %>% 
  dplyr::select(System, Species, HostGroup, ContactType, BehaviourClass, 
                Response,
                N, NID, Area, NYears, Estimate, SE, t, R2, LinearLower, LinearUpper) %>% 
  write.csv("Output/SuppTableReplicates.csv", row.names = F)

# Making saturation meta-analysis equivalent ####

SaturationEstimates <- readRDS("Output/SaturationEstimates.rds") # Made in 03c_
# SaturationEstimates <- readRDS("Output/GAMInflectionDF.rds") %>% dplyr::select(-Portion) %>% mutate(R2 = Estimate) %>% mutate(Portion = rep(c("First", "Last"), 60)) # Produced in 03b_

SaturationMetaDF <- 
  SummaryDF %>% 
  mutate_at("SystemReplicate", ~str_replace(.x, "\n", "_")) %>% 
  left_join(SaturationEstimates, by = c("SystemReplicate", "Response")) %>% 
  dplyr::select(SystemReplicate, Species, Behaviour, HostGroup, 
                N:NYears, 
                Response,
                Estimate:Portion) %>% 
  mutate(Sig = as.numeric((LinearLower*LinearUpper) > 0 ))

SaturationMetaDF %>% # MetaDF is the data frame for running the meta-analyses, not their outputs
  write.csv("Output/SaturationMetaDF.csv", row.names = F)

