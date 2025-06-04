
# 00_Phylopic Setup ####

library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce)
library("rnaturalearth"); library("rnaturalearthdata"); library(rotl); library(ape); library(ggtree)
library(fs); library(ggimage)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

World <- ne_countries(scale = "medium", returnclass = "sf")

World %<>% filter(!continent == "Antarctica")

SummaryTable <- 
  gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1fTcR8cefcihKxO277VNXci_SKDPXEbV5exg6HVg42ks/edit#gid=1894078012") %>% 
  filter(!is.na(RName))#, (Initial == "Y" | RName %in% c("MorelleBoars")))

SummaryTable %>% 
  dplyr::select(RName, Species) %>% 
  unique %>% data.frame %>% 
  write.csv("SystemList.csv", row.names = F)

SummaryTable %<>% 
  mutate_at("HostGroup", 
            ~case_when(.x %in% c("Ungulate", "Elephant") ~ "LargeHerbivore",
                       .x %in% c("Rodent") ~ "SmallMammal",
                       .x %in% c("Fish", "Cetacean", "Shark") ~ "Aquatic",
                       .x %in% c("Insect", "Reptile") ~ "Ectotherm",
                       TRUE ~ .x))

Images <- "Phylopics_Species" %>% dir_ls(regex = ".png") %>% sort 

# SummaryTable %<>% mutate(ImageFile = paste0("Phylopics_System/", RName, ".png"))

SummaryTable %<>% mutate(ImageFile = paste0("Phylopics_Species/", Species %>% str_replace(" ", "_"), ".png"))

AlberTree <- # Species level
  SummaryTable %>%
  pull(Species) %>% sort %>% 
  tnrs_match_names(do_approximate_matching = T ) %>% 
  data.frame %>% 
  pull(ott_id) %>% 
  tol_induced_subtree(ott_ids = .)

AlberTree$tip.label %<>% strip_ott_ids

AlberTree %<>% 
  compute.brlen(method = "Grafen") %>% 
  ladderize

AlberTree$tip.label %<>% str_replace("_", " ") %>% str_replace("maniculatus", "leucopus")

Delta <- "\U2022"

AlberTree$tip.label %>% #str_replace(" ", "_") %>% 
  map_chr(~SummaryTable %>% filter(Species == .x) %>% pull(HostGroup) %>% as.character %>% unique) ->
  Groups

groupInfo <- split(AlberTree$tip.label, factor(Groups, levels = unique(Groups)))

AlberTree <- groupOTU(AlberTree, groupInfo)

SummaryTable %>% 
  dplyr::select(RName, Species) %>% mutate_at("Species", ~factor(.x, levels = AlberTree$tip.label)) %>%  
  count(Species) %>%
  mutate(Deltas = map_chr(n, ~paste(rep(Delta, .x), collapse = " "))) ->
  SpeciesCounts

# SpeciesCounts <- SpeciesCounts[rep(1:nrow(SpeciesCounts), SpeciesCounts$n),]

SpeciesCounts %<>% mutate(Y = 1:n())

SpeciesCounts %<>% 
  left_join(SummaryTable %>% 
              group_by(Species) %>% summarise_at(c("System", "ImageFile", "HostGroup"), ~list(unique(.x))) %>% 
              dplyr::select(Species, System, ImageFile, HostGroup)) %>% 
  unnest(System, ImageFile, HostGroup) %>% 
  group_by(Species) %>% mutate(X = (1:n())/10 + 0.9) %>% 
  data.frame

# Species-level ####

if(1){
  
  SummaryTable %>% 
    dplyr::select(RName, Species) %>% 
    unique %>% unnest(Species) %>%
    count(Species) %>%
    mutate(n = ifelse(n == 1, 0, n)) %>%
    mutate(Deltas = map_chr(n, ~paste(rep(Delta, .x), collapse = " "))) ->
    
    SpeciesCounts
  
  # MetaDF <- read.csv("Output/MetaDF.csv")
  # 
  # MetaDF %<>% mutate_at("Behaviour", ~factor(.x, levels = rev(c("Direct", "Proximity", "Indirect"))))
  
  # SpeciesCounts %<>% 
  #   left_join(MetaDF %>% filter(Behaviour == "Indirect") %>% dplyr::select(Species, Behaviour)) %>% 
  #   mutate(Species2 = ifelse(!is.na(Behaviour), paste0("**", Species, "**"), Species)) %>% 
  #   mutate(Name = paste(Species2, Deltas) %>% str_trim)
  
  SpeciesCounts %<>% mutate(Name = paste(Species, Deltas) %>% str_trim)
  
  SpeciesCountVector <- SpeciesCounts$Name
  names(SpeciesCountVector) <- SpeciesCounts$Species
  
  SpeciesCountVector <- SpeciesCountVector[AlberTree$tip.label]
  
  PhyloImageVector <- paste0("Phylopics_Species/", 
                             SpeciesCountVector %>% 
                               names %>% 
                               str_replace(" ", "_"), ".png")
  
  AlberTree$tip.label <- SpeciesCountVector %>% str_replace_all("_", " ")
  
}
