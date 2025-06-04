
# 06_Capricorn Figures ####

{
  
  library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
  library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce)
  library("rnaturalearth"); library("rnaturalearthdata"); library(rotl); library(ape); library(ggtree)
  
  theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))
  
  source("R/03b_Summarising GAMMs.R")
  
  source("R/04_Combining Summary Data.R")
  
  source("R/00_Phylopic Setup.R")
  
  dir_create("Figures")
  
}

# Figure 1 is a schematic

# Figure 2: Infographic #### 
# (Put together in post)

# Panel 2A: Map ####

LongSecVar <- seq(from = -180, to = 180, length.out = round(nrow(SummaryTable)/2)) %>% rep(each = 2)
LatSecVar <- rep(c(100, -70), round(nrow(SummaryTable)/2))

Systems <- 
  SummaryTable %>% arrange(SiteLong) %>% 
  pull(System)

OrderedTops <- c("Palmyra sharks", "Golden-crowned sparrow",
                 "Prairie dogs", "Mountain Lake mice", 
                 "Bottlenose whales",
                 "Soay sheep", "Rum deer",
                 "Moray dolphins", "Firth of Tay dolphins",
                 "Kielder voles", "Wytham tits",
                 "Wytham wood mice", "European boar", 
                 "Ein Gedi hyraxes", "Vulturine guineafowl",
                 "Reticulated giraffes", "Kenyan elephants",
                 "Shark Bay dolphins") 

OrderedBottoms <- c("Moorea sharks", "Acorn woodpeckers",
                    "Desert tortoises", "Potomac dolphins", 
                    "Falkirk wood mice",
                    "Cornish cattle", 
                    "Wild crickets", "Liverpool wood mice", 
                    "Woodchester badgers", "Wytham badgers", 
                    "Silwood wood mice", "Kalahari meerkats",
                    "African dogs", 
                    "Multimammate mice", 
                    "Serengeti lions",
                    "Chagos sharks", "Sleepy lizards", 
                    "Water dragons")

rbind(OrderedTops, OrderedBottoms) %>% c -> OrderedAll 

SummaryTable <- 
  SummaryTable %>% 
  mutate(System = fct_relevel(System, OrderedAll)) %>% 
  arrange(System) %>% 
  mutate(LongSec = LongSecVar, LatSec = LatSecVar)

LongSummaryTable <- 
  SummaryTable %>% 
  pivot_longer(cols = c(SiteLong, LongSec), 
               names_to = "Begin",
               values_to = "SiteLong") %>% 
  mutate(Begin = "Begin") %>% dplyr::select(-SiteLat) %>% 
  bind_cols(SummaryTable %>% 
              pivot_longer(cols = c(SiteLat, LatSec), 
                           names_to = "End",
                           values_to = "SiteLat") %>% 
              mutate(Begin = "End") %>% dplyr::select(SiteLat))

(FullMap <-
    SummaryTable %>% 
    ggplot(aes(x = SiteLong, y = SiteLat)) +
    geom_sf(data = World, inherit.aes = F, fill = "grey", 
            colour = "grey", size = 2) +
    geom_path(data = LongSummaryTable, aes(group = System)) +
    geom_text(aes(x = LongSec, y = LatSec +20*(as.numeric(LatSec > 0 ) - 0.5)*2, 
                  label = System, 
                  hjust = ifelse(LatSec < 0, 1, 0)), 
              size = 6,
              angle = 90) + 
    geom_image(aes(x = LongSec, y = LatSec + 10*(as.numeric(LatSec > 0 ) - 0.5)*2,
                   image = ImageFile, colour = HostGroup), #by = "height",
               asp = 1) +
    coord_sf(datum = sf::st_crs(3857), ylim = c(-160, 200)) +
    labs(x = NULL, y = NULL) +
    scale_colour_brewer(palette = "Spectral") + 
    scale_fill_brewer(palette = "Spectral", 
                      labels = c("Aquatic", "Bird", "Carnivore", "Ectotherm", "Large herbivore", "Small mammal")) + 
    theme_void() + 
    theme(legend.position = "none") +
    NULL)

ggsave("Figures/PhyloMap.jpeg", units = "mm", height = 300, width = 300, dpi = 600)

# Panel 2B: Phylogeny ####

(TreePlot <- 
   AlberTree %>% 
   ggtree() +
   geom_image(data = . %>% filter(between(x, 0.975, 1.05)),
              aes(x = x + 0.11, colour = group, fill = group), by = "height", height = 0.03,
              image = PhyloImageVector) +
   geom_tiplab(aes(x = x + 0.23), 
               size = 6.5,
               colour = "black")  +
   scale_x_continuous(limits = c(0, 3.5))  + # Best for /?
   labs(colour = "Animal group", fill = "Animal group") +
   theme(legend.position = c(0, 1), legend.justification = c(0, 1)) +
   scale_colour_brewer(palette = "Spectral", 
                       labels = c("Aquatic", "Bird", "Carnivore", "Ectotherm", "Large herbivore", "Small mammal")))

ggsave("Figures/PhyloTree.jpeg", units = "mm", height = 275, width = 275, dpi = 600)

# Panel 1C: Density Schematic ###

source("R/06b_Density Schematic.R")

# Figure 3: Meta-analysis Outputs ####

# Panel 3A: Effect estimates ####

MetaDF <- read.csv("Output/MetaDFOutput.csv")

# "Output/MetaAnalyses" %>% dir_ls(regex = "Strength_Meta") %>% 
#   map(read.csv) %>% bind_rows#(.id = "Behaviour")

MetaDF %<>% 
  mutate_at("Behaviour", ~str_replace(.x, "Proximity", "Social") %>% str_replace("Indirect", "Spatial")) %>% 
  mutate_at("Behaviour", 
            ~factor(.x, levels = c("Spatial", "Social")))

MetaDF <- 
  SummaryTable %>% dplyr::select(RName, ImageFile) %>% 
  left_join(MetaDF, ., by = c("System" = "RName"))

(Panel3A <- 
    MetaDF %>% 
    filter(Response == "Strength") %>% 
    mutate(Sig = as.numeric((LinearLower*LinearUpper) > 0 )) %>%
    arrange(Estimate) %>% 
    mutate_at("SystemReplicate", ~factor(.x, levels = unique(.x))) %>% 
    ggplot(aes(SystemReplicate, Estimate)) +
    geom_hline(yintercept = 0, lty = 2, colour = "light grey") +
    geom_errorbar(aes(colour = Behaviour,
                      ymin = LinearLower, ymax = LinearUpper, alpha = as.factor(Sig)), 
                  width = 0.3)  +
    geom_point(colour = "white", size = 3.5, alpha = 0.8) +
    geom_point(aes(colour = Behaviour, shape = Behaviour), size = 2.5) +
    # geom_image(aes(image = ImageFile, y = LinearUpper + 0.1)) +
    coord_flip() + 
    theme(axis.line.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank()) +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1)) +
    labs(y = "Linear density effect estimate", shape = "Contact type") +
    scale_alpha_manual(values = c(0.3, 1), guide = "none") +
    scale_shape_manual(name = "Network type", limits = c("Spatial", "Social"), 
                       values = c(17, 16)) +
    scale_colour_manual(name = "Network type", limits = c("Spatial", "Social"),
                        values = c(AlberColours[[2]], AlberColours[[1]])) +
    NULL)

# Panel 3B: Spatial versus Social ####

# LinearPredDF <- readRDS("Output/LinearMetaPredictions.rds")

# "Output/MetaAnalyses" %>% dir_ls(regex = "Strength_LinearMeta") %>%
#   map(readRDS) %>% bind_rows(.id = "Behaviour")

LinearPredDF <- 
  "Output/MetaAnalyses" %>% dir_ls(regex = "Strength_LinearMetaPredictions") %>% map(readRDS) %>% 
  map(~.x %>% slice(n())) %>% 
  bind_rows(.id = "Rep") %>% 
  mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_LinearMetaPredictions.rds")) %>% 
  filter(Rep != "Indirect_Associations") %>% 
  separate(Rep, sep = "_", into = c("Behaviour", "Response")) %>% 
  mutate_at("Behaviour", ~str_replace(.x, "Proximity", "Social") %>% 
              str_replace("Indirect", "Spatial")) %>% 
  mutate_at("Behaviour", ~factor(.x, levels = levels(MetaDF$Behaviour)))

MetaDF <- read.csv("Output/MetaDFOutput.csv")

(Panel3B <- 
    LinearPredDF %>% 
    mutate_at("Behaviour", ~str_replace(.x, "Proximity", "Social") %>% str_replace("Indirect", "Spatial")) %>% 
    mutate_at("Behaviour", ~factor(.x, levels = c("Spatial", "Social"))) %>% 
    # filter(NYears == last(NYears)) %>% 
    ggplot() +
    geom_hline(yintercept = 0,linetype = 2) +
    geom_jitter(data = MetaDF %>% filter(Response == "Strength") %>% 
                  mutate_at("Behaviour", ~str_replace(.x, "Proximity", "Social") %>% str_replace("Indirect", "Spatial")) %>% 
                  mutate_at("Behaviour", ~factor(.x, levels = c("Spatial", "Social"))), 
                aes(x = Behaviour, colour = Behaviour, y = R, size = 1/vi), width = 0.1, alpha = 0.25) +
    geom_segment(aes(x = Behaviour, xend = Behaviour,
                     y = Lower, yend = Upper)) +
    geom_point(aes(x = Behaviour, y = Fit, group = Behaviour, colour = Behaviour), 
               size = 5, shape = 16,
               colour = "black") +
    guides(size = "none") +
    labs(x = "Network type") +
    labs(y = "Density effect r value") +
    scale_size_continuous(range = c(1,3)) +
    scale_colour_manual(values = c(AlberColours[[2]], AlberColours[[1]])) +
    theme(legend.position = "none"))

# Panel 3C: Saturation ####

# SaturationPredDF <- readRDS("Output/SaturationMetaPredictions.rds")

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
  mutate_at("Behaviour", ~str_replace(.x, "Proximity", "Social") %>% str_replace("Indirect", "Spatial")) %>% 
  mutate_at("Behaviour", ~factor(.x, levels = c("Spatial", "Social"))) %>% 
  filter(Response != "Degree")

SaturationMetaDF <- read.csv("Output/SaturationMetaDFOutput.csv")

(Panel3C <- 
    SaturationPredDF %>% 
    mutate(X = case_when(Behaviour == "Spatial" & Portion == "First" ~ 0.5,
                         Behaviour == "Spatial" & Portion == "Last" ~ 1,
                         Behaviour == "Social" & Portion == "First" ~ 2,
                         Behaviour == "Social" & Portion == "Last" ~ 2.5)) %>% 
    # mutate_at("Behaviour", ~str_replace(.x, "Proximity", "Direct")) %>% 
    filter(NYears == last(NYears)) %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_jitter(#position = position_dodge(w = 0.3), 
      data = SaturationMetaDF %>% 
        mutate(X = case_when(Behaviour == "Indirect" & Portion == "First" ~ 0.5,
                             Behaviour == "Indirect" & Portion == "Last" ~ 1,
                             Behaviour == "Proximity" & Portion == "First" ~ 2,
                             Behaviour == "Proximity" & Portion == "Last" ~ 2.5)) %>% 
        mutate_at("Behaviour", ~str_replace(.x, "Proximity", "Direct")), 
      aes(x = X, colour = Portion, y = R, size = 1/vi), 
      width = 0.1,
      alpha = 0.25) +
    geom_segment(#position = position_dodge(w = 0.3), 
      aes(x = X, xend = X, group = Portion,
          y = Lower, yend = Upper)) +
    geom_point(#position = position_dodge(w = 0.3), 
      aes(x = X, y = Fit, group = Portion, colour = Portion), size = 5, shape = 16,
      colour = "black") +
    guides(size = "none") +
    scale_size_continuous(range = c(1,3)) +
    labs(y = "Density effect r value") +
    scale_x_continuous(
      limits = c(0, 3),
      breaks = c(0.75, 2.25),
      labels = c("Spatial", "Social"), name = "Network type") +
    scale_colour_manual(values = c(AlberColours[[3]], AlberColours[["Pink"]])))

(Panel3A/(Panel3B + Panel3C)) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure3.jpeg", units = "mm", height = 220, width = 150)

# Figure 4: GAMM Smooths ####

Response <- "Strength"

PlotNodeData <- 
  NodeData %>%
  left_join(SummaryTable %>% rename(SystemLabel = System, System = RName), 
            by = c("System")) %>% mutate_at(c("System", "SystemLabel"), as.character) %>% 
  #   mutate_at("SystemReplicate", ~str_replace_all(.x, c(System = SystemLabel)))
  mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
  # mutate(SystemReplicate = paste0(SystemLabel, "\n", Behaviour)) %>% 
  mutate(SystemReplicate = paste0(System, "\n", Behaviour, "\n", Response)) %>%
  group_by(SystemReplicate, Site) %>%
  mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) %>%
  mutate_at("SystemReplicate", ~factor(.x, levels = levels(PredictedData$SystemReplicate)))

# PosterPlotData %<>%
#   mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder))

PredictedData %<>% # Cutting down to only a middling K
  dplyr::select(SystemReplicate, matches("K4$"), Density.Annual) %>%
  rename_all(~str_remove(.x, "[.]K4"))

PosterPlotData <-
  PredictedData %>% filter(str_detect(SystemReplicate, "Strength")) %>% 
  left_join(SummaryDF[,c("SystemReplicate", "Species")] %>% 
              mutate_at("SystemReplicate", ~paste0(.x, "\nStrength")) %>% 
              unique, by = c("SystemReplicate")) %>% 
  # bind_rows(.id = "SystemReplicate") %>%
  # pivot_longer(contains("Fit")) %>% arrange(desc(name)) %>% filter(name == "Fit.K4") %>%
  # mutate_at("name", ~factor(.x, levels = unique(name))) %>%
  mutate_at("Species", ~factor(.x, levels = str_replace(AlberTree$tip.label %>% str_remove_all(" •"), "_", " "))) %>% 
  arrange(Species, SystemReplicate)

SystemOrder <- PosterPlotData$SystemReplicate %>% unique

# if(file.exists("SavedSummaryTable.csv")) SummaryTable <- read.csv("SavedSummaryTable.csv")
# 
# SummaryTable <- 
#   gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1fTcR8cefcihKxO277VNXci_SKDPXEbV5exg6HVg42ks/edit#gid=1894078012") %>% 
#   filter(!is.na(RName), Initial == "Y")
# 
# SummaryTable %>% write.csv("SavedSummaryTable.csv", row.names = F)

# source("00_Phylopic Setup.R")

ImageDF <-
  InterpretDF %>% 
  # mutate(System = str_split(SystemReplicate, "\n") %>% map_chr(1)) %>% 
  left_join(SummaryTable[,c("RName", "HostGroup", "ImageFile")], by = c("System" = "RName"))# %>% 
# left_join(NodeData %>%
#             # mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
#             group_by(SystemReplicate, Site) %>%
#             # mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
#             mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) %>% 
#             ungroup %>% group_by(SystemReplicate) %>% 
#             summarise(Associations = max(Associations),
#                       MaxDensity.Annual = max(Density.Annual))
# )

ContactTranslate <- c("Trapping" = "-trapping",
                      "Sharing" = " sharing",
                      "Proximity" = " proximity",
                      "Collar" = " collar",
                      "Fight" = "Fighting",
                      "Mate" = "Mating")

PredictedData %<>%
  mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
  filter(!is.na(SystemReplicate))

PredictedData %>% 
  mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
  ggplot(aes(Density.Annual, Fit)) +
  geom_point(data = PlotNodeData,
             aes(Density.Annual, Strength), inherit.aes = F, alpha = 0.1) +
  geom_line(size = 2, colour = "grey") + 
  geom_line(data = PredictedData %>% filter(Slope == "P") %>% 
              mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)), 
            aes(group = Group),
            size = 2, colour = AlberColours[[2]]) +
  geom_line(data = PredictedData %>% filter(Slope == "N") %>% 
              mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)), 
            aes(group = Group),
            size = 2, colour = AlberColours[[1]]) +
  facet_wrap(~SystemReplicate, scales = "free", nrow = 8) +
  theme(axis.text = element_blank()) +
  scale_linetype_manual(values = c(2, 1), guide = "none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  geom_text(data = InterpretDF %>% dplyr::select(-Label) %>%
              left_join(SummaryTable[,c("RName", "System")] %>% rename(Label = System),
                        by = c("System" = "RName")) %>%
              mutate_at("SystemReplicate", ~paste0(.x, "\nStrength")) %>% 
              mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
              filter(!is.na(SystemReplicate)),
            hjust = 0, size = 3,
            alpha = 0.6,
            # vjust = 1,
            aes(label = Label, y = Fit, x = MinDensity.Annual)) +
  geom_text(data = InterpretDF %>%
              mutate_at("SystemReplicate", ~paste0(.x, "\nStrength")) %>% 
              mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
              filter(!is.na(SystemReplicate)),
            hjust = 0, size = 3,
            alpha = 0.6,
            aes(label = str_replace_all(Behaviour, ContactTranslate),
                y = Strength2,
                x = MinDensity.Annual)) + # %>% substr(1, 4))) +
  labs(x = "Density", y = Response) +
  geom_image(data = ImageDF %>%
               # mutate_at("SystemReplicate", ~str_replace_all(.x, "_", "\n")) %>%
               mutate_at("SystemReplicate", ~paste0(.x, "\nStrength")) %>% 
               mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
               filter(!is.na(SystemReplicate)),
             hjust = 0.25, #vjust = 1,
             asp = 1.25,
             size = 0.15,
             # colour = "white", label.colour = "black",
             aes(image = ImageFile, 
                 colour = HostGroup,
                 y = Fit*0.9, x = MaxDensity.Annual)) +
  scale_colour_brewer(palette = "Spectral", 
                      labels = c("Aquatic", "Bird", "Carnivore", "Ectotherm", "Large herbivore", "Small mammal")) +
  theme(legend.position = "top") + labs(colour = "Animal group") +
  guides(colour = guide_legend(# reverse = F,
    # direction = "horizontal",
    # title.position = "left",
    # title.vjust = 0.25, title.hjust = 0.5,
    # label.position = "top",
    # title.theme = element_text(size = 15),
    # label.hjust = 0.5,
    # label.vjust = 1.5,
    # label.theme = element_text(angle = 0), 
    nrow = 1)) +
  NULL

ggsave("Figures/Figure4_Strength.jpeg", 
       width = 350, height = 350, units = "mm", dpi = 600)

# ggsave("Figures/Figure4_Strength.jpeg", 
#        width = 350*8.27/11.69, height = 350, units = "mm", dpi = 600)

# ggsave("Figures/Figure4Poster.jpeg", 
#        width = 350*8.27/11.69, height = 350, units = "mm", dpi = 600)

# ~~~~ Supplementary Figures ####

# Supplementary Table: Effect sizes from meta-analysis ####

ModelList <- readRDS("Output/MetaAnalyses/FullLinearMetaModelList.rds")

MetaEstimates <- read.csv(paste0("Output/StrengthMetaAnalysis/", FocalRep, "_Table SX.csv"))

(ModelList$`Behaviour + NID + Area + NYears` %>% 
  summary)[1:7][-c(1, 3, 4)] %>% 
  data.frame %>% 
  rownames_to_column("Variable") %>% 
  rename(Estimate = beta, P = pval, CI_Lower = `ci.lb`, CI_Upper = `ci.ub`) %>% 
  mutate_at(c(2, 4, 5), ~round(.x, 3)) %>% 
  write.csv("Output/SuppTable3.csv", row.names = F)
  
# Supplementary Figure 1: Radius density ####

RadiusData <- readRDS("Intermediate/RadiusData.rds")

RadiusData %>% 
  group_by(System, Site) %>% 
  mutate_at("Density.Annual", ~c(scale(.x))) %>% 
  ggplot(aes(Density.Annual, RadCount)) + 
  geom_point(aes(alpha = SystemSizeAlpha), colour = AlberColours[["Blue"]]) + 
  geom_smooth() +
  facet_wrap(~System, scales = "free") +
  ggpubr::stat_cor(aes(y = RadCount*1.2), 
                   method = "spearman") +
  # ggpubr::stat_cor() +
  theme(legend.position = "none") +
  NULL

ggsave("Figures/RadiusComparison.jpeg", 
       width = 350, height = 350, units = "mm", dpi = 600)

  
# Supplementary Figure 3: Heteroskedasticity stuff ####


library(lmtest)

BPValues <- 
  LMList %>% 
  map("Strength") %>% 
  map(bptest) %>% 
  map(1) %>% as.data.frame %>% 
  pivot_longer(1:ncol(.)) %>% 
  separate(name, sep = "_", into = c("System", "ContactType")) %>% 
  mutate(SpatialBehaviour = ifelse(ContactType == "HRO", "Spatial", "Social")) %>% 
  rename(BP = value)

BPValues %>% saveRDS("Intermediate/BPValues.rds")

BPValues %>% 
  SinaGraph("SpatialBehaviour", "BP")

bptest(LMList[[1]]$Strength)[[1]]

{
  
  library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
  library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce)
  library("rnaturalearth"); library("rnaturalearthdata"); library(rotl); library(ape); library(ggtree)
  
  theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))
  
  source("R/04_Combining Summary Data.R")
  
  source("R/00_Phylopic Setup.R")
  
  # MetaDF ####
  
  MetaDF <- read.csv("Output/MetaDFOutput.csv")
  
  # "Output/MetaAnalyses" %>% dir_ls(regex = "Strength_Meta") %>% 
  #   map(read.csv) %>% bind_rows#(.id = "Behaviour")
  
  MetaDF %<>% mutate_at("Behaviour", ~factor(.x, levels = rev(c("Indirect", "Proximity"))))
  
  MetaDF <- SummaryTable %>% dplyr::select(RName, ImageFile) %>% left_join(MetaDF, ., by = c("System" = "RName"))
  
}

MetaDF %<>% left_join(BPValues, by = c("System", "ContactType"))

MetaDF %>% 
  ggplot(aes(BP, Estimate)) + geom_point() +
  geom_smooth(method = lm)

MetaDF %>% 
  ggplot(aes(Estimate, BP)) + geom_point() +
  geom_smooth(method = lm, aes(colour = SpatialBehaviour))


BPValues %>% 
  SinaGraph("SpatialBehaviour", "BP") + scale_y_log10() +
  labs(x = "Contact type") +
  labs(colour = "Contact\ntype") +
  
  (
    
    MetaDF %>% 
      ggplot(aes(Estimate, BP)) + geom_point() + theme(legend.position = "none") +
      geom_smooth(method = lm, aes(colour = SpatialBehaviour)) +
      scale_y_log10() +
      labs(colour = "Contact\ntype")
    
  ) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Supplementary_BPTests.jpeg", 
       width = 250, height = 150, units = "mm", dpi = 600)


# Supplementary Figure X: scaled linear ####

NodeData %>%
  group_by(System, Behaviour, Site) %>%
  mutate_at(c("Strength", "Density.Annual"), ~c(scale(.x))) %>%
  ggplot(aes(Density.Annual, Strength)) +
  geom_point(aes(alpha = SystemSizeAlpha),
             colour = AlberColours[[1]]) +
  geom_smooth(colour = "black", fill = NA, method = lm) +
  facet_wrap(~System + Behaviour, scales = "free", ncol = 8) +
  # theme_void() +
  theme(axis.text = element_blank()) +
  # theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_x_continuous(breaks = c(-10:10)) +
  scale_y_continuous(breaks = c(-10:10)) +
  theme(legend.position = "none") +
  labs(x = expression(Density/km^{2})) +
  ggpubr::stat_cor(aes(y = Strength*1.5))

ggsave("Figures/ScaledLinear.jpeg",
       units = "mm", height = 350, width = 350, dpi = 600)

# Supplementary Figure X: Log-log linear ####

NodeData %>%
  group_by(System, Behaviour, Site) %>%
  mutate_at(c("Strength", "Density.Annual"), ~c(log10(.x+1))) %>%
  ggplot(aes(Density.Annual, Strength)) +
  geom_point(aes(alpha = SystemSizeAlpha),
             colour = AlberColours[[1]]) +
  geom_smooth(colour = "black", fill = NA, method = lm) +
  facet_wrap(~System + Behaviour, scales = "free", ncol = 8) +
  # theme_void() +
  theme(axis.text = element_blank()) +
  # theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_x_continuous(breaks = c(-10:10)) +
  scale_y_continuous(breaks = c(-10:10)) +
  theme(legend.position = "none") +
  labs(y = "log(Strength+1)", 
       x = expression(log(Density/km^{2}+1))) +
  ggpubr::stat_cor(aes(y = Strength*1.25))

ggsave("Figures/LogLogLinear.jpeg",
       units = "mm", height = 350, width = 350, dpi = 600)



# Smooths without points ####

source("03b_Summarising GAMMs.R")

PredictedData %>% 
  ggplot(aes(Density.Annual, Fit)) +
  geom_line(size = 2, colour = "grey") + 
  geom_line(data = PredictedData %>% filter(Slope == "P"), 
            aes(group = Group),
            size = 2, colour = AlberColours[[2]]) +
  geom_line(data = PredictedData %>% filter(Slope == "N"), 
            aes(group = Group),
            size = 2, colour = AlberColours[[1]]) +
  facet_wrap(~SystemReplicate, scales = "free", nrow = 10) +
  theme(axis.text = element_blank()) +
  scale_linetype_manual(values = c(2, 1), guide = "none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  geom_text(data = InterpretDF, 
            hjust = 0, size = 3, vjust = 1,
            aes(label = Label, y = Strength3, x = MinDensity.Annual)) +
  labs(x = "Density", y = Response)

ggsave("Figures/SupplementaryFigure2.jpeg", 
       width = 350*8.27/11.69, height = 350, units = "mm")

# Table of sample sizes ####



# Figure of meta-analytical estimates ####

MetaDF
Proximity_StrengthLinearMetaModelList <- readRDS("Output/MetaAnalyses/Proximity_StrengthLinearMetaModelList.rds")
Indirect_StrengthLinearMetaModelList <- readRDS("Output/MetaAnalyses/Indirect_StrengthLinearMetaModelList.rds")

Proximity_StrengthLinearMetaModelList

# Smooth figures for presentations etc ####

PredictedData %>% 
  # group_by(SystemReplicate) %>% mutate_at(c("Density.Annual", "Fit"), ~scales::rescale(.x, to = c(0, 1))) %>% 
  mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
  ggplot(aes(Density.Annual, Fit, group = SystemReplicate)) +
  # geom_point(data = PlotNodeData,
  #            aes(Density.Annual, Associations), inherit.aes = F, alpha = 0.1) +
  geom_line(size = 2, colour = "grey") + 
  geom_line(data = PredictedData %>% 
              # group_by(SystemReplicate) %>% mutate_at(c("Density.Annual", "Fit"), ~scales::rescale(.x, to = c(0, 1))) %>% 
              filter(Slope == "P") %>% 
              mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)), 
            alpha = 0.5,
            aes(group = paste0(Group, SystemReplicate)),
            size = 2, colour = AlberColours[[2]]) +
  geom_line(data = PredictedData %>% 
              # group_by(SystemReplicate) %>% mutate_at(c("Density.Annual", "Fit"), ~scales::rescale(.x, to = c(0, 1))) %>% 
              filter(Slope == "N") %>% 
              mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)), 
            aes(group = paste0(Group, SystemReplicate)),
            alpha = 0.5,
            size = 2, colour = AlberColours[[1]]) +
  theme(axis.text = element_blank()) +
  scale_linetype_manual(values = c(2, 1), guide = "none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  labs(x = "Density", y = Response)  +
  NULL

# Log-log linear fits ####

NodeData %>% 
  # group_by(SystemReplicate) %>% mutate_at(c("Density.Annual", "Fit"), ~scales::rescale(.x, to = c(0, 1))) %>% 
  mutate_at("SystemReplicate", ~factor(.x, levels = SystemOrder)) %>% 
  ggplot(aes(Density.Annual, Strength, group = SystemReplicate)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = lm) 
  theme(axis.text = element_blank()) +
  scale_linetype_manual(values = c(2, 1), guide = "none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  labs(x = "Density", y = Response)  +
  NULL


# NodeData <- readRDS("Intermediate/NodeData.rds")

NodeData %>%
  group_by(System, Behaviour, Site) %>%
  mutate_at(c("Strength", "Density.Annual"), ~c(log10(.x+1))) %>%
  ggplot(aes(Density.Annual, Strength)) +
  geom_point(aes(alpha = SystemSizeAlpha),
             colour = AlberColours[[1]]) +
  geom_smooth(colour = "black", fill = NA, method = lm) +
  facet_wrap(~System + Behaviour, scales = "free", ncol = 8) +
  theme_void() +
  theme(axis.text = element_blank()) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  scale_x_continuous(breaks = c(-10:10)) +
  scale_y_continuous(breaks = c(-10:10)) +
  theme(legend.position = "none") +
  labs(x = expression(Density/km^{2}))

ggsave("Figures/LogLogLinear.jpeg",
       units = "mm", height = 350, width = 350)

# NodeData %>% 
#   group_by(System, Behaviour, Site) %>%
#   mutate_at(c("Associations", "Density.Annual"), ~c(scale(.x))) %>%
#   ggplot(aes(Density.Annual, Associations)) + 
#   # geom_point(aes(alpha = SystemSizeAlpha), 
#   #            colour = AlberColours[[1]]) +
#   geom_smooth(colour = "black", fill = NA) +
#   facet_wrap(~System + Behaviour, scales = "free", ncol = 8) +
#   theme_void() +
#   theme(axis.text = element_blank()) +
#   theme(strip.background = element_blank(), strip.text = element_blank()) +
#   scale_x_continuous(breaks = c(-10:10)) +
#   scale_y_continuous(breaks = c(-10:10)) +
#   theme(legend.position = "none") +
#   labs(x = expression(Density/km^{2}))
# 
# ggsave("Figures/ExampleFigure2.jpeg", 
#        units = "mm", height = 200, width = 350)

# Figure 3 ####

# InterpretDF %>%
SummaryDF %>%
  # left_join(MetaDF) %>%
  ggMMplot("Behaviour", "Interpretation", Just = T) +
  labs(fill = "Shape") +
  lims(y = c(0, 1.2), x = c(0, 1.05)) +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(NA, 1.25)) +
  coord_fixed() +
  
  SummaryDF %>% ggMMplot("HostGroup", "Interpretation", Just = T) +
  labs(fill = "Shape") + 
  scale_fill_brewer(palette = "Spectral") +
  lims(y = c(0, 1.2), x = c(0, 1.05)) + coord_fixed() +
  
  plot_layout(guides = "collect")

ggsave("Figures/SupplementaryFigure3.jpeg", units = "mm", height = 160, width = 320)


SummaryDF %>% ggMMplot("Indirect", "Interpretation"
) +
  labs(fill = "Shape") + 
  coord_fixed()


# Showing multiple K smooths ####

FullPredictedData <- read.csv("Output/FullPredictedSmooths.csv")

PredictedData %>% 
  bind_rows(.id = "SystemReplicate") %>% 
  pivot_longer(contains("Fit")) %>% arrange(desc(name)) %>% 
  mutate_at("name", ~factor(.x, levels = unique(name))) %>% 
  ggplot(aes(Density.Annual, value)) +
  geom_point(data = NodeData %>%
               group_by(SystemReplicate, Site) %>%
               mutate_at(c("Density.Annual", Responses), ~c(scale(.x))) ,
             aes(Density.Annual, Associations), inherit.aes = F, alpha = 0.1) +
  geom_line(aes(colour = name),
            size = 2) +
  facet_wrap(~SystemReplicate, scales = "free", nrow = 10) +
  scale_colour_discrete_sequential(palette = AlberPalettes[[1]]) +
  theme(legend.position = "none") +
  theme(strip.text = element_blank())

ggsave("Figures/MultipleFits.jpeg", 
       width = 350*8.27/11.69, height = 350, 
       units = "mm", dpi = 600)

# Misc ####

MetaDF %>% 
  mutate(System = str_split(SystemReplicate, "_") %>% map_chr(1)) %>% 
  ggplot(aes(Behaviour, colour = Behaviour, Estimate)) +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]], AlberColours[[3]])) +
  geom_line(aes(group = System)) +
  geom_point(size = 3, colour = "white") +
  geom_point(size = 2) +
  NULL  

# Year Sampling Data ####

YearDF <- 
  NodeData %>% 
  group_by(System) %>% 
  summarise(MinYear = min(Year),
            MaxYear = max(Year)) %>% 
  mutate_at(2:3, as.numeric) %>% 
  na.omit %>% data.frame %>% 
  left_join(SummaryTable %>% 
              dplyr::select(SystemLabel = System, 
                            HostGroup,
                            System = RName, Species)) %>%
  mutate_at("Species", ~str_replace(.x, " ", "_"))

YearDF %<>% 
  arrange(MinYear) %>% 
  mutate_at(c("System", "SystemLabel"), ~fct_relevel(.x, rev(unique(.x)))) %>% 
  mutate_at("HostGroup", 
            ~case_when(.x %in% c("Ungulate", "Elephant") ~ "LargeHerbivore",
                       .x %in% c("Rodent") ~ "SmallMammal",
                       .x %in% c("Fish", "Cetacean", "Shark") ~ "Aquatic",
                       .x %in% c("Insect", "Reptile") ~ "Ectotherm",
                       TRUE ~ .x))

YearPhylos <-
  YearDF %>% 
  mutate(ImageFile = paste0("Phylopics_Species/", Species, ".png")) %>% 
  mutate_at("MinYear", ~.x - 3) %>% 
  filter(MinYear > 0)

DecadeDF <- 
  data.frame(Y = 1:6*10 + 1960)

YearDF %>% 
  filter(MinYear > 1) %>% 
  pivot_longer(2:3) %>% 
  ggplot(aes(SystemLabel, value)) +
  geom_vline(aes(xintercept = SystemLabel),
             alpha = 0.1) +
  geom_hline(yintercept = DecadeDF$Y,
             alpha = 0.1) +
  geom_image(data = YearPhylos,
             aes(image = ImageFile,
                 colour = HostGroup,
                 y = MinYear)) +
  geom_line(aes(group = System,
                colour = HostGroup), 
            # colour = "dark grey",
            linewidth = 2) +
  geom_point(size = 5,
             aes(group = System, 
                 colour = HostGroup)) +
  labs(y = "Year", x = NULL, colour = "Host group") + 
  theme(legend.position = c(0.1, 0.1),
        legend.background = element_rect(fill = "white"),
        legend.justification = c(0, 0)) +
  coord_flip() +
  scale_colour_brewer(palette = "Spectral", 
                      labels = c("Aquatic", "Bird", "Carnivore", "Ectotherm", "Large herbivore", "Small mammal")) + 
  scale_fill_brewer(palette = "Spectral", 
                    labels = c("Aquatic", "Bird", "Carnivore", "Ectotherm", "Large herbivore", "Small mammal")) 

ggsave("Figures/TimeFigure.jpeg", units = "mm", 
       width = 200, height = 200)

ggsave("Figures/TimeFigure2.jpeg", units = "mm", 
       width = 250, height = 250)

# Looking at distributions ####

MetaDF <- read.csv("MetaDF.csv")

MetaDF %<>% mutate_at("Behaviour", ~factor(.x, levels = rev(c("Direct", "Proximity", "Indirect"))))

NodeData %<>% dplyr::select(-Behaviour) %>% left_join(MetaDF[,c("SystemReplicate", "Behaviour")])

NodeData %>% ggplot(aes(Density.Annual, Associations)) + 
  geom_density_2d() +
  facet_wrap(~Behaviour, scales = "free")

NodeData %>% ggplot(aes(Density.Annual, Associations)) + 
  geom_density_2d(aes(colour = SystemReplicate)) +
  facet_wrap(~Behaviour, scales = "free")

NodeData %>% 
  group_by(SystemReplicate, Site) %>% 
  mutate_at(c("Density.Annual", "Associations"), ~c(scale(.x))) %>% 
  mutate_at(c("Density.Annual", "Associations"), ~c(scale(rank(.x)))) %>% 
  # ggplot(aes(Density.Annual, colour = SpatialBehaviour, group = SystemReplicate)) + 
  ggplot(aes(Density.Annual, colour = Behaviour)) + 
  geom_density()

NodeData %>% 
  group_by(SystemReplicate, Site) %>% 
  # mutate_at(c("Density.Annual", "Associations"), ~c(scale(.x))) %>% 
  mutate_at(c("Density.Annual", "Associations"), ~c(scale(rank(.x)))) %>% 
  # ggplot(aes(Associations, colour = Behaviour)) + 
  ggplot(aes(Associations, colour = Behaviour, group = SystemReplicate)) + 
  geom_density()

(FocalSystems <- "Output/GAMM" %>% list.files %>% str_remove(".rds$"))

FocalSystems %>% 
  map(function(FocalSystem){
    
    NodeData %>% 
      filter(SystemReplicate == FocalSystem) %>% 
      ggplot(aes(Density.Annual, Associations)) + 
      geom_density_2d(aes(colour = SystemReplicate)) +
      theme(legend.position = "none")
    
  }) %>% ArrangeCowplot()

FocalSystem <- FocalSystems[1]

RankEffects <- 
  FocalSystems %>% 
  map(function(FocalSystem){
    
    LM <- 
      NodeData %>% 
      filter(SystemReplicate == FocalSystem) %>% 
      mutate_at(c("Density.Annual", "Associations"), ~c(scale(rank(.x)))) %>% 
      lm(data = ., Associations ~ Density.Annual)
    
    LM %>% 
      confint %>% data.frame %>% dplyr::select(Lower = 1, Upper = 2) %>% 
      # rownames_to_column %>% 
      slice(2) %>% 
      bind_cols(Density.Annual = data.frame(coef(LM))[2,], .)
    
  }) %>% bind_rows() %>% 
  mutate(SystemReplicate = FocalSystems) %>% 
  left_join(MetaDF) %>% 
  arrange(Density.Annual) %>% mutate_at("SystemReplicate", ~factor(.x, levels = unique(.x))) %>% 
  ggplot(aes(SystemReplicate, Density.Annual)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, colour = Behaviour)) +
  geom_point(aes(colour = Behaviour)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

OriginalEffects <- 
  FocalSystems %>% 
  map(function(FocalSystem){
    
    LM <- 
      NodeData %>% 
      filter(SystemReplicate == FocalSystem) %>% 
      mutate_at(c("Density.Annual", "Associations"), ~c(scale((.x)))) %>% 
      lm(data = ., Associations ~ Density.Annual)
    
    LM %>% 
      confint %>% data.frame %>% dplyr::select(Lower = 1, Upper = 2) %>% 
      # rownames_to_column %>% 
      slice(2) %>% 
      bind_cols(Density.Annual = data.frame(coef(LM))[2,], .)
    
  }) %>% bind_rows() %>% 
  mutate(SystemReplicate = FocalSystems) %>% 
  left_join(MetaDF) %>% 
  arrange(Density.Annual) %>% mutate_at("SystemReplicate", ~factor(.x, levels = unique(.x))) %>% 
  ggplot(aes(SystemReplicate, Density.Annual)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, colour = Behaviour)) +
  geom_point(aes(colour = Behaviour)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

(OriginalEffects|RankEffects) + 
  plot_layout(guides = "collect")

# Plotting both ####

NodeData %<>% 
  group_by(SystemReplicate, Site) %>% 
  summarise(XCentre = mean(X, na.rm = T), YCentre = mean(Y, na.rm = T)) %>% 
  left_join(NodeData, .)

NodeData$CentreDistance <- 
  NodeData %>% 
  Pythagoreg(c("X", "XCentre"), c("Y", "YCentre"))

OriginalvRankDF <- 
  FocalSystems %>% 
  map(function(FocalSystem){
    
    LM1 <- 
      NodeData %>% 
      filter(SystemReplicate == FocalSystem) %>% 
      mutate_at(c("Density.Annual", "Associations", "CentreDistance"), ~c(scale((.x)))) %>% 
      lm(data = ., Associations ~ Density.Annual + CentreDistance)
    
    LM2 <- 
      NodeData %>% 
      filter(SystemReplicate == FocalSystem) %>% 
      mutate_at(c("Density.Annual", "Associations", "CentreDistance"), ~c(scale(rank(.x)))) %>% 
      lm(data = ., Associations ~ Density.Annual + CentreDistance)
    
    LM1 %>% 
      confint %>% data.frame %>% dplyr::select(Lower = 1, Upper = 2) %>% 
      # rownames_to_column %>% 
      slice(2) %>% 
      bind_cols(Density.Annual = coef(LM1)[2], .) %>% 
      bind_cols(LM2 %>% 
                  confint %>% data.frame %>% dplyr::select(Lower = 1, Upper = 2) %>% 
                  # rownames_to_column %>% 
                  slice(2) %>% 
                  bind_cols(Density.Annual = coef(LM2)[2], .) %>% 
                  rename_all(~paste0(.x, ".Rank")))
    
  })

OriginalvRankDF %>% bind_rows() %>% mutate(SystemReplicate = FocalSystems) %>% 
  left_join(MetaDF) %>% 
  ggplot(aes(Density.Annual, Density.Annual.Rank)) +
  geom_errorbar(aes(ymin = Lower.Rank, ymax = Upper.Rank, colour = Behaviour), 
                alpha = 0.3,
                width = 0) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, colour = Behaviour),                 
                 alpha = 0.3,
                 width = 0.1) +
  geom_point(aes(colour = Behaviour), size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_cor()

OriginalvRankDF %>% bind_rows() %>% mutate(SystemReplicate = FocalSystems) %>% 
  left_join(MetaDF) %>% 
  arrange(Density.Annual.Rank) %>% mutate_at("SystemReplicate", ~factor(.x, levels = unique(.x))) %>% 
  ggplot(aes(SystemReplicate, Density.Annual.Rank)) +
  geom_errorbar(aes(ymin = Lower.Rank, ymax = Upper.Rank, colour = Behaviour)) +
  geom_point(aes(colour = Behaviour)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

OriginalvRankDF %>% bind_rows() %>% mutate(SystemReplicate = FocalSystems) %>% 
  left_join(MetaDF) %>% 
  SinaGraph("Behaviour", "Density.Annual.Rank")

# Splitting inflection point ####

FocalSystem <- FocalSystems[1]

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
    
    ReturnDF %<>% rename_at(vars(Density.Annual:SE), ~paste0(.x, "_Last"))
    
    # ReturnDF %<>% mutate(N.First = nrow(TestDF1), N.Middle = nrow(TestDF2), N.Last = nrow(TestDF3))
    
    ReturnDF
    
    # list(TestDF1, TestDF2, TestDF3) %>% return
    
  })

SplitPoints %<>% 
  bind_rows() %>% 
  mutate(SystemReplicate = FocalSystems %>% setdiff("KenyanElephants_GOG")) %>% 
  left_join(MetaDF[,c("SystemReplicate", "Behaviour")])

SplitPoints %>% pivot_longer(cols = c(Density.Annual_First, Density.Annual_Middle, Density.Annual_Last)) %>% 
  mutate_at("Behaviour", ~str_replace(.x, "Direct", "Proximity")) %>% 
  lm(data = ., value ~ name + Behaviour) %>% summary

SplitPoints %>% pivot_longer(cols = c(Density.Annual_First, SE_First, Density.Annual_Middle, SE_Middle, Density.Annual_Last, SE_Last)) %>% 
  separate(name, sep = "_", into = c("Effect", "Portion")) %>% 
  pivot_wider(names_from = Effect, values_from = value) %>% data.frame()

SplitPoints %>% pivot_longer(cols = c(Density.Annual_First, SE_First, Density.Annual_Middle, SE_Middle, Density.Annual_Last, SE_Last)) %>% 
  separate(name, sep = "_", into = c("Effect", "Portion")) %>% 
  pivot_wider(names_from = Effect, values_from = value) %>% data.frame() %>% 
  ggplot(aes(SystemReplicate, Density.Annual, colour = Portion)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, shape = Behaviour), 
                alpha = 0.3,
                width = 0) +
  geom_point(aes(colour = Portion, shape = Behaviour), size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

SplitPoints %>% 
  ggplot(aes(Density.Annual.x, Density.Annual.y)) +
  geom_errorbar(aes(ymin = Density.Annual.y - SE.y, ymax = Density.Annual.y + SE.y, colour = Behaviour), 
                alpha = 0.3,
                width = 0) +
  geom_errorbarh(aes(xmin = Density.Annual.x - SE.x, xmax = Density.Annual.x + SE.x, colour = Behaviour),                 
                 alpha = 0.3) +
  geom_point(aes(colour = Behaviour), size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  NULL



