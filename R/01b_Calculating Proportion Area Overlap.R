
# X_PropAreaOverlap.R ####

library(dplyr);library(sf);library(ggplot2);library(patchwork);library(ggpubr); library(tidyverse)
library(magrittr)

# loading the data ####

spat_df <- readRDS("Intermediate/AllObservationDF.rds")

spat_df$`Datasets/FieldingCows` %<>% mutate(Year = Site)

spat_df %<>% 
  map(function(a){
    
    if(any(str_detect(colnames(a), "Site"))){
      
      a %>% mutate(Rep = paste0(Site, "_", Year)) %>% return
      
    }else{
      
      a %>% rename(Rep = Year) %>% rename
      
    }
    
  })

# Function to calculate centroids and bounding box area

calculate_centroid_area <- function(df) {
  
  census_centroid <- df %>%
    dplyr::group_by(Id, Rep) %>%
    dplyr::summarise(X = mean(X, na.rm = TRUE), Y = mean(Y, na.rm = TRUE))
  
  census_area <- st_as_sf(census_centroid, coords = c("X", "Y")) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_area()
  
  centroid_radius <- as.numeric(sqrt(census_area / (20 * pi)))
  
  list(centroid = census_centroid, radius = centroid_radius)
  
}

# Function to calculate points within buffer radius
calculate_point_radius <- function(df, year, buffer_dist) {
  census_year <- df %>% filter(Rep == year)
  census_year_sf <- st_as_sf(census_year, coords = c("X", "Y"))
  
  circle <- st_buffer(census_year_sf, dist = buffer_dist)
  point_count <- as.data.frame(st_intersects(circle, census_year_sf) |> lengths()) %>%
    dplyr::rename(RadCount = 1) %>%
    dplyr::mutate(RadCount = RadCount - 1, Id = census_year$Id, Rep = year) %>%
    dplyr::select(Id, RadCount, Rep)
  
  return(point_count)
}

# Function to process each species
process_species <- function(df) {
  
  df <- df %>% dplyr::group_by(Rep, Id) %>% dplyr::arrange(Rep, Id)
  
  centroid_info <- calculate_centroid_area(df)
  
  # years <- seq(min(df$Year), max(df$Year), 1)
  
  years <- df$Rep %>% unique %>% sort
  
  points_count <- bind_rows(lapply(years, calculate_point_radius, df = centroid_info$centroid, buffer_dist = centroid_info$radius))
  
  return(points_count)
  
}

# Function to generate plots
generate_plot <- function(localdensity, title) {
  
  ggplot(localdensity, aes(x = RadCount, y = Density.Annual)) +
    geom_point() +
    geom_smooth() +
    xlab("1/20th total area - based \nLocal Density") +
    ylab("KD-based \nLocal Density") +
    ggtitle(title) +
    stat_cor(method = "pearson", label.x.npc=0.75, label.y.npc=0.025, aes(label = ..r.label..))
  
}

spat_df %<>% map(process_species)

spat_df %>% 
  saveRDS("Intermediate/PropAreaOverlap_list.rds")

# Trying my way ####

NodeData <- readRDS("Intermediate/NodeData.rds")

RadiusData <- 
  NodeData %>% 
  dplyr::select(Id, Year, Density.Annual, System, Site, X, Y, 
                SystemSizeAlpha) %>% 
  # group_by(Id, Year, Density.Annual, System, Site, X, Y) %>% 
  # summarise_at("SystemSizeAlpha", max) %>% 
  # dplyr::select(Id, Year, Density.Annual, System, Site, X, Y) %>% nrow
  unique %>% 
  split.data.frame(.$System) %>% 
  map(function(a){
    
    if(!all(a$Site == "Site")){
      
      a %>% mutate(Rep = paste0(Site, "_", Year)) %>% return
      
    }else{
      
      a %>% rename(Rep = Year) %>% rename
      
    }
    
  })

RadiusData <- 
  RadiusData %>% 
  map(process_species) %>% 
  bind_rows(.id = "System") %>% 
  left_join(RadiusData %>% bind_rows(.id = "System"), by = c("Id", "Rep", "System"))

library(cowplot)

theme_set(theme_cowplot())

RadiusData %>% saveRDS("Intermediate/RadiusData.rds")

RadiusData %>% 
  group_by(System, Site) %>% 
  mutate_at("Density.Annual", ~c(scale(.x))) %>% 
  ggplot(aes(Density.Annual, RadCount)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~System, scales = "free") +
  ggpubr::stat_cor()


