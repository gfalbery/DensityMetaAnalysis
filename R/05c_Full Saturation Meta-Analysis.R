
# 05c_Full Saturation Meta-Analysis.R ####

{
  
  library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
  library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)
  library(metafor); library(ggplot2); library(MuMIn); library(tidyr)
  library(rotl); library(ape); library(corrplot); library(emmeans); library(patchwork)
  
  source("R/00_DanFunctions.R")
  
  source("R/04_Combining Summary Data.R")
  
  MetaDF <- read.csv("Output/SaturationMetaDF.csv")
  
  MetaDF %<>% na.omit
  
  MetaDF %<>% 
    mutate(rArea = Area) %>% 
    mutate_at("Area", log10)
  
  MetaDF %<>% 
    mutate(rBehaviour = Behaviour) %>% 
    mutate(Behaviour = factor(str_replace_all(Behaviour, c("Direct" = "Proximity")), levels = c("Proximity", "Indirect")))
  
}

## extract systems

MetaDF %>% 
  separate(SystemReplicate, sep = "_", into = c("System", "ContactType"), remove = F) %>% 
  dplyr::select(c("System", "ContactType")) ->
  MetaDF[,c("System", "ContactType")]

# Sorting Phylogeny ####

Species <- MetaDF$Species %>% unique

SummaryTable %>%
  pull(Species) %>% sort %>% 
  tnrs_match_names(do_approximate_matching = T ) %>% 
  data.frame %>% 
  pull(ott_id) %>% 
  tol_induced_subtree(ott_ids = .)

Phylogeny <- Species %>% 
  tnrs_match_names(do_approximate_matching = T) %>% 
  data.frame

MetaDF %<>% 
  mutate(search_string = Species %>% tolower) %>% 
  left_join(Phylogeny, by = "search_string")

MetaDF[,c("approximate_match", "is_synonym", "flags", "number_matches")] <- NULL

tree = tol_induced_subtree(ott_ids = Phylogeny$ott_id)

tree$tip.label = strip_ott_ids(tree$tip.label) # remove ott information from the tips

is.binary(tree) # check binary

tree = compute.brlen(tree,method = "Grafen") # assign branch lengths
is.ultrametric(tree)

tree = ladderize(tree)

saveRDS(tree, "Intermediate/AlberTree.rds")

cmatrix = vcv.phylo(tree,cor = T)

mean(cmatrix)
median(cmatrix)
cmat2 = cmatrix
diag(cmat2) = NA
cmat2 = cmat2[lower.tri(cmat2)]
mean(cmat2)
median(cmat2)

MetaDF %<>% mutate(Tip = str_replace_all(unique_name, " ", "_"))

MetaDF %<>% 
  mutate(Obs = 1:n()) %>% 
  mutate(Study = as.factor(SystemReplicate)) %>% 
  mutate(Species = Tip)

MetaDF %<>% mutate(R = Estimate*sqrt(R2)/abs(Estimate))

MetaDF %<>% escalc(data = ., ri = R, ni = N, measure = "ZCOR")

MetaDF %<>% 
  summary(transf = transf.ztor) %>% 
  dplyr::select(RLower = ci.lb, RUpper = ci.ub) %>% 
  bind_cols(MetaDF, .) %>% data.frame

MetaDF %<>% 
  mutate(Rep = Response)

# MetaDF %<>% filter(Response != "Degree")

# Adding Replicates ####

ModelReps <- MetaDF$Rep %>% unique
FocalRep <- ModelReps[2]

MetaDF %<>% 
  arrange(HostGroup) %>% 
  mutate_at("HostGroup", ~factor(.x, levels = sort(unique(.x))))

dir_create("Output/SaturationMetaAnalyses")

for(FocalRep in ModelReps){
  
  print(FocalRep)
  
  TestDF <- MetaDF %>% filter(Rep == FocalRep)
  
  ## phylogenetic meta for variation
  
  # rem = rma.mv(yi = yi,V = vi,
  #              random = list(~1|Study/Obs, ~1|Tip, ~1|Species),
  #              R = list(Tip = cmatrix[TestDF$Tip, TestDF$Tip]),
  #              control = list(optimizer = "optim", optmethod = "BFGS"),
  #              method = "REML", mods = ~1, data = TestDF)
  # 
  # ## back transform
  # 
  # transf.ztor(predict(rem)$pred)
  # 
  # r2f(rem) ## works
  # 
  # summary(rem)
  # 
  # ## I2
  # rem_i2 = i2(rem)
  # round(rem_i2$I2/100,3)
  # round(rem_i2$allI2/100,3)
  # 
  # ## other Is
  # round(rem_i2$allI2[1]/sum(rem_i2$allI2),2)
  # round(rem_i2$allI2[2]/sum(rem_i2$allI2),2) ## h2
  # round(rem_i2$allI2[3]/sum(rem_i2$allI2),2)
  
  # Explanatory models ####
  
  # ExplanatoryList <- list(c("Behaviour", "HostGroup", "NID", "Area", "NYears"),
  #                         c("Behaviour", "NID", "Area", "NYears"),
  #                         c("HostGroup", "NID", "Area", "NYears"))
  
  ExplanatoryList <- list(c("Portion", "Behaviour", "Portion:Behaviour", "NID", "Area", "NYears"),
                          c("Portion", "Behaviour", "NID", "Area", "NYears"),
                          c("HostGroup", "NID", "Area", "NYears"),
                          c("NID", "Area", "NYears"))[1:2]
  
  ModelList <- list()
  
  MetaEstimates <- 
    ExplanatoryList %>% 
    map(function(Covar){
      
      print(Covar)
      
      Formula <- paste("~", paste0(Covar, collapse = " + ")) %>% as.formula
      
      Model <- rma.mv(yi = yi, V = vi,
                      random = list(~1|Study/Obs, ~1|Tip,
                                    ~1|Species),
                      R = list(Tip = cmatrix),
                      control = list(optimizer = "optim", optmethod = "BFGS"),
                      method = "ML", 
                      mods = Formula,
                      data = TestDF %>% filter(yi > -0.5))
      
      ModelList[[paste0(Covar, collapse = " + ")]] <<- Model
      
      REMLMod = update(Model, method = "REML")
      
      data.frame(Predictor = paste0(Covar, collapse = " + "),
                 K = length(coef(Model)),
                 AICc = AICc(Model),
                 R2 = r2f(REMLMod))
      
    }) %>% bind_rows(.id = "X")
  
  Formula <- paste("~", paste0(ExplanatoryList[[2]], collapse = " + ")) %>% as.formula
  
  Model <- rma.mv(yi = yi, V = vi,
                  random = list(~1|Study/Obs, ~1|Tip, 
                                ~1|Species),
                  R = list(Tip = cmatrix),
                  control = list(optimizer = "optim", optmethod = "BFGS"),
                  method = "ML", 
                  mods = Formula,
                  data = TestDF)
  
  # REMLMod = update(Model, method = "REML")
  REMLMod = Model
  
  (MetaEstimates %<>% 
      bind_rows(
        data.frame(Predictor = paste0(paste0(ExplanatoryList[[2]], collapse = " + "), " - Phylo"),
                   K = length(coef(Model)),
                   AICc = AICc(Model),
                   R2 = r2f(REMLMod))) %>% 
      arrange(AICc) %>% 
      mutate(Delta = round(AICc - min(AICc), 2)) %>% 
      mutate(wi = round(Weights(AICc), 2)))
  
  MetaEstimates %>% write.csv(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_Table SX.csv"),
                              row.names = F)
  
  Model %>% saveRDS(paste0("Output/SaturationMetaAnalyses/", FocalRep, "LinearMetaModel.rds"))
  
  Top <- Model #List[[MetaEstimates[1, "Predictor"]]]
  
  msum = summary(Top)
  
  Top %>% summary %>% extract(c("beta", "se", "zval", "pval", "ci.lb", "ci.ub")) %>% 
    data.frame %>% 
    mutate_all(~round(.x, 2)) %>% rename_all(~.x %>% str_remove("val") %>% CamelConvert) %>% 
    rownames_to_column("Coef") %>% 
    write.csv(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_Table SY.csv"))
  
  # Predicting from the model ####
  
  FullPredDF <- 
    MakePredictDF(MetaDF[,c(ExplanatoryList[[2]])], 
                  HoldNumeric = c("NID", "Area")) %>% 
    expand.grid()
  
  ModelMatrix <- model.matrix(Formula, data = FullPredDF)[,-c(1)]
  
  FittedValues  <- predict(Top, 
                           newmods = ModelMatrix, 
                           transf = transf.ztor)
  
  FullPredDF <- FittedValues %>% data.frame %>% 
    dplyr::select(Fit = pred, Lower = 2, Upper = 3) %>% 
    bind_cols(FullPredDF, .)
  
  FullPredDF %>% saveRDS(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_LinearMetaPredictions.rds"))
  MetaDF %>% write.csv(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_MetaDFOutput.csv"),
                       row.names = F)
  
}

# Quick Check ####

if(0){
  
  MetaOutputs <- 
    "Output/SaturationMetaAnalyses" %>% dir_ls(regex = "LinearMetaPredictions") %>% map(readRDS) %>% 
    map(~.x %>% group_by(Portion) %>% slice(n())) %>% 
    bind_rows(.id = "Rep") %>% 
    mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_LinearMetaPredictions.rds")) %>% 
    filter(Rep != "Indirect_Associations") %>% 
    mutate_at("Rep", ~factor(.x, levels = c("Proximity_Associations", paste0(rep(c("Proximity", "Indirect")), 
                                                                             "_",
                                                                             rep(c("Degree", "Strength"), each = 2))))) %>% 
    separate(Rep, sep = "_", into = c("Behaviour", "Response"), remove = F) %>%
    mutate_at("Behaviour", ~factor(.x, levels = levels(MetaDF$Behaviour)))
  
  MetaDF %>% 
    filter(Rep != "Indirect_Associations") %>% filter(Estimate < 2.5, Estimate > -2) %>% 
    mutate_at("Rep", ~factor(.x, levels = c("Proximity_Associations", paste0(rep(c("Proximity", "Indirect")), 
                                                                             "_",
                                                                             rep(c("Degree", "Strength"), each = 2))))) %>% 
    # mutate_at("Estimate", rank) %>% 
    ggplot(aes(Rep, Estimate, colour = Portion)) + 
    # geom_point(position = position_jitter(w = 0.5)) +
    # geom_sina(alpha = 0.45) +
    geom_boxplot() +
    # geom_label(aes(label = SystemReplicate)) +
    geom_errorbar(data = MetaOutputs, aes(y = Fit, group = Portion, ymin = Lower, ymax = Upper),
                  position = position_dodge(w = 0.85),
                  colour ="black", width = 0.3) +
    geom_point(data = MetaOutputs, aes(y = Fit, group = Portion), 
               position = position_dodge(w = 0.85),
               colour ="black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_manual(values = c(AlberColours[[2]], AlberColours[[1]]))
  
  # Comparing model effects ###
  
  "Output/SaturationMetaAnalyses/" %>% dir_ls(regex = "LinearMetaModel.rds") %>% 
    map(readRDS)
  
  "Output/SaturationMetaAnalyses/" %>% dir_ls(regex = "Table SX") %>% 
    map(read.csv) %>% 
    bind_rows(.id = "Rep") %>% 
    mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_Table SX.csv")) %>% 
    filter(Rep != "Indirect_Associations") %>% 
    # filter(Delta < 2) %>% 
    filter(X == 1)
  
  # Looking at full outputs ####
  
  "Output/SaturationMetaAnalyses/" %>% dir_ls(regex = "Table SY") %>% 
    map(read.csv) %>% 
    bind_rows(.id = "Rep") %>% 
    mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_Table SY.csv")) %>% 
    filter(Rep != "Indirect_Associations") %>% 
    separate(Rep, sep = "_", into = c("Behaviour", "Response"), remove = F) %>% 
    mutate_at("Behaviour", ~factor(.x, levels = levels(MetaDF$Behaviour))) %>% 
    mutate_at("Rep", ~factor(.x, levels = c("Proximity_Associations", paste0(rep(c("Proximity", "Indirect")), 
                                                                             "_",
                                                                             rep(c("Degree", "Strength"), each = 2))))) %>% 
    filter(Coef == "PortionLast") %>% 
    # filter(Delta < 2) %>% 
    ggplot(aes(Response, Beta, colour = Behaviour)) + 
    geom_hline(yintercept = 0, lty = 2) +
    geom_errorbar(aes(#y = Fit, group = Behaviour, 
      ymin = Ci.lb, ymax = Ci.ub),
      position = position_dodge(w = 0.85), 
      width = 0.3) +
    geom_point(aes(#y = Fit, 
      group = Behaviour), 
      position = position_dodge(w = 0.85),
      colour ="black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Now testing and outputting first and last separately ####
  
  MetaDF %<>% mutate(Rep2 = paste0(Rep, "_", Portion))
  
  ModelReps <- MetaDF$Rep2 %>% unique()
  
  for(FocalRep in ModelReps){
    
    print(FocalRep)
    
    TestDF <- MetaDF %>% filter(Rep2 == FocalRep)
    
    ## phylogenetic meta for variation
    
    rem = rma.mv(yi = yi,V = vi,
                 random = list(~1|Study/Obs, ~1|Tip, ~1|Species),
                 R = list(Tip = cmatrix[TestDF$Tip, TestDF$Tip]),
                 control = list(optimizer = "optim", optmethod = "BFGS"),
                 method = "REML", mods = ~1, data = TestDF)
    
    ## back transform
    
    transf.ztor(predict(rem)$pred)
    
    r2f(rem) ## works
    
    summary(rem)
    
    ## I2
    rem_i2 = i2(rem)
    round(rem_i2$I2/100,3)
    round(rem_i2$allI2/100,3)
    
    ## other Is
    round(rem_i2$allI2[1]/sum(rem_i2$allI2),2)
    round(rem_i2$allI2[2]/sum(rem_i2$allI2),2) ## h2
    round(rem_i2$allI2[3]/sum(rem_i2$allI2),2)
    
    # Explanatory models ####
    
    # ExplanatoryList <- list(c("Behaviour", "HostGroup", "NID", "Area", "NYears"),
    #                         c("Behaviour", "NID", "Area", "NYears"),
    #                         c("HostGroup", "NID", "Area", "NYears"))
    
    ExplanatoryList <- list(c("HostGroup", "NID", "Area", "NYears"),
                            c("NID", "Area", "NYears"))
    
    ModelList <- list()
    
    MetaEstimates <- 
      ExplanatoryList %>% 
      map(function(Covar){
        
        print(Covar)
        
        Formula <- paste("~", paste0(Covar, collapse = " + ")) %>% as.formula
        
        Model <- rma.mv(yi = yi, V = vi,
                        random = list(~1|Study/Obs, ~1|Tip,
                                      ~1|Species),
                        R = list(Tip = cmatrix),
                        control = list(optimizer = "optim", optmethod = "BFGS"),
                        method = "ML", 
                        mods = Formula,
                        data = TestDF %>% filter(yi > -0.5))
        
        ModelList[[paste0(Covar, collapse = " + ")]] <<- Model
        
        REMLMod = update(Model, method = "REML")
        
        data.frame(Predictor = paste0(Covar, collapse = " + "),
                   K = length(coef(Model)),
                   AICc = AICc(Model),
                   R2 = r2f(REMLMod))
        
      }) %>% bind_rows(.id = "X")
    
    Formula <- paste("~", paste0(ExplanatoryList[[1]], collapse = " + ")) %>% as.formula
    
    Model <- rma.mv(yi = yi, V = vi,
                    random = list(~1|Study/Obs, ~1|Tip, 
                                  ~1|Species),
                    R = list(Tip = cmatrix),
                    control = list(optimizer = "optim", optmethod = "BFGS"),
                    method = "ML", 
                    mods = Formula,
                    data = TestDF)
    
    REMLMod = update(Model, method = "REML")
    
    (MetaEstimates %<>% 
        bind_rows(
          data.frame(Predictor = paste0(paste0(ExplanatoryList[[2]], collapse = " + "), " - Phylo"),
                     K = length(coef(Model)),
                     AICc = AICc(Model),
                     R2 = r2f(REMLMod))) %>% 
        arrange(AICc) %>% 
        mutate(Delta = round(AICc - min(AICc), 2)) %>% 
        mutate(wi = round(Weights(AICc), 2)))
    
    MetaEstimates %>% write.csv(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_Table SX.csv"),
                                row.names = F)
    
    Model %>% saveRDS(paste0("Output/SaturationMetaAnalyses/", FocalRep, "LinearMetaModel.rds"))
    
    Top <- Model #List[[MetaEstimates[1, "Predictor"]]]
    
    msum = summary(Top)
    
    Top %>% summary %>% extract(c("beta", "se", "zval", "pval", "ci.lb", "ci.ub")) %>% 
      data.frame %>% 
      mutate_all(~round(.x, 2)) %>% rename_all(~.x %>% str_remove("val") %>% CamelConvert) %>% 
      rownames_to_column("Coef") %>% 
      write.csv(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_Table SY.csv"))
    
    # Predicting from the model ####
    
    FullPredDF <- 
      MakePredictDF(MetaDF[,c(ExplanatoryList[[1]])], 
                    HoldNumeric = c("NID", "Area")) %>% 
      expand.grid()
    
    ModelMatrix <- model.matrix(Formula, data = FullPredDF)#[,-c(1)]
    
    FittedValues  <- predict(Top, 
                             newmods = ModelMatrix, 
                             transf = transf.ztor)
    
    FullPredDF <- FittedValues %>% data.frame %>% 
      dplyr::select(Fit = pred, Lower = 2, Upper = 3) %>% 
      bind_cols(FullPredDF, .)
    
    FullPredDF %>% saveRDS(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_LinearMetaPredictions.rds"))
    MetaDF %>% write.csv(paste0("Output/SaturationMetaAnalyses/", FocalRep, "_MetaDFOutput.csv"),
                         row.names = F)
    
  }
  
  "Output/SaturationMetaAnalyses/" %>% dir_ls(regex = "Table SY") %>% 
    map(read.csv) %>% 
    bind_rows(.id = "Rep") %>% filter(str_detect(Rep, "_First|_Last")) %>% 
    mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_Table SY.csv")) %>% 
    filter(Rep != "Indirect_Associations") %>% 
    separate(Rep, sep = "_", into = c("Behaviour", "Response", "Portion"), remove = F) %>% 
    mutate_at("Behaviour", ~factor(.x, levels = levels(MetaDF$Behaviour))) %>% 
    mutate_at("Rep", ~factor(.x, levels = c("Proximity_Associations", paste0(rep(c("Proximity", "Indirect")), 
                                                                             "_",
                                                                             rep(c("Degree", "Strength"), each = 2))))) %>% 
    filter(Coef == "PortionLast") %>% 
    # filter(Delta < 2) %>% 
    ggplot(aes(Response, Beta, colour = Behaviour)) + 
    geom_hline(yintercept = 0, lty = 2) +
    geom_errorbar(aes(#y = Fit, group = Behaviour, 
      ymin = Ci.lb, ymax = Ci.ub),
      position = position_dodge(w = 0.85), 
      width = 0.3) +
    geom_point(aes(#y = Fit, 
      group = Behaviour), 
      position = position_dodge(w = 0.85),
      colour ="black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Old code below #####
  
  ## phylogenetic meta for variation
  
  rem = rma.mv(yi = yi,V = vi,
               random = list(~1|Study/Obs, ~1|Tip, ~1|Species),
               R = list(Tip = cmatrix),
               method = "REML", mods = ~1, data = MetaDF)
  
  ## back transform
  
  transf.ztor(predict(rem)$pred)
  
  r2f(rem) ## works
  
  summary(rem)
  
  ## I2
  rem_i2 = i2(rem)
  round(rem_i2$I2/100,3)
  round(rem_i2$allI2/100,3)
  
  ## other Is
  round(rem_i2$allI2[1]/sum(rem_i2$allI2),2)
  round(rem_i2$allI2[2]/sum(rem_i2$allI2),2) ## h2
  round(rem_i2$allI2[3]/sum(rem_i2$allI2),2)
  
  # Explanatory models ####
  
  ExplanatoryList <- list(c("Portion", "Behaviour", "HostGroup", "NID", "Area", "NYears"),
                          c("Portion", "Behaviour", "NID", "Area", "NYears"),
                          c("Portion:Behaviour", "NID", "Area", "NYears"),
                          c("Portion*Behaviour", "NID", "Area", "NYears"),
                          c("Portion", "HostGroup", "NID", "Area", "NYears"))
  
  ModelList <- list()
  
  (MetaEstimates <- 
      ExplanatoryList %>% 
      map(function(Covar){
        
        print(Covar)
        
        Formula <- paste("~", paste0(Covar, collapse = " + ")) %>% as.formula
        
        Model <- rma.mv(yi = yi, V = vi,
                        random = list(~1|Study/Obs, ~1|Tip, 
                                      ~1|Species),
                        R = list(Tip = cmatrix),
                        method = "ML", 
                        mods = Formula,
                        data = MetaDF)
        
        ModelList[[paste0(Covar, collapse = " + ")]] <<- Model
        
        REMLMod = update(Model, method = "REML")
        
        data.frame(Predictor = paste0(Covar, collapse = " + "),
                   K = length(coef(Model)),
                   AICc = AICc(Model),
                   R2 = r2f(REMLMod))
        
      }) %>% bind_rows())
  
  Formula <- paste("~", paste0(ExplanatoryList[[2]], collapse = " + ")) %>% as.formula
  
  Model <- rma.mv(yi = yi, V = vi,
                  random = list(~1|Study/Obs, #~1|Tip, 
                                ~1|Species),
                  # R = list(Tip = cmatrix),
                  method = "ML", 
                  mods = Formula,
                  data = MetaDF)
  
  REMLMod = update(Model, method = "REML")
  
  (MetaEstimates %<>% 
      bind_rows(
        data.frame(Predictor = paste0(paste0(ExplanatoryList[[2]], collapse = " + "), " - Phylo"),
                   K = length(coef(Model)),
                   AICc = AICc(Model),
                   R2 = r2f(REMLMod))) %>% 
      arrange(AICc) %>% 
      mutate(Delta = round(AICc - min(AICc), 2)) %>% 
      mutate(wi = round(Weights(AICc), 2)))
  
  MetaEstimates %>% write.csv("Figures/Table SZ.csv")
  
  Model %>% saveRDS("Output/SaturationMetaModel.rds")
  
  Top <- Model #List[[MetaEstimates[1, "Predictor"]]]
  
  Top %>% summary %>% extract(c("beta", "se", "zval", "pval")) %>% 
    data.frame %>% 
    mutate_all(~round(.x, 2)) %>% rename_all(~.x %>% str_remove("val") %>% CamelConvert) %>% 
    rownames_to_column("Coef") %>% 
    write.csv("Figures/Table SZZ.csv")
  
  # Predicting from the model ####
  
  FullPredDF <- 
    MakePredictDF(MetaDF[,c(ExplanatoryList[[2]])], 
                  HoldNumeric = c("NID", "Area")) %>% 
    expand.grid()
  
  ModelMatrix <- model.matrix(Formula, data = FullPredDF)[,-1]
  
  FittedValues  <- predict(Top, newmods = ModelMatrix, 
                           transf = transf.ztor)
  
  FullPredDF <- FittedValues %>% data.frame %>% 
    dplyr::select(Fit = pred, Lower = 2, Upper = 3) %>% 
    bind_cols(FullPredDF, .)
  
  FullPredDF %>% saveRDS("Output/SaturationMetaPredictions.rds")
  MetaDF %>% write.csv("Output/SaturationMetaDFOutput.csv")
  
  # Now trying with change in slope?? ####
  
  # Plotting ####
  
  theme_set(theme_cowplot())
  
  FullPredDF %>% 
    filter(NYears == last(NYears), Behaviour == "Indirect") %>% 
    ggplot() +
    geom_hline(yintercept = 0,linetype = 2) +
    geom_jitter(data = MetaDF, aes(x = Portion, colour = Portion, y = R, size = 1/vi), width = 0.1, alpha = 0.25) +
    geom_segment(aes(x = Portion, xend = Portion,
                     y = Lower, yend = Upper)) +
    geom_point(aes(x = Portion, y = Fit, colour = Portion), size = 5, shape = 15) +
    guides(size = "none") +
    scale_size_continuous(range = c(1,3)) +
    theme(legend.position = "none")
  
  FullPredDF %>% 
    filter(NYears == last(NYears)) %>% 
    # mutate_at("Portion", ~paste0(.x, " 50%")) %>% 
    ggplot() +
    geom_hline(yintercept = 0,linetype = 2) +
    geom_jitter(data = MetaDF, 
                position = position_dodge(w = 0.75),
                aes(x = Portion, colour = Behaviour, y = R, size = 1/vi), alpha = 0.25) +
    geom_segment(#position = position_dodge(w = 0.5),
      aes(x = c(0.8, 1.8, 1.2, 2.2), xend = c(0.8, 1.8, 1.2, 2.2), group = Behaviour,
          y = Lower, yend = Upper)) +
    geom_point(position = position_dodge(w = 0.75),
               aes(x = Portion, y = Fit, colour = Behaviour), size = 5, shape = 15) +
    guides(size = "none") +
    scale_size_continuous(range = c(1,3)) +
    scale_x_discrete(labels = c("First 50%", "Last 50%"))
  
  ## repeat for Nyear
  newgrid = data.frame(expand.grid(NID = mean(data$NID),
                                   Area = mean(data$Area),
                                   NYears = seq(min(data$NYears),max(data$NYears),l = 50),
                                   Behaviour = levels(data$Behaviour))) 
  
  ## model matrix sans intercept
  predgrid = model.matrix(f,data = newgrid)[,-1]
  
  ## predict
  pred = predict(top,newmods = predgrid,transf = transf.ztor)
  
  ## attach predictions
  pred = data.frame(pred,newgrid)
  
  ## trim
  pred = pred[pred$Behaviour == "Proximity",]
  
  ## save
  ypred = pred
  
  ## plot
  p2 = ggplot()+
    geom_hline(yintercept = 0,linetype = 2)+
    geom_ribbon(data = ypred,aes(x = NYears,ymin = ci.lb,ymax = ci.ub),fill = "grey80")+
    geom_line(data = ypred,aes(x = NYears,y = pred),size = 2)+
    geom_jitter(data = data,aes(x = NYears,y = r,size = 1/vi),width = 0.1,alpha = 0.25)+
    th+
    guides(size = "none")+
    scale_size_continuous(range = c(1,3))+
    ylim(-0.2,0.8)
  
  ## combine plots
  p1+p2
  
  ## make pretty
  g1 = p1+
    ylab(expression(paste("correlation coefficient (",italic(r),")")))
  
  ## next
  g2 = p2+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  ## combine
  g1+g2+plot_layout(widths = c(1,0.75))
  
}