
# 05d_Strength Meta-Analysis.R ####

{
  
  library(tidyverse); library(ggregplot); library(INLA); library(magrittr); library(cowplot); library(colorspace)
  library(patchwork); library(glue); library(ggpointdensity); library(broom); library(ggforce); library(fs)
  library(metafor); library(ggplot2); library(MuMIn); library(tidyr)
  library(rotl); library(ape); library(corrplot); library(emmeans); library(patchwork)
  
  source("R/00_DanFunctions.R")
  
  source("R/04_Combining Summary Data.R")
  
  MetaDF <- read.csv("Output/MetaDF.csv")
  
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
  # mutate(SpatialBehaviour = c("NonSpatial", "Spatial")[as.numeric(ContactType == "HRO") + 1]) %>% 
  # mutate(Rep = paste0(Behaviour, "_", Response)) %>% 
  mutate(Rep = Response)

# MetaDF %<>% filter(Response != "Degree")

# Adding Replicates ####

ModelReps <- MetaDF$Rep %>% unique
FocalRep <- ModelReps[1]

dir_create("Output/StrengthMetaAnalysis")

# for(FocalRep in ModelReps){

print(FocalRep)

TestDF <- MetaDF %>% filter(Rep == FocalRep)

## phylogenetic meta for variation

# rem = rma.mv(yi = yi,V = vi,
#              random = list(~1|Study/Obs, ~1|Tip, ~1|Species),
#              R = list(Tip = cmatrix[TestDF$Tip, TestDF$Tip]),
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

ExplanatoryList <- list(c("Behaviour", "HostGroup", "NID", "Area", "NYears"),
                        c("Behaviour", "NID", "Area", "NYears"),
                        c("HostGroup", "NID", "Area", "NYears"))

# ExplanatoryList <- list(c("HostGroup", "NID", "Area", "NYears"),
#                         c("NID", "Area", "NYears"))

ModelList <- list()

MetaEstimates <- 
  ExplanatoryList %>% 
  map(function(Covar){
    
    print(Covar)
    
    Formula <- paste("~", paste0(Covar, collapse = " + ")) %>% as.formula
    
    Model <- rma.mv(yi = yi, V = vi,
                    random = list(~1|Study/Obs, ~1|Tip, 
                                  ~1|Species),
                    # R = list(Tip = cmatrix),
                    method = "ML", 
                    control = list(optimizer = "optim", optmethod = "BFGS"),
                    mods = Formula,
                    data = TestDF)
    
    ModelList[[paste0(Covar, collapse = " + ")]] <<- Model
    
    REMLMod = update(Model, method = "REML")
    
    data.frame(Predictor = paste0(Covar, collapse = " + "),
               K = length(coef(Model)),
               AICc = AICc(Model),
               R2 = r2f(REMLMod))
    
  }) %>% bind_rows()

Formula <- paste("~", paste0(ExplanatoryList[[2]], collapse = " + ")) %>% as.formula

Model <- rma.mv(yi = yi, V = vi,
                random = list(~1|Study/Obs, #~1|Tip, 
                              ~1|Species),
                # R = list(Tip = cmatrix),
                method = "ML", 
                control = list(optimizer = "optim", optmethod = "BFGS"),
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

MetaEstimates %>% write.csv(paste0("Output/StrengthMetaAnalysis/", FocalRep, "_Table SX.csv"))

Model %>% saveRDS(paste0("Output/StrengthMetaAnalysis/", FocalRep, "LinearMetaModel.rds"))

Top <- Model #List[[MetaEstimates[1, "Predictor"]]]

msum = summary(Top)

Top %>% summary %>% extract(c("beta", "se", "zval", "pval")) %>% 
  data.frame %>% 
  mutate_all(~round(.x, 2)) %>% rename_all(~.x %>% str_remove("val") %>% CamelConvert) %>% 
  rownames_to_column("Coef") %>% 
  write.csv(paste0("Output/StrengthMetaAnalysis/", FocalRep, "_Table SY.csv"))

# Predicting from the model ####

FullPredDF <- 
  MakePredictDF(MetaDF[,c(ExplanatoryList[[2]])], 
                HoldNumeric = c("NID", "Area")) %>% 
  expand.grid()

ModelMatrix <- model.matrix(Formula, data = FullPredDF)[,-1]

FittedValues  <- predict(Top, 
                         newmods = ModelMatrix, 
                         transf = transf.ztor)

FullPredDF <- FittedValues %>% data.frame %>% 
  dplyr::select(Fit = pred, Lower = 2, Upper = 3) %>% 
  bind_cols(FullPredDF, .)

FullPredDF %>% saveRDS(paste0("Output/StrengthMetaAnalysis/", FocalRep, "_LinearMetaPredictions.rds"))
MetaDF %>% write.csv(paste0("Output/StrengthMetaAnalysis/", FocalRep, "_MetaDFOutput.csv"))

# }

# Quick Check ####

MetaOutputs <- 
  "Output/StrengthMetaAnalysis" %>% dir_ls(regex = "LinearMetaPredictions") %>% map(readRDS) %>% 
  map(~.x %>% slice(n())) %>% 
  bind_rows(.id = "Rep") %>% 
  mutate_at("Rep", ~.x %>% str_split("/") %>% map_chr(last) %>% str_remove("_LinearMetaPredictions.rds")) %>% 
  filter(Rep != "Indirect_Associations") %>% 
  separate(Rep, sep = "_", into = c("Behaviour", "Response")) %>% 
  mutate_at("Behaviour", ~factor(.x, levels = levels(MetaDF$Behaviour)))

MetaDF %>% 
  filter(Rep != "Indirect_Associations") %>% 
  ggplot(aes(Response, Estimate, colour = Behaviour)) + 
  # geom_point(position = position_jitter(w = 0.5)) +
  geom_sina(alpha = 0.45) +
  geom_errorbar(data = MetaOutputs, aes(y = Fit, group = Behaviour, ymin = Lower, ymax = Upper),
                position = position_dodge(w = 0.85),
                colour ="black", width = 0.3) +
  geom_point(data = MetaOutputs, aes(y = Fit, group = Behaviour), 
             position = position_dodge(w = 0.85),
             colour ="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(values = c(AlberColours[[2]], AlberColours[[1]]))



