anno <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B1_2.csv")

sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/breast_qcd/B1_2_qcd_artefact_yes.rds")

CD <- as.data.frame(colData(sce))

CD <- CD %>% 
  select(-Region) %>%
  left_join(anno)
head(CD)

sce$Region <- CD$Region
sce <- sce[, sce$Region != "Artefact_Fold_exclude"]

saveRDS(sce, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/breast_qcd/B1_2_qcd.rds")
