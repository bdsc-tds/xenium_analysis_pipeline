library(readxl)
old_anno <- read_excel("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/Visium_combined_annotations_CS_ED.xlsx")
lung_anno <- old_anno %>%
  filter(Section %in% c("L1_2", "L1_4", "L2_2", "L4_2")) %>%
  select(Section, Barcode, Region)

lung_12 <- lung_anno %>% filter(Section == "L1_2") %>% select(-Section)
write.csv(lung_12, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/L1_2.csv", row.names = FALSE)
lung_12 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/L1_2.csv")
head(lung_12)
#              Barcode            Region             
# 1 AACACTTGGCAAGGAA-1 Tumor_Stroma_mix   
# 2 AACAGGATTCATAGTT-1 Intratumoral_Stroma
# 3 AACAGGTTATTGCACC-1 Intratumoral_Stroma
# 4 AACCACTGCCATAGCC-1 Intratumoral_Vessel
# 5 AACCGCCAGACTACTT-1 Tumor_pure         
# 6 AACCTACTGTAACTCA-1 Tumor_Stroma_mix   

lung_14 <- lung_anno %>% filter(Section == "L1_4") %>% select(-Section)
write.csv(lung_14, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/L1_4.csv", row.names = FALSE)

lung_22 <- lung_anno %>% filter(Section == "L2_2") %>% select(-Section)
write.csv(lung_22, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/L2_2.csv", row.names = FALSE)

lung_32 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/L3_2.csv")
lung_32 <- lung_32 %>% 
  rename(Region = KvL)
head(lung_32)
#              Barcode              Region
# 1 AACACTTGGCAAGGAA-1     Immune_Cell_mix
# 2 AACAGGATTCATAGTT-1     Immune_Cell_mix
# 3 AACAGGTTATTGCACC-1   Most_likely_Tumor
# 4 AACAGGTTCACCGAAG-1         Lymphocytes
# 5 AACAGTCAGGCTCCGC-1 Intratumoral_Stroma
# 6 AACATACTCATATGCG-1     Immune_Cell_mix
lung_32 <- lung_32 %>%
  mutate(Region = ifelse(Region == "Artefact_Fold_exclude", "Exclude", Region))
lung_32 <- lung_32 %>%
  mutate(Region = ifelse(Region == "", "Exclude", Region))
write.csv(lung_32, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/L3_2.csv", row.names = FALSE)

lung_42 <- lung_anno %>% filter(Section == "L4_2") %>% select(-Section)
write.csv(lung_42, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/L4_2.csv", row.names = FALSE)


# Breast ------------------------------------------------------------------
breast_14 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B1_4.csv")
breast_14 <- breast_14 %>% rename(Region = CS)
breast_14 <- breast_14 %>%
  mutate(Region = ifelse(Region == "Artefact_Fold_exclude", "Exclude", Region))
write.csv(breast_14, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B1_4.csv", row.names = FALSE)

breast_22 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B2_2.csv")
breast_22 <- breast_22 %>% rename(Region = CS)
write.csv(breast_22, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B2_2.csv", row.names = FALSE)

breast_32 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B3_2.csv")
breast_32 <- breast_32 %>% rename(Region = CS)
breast_32 <- breast_32 %>%
  mutate(Region = ifelse(Region == "Fold_Exclude", "Exclude", Region))
write.csv(breast_32, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B3_2.csv", row.names = FALSE)

breast_42 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B4_2.csv")
breast_42 <- breast_42 %>% rename(Region = CS)
breast_42 <- breast_42 %>%
  mutate(Region = ifelse(Region == "", "Exclude", Region))
write.csv(breast_42, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/B4_2.csv", row.names = FALSE)


# DLBCL -----------------------------------------------------------------
dlbcl_1 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_1.csv")
dlbcl_1 <- dlbcl_1 %>% rename(Region = CS)
dlbcl_1 <- dlbcl_1 %>%
  mutate(Region = ifelse(Region == "Empty", "Exclude", Region))
write.csv(dlbcl_1, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_1.csv", row.names = FALSE)

dlbcl_2 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_2.csv")
dlbcl_2 <- dlbcl_2 %>% rename(Region = CS)
dlbcl_2 <- dlbcl_2 %>%
  mutate(Region = ifelse(Region == "Empty", "Exclude", Region))
write.csv(dlbcl_2, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_2.csv", row.names = FALSE)

dlbcl_3 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_3.csv")
dlbcl_3 <- dlbcl_3 %>% rename(Region = CS)
write.csv(dlbcl_3, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_3.csv", row.names = FALSE)

dlbcl_4 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_4.csv")
dlbcl_4 <- dlbcl_4 %>% rename(Region = CS)
dlbcl_4 <- dlbcl_4 %>%
  mutate(Region = ifelse(Region == "Empty", "Exclude", Region))
write.csv(dlbcl_4, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_4.csv", row.names = FALSE)

dlbcl_5 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_5.csv")
dlbcl_5 <- dlbcl_5 %>% rename(Region = CS)
dlbcl_5 <- dlbcl_5 %>%
  mutate(Region = ifelse(Region == "Empty", "Exclude", Region))
write.csv(dlbcl_5, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_5.csv", row.names = FALSE)

dlbcl_6 <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_6.csv")
dlbcl_6 <- dlbcl_6 %>% rename(Region = CS)
dlbcl_6 <- dlbcl_6 %>%
  mutate(Region = ifelse(Region == "Empty", "Exclude", Region))
write.csv(dlbcl_6, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_6.csv", row.names = FALSE)









