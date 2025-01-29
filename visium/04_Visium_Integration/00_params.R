if(disease == "breast"){
  sample_names <- c("B1_2_OPHI", "B1_4_OPHI", "B2_2_1256", "B3_2_1GVR", "B4_2_1FHZ")
  save_names   <- c("B1_2",      "B1_4",      "B2_2",      "B3_2",      "B4_2"     )
  geo_names    <- c("B1-1",      "B1-3",      "B2-1",      "B3-1",      "B4-1"     )
  nsamples = 5
}else if(disease == "lung"){
  sample_names <- c("L1_2_0PSV", "L1_4_0PSV", "L2_2_0WMU", "L3_2_1G73", "L4_2_1GA2")
  save_names   <- c("L1_2",      "L1_4",      "L2_2",      "L3_2",      "L4_2"     )
  geo_names    <- c("L1-1",      "None",      "L2-1",      "L3-3",      "L4-1"     )
  nsamples = 5
}else if(disease == "dlbcl"){
  sample_names <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
  save_names   <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
  nsamples = 6
}