plan = drake_plan(
  
  # load data -------------------------------------------
  glycans = fread(here("input_data/glycans_corrected_20190409.csv")),
  glycans_raw = fread(here("input_data/glycans_raw.csv")),  
  # IPs
  ips = fread(here("input_data/immunopheno.corrected.csv")),
  #Popante output
  ip_igx_univar = fread(here("input_data/glycans_20190409_immunopheno_corrected_cleanedSimon.tsv"), header = TRUE) %>%
    set_names(names(.) %>% tolower(.) %>% str_replace_all(., c("\\^"="", "[\\.]+"="_"))),
  # annotations
  ip_anno = fread(here("input_data/data_annotations/all_immunophenotypes_annotation_av.csv")) %>%
    set_names(names(.) %>% tolower(.) %>% str_replace_all(., c( "\\.$"="", "[\\.]+"="_", " "="_"))) %>%
    mutate(composite_lin_source = ifelse((source=="Lin"|source=="MFI"), source, lineage)),
  IgA_anno_names_raw = readxl::read_xlsx(here("input_data/data_annotations/IgA_explanatory_overview.xlsx"), sheet = "raw_names"),
  IgA_anno_names_derived = readxl::read_xlsx(here("input_data/data_annotations/IgA_explanatory_overview.xlsx"), sheet = "derived_names"),
  twin_fam = fread(here("input_data/data_annotations/TwinDetails_110119.csv")),
  
  # pre-process -------------------------------------------
  
  # QC filter + FDR
  ip_igx_univar_good = ip_igx_univar %>%
    # dplyr::select(predictor) %>% 
    left_join(dplyr::select(ip_anno, set_name, robust_mario_qc), by=c("predictor"="set_name")) %>%
    dplyr::filter(robust_mario_qc == "Good") %>%
    dplyr::select(-robust_mario_qc) %>% 
    mutate(pv_adj_global = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue),
  
  # annotate with subgroups of IgX and IP
  ip_igx_univar_good_anno = ip_igx_univar_good %>%
    left_join(ip_anno %>% 
                dplyr::select(set_name, subset_name, composite_lin_source), 
              by=c("predictor"="set_name")) %>% 
    mutate(IgX = ifelse(str_detect(response, pattern = "IgG"), "IgG", "IgA")) %>% 
    dplyr::select(response, IgX, subset_name, composite_lin_source, everything()),
  
  # data exploration -------------------------------------------
  
  #qq plot
  qq = make_univar_qq(v_pval = ip_igx_univar_good$pvalue, fig_path = file_out("results/figures/univar_qqplot.png")),
  
  # glycan marginal distributions pdf
  glycan_marginal_pdf = make_glycan_marginal_plots(glycans_raw, fig_path = file_out("results/figures/glycan_marginal_plots.pdf")),
  
  # top association scatter plots pdf
  top_assoc_scatter_pdf = make_top_assoc_scatters(ntop = 100, ip_igx_univar_good, glycans, ips, file_path = file_out("results/figures/top_assoc_scatters.pdf")),
  
  # number of significant IgG, IgA and IP ~ FDR threshold
  n_associated_m = get_fdr_nfeat_associated_univar(fdr_thresholds = seq(0,1,0.01), ip_igx_univar_good = ip_igx_univar_good),
  
  # enrichment -------------------------------------------
  
  target_indices = unique(floor(10^seq(0,log10(nrow(ip_igx_univar_good)), length.out = 1e4))),
  
  # cumulative mean number of IgG among top associated IgX
  running_IgG_perc =  cummean(str_detect(ip_igx_univar_good$response, pattern = "IgG")),
  
  # enrichment of lin_sources ~ FDR
  lin_source_props_df = get_lin_source_props(target_indices,
                                             ip_igx_univar_good,
                                             ip_anno),
  # enrichment of lin_sources ~ FDR+IgX
  lin_source_IgX_props_df = get_lin_source_props_byIgX(target_indices,
                                             ip_igx_univar_good_anno),
  
  # WGCNA  -------------------------------------------
  
  
  glycan_cor = cor(x=glycans[,-(1:5)], use = "pairwise.complete.obs", method = "pearson"),
  
  # generate report -------------------------------------------
  
  report = target(
    command = rmarkdown::render(
      knitr_in("explore_univariate_igx.Rmd"),
      output_file = file_out(here::here("results/report.html")),
      quiet = TRUE
    ),
    trigger = trigger(change = digest(file.info(file_in("results/figures/univar_qqplot.png"))))  # problem: seems to interfere with automatic dependency detection in Rmd
  )
)
