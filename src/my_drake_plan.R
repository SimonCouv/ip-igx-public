plan = drake_plan(
  
  # load data -------------------------------------------
  glycans = fread(here("input_data/glycans_corrected_20190409.csv")),
  glycans_raw = fread(here("input_data/glycans_raw.csv")),  
  # IPs
  ips = fread(here("input_data/immunopheno.corrected.csv")),
  ips_raw = fread(here("input_data/immunopheno.raw.csv")),
  #Popante output
  ip_igx_univar = fread(here("input_data/glycans_20190409_immunopheno_corrected_cleanedSimon.tsv"), header = TRUE) %>%
    set_names(names(.) %>% tolower(.) %>% str_replace_all(., c("\\^"="", "[\\.]+"="_"))),
  # annotations
  ip_anno = fread(here("input_data/data_annotations/all_immunophenotypes_annotation_av.csv")) %>%
    mutate_all(na_if, y="") %>%   # set all empty string fields to NA 
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
  
  # missingness filter
  
  glycans_mf = glycans[rowMeans(is.na(glycans)) < 0.6, colMeans(is.na(glycans)) < 0.6, with=FALSE],
  glycans_raw_mf = glycans[rowMeans(is.na(glycans_raw)) < 0.6, colMeans(is.na(glycans_raw)) < 0.6, with=FALSE],
  
  # IP missingness
  ip_missing = ips[, lapply(.SD, function(x) mean(is.na(x))), .SDcols=-(FID:IID)] %>%
    data.table::transpose(.) %>% 
    .[,.(set_name = names(ips)[-(1:2)],
         missing_prop = V1)
      ] %>% 
    .[order(-missing_prop)],
  
  ip_raw_missing = ips_raw[, lapply(.SD, function(x) mean(is.na(x))), .SDcols=-(FID:IID)] %>% 
    data.table::transpose(.) %>% 
    .[,.(set_name = names(ips_raw)[-(1:2)],
         missing_prop = V1)
      ] %>% 
    .[order(-missing_prop)],
  
  # IP mean per-feature relative counts
  ip_raw_mean = ips_raw[, lapply(.SD, mean, na.rm=TRUE), .SDcols = -(FID:IID)] %>%
    data.table::transpose(.) %>%
    .[,.(set_name = colnames(ips_raw)[-(1:2)],
         mean_signal = V1)
      ],
  
  # data exploration -------------------------------------------
  
  #qq plot
  qq = make_univar_qq(v_pval = ip_igx_univar_good$pvalue, fig_path = file_out("results/figures/univar_qqplot.png")),
  
  # glycan marginal distributions pdf
  raw_glycan_marginal_pdf = make_marginal_plots(glycans_raw_mf[,-(1:5)], fig_path = file_out("results/figures/raw_glycan_marginal_plots.pdf")),  # important as input to lmer models
  glycan_marginal_pdf = make_marginal_plots(glycans_mf[,-(1:5)], fig_path = file_out("results/figures/glycan_marginal_plots.pdf")),              # important as input to a.o. association scatters
  ip_marginal_pdf = make_marginal_plots(ips[,3:100], fig_path = file_out("results/figures/ips_marginal_plots.pdf")),              # important as input to a.o. association scatters
  
  # top association scatter plots pdf
  top_assoc_scatter_pdf = make_top_assoc_scatters(ntop = 100, ip_igx_univar_good, glycans_mf, ips, file_path = file_out("results/figures/top_assoc_scatters.pdf")),
  
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
  
  # Correlation structures  -------------------------------------------
  
  glycans_fam_adj = glycans_lmer_adjust(glycans_raw_mf),
  glycan_residuals_pdf = make_marginal_plots(glycans_fam_adj, fig_path = file_out("results/figures/glycan_res_marginal_plots.pdf")),             # diagnostic of LME models
  
  
  # WGCNA  -------------------------------------------
  
  # 1. scale-free topology: power parameter
  
  scale_free_powers = c(c(1:10), seq(from = 12, to=20, by=2)),
  scale_free_common_pars = list(verbose= 0,  moreNetworkConcepts = TRUE, powerVector = scale_free_powers),
  scale_free_iter_pars = expand.grid(networktype = c('signed', 'signed hybrid', 'unsigned'),
                                      corFnc = c('cor', 'bicor'),
                                      stringsAsFactors = FALSE),
  derived_glycan_pos = 108:151,
  
  glycans_scale_free = optimize_network_pars(data=glycans_fam_adj, common_pars = scale_free_common_pars, pars = scale_free_iter_pars, powers = scale_free_powers, plt_title = "glycans"),
  glycans_scale_free_noderiv = optimize_network_pars(data=glycans_fam_adj[, -derived_glycan_pos], scale_free_common_pars, scale_free_iter_pars, scale_free_powers, plt_title = "glycans"),
  
  # 2. parameter search ~ modularity
  
  glycan_cross_iter_pars = list(figdir = here(figdir, "WGCNA_pars/glycans/")),
  
  glycan_network_pars = read_tsv(here(cluster_pars_dir,"glycan_pars.tsv")) %>%
    dplyr::select(-aim),
  glycan_cut_pars = expand.grid(deep =c(4, 1),
                                 minModuleSize=c(3,5,10,20)),
  glycan_modality = ifelse(str_detect(names(glycans_fam_adj), pattern = "IgG"),"IgG","IgA"),
  
  glycan_module_stats = wgcna_parameter_search(data=glycans_fam_adj,
                                               glycan_network_pars, 
                                               glycan_cut_pars,
                                               glycan_cross_iter_pars, 
                                               create_plot = TRUE,
                                               feat_modality=glycan_modality),
  glycan_module_stats_noderiv = wgcna_parameter_search(data=glycans_fam_adj[,-derived_glycan_pos],
                                                       glycan_network_pars, 
                                                       glycan_cut_pars,
                                                       glycan_cross_iter_pars, 
                                                       create_plot = TRUE,
                                                       feat_modality=glycan_modality[-derived_glycan_pos]),
  
  network_module_heatmaps = make_module_heatmap_pdf(plots = glycan_module_stats$plots,
                                                    file_path = file_out(here(figdir, "WGCNA_pars/glycans/overview.pdf")),
                                                    pointsize = 2,
                                                    width = 10,
                                                    height=10,
                                                    network_pars = glycan_network_pars),
  
  network_module_heatmaps_noderiv = make_module_heatmap_pdf(plots = glycan_module_stats_noderiv$plots,
                                                    file_path = file_out(here(figdir, "WGCNA_pars/glycans/overview_noderiv.pdf")),
                                                    pointsize = 2,
                                                    width = 10,
                                                    height=10,
                                                    network_pars = glycan_network_pars),
  
  modularity_overview = write_csv(glycan_module_stats$modularity_overview, path = file_out(here(resultdir, "modularity_overview.csv"))),
  modularity_overview_noderiv = write_csv(glycan_module_stats_noderiv$modularity_overview, path = file_out(here(resultdir, "modularity_overview_noderiv.csv"))),
  
  # generate report -------------------------------------------
  
  # rmd_clean_name = clean_rmd_for_render(rmd = file_in("explore_univariate_igx.Rmd"), out_name = file_out("explore_univariate_igx_clean.Rmd"), "View\\("),
  
  report = target(
    command = rmarkdown::render(
      knitr_in("explore_univariate_igx.Rmd"),
      output_file = file_out(here::here("results/report.html")),
      quiet = TRUE
    ),
    trigger = trigger(change = digest(file.info(file_in("results/figures/univar_qqplot.png"))))  # problem: seems to interfere with automatic dependency detection in Rmd
  )
)
