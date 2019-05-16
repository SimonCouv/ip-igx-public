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
  ip_batch = fread(here("/input_data/immuno_poppante_AGE_BATCH.covar")),
  
  # pre-process -------------------------------------------
  
  # intersect 
  
  overlap_samples_raw = intersect(glycans_raw$IID, ips_raw$IID),
  overlap_samples = intersect(glycans$IID, ips$IID),
  
  glycans_o = glycans[glycans$IID %in% overlap_samples,],
  glycans_raw_o = glycans_raw[glycans_raw$IID %in% overlap_samples_raw,],
  
  ips_o = ips[ips$IID %in% overlap_samples,],
  ips_raw_o = ips_raw[ips_raw$IID %in% overlap_samples_raw],
  
  # glycan missingness filter
  
  glycans_omf = glycans_o[, colMeans(is.na(glycans_o)) < 0.2, with=FALSE] %>% .[rowMeans(is.na(.)) < 0.6, ],
  glycans_raw_omf = glycans_raw_o[, colMeans(is.na(glycans_raw_o)) < 0.2, with=FALSE] %>% .[rowMeans(is.na(.)) < 0.6, ],
  
  # IP missingness
  ip_missing = ips_o[, lapply(.SD, function(x) mean(is.na(x))), .SDcols=-(FID:IID)] %>%
    data.table::transpose(.) %>% 
    .[,.(set_name = names(ips_o)[-(1:2)],
         missing_prop = V1)
      ] %>% 
    .[order(-missing_prop)],
  
  ip_raw_missing = ips_raw_o[, lapply(.SD, function(x) mean(is.na(x))), .SDcols=-(FID:IID)] %>% 
    data.table::transpose(.) %>% 
    .[,.(set_name = names(ips_raw_o)[-(1:2)],
         missing_prop = V1)
      ] %>% 
    .[order(-missing_prop)],
  
  # IP missingness filter
  ips_omf = ips_o[, colMeans(is.na(ips_o)) < 0.2, with=FALSE] %>% .[rowMeans(is.na(.)) < 0.6, ],
  ips_raw_omf = ips_raw_o[, colMeans(is.na(ips_raw_o)) < 0.2, with=FALSE] %>% .[rowMeans(is.na(.)) < 0.6, ],
  
  # IP mean per-feature relative counts
  ip_raw_mean = ips_raw_o[, lapply(.SD, mean, na.rm=TRUE), .SDcols = -(FID:IID)] %>%
    data.table::transpose(.) %>%
    .[,.(set_name = colnames(ips_raw_o)[-(1:2)],
         mean_signal = V1)
      ],
  ip_mean = ips_o[, lapply(.SD, mean, na.rm=TRUE), .SDcols = -(FID:IID)] %>%
    data.table::transpose(.) %>%
    .[,.(set_name = colnames(ips_o)[-(1:2)],
         mean_signal = V1)
      ],
  
  # intersect filtered datasets
  
  overlap_omf_samples_raw = intersect(glycans_raw_omf$IID, ips_raw_omf$IID),
  overlap_omf_samples = intersect(glycans_omf$IID, ips_omf$IID),
  
  glycans_omfo = glycans_omf[glycans_omf$IID %in% overlap_omf_samples,],
  glycans_raw_omfo = glycans_raw_omf[glycans_raw_omf$IID %in% overlap_omf_samples_raw,],
  
  ips_omfo = ips_omf[ips_omf$IID %in% overlap_omf_samples,],
  ips_raw_omfo = ips_raw_omf[ips_raw_omf$IID %in% overlap_omf_samples_raw],
  
  # quantile normalize
  glycans_raw_omfo_qn = cbind(
    glycans_raw_omfo[,1:5],
    normalize.quantiles(as.matrix(glycans_raw_omfo[,-(1:5)]))
  ) %>% set_names(names(glycans_raw_omfo)),
  
  
  # data exploration -------------------------------------------
  # glycan marginal distributions pdf
  raw_glycan_marginal_pdf = make_marginal_plots(glycans_raw_omf[,-(1:5)], fig_path = file_out("results/figures/raw_glycan_marginal_plots.pdf")),  # important as input to lmer models
  glycan_marginal_pdf = make_marginal_plots(glycans_omf[,-(1:5)], fig_path = file_out("results/figures/glycan_marginal_plots.pdf")),              # important as input to a.o. association scatters
  
  # IP marginal distributions pdf
  raw_ip_marginal_pdf = make_marginal_plots(ips_raw_omf[,3:100], fig_path = file_out("results/figures/raw_ips_marginal_plots_top100.pdf")),              # important as input to a.o. association scatters
  ip_marginal_pdf = make_marginal_plots(ips_omf[,3:100], fig_path = file_out("results/figures/ips_marginal_plots_top100.pdf")),              # important as input to a.o. association scatters
  
  # batch effects -------------------------------------------
  # batch effects: riPCA-imputed glycans
  # find optimal number of PCs to retain in the imputation, using CV
  glycans_ncomp = future_map(.x = list("raw" = glycans_raw_omfo[,-(1:5)],
                                "corr" = glycans_omfo[,-(1:5)], 
                                "qn" = glycans_raw_omfo_qn[,-(1:5)],
                                "famadj" = glycans_fam_adj,
                                "qn_famadj" = glycans_qn_fam_adj),
                      .f = estim_ncpPCA,
                      ncp.max = 15, 
                      ncp.min = 0, 
                      scale=TRUE,
                      verbose=TRUE),
  
  # impute using optimal number of PCs
  glycans_imp = future_map2(
    .x = map(
      .x = list("raw" = glycans_raw_omfo[,-(1:5)],
           "corr" = glycans_omfo[,-(1:5)], 
           "qn" = glycans_raw_omfo_qn[,-(1:5)],
           "famadj" = glycans_fam_adj,
           "qn_famadj" = glycans_qn_fam_adj),
      ~as.data.frame(.x)),
    .y = map(glycans_ncomp, "ncp"),
    .f = function(x,y,...){imputePCA(X=x,ncp=y,...)},
    scale=TRUE,
    method="Regularized"
  ),
  
  # standard PCA on imputed datasets, with internal scaling
  glycans_pca = future_map2(
    .x = map2(
      .x = glycans_imp,
      .y = list(glycans_raw_omfo,
              glycans_omfo,
              glycans_raw_omfo_qn,
              glycans_raw_omfo,
              glycans_raw_omfo),
      .f = function(x,y){cbind(y[, 'Plate_NO'], x$completeObs)}
    ),
    .y = map(glycans_ncomp, "ncp"),
    .f =  function(x,y,...){FactoMineR::PCA(X = x, ncp = y, ...)},
    quali.sup = 1,
    graph = FALSE
  ),
  
  # formal test (Reese2013)
  glycans_gPCA = future_map2(
    .x = glycans_imp,
    .y = list(glycans_raw_omfo,
              glycans_omfo,
              glycans_raw_omfo_qn,
              glycans_raw_omfo,
              glycans_raw_omfo),
    .f = function(x,y,...){
      gPCA.batchdetect(x=x$completeObs, 
                       batch=y$Plate_NO)
    },
    center = FALSE,
    scaleY = TRUE,
    nperm=20000,
    seed=5
  ),
  
  # Univariate analysis  -------------------------------------------
  
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
  
  #qq plot
  qq = make_univar_qq(v_pval = ip_igx_univar_good$pvalue, fig_path = file_out("results/figures/univar_qqplot.png")),
  
  # top association scatter plots pdf
  top_assoc_scatter_pdf = make_top_assoc_scatters(ntop = 100, ip_igx_univar_good, glycans_omf, ips_omf, file_path = file_out("results/figures/top_assoc_scatters.pdf")),
  
  # number of significant IgG, IgA and IP ~ FDR threshold
  n_associated_m = get_fdr_nfeat_associated_univar(fdr_thresholds = seq(0,1,0.01), ip_igx_univar_good = ip_igx_univar_good),
  
  # enrichment analyses:
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
  
  glycans_fam_adj = glycans_lmer_adjust(glycans_raw_omfo, form=glycan ~ (1|FID) + Plate_NO + Age),
  glycans_qn_fam_adj = glycans_lmer_adjust(glycans_raw_omfo_qn, form=glycan ~ (1|FID) + Age),  # don't adjust for plate, as no longer strong batch effect after QN (see gPCA)
  glycan_residuals_pdf = make_marginal_plots(glycans_fam_adj, fig_path = file_out("results/figures/glycan_res_marginal_plots.pdf")),             # diagnostic of LME models
  
  
  # WGCNA  -------------------------------------------
  
  # 1. scale-free topology: power parameter
  
  scale_free_powers = c(c(1:10), seq(from = 12, to=20, by=2)),
  scale_free_common_pars = list(verbose= 0,  moreNetworkConcepts = TRUE, powerVector = scale_free_powers),
  scale_free_iter_pars = expand.grid(networktype = c('signed', 'signed hybrid', 'unsigned'),
                                     corFnc = c('cor', 'bicor'),
                                     stringsAsFactors = FALSE),
  derived_glycan_pos = str_detect(names(glycans_fam_adj), pattern = "_"),
  
  glycans_scale_free = optimize_network_pars(data=glycans_qn_fam_adj, common_pars = scale_free_common_pars, pars = scale_free_iter_pars, powers = scale_free_powers, plt_title = "glycans"),
  glycans_scale_free_noderiv = optimize_network_pars(data=glycans_qn_fam_adj[, -!derived_glycan_pos], scale_free_common_pars, scale_free_iter_pars, scale_free_powers, plt_title = "glycans"),
  
  # 2. parameter search ~ modularity
  
  glycan_cross_iter_pars = list(figdir = "results/figures/WGCNA_pars/glycans/"),
  
  glycan_network_pars = read_tsv("WGCNA_parameters/glycan_pars.tsv") %>%
    dplyr::select(-aim),
  glycan_cut_pars = expand.grid(deep =c(4, 1),
                                minModuleSize=c(3,5,10,20)),
  glycan_modality = ifelse(str_detect(names(glycans_fam_adj), pattern = "IgG"),"IgG","IgA"),
  
  glycan_module_stats = wgcna_parameter_search(data=glycans_qn_fam_adj,
                                               glycan_network_pars, 
                                               glycan_cut_pars,
                                               glycan_cross_iter_pars, 
                                               create_plot = TRUE,
                                               feat_modality=glycan_modality),
  glycan_module_stats_noderiv = wgcna_parameter_search(data=glycans_qn_fam_adj[,!derived_glycan_pos],
                                                       glycan_network_pars, 
                                                       glycan_cut_pars,
                                                       glycan_cross_iter_pars, 
                                                       create_plot = TRUE,
                                                       feat_modality=glycan_modality[!derived_glycan_pos]),
  
  network_module_heatmaps = make_module_heatmap_pdf(plots = glycan_module_stats$plots,
                                                    file_path = file_out("results/figures/WGCNA_pars/glycans/overview.pdf"),
                                                    pointsize = 2,
                                                    width = 10,
                                                    height=10,
                                                    network_pars = glycan_network_pars),
  
  network_module_heatmaps_noderiv = make_module_heatmap_pdf(plots = glycan_module_stats_noderiv$plots,
                                                            file_path = file_out("results/figures/WGCNA_pars/glycans/overview_noderiv.pdf"),
                                                            pointsize = 2,
                                                            width = 10,
                                                            height=10,
                                                            network_pars = glycan_network_pars),
  
  modularity_overview = write_csv(glycan_module_stats$modularity_overview, path = file_out("results/modularity_overview.csv")),
  modularity_overview_noderiv = write_csv(glycan_module_stats_noderiv$modularity_overview, path = file_out("results/modularity_overview_noderiv.csv")),
  
  # generate report -------------------------------------------
  
  # rmd_clean_name = clean_rmd_for_render(rmd = file_in("explore_univariate_igx.Rmd"), out_name = file_out("explore_univariate_igx_clean.Rmd"), "View\\("),
  
  report = target(
    command = rmarkdown::render(
      knitr_in("EDA.Rmd"),
      output_file = file_out(here::here("results/report.html")),
      quiet = TRUE
    ),
    trigger = trigger(change = digest(file.info(file_in("results/figures/univar_qqplot.png"))))  # problem: seems to interfere with automatic dependency detection in Rmd
  )
)
