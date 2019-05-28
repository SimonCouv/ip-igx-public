plan = drake_plan(
  
  # load data -------------------------------------------
  glycans = fread(file_in("input_data/glycans_corrected_20190409.csv")),
  glycans_raw = fread(file_in("input_data/glycans_raw.csv")),  
  # IPs
  ips = fread(file_in("input_data/immunopheno.corrected.csv")),
  ips_raw = fread(file_in("input_data/immunopheno.raw.csv")),
  #Popante output
  ip_igx_univar = fread(file_in("input_data/glycans_20190409_immunopheno_corrected_cleanedSimon.tsv"), header = TRUE) %>%
    set_names(names(.) %>% tolower(.) %>% str_replace_all(., c("\\^"="", "[\\.]+"="_"))),
  # annotations
  ip_anno = fread(file_in("input_data/data_annotations/all_immunophenotypes_annotation_av.csv")) %>%
    mutate_all(na_if, y="") %>%   # set all empty string fields to NA 
    set_names(names(.) %>% tolower(.) %>% str_replace_all(., c( "\\.$"="", "[\\.]+"="_", " "="_"))) %>%
    mutate(composite_lin_source = ifelse((source=="Lin"|source=="MFI"), source, lineage)),
  ip_anno_dt = as.data.table(ip_anno),
  IgA_anno_names_raw = readxl::read_xlsx(file_in("input_data/data_annotations/IgA_explanatory_overview.xlsx"), sheet = "raw_names"),
  IgA_anno_names_derived = readxl::read_xlsx(file_in("input_data/data_annotations/IgA_explanatory_overview.xlsx"), sheet = "derived_names"),
  twin_fam = fread(file_in("input_data/data_annotations/TwinDetails_110119.csv")),
  ip_batch = fread(file_in("input_data/immuno_poppante_AGE_BATCH.covar")),
  
  # pre-process (omfo) -------------------------------------------
  
  # intersect 
  
  overlap_samples_raw = intersect(glycans_raw$IID, ips_raw$IID),
  overlap_samples = intersect(glycans$IID, ips$IID),
  
  glycans_o = glycans[glycans$IID %in% overlap_samples,],
  glycans_raw_o = glycans_raw[glycans_raw$IID %in% overlap_samples_raw,],
  
  ips_o = ips[ips$IID %in% overlap_samples,],
  ips_raw_o = ips_raw[ips_raw$IID %in% overlap_samples_raw],
  
  # cap CSF to 100%, and set P3 CSFs for IID 53612 to NA
  ips_o_capped = cap_ips(ips_o, ip_anno_dt),
  ips_raw_o_capped = cap_ips(ips_raw_o, ip_anno_dt),
  
  # glycan missingness filter
  
  glycans_omf = glycans_o[, colMeans(is.na(glycans_o)) < 0.2, with=FALSE] %>%
    .[rowMeans(is.na(.)) < 0.6, ],
  glycans_raw_omf = glycans_raw_o[, colMeans(is.na(glycans_raw_o)) < 0.2, with=FALSE] %>%
    .[rowMeans(is.na(.)) < 0.6, ],
  
  # IP missingness
  ip_missing = ips_o_capped[, lapply(.SD, function(x) mean(is.na(x))), .SDcols=-(FID:IID)] %>%
    data.table::transpose(.) %>% 
    .[,.(set_name = names(ips_o_capped)[-(1:2)],
         missing_prop = V1)
      ] %>% 
    .[order(-missing_prop)],
  
  ip_raw_missing = ips_raw_o_capped[, lapply(.SD, function(x) mean(is.na(x))), .SDcols=-(FID:IID)] %>% 
    data.table::transpose(.) %>% 
    .[,.(set_name = names(ips_raw_o_capped)[-(1:2)],
         missing_prop = V1)
      ] %>% 
    .[order(-missing_prop)],
  
  # IP: filter for zero-variance, col(IP)-wise missingness, final_qc_max, old_p6, then row(sample)-wise missingness
  ips_omf = filter_ips(
    ips=ips_o_capped,
    ip_anno_dt=ip_anno_dt,
    row_miss_threshold = 0.6,
    col_miss_threshold = 0.2
  ),
  ips_raw_omf = filter_ips(
    ips=ips_raw_o_capped,
    ip_anno_dt=ip_anno_dt,
    row_miss_threshold = 0.6,
    col_miss_threshold = 0.2
  ),
  
  # intersect filtered datasets
  
  overlap_omf_samples_raw = intersect(glycans_raw_omf$IID, ips_raw_omf$IID),
  overlap_omf_samples = intersect(glycans_omf$IID, ips_omf$IID),
  
  glycans_omfo = glycans_omf[glycans_omf$IID %in% overlap_omf_samples,],
  glycans_raw_omfo = glycans_raw_omf[glycans_raw_omf$IID %in% overlap_omf_samples_raw,],
  
  ips_omfo = ips_omf[ips_omf$IID %in% overlap_omf_samples,],
  ips_raw_omfo = ips_raw_omf[ips_raw_omf$IID %in% overlap_omf_samples_raw],
  
  #  transform features -------------------------------------------
  #  (min-max) -> arcsin(sqrt) -> scale -> QN
  
  glycans_raw_omfo_trans =  cbind(
    glycans_raw_omfo[,1:5],
    transform_features(glycans_raw_omfo[,-(1:5)]/100)
  ) %>% set_names(names(glycans_raw_omfo)),
  
  SPEL_omfo = ip_anno_dt[(set_name %in% names(ips_raw_omfo)) & source == "MFI",
                         set_name],
  ips_raw_omfo_mm = ips_raw_omfo[,-(1:2)][    # min_max_scale the SPEL columns only
    ,(SPEL_omfo) := lapply(.SD, min_max_scale), .SDcols = SPEL_omfo
    ],
  ips_raw_omfo_trans = cbind(
    ips_raw_omfo[,1:2],
    transform_features(ips_raw_omfo_mm)
  ) %>% set_names(names(ips_raw_omfo)),
  
  # Relatedness adjustment -------------------------------------------
  
  # glycans
  glycans_fam_adj = lmer_adjust(
    raw = glycans_raw_omfo,
    factor_cols = c("FID","IID","Plate_NO","Sex"),
    covar_pos = 1:5,
    form=feat ~ (1|FID) + Plate_NO + Age
  ),
  
  glycans_qn_fam_adj = lmer_adjust(
    raw = glycans_raw_omfo_qn,
    factor_cols = c("FID","IID", "Sex"),
    covar_pos = 1:5,
    form=feat ~ (1|FID) + Age # don't adjust for plate, as no longer strong batch effect after QN (see gPCA)
  ),  
  
  glycan_residuals_pdf = make_marginal_plots(glycans_fam_adj, fig_path = file_out("results/figures/glycan_res_marginal_plots.pdf")),             # diagnostic of LME models
  
  # test within-target parallelism
  par_dummy = call_future(),
  
  
  # data exploration -------------------------------------------
  # IP mean per-feature relative counts
  ip_raw_mean = ips_raw_omfo[, lapply(.SD, mean, na.rm=TRUE), .SDcols = -(FID:IID)] %>%
    data.table::transpose(.) %>%
    .[,.(set_name = colnames(ips_raw_omfo)[-(1:2)],
         mean_signal = V1)
      ],
  ip_mean = ips_omfo[, lapply(.SD, mean, na.rm=TRUE), .SDcols = -(FID:IID)] %>%
    data.table::transpose(.) %>%
    .[,.(set_name = colnames(ips_omfo)[-(1:2)],
         mean_signal = V1)
      ],
   
  # glycan marginal distributions pdf
  raw_glycan_marginal_pdf = make_marginal_plots(glycans_raw_omf[,-(1:5)], fig_path = file_out("results/figures/raw_glycan_marginal_plots.pdf")),  # important as input to lmer models
  glycan_marginal_pdf = make_marginal_plots(glycans_omf[,-(1:5)], fig_path = file_out("results/figures/glycan_marginal_plots.pdf")),              # important as input to a.o. association scatters
  
  # IP marginal distributions pdf
  raw_ip_marginal_pdf = make_marginal_plots(ips_raw_omf[,3:100], fig_path = file_out("results/figures/raw_ips_marginal_plots_top100.pdf")),              # important as input to a.o. association scatters
  ip_marginal_pdf = make_marginal_plots(ips_omf[,3:100], fig_path = file_out("results/figures/ips_marginal_plots_top100.pdf")),              # important as input to a.o. association scatters
  
  # batch effects - glycans -------------------------------------------
  # batch effects: riPCA-imputed glycans
  # find optimal number of PCs to retain in the imputation, using CV, on scaled data
  glycans_scaled_ncomp = future_map(.x = list(
    "raw" = scale(glycans_raw_omfo[,-(1:5)]),
    "corr" = scale(glycans_omfo[,-(1:5)]), 
    "qn" = glycans_raw_omfo_qn[,-(1:5)],
    "famadj" = scale(glycans_fam_adj),
    "qn_famadj" = glycans_qn_fam_adj),         # scale input
  .f = estim_ncpPCA,
  ncp.max = 15, 
  ncp.min = 0, 
  scale=FALSE,  #because inputs are already scaled
  verbose=TRUE),
  
  # impute using optimal number of PCs, on scaled data
  glycans_scaled_riPCA_imp = future_map2(
    .x = map(
      .x =  list(
        "raw" = scale(glycans_raw_omfo[,-(1:5)]),
        "corr" = scale(glycans_omfo[,-(1:5)]), 
        "qn" = glycans_raw_omfo_qn[,-(1:5)],
        "famadj" = scale(glycans_fam_adj),
        "qn_famadj" = glycans_qn_fam_adj),
      ~as.data.frame(scale(.x))),                    #scale input IF not already scaled
    .y = map(glycans_scaled_ncomp, "ncp"),
    .f = function(x,y,...){imputePCA(X=x,ncp=y,...)},
    scale=FALSE,  #because inputs are already scaled (required because FactoMineR::PCA does not provide internal scaling)
    method="Regularized"
  ),
  
  # mean impute
  glycans_mean_imp = map(
    .x = list(
      "raw" = glycans_raw_omfo[,-(1:5)],
      "corr" = glycans_omfo[,-(1:5)], 
      "qn" = glycans_raw_omfo_qn[,-(1:5)],
      "famadj" = glycans_fam_adj,
      "qn_famadj" = glycans_qn_fam_adj),
    .f = ~map(.x, function(x){x[is.na(x)]<- mean(x, na.rm=TRUE); x}) %>% 
      as_tibble()
  ),
  
  # compare feature-wise distributions between input and imputed data points
  p_riPCA_input_impute_comparison = compare_input_imputed(
    raw_df_l = list(
      "raw" = scale(glycans_raw_omfo[,-(1:5)]),
      "corr" = scale(glycans_omfo[,-(1:5)]), 
      "qn" = glycans_raw_omfo_qn[,-(1:5)],
      "famadj" = scale(glycans_fam_adj),
      "qn_famadj" = glycans_qn_fam_adj),
    imputed_df_l = map(glycans_scaled_riPCA_imp, "completeObs")
  ),
  
  p_mean_input_impute_comparison = compare_input_imputed(
    raw_df_l = list("raw" = glycans_raw_omfo[,-(1:5)],
              "corr"= glycans_omfo[,-(1:5)],
              "qn"= glycans_raw_omfo_qn[,-(1:5)],
              "famadj"= glycans_fam_adj,
              "qn_famadj" =  glycans_qn_fam_adj),
    imputed_df_l = glycans_mean_imp
  ),
  
  # PCA analysis, for PC boxplots per plate
  # standard PCA on riPCA-imputed datasets, on scaled data
  glycans_riPCAimp_pca = future_map2(
    .x = map2(
      .x = glycans_scaled_riPCA_imp,
      .y = list(glycans_raw_omfo,
                glycans_omfo,
                glycans_raw_omfo_qn,
                glycans_raw_omfo,
                glycans_raw_omfo),
      .f = function(x,y){cbind(y[, 'Plate_NO'], x$completeObs)}
    ),
    .y = map(glycans_scaled_ncomp, "ncp"),
    .f =  function(x,y,...){FactoMineR::PCA(X = x, ncp = y, ...)},
    quali.sup = 1,
    graph = FALSE
  ),
  
  # standard PCA on mean-imputed datasets, with prior scaling
  glycans_meanimp_pca= future_map(
    .x = map2(
      .x = map(glycans_mean_imp, ~scale(.x)),
      .y = list(glycans_raw_omfo,
                glycans_omfo,
                glycans_raw_omfo_qn,
                glycans_raw_omfo,
                glycans_raw_omfo),
      .f = function(x,y){cbind(y[, 'Plate_NO'], x)}
    ),
    .f =  function(x,y,...){FactoMineR::PCA(X = x, ...)},
    quali.sup = 1,
    graph = FALSE,
    ncp=15
  ),
  
  # formal test (Reese2013)
  glycans_gPCA_meanimp = future_map2(
    .x = map(glycans_mean_imp, ~scale(.x)),  #scale inputs to give same weight
    .y = list(glycans_raw_omfo,
              glycans_omfo,
              glycans_raw_omfo_qn,
              glycans_raw_omfo,
              glycans_raw_omfo),
    .f = function(x,y,...){
      gPCA.batchdetect(x=x, batch=y$Plate_NO)
    },
    center = TRUE,  # data is centered [counter-intuitive parameter name, see ?gPCA.batchdetect]
    scaleY = TRUE,  # adjust for unequal batch sizes
    nperm=20000,
    seed=5
  ),
  # formal test (Reese2013)
  glycans_gPCA_riPCAimp = future_map2(
    .x = glycans_scaled_riPCA_imp,  #already scaled prior to imputation
    .y = list(glycans_raw_omfo,
              glycans_omfo,
              glycans_raw_omfo_qn,
              glycans_raw_omfo,
              glycans_raw_omfo),
    .f = function(x,y,...){
      gPCA.batchdetect(x=x$completeObs, batch=y$Plate_NO)
    },
    center = TRUE,  # data is centered [counter-intuitive parameter name, see ?gPCA.batchdetect]
    scaleY = TRUE,  # adjust for unequal batch sizes
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
  
  glycan_network_pars = read_tsv(file_in("WGCNA_parameters/glycan_pars_qn_famadj.tsv")) %>%
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
  
  glycan_modularity_overview = write_csv(glycan_module_stats$modularity_overview, path = file_out("results/modularity_overview.csv")),
  glycan_modularity_overview_noderiv = write_csv(glycan_module_stats_noderiv$modularity_overview, path = file_out("results/modularity_overview_noderiv.csv")),
  
  # Immunophenotypes
  # ips_scale_free = optimize_network_pars(data=glycans_qn_fam_adj, common_pars = scale_free_common_pars, pars = scale_free_iter_pars, powers = scale_free_powers, plt_title = "glycans"),
  
  
  
  # generate report -------------------------------------------
  
  # rmd_clean_name = clean_rmd_for_render(rmd = file_in("explore_univariate_igx.Rmd"), out_name = file_out("explore_univariate_igx_clean.Rmd"), "View\\("),
  
  # report = target(
  #   command = rmarkdown::render(
  #     knitr_in("EDA.Rmd"),
  #     output_file = file_out(here::here("results/report.html")),
  #     quiet = TRUE
  #   ),
  #   trigger = trigger(change = digest(file.info(file_in("results/figures/univar_qqplot.png"))))  # problem: seems to interfere with automatic dependency detection in Rmd
  # )
)
