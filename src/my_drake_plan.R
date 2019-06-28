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
    mutate(lineage=recode(lineage, `04-Aug`="4-8"),
           composite_lin_source = ifelse((source=="Lin"|source=="MFI"), source, lineage)),
  ip_anno_dt = as.data.table(ip_anno),
  IgA_anno_names_raw = readxl::read_xlsx(file_in("input_data/data_annotations/IgA_explanatory_overview.xlsx"), sheet = "raw_names"),
  IgA_anno_names_derived = readxl::read_xlsx(file_in("input_data/data_annotations/IgA_explanatory_overview.xlsx"), sheet = "derived_names"),
  twin_fam = fread(file_in("input_data/data_annotations/TwinDetails_110119.csv")),
  ip_batch = fread(file_in("input_data/immuno_poppante_AGE_BATCH.covar")),
  
  # pre-process (omfo) -------------------------------------------
  
  # intersect 
  
  overlap_samples_raw = setdiff(intersect(glycans_raw$IID, ips_raw$IID), 
                                c(99611, 99612)),   # see EDA.Rmd, Section batch effects -> IPs -> Experimental design
  overlap_samples = setdiff(intersect(glycans$IID, ips$IID),
                            c(99611, 99612)),   # see EDA.Rmd, Section batch effects -> IPs -> Experimental design
  
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
  
  #for descriptive stats
  ip_batch_omfo = ip_batch[ip_batch$IID %in% overlap_omf_samples_raw,],
  
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
  glycans_trans_fam_adj = lmer_adjust(
    raw = glycans_raw_omfo_trans,
    factor_cols = c("FID","IID", "Sex"),
    covar_pos = 1:5,
    form=feat ~ (1|FID) + Age # don't adjust for plate, as no longer strong batch effect after QN (see gPCA)
  ),

  ips_trans_fam_adj_l = lmer_adjust_list(
    raw = ip_batch[ips_raw_omfo_trans, on=.(FID, IID),nomatch=0],
    factor_cols = c("batch", "FID", "IID"),
    covar_pos = 1:4,
    form=feat ~ (1|FID) + age),  # don't adjust for plate, as no longer strong batch effect after QN (see gPCA)
  
  ips_trans_fam_adj = fam_adj_collect(
    fam_adj_l=ips_trans_fam_adj_l,
    raw=ip_batch[ips_raw_omfo_trans, on=.(FID, IID),nomatch=0]),
  
  ips_trans_fam_adj_tsv = write_tsv(
    ips_trans_fam_adj,
    path = file_out("data/ips_trans_fam_adj.tsv"),
    col_names = TRUE),
  
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
    "qn" = glycans_raw_omfo_trans[,-(1:5)],
    # "famadj" = scale(glycans_fam_adj),
    "qn_famadj" = glycans_trans_fam_adj),         # scale input
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
        "qn" = glycans_raw_omfo_trans[,-(1:5)],
        # "famadj" = scale(glycans_fam_adj),
        "qn_famadj" = glycans_trans_fam_adj),
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
      "qn" = glycans_raw_omfo_trans[,-(1:5)],
      # "famadj" = glycans_fam_adj,
      "qn_famadj" = glycans_trans_fam_adj),
    .f = ~map(.x, function(x){x[is.na(x)]<- mean(x, na.rm=TRUE); x}) %>% 
      as_tibble()
  ),
  
  # compare feature-wise distributions between input and imputed data points
  p_riPCA_input_impute_comparison = compare_input_imputed(
    raw_df_l = list(
      "raw" = scale(glycans_raw_omfo[,-(1:5)]),
      "corr" = scale(glycans_omfo[,-(1:5)]), 
      "qn" = glycans_raw_omfo_trans[,-(1:5)],
      # "famadj" = scale(glycans_fam_adj),
      "qn_famadj" = glycans_trans_fam_adj),
    imputed_df_l = map(glycans_scaled_riPCA_imp, "completeObs")
  ),
  
  p_mean_input_impute_comparison = compare_input_imputed(
    raw_df_l = list("raw" = glycans_raw_omfo[,-(1:5)],
              "corr"= glycans_omfo[,-(1:5)],
              "qn"= glycans_raw_omfo_trans[,-(1:5)],
              # "famadj"= glycans_fam_adj,
              "qn_famadj" =  glycans_trans_fam_adj),
    imputed_df_l = glycans_mean_imp
  ),
  
  # PCA analysis, for PC boxplots per plate
  # standard PCA on riPCA-imputed datasets, on scaled data
  glycans_riPCAimp_pca = future_map2(
    .x = map2(
      .x = glycans_scaled_riPCA_imp,
      .y = list(glycans_raw_omfo,
                glycans_omfo,
                glycans_raw_omfo_trans,
                # glycans_raw_omfo,
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
                glycans_raw_omfo_trans,
                # glycans_raw_omfo,
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
              glycans_raw_omfo_trans,
              # glycans_raw_omfo,
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
              glycans_raw_omfo_trans,
              # glycans_raw_omfo,
              glycans_raw_omfo),
    .f = function(x,y,...){
      gPCA.batchdetect(x=x$completeObs, batch=y$Plate_NO)
    },
    center = TRUE,  # data is centered [counter-intuitive parameter name, see ?gPCA.batchdetect]
    scaleY = TRUE,  # adjust for unequal batch sizes
    nperm=20000,
    seed=5
  ),
  
  
  # batch effects: IPs ---------------------------------------------------------
  # only use mean-impute

  # mean-impute 
  # ips_scaled_imp = future_map(
  #   .x =  list("raw" = ips_raw_omfo[,-(1:2)],
  #          "corr" = ips_omfo[,-(1:2)],
  #          "qn" = ips_raw_omfo_qn[,-(1:2)],
  #          "famadj" = ips_fam_adj,
  #          "qn_famadj" = ips_qn_fam_adj),
  #   .f = ~map(.x, function(x){x[is.na(x)]<- mean(x, na.rm=TRUE); x}) %>% 
  #     as_tibble()
  # ),
  
  # # impute using optimal number of PCs, without prior scaling
  # glycans_imp = future_map2(
  #   .x = map(
  #     .x = list("raw" = glycans_raw_omfo[,-(1:5)],
  #          "corr" = glycans_omfo[,-(1:5)], 
  #          "qn" = glycans_raw_omfo_qn[,-(1:5)],
  #          "famadj" = glycans_fam_adj,
  #          "qn_famadj" = glycans_trans_fam_adj),
  #     ~as.data.frame(.x)),                    # do NOT scale input
  #   .y = map(glycans_scaled_ncomp, "ncp"),
  #   .f = function(x,y,...){imputePCA(X=x,ncp=y,...)},
  #   scale=TRUE,  # use internal scaling
  #   method="Regularized"
  # ),
  # 
  # # standard PCA on imputed datasets, on scaled data
  # glycans_scaled_pca = future_map2(
  #   .x = map2(
  #     .x = glycans_scaled_riPCA_imp,
  #     .y = list(glycans_raw_omfo,
  #             glycans_omfo,
  #             glycans_raw_omfo_qn,
  #             glycans_raw_omfo,
  #             glycans_raw_omfo),
  #     .f = function(x,y){cbind(y[, 'Plate_NO'], x$completeObs)}
  #   ),
  #   .y = map(glycans_scaled_ncomp, "ncp"),
  #   .f =  function(x,y,...){FactoMineR::PCA(X = x, ncp = y, ...)},
  #   quali.sup = 1,
  #   graph = FALSE
  # ),
  # 
  # # formal test (Reese2013)
  # glycans_gPCA = future_map2(
  #   .x = glycans_scaled_riPCA_imp,
  #   .y = list(glycans_raw_omfo,
  #             glycans_omfo,
  #             glycans_raw_omfo_qn,
  #             glycans_raw_omfo,
  #             glycans_raw_omfo),
  #   .f = function(x,y,...){
  #     gPCA.batchdetect(x=x$completeObs, 
  #                      batch=y$Plate_NO)
  #   },
  #   center = TRUE,  # data is centered [counter-intuitive parameter name, see ?gPCA.batchdetect]
  #   scaleY = TRUE,  # adjust for unequal batch sizes
  #   nperm=20000,
  #   seed=5
  # ),
  

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
  
  
  
  
  # WGCNA glycans -------------------------------------------
  
  # 1. scale-free topology: power parameter
  
  scale_free_powers = c(c(1:10), seq(from = 12, to=20, by=2)),
  scale_free_common_pars = list(verbose= 0,  moreNetworkConcepts = TRUE, powerVector = scale_free_powers),
  scale_free_iter_pars = expand.grid(networktype = c('signed', 'signed hybrid', 'unsigned'),
                                     corFnc = c('cor', 'bicor'),
                                     stringsAsFactors = FALSE) %>% 
    as_tibble() %>% 
    left_join(tibble(corFnc = c('cor', 'cor','bicor'),
                         corOptions = list(list(use="p"),
                                           list(use="p", method="spearman"),
                                           list(use="p"))),
              by="corFnc"),
  
  derived_glycan_pos = str_detect(names(glycans_trans_fam_adj), pattern = "_"),
  glycans_scale_free = optimize_network_pars(data=glycans_trans_fam_adj, common_pars = scale_free_common_pars, pars = scale_free_iter_pars, plt_title = "glycans"),
  glycans_scale_free_noderiv = optimize_network_pars(data=glycans_trans_fam_adj[, -!derived_glycan_pos],common_pars = scale_free_common_pars, pars =scale_free_iter_pars, plt_title = "glycans"),
  
  # 2. parameter search ~ modularity
  
  # parameters 
  
  glycan_network_pars = read_tsv(file_in("WGCNA_parameters/glycan_pars_trans.tsv")) %>%
    dplyr::select(-aim) %>% 
    group_by(corfnc, networktype, method) %>% 
    summarise(powers=list(softpower)),
  glycan_network_pars_noderiv = read_tsv(file_in("WGCNA_parameters/glycan_pars_trans_noderiv.tsv")) %>%
    dplyr::select(-aim) %>% 
    group_by(corfnc, networktype, method) %>% 
    summarise(powers=list(softpower)),
  glycan_modality = ifelse(str_detect(names(glycans_trans_fam_adj), pattern = "IgG"),"IgG","IgA"),
  
  # module-parameter search
  
  # impute for calculation of PC1s
  # 
  # impute in sample space (via t()), so as not to induce strong correlations
  # between features based only on correlation on a small subset of samples 
  # (given that currently up to 20% feature missingness is allowed)
  # 
  # use pre-imputed df, so imputation (in sample space) takes all features into account, 
  # not only the features in the modules (as is the case for internal knn.impute 
  # in WGCNA::moduleEigengenes). This assures that PC1 can always be used, 
  # therefore the interpretation is the same for the representatives of all 
  # modules (NOT: some representatives are PC1s, some are hubGenes)
  glycans_trans_fam_adj_kkn_imp = t(impute::impute.knn(t(as.matrix(glycans_trans_fam_adj)),
                                                   rowmax = 0.6)$data),
  
  p_knn_input_impute_comparison = compare_input_imputed(
    raw_df_l = list(glycans_trans_fam_adj),
    imputed_df_l = list(glycans_trans_fam_adj_kkn_imp)
  ),
  
  glycan_modules = pmap(
    .l = glycan_network_pars,
    .f = wgcna_parameter_search,
    data=glycans_trans_fam_adj,
    cut_pars=expand.grid(deep =c(4, 1),
                         minModuleSize=c(3,5,10,20)),
    create_plot = TRUE,
    feat_modality=glycan_modality
  ),
  
  glycan_modules_df = wrangle_glycan_module_df(
    modules = glycan_modules,
    network_pars = glycan_network_pars,
    dat = glycans_trans_fam_adj_kkn_imp
  ),
  
  glycan_modules_noderiv = pmap(
    .l = glycan_network_pars_noderiv,
    .f = wgcna_parameter_search,
    data=glycans_trans_fam_adj[,!derived_glycan_pos],
    cut_pars=expand.grid(deep =c(4, 1),
                         minModuleSize=c(3,5,10,20)),
    create_plot = TRUE,
    feat_modality=glycan_modality[!derived_glycan_pos]
  ),
  
  glycan_modules_df_noderiv = wrangle_glycan_module_df(
    modules = glycan_modules_noderiv,
    network_pars = glycan_network_pars_noderiv,
    dat = glycans_trans_fam_adj_kkn_imp[,!derived_glycan_pos]
  ),
  
  glycans_mod_varExpl_plots = modularity_varExpl_plots(glycan_modules_df_noderiv, 
                                                        varExpl_trh = 0.5,
                                                        ncol=1),  # for report

  glycan_select_scores = glycan_modules_df_noderiv %>% 
    unnest(cut_stats) %>% 
    unnest(module_stats)  %>% 
    dplyr::filter(
      corfnc=="bicor",
      networktype=="signed_hybrid",
      powers==6,
      cut==3
    ) %>% 
    dplyr::select(colors, score) %>%
    spread(colors,score) %>% 
    map(unlist) %>% 
    as_tibble(),
  
  glycan_select_scores_tsv = write_tsv(
    glycan_select_scores,
    path = file_out("data/glycan_select_scores.tsv")
  ),
  
  glycan_select_module_anno = read_tsv(
    file_in("WGCNA_parameters/glycan_modules_anno_bicor_signed-hybrid_power6_cut3.tsv")
  ),
  

  
  # glycans modules - IPs -------------------------------------------
  glycan_modules_cor_ips = cor_p_bh(
    x = glycan_select_scores,
    y = ips_trans_fam_adj,
    method="pearson",
    use="pairwise.complete.obs"
  ),
  
  # annotate with subgroups of IgX and IP
  glycan_modules_cor_ips_anno = glycan_modules_cor_ips %>%
    left_join(ip_anno %>% 
                dplyr::select(set_name, subset_name, composite_lin_source), 
              by=c("IP"="set_name")) %>% 
    left_join(glycan_select_module_anno, by="glycan_module") %>% 
    dplyr::select(module_tag,n, igx, site, chain, subset_name, 
                  composite_lin_source, correlation, p.val, p.val_bh, p.val_storey, 
                  everything()),
  
  glycan_modules_cor_ips_anno_tsv =  write_tsv(
    glycan_modules_cor_ips_anno,
    path = file_out("data/glycan_modules_cor_ips_anno.tsv")
  ),
  
  #qq plot
  igx_mod_ips_qq = make_univar_qq(v_pval = glycan_modules_cor_ips$p.val, fig_path = file_out("results/figures/igx_mods_ips_qqplot.png")),

  # top association scatter plots pdf
  # top_assoc_scatter_pdf = make_top_assoc_scatters(ntop = 100, glycan_modules_cor_ips_anno, glycan_select_scores, ips_trans_fam_adj, file_path = file_out("results/figures/igx_modules_ips_top_assoc_scatters.pdf")),

  # number of significant IgG, IgA and IP ~ FDR threshold
  igx_mod_ips_n_associated_m = get_fdr_nfeat_associated_univar_igx_mods(
    fdr_thresholds = seq(0,1,0.0005),
    data = glycan_modules_cor_ips_anno,
    p_adj="p.val_storey"
  ),

  # enrichment analyses:
  igx_mod_ips_target_indices = unique(floor(10^seq(0,log10(nrow(glycan_modules_cor_ips)), length.out = 1e4))),

  # cumulative mean number of IgG among top associated IgX
  igx_mod_ips_running_IgG_perc =  cummean(str_detect(glycan_modules_cor_ips$glycan_module, pattern = "G")),

  # enrichment of lin_sources ~ FDR
  # igx_mod_ips_lin_source_props_df = get_lin_source_props_igx_mod_ips(
  #   target_indices=glycan_modules_cor_ips_anno,
  #   data=glycan_modules_cor_ips_anno,
  #   ip_anno=ip_anno),
  
  # enrichment of lin_sources ~ FDR+IgX
  # igx_mod_ips_lin_source_IgX_props_df = get_lin_source_props_byIgX_igx_mod_ips(
  #   target_indices=igx_mod_ips_target_indices,
  #   data=glycan_modules_cor_ips_anno,
  #   igx_var = igx),
  
  # WGCNA IPs -------------------------------------------
  
  #from DTR cluster
  ip_module_stats = readRDS(file_in("data/ip_module_stats_combined.RDS")),
  
  ips_trans_fam_adj_knn_imp = t(impute::impute.knn(t(as.matrix(ips_trans_fam_adj)),
                                                   rowmax = 0.6)$data),
  
  ip_modules_df = wrangle_ip_module_df(
    ip_module_stats_df = ip_module_stats,
    dat = ips_trans_fam_adj_knn_imp
  ),
  
  IP_mod_varExpl_plots = modularity_varExpl_plots(ip_modules_df, 
                                                  varExpl_trh = 0.5,
                                                  ncol=1),
  
  ip_select_modules = ip_modules_df %>% 
    unnest(cut_stats) %>% 
    unnest(module_stats)  %>% 
    dplyr::filter(
      corfnc=="spear",
      networktype=="signed_hybrid",
      powers==5,
      cut==8
    ),
  
  ip_select_scores = ip_select_modules %>% 
    dplyr::select(colors, score) %>%
    spread(colors,score) %>% 
    map(unlist) %>% 
    as_tibble(),
  
  ip_select_scores_tsv = write_tsv(
    ip_select_scores,
    path = file_out("data/ip_select_scores.tsv")
  ),
  
  ip_select_PC1_cor_tmp = ip_modules_df %>% 
    unnest(cut_stats) %>% 
    dplyr::filter(
      corfnc=="spear",
      networktype=="signed_hybrid",
      powers==5,
      cut==8
    ) %>% 
    unnest(feature_stats) %>% 
    dplyr::filter(colors!="grey") %>% 
    dplyr::select(feat, colors) %>%
    mutate(
      PC1_cor = calc_PC1_cor(IP=feat, 
                             color=colors,
                             ip_select_scores=ip_select_scores,
                             ips_trans_fam_adj=ips_trans_fam_adj,
                             use="p",
                             method="pearson")
    ) %>% 
    left_join(ip_anno %>% 
                dplyr::select(set_name, subset_name, composite_lin_source), 
              by=c("feat"="set_name")) %>% 
    nest(-colors, .key = "feature_data")%>% 
    dplyr::filter(colors!="grey"),
  
  # split up long calculation as multiple targets, see drake manual, section large plans
  # more efficient than using within-target parallelism
  core_CD_df_sliced = target(
    parse_ip_subsets_helper(
      x=ip_select_PC1_cor_tmp,
      ips_trans_fam_adj = ips_trans_fam_adj),
    transform=split(ip_select_PC1_cor_tmp, slices=6)
  ),
  
  core_CD_df = target(
    c(core_CD_df_sliced),
    transform=combine(core_CD_df_sliced)
  ),
  
  ip_select_PC1_cor= ip_select_PC1_cor_tmp %>%
    mutate(
      core_CD_df = core_CD_df,
      feature_core_data = map2(
        feature_data,
        core_CD_df,
        ~left_join(.x,
                   .y %>% 
                     unnest(.sep="_") %>%
                     dplyr::select(-CD_table_equiv_feat) %>%
                     set_names(str_remove(names(.), pattern="CD_table_")),
                   # dplyr::rename("feat"="CD_table_feat"),
                   by=c("feat", "composite_lin_source"))  #mind that features in core_CD_df are subset of those in feature_data
      )
    ) %>%
    dplyr::select(-one_of(c("feature_data", "core_CD_df"))),
  
  ip_select_module_anno =ip_modules_df %>% 
    unnest(cut_stats) %>% 
    dplyr::filter(
      corfnc=="spear",
      networktype=="signed_hybrid",
      powers==5,
      cut==8
    ) %>% 
    unnest(feature_stats) %>% 
    left_join(ip_anno %>% 
                dplyr::select(set_name, subset_name, composite_lin_source), 
              by=c("feat"="set_name")) %>% 
    count(colors, composite_lin_source) %>% 
    group_by(colors) %>% 
    mutate(module_size = sum(n),
           n_source = n_distinct(composite_lin_source),
           prop = n/module_size, 
           source_label = sprintf("%s(%s%%)", composite_lin_source, round(prop*100))) %>% 
    dplyr::filter(prop>0.05 & colors!="grey") %>% 
    arrange(desc(prop)) %>% 
    summarise(label_perc = paste0(source_label, collapse = "_"),
              label = paste0(composite_lin_source, collapse = "_"),
              label_major = composite_lin_source[which.max(prop)],
              module_size=unique(module_size),
              max_prop = max(prop),
              n_source = n_distinct(composite_lin_source)
    ) %>% 
    arrange(desc(module_size)) %>% 
    mutate(module_tag = paste0("I", 1:nrow(.))) %>% 
    left_join(dplyr::select(ip_select_modules,colors, varExplained), ., by="colors") %>% 
    left_join(ip_select_PC1_cor, by="colors"),
  
  # glycans modules - WGCNA modules -------------------------------------------
  glycan_modules_cor_ips_modules = cor_p_bh(
    x = glycan_select_scores,
    y = ip_select_scores,
    xname="glycan_module",
    yname="IP_module",
    method="pearson",
    use="pairwise.complete.obs"
  ),
  
  glycan_modules_cor_ips_modules_anno = glycan_modules_cor_ips_modules %>%
    left_join(dplyr::select(ip_select_module_anno, colors, 
                            varExplained_IP=varExplained, label_perc, label_major, 
                            module_size, module_tag, IP_feature_core_data=feature_core_data),
              by=c("IP_module"="colors")) %>% 
    left_join(glycan_select_module_anno, by="glycan_module", suffix=c("_IP", "_glycan")) %>% 
    dplyr::select(module_tag_glycan, igx, site, fucose, sialic_acid, hexose,
                  module_tag_IP, label_perc, IP_feature_core_data, p.val_storey,
                  correlation, p.val, p.val_bh, module_size_glycan=n, module_size_IP=module_size,
                  everything()),
  
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
  result_report = target(
    command = rmarkdown::render(
      knitr_in("cor_WGCNA.Rmd"),
      output_file = file_out(here::here("results/cor_WGCNA.html")),
      quiet = TRUE
    )
  )
  
  
)
