plan_slurm = drake_plan(
  
  # see template file specified in make_slurm.R 
  # for valid list names for the resources argument
  # memory -> --mem-per-cpu
  # ncpus -> --cpus-per-task (and ntasks=1, so sets the #cpus for the target)
  
  ips_trans_fam_adj = target(
    fread(file_in("data/ips_trans_fam_adj.tsv")),
    resources=list(ncpus=1,  memory=1024*10),
    hpc = FALSE
  ),
  
  ip_anno = target(
    fread(file_in("data/data_annotations/all_immunophenotypes_annotation_av.csv")) %>%
      mutate_all(na_if, y="") %>%   # set all empty string fields to NA 
      set_names(names(.) %>% tolower(.) %>% str_replace_all(., c( "\\.$"="", "[\\.]+"="_", " "="_"))) %>%
      mutate(lineage=recode(lineage, `04-Aug`="4-8"),
             composite_lin_source = ifelse((source=="Lin"|source=="MFI"), source, lineage)),
    resources=list(ncpus=1,  memory=1024*10),
    hpc = FALSE
  ),
  
  # correlation matrices -------------------------------------------
  
  # dummy = target({print("echoing SGE");system("echo $SGE_CELL")}, resources=list(ncpus=1, memory=100)),
  
  ips_pearson = target(fread("data/ips_slurm_blockwise_pearson.tsv"),resources = list(ncpus=2, memory=1024*40)),
  ips_bicor = target(fread("data/ips_slurm_blockwise_bicor.tsv"),resources = list(ncpus=2, memory=1024*40)),
  ips_spear = target(fread("data/ips_slurm_blockwise_spearman.tsv"),resources = list(ncpus=2, memory=1024*40)),

  p_ips_cor_density = target(
    get_ip_cor_density(corr_like = corr_like,
                       ips_trans_fam_adj = ips_trans_fam_adj,
                       ip_anno = ip_anno),
    transform= map(corr_like=c(ips_pearson, ips_bicor, ips_spear)),
    resources=list(ncpus=1,  memory=1024*210),
  ),

  # # WGCNA IPs -------------------------------------------
  # 
  # # 1. scale-free topology: power parameter
  # 
  # similarity = target(
  #   similarities_from_cor(corr_like, networktype),  # list of lists: level1= corFnc, level2=networktype
  #   transform = cross(networktype=c('signed', 'signed_hybrid', 'unsigned'), corr_like=c(ips_pearson, ips_bicor, ips_spear)),
  #   resources = list(ncpus=1, memory=1024*70)
  # ),
  # 
  # ips_scale_free = target(
  #   pickSoftThreshold.fromSimilarity(similarity=as.matrix(similarity),
  #                                    powerVector = c(c(1:10), seq(from = 12, to=20, by=2)),
  #                                    verbose = 0,
  #                                    moreNetworkConcepts = TRUE) %>%
  #     .$fitIndices,
  #   transform = map(similarity,
  #                   .id = c(networktype, corr_like)),
  #   resources = list(ncpus=1, memory=1024*70)
  # ),
  # 
  # sft_df = target(
  #   bind_rows_plan(ips_scale_free, drakeplan=ignore(plan_slurm), vars=c("networktype", "corr_like")),
  #   transform = combine(ips_scale_free),
  #   hpc=FALSE
  # ),
  # 
  # sf_plots = target(
  #   sf_plot(sft_df),
  #   hpc=FALSE
  # ),
  # 
  # sf_plots_rds = target(
  #   saveRDS(sf_plots, file = file_out("results/figures/WGCNA_pars/IPs/scale_free/sf_plots.RDS")),
  #   hpc=FALSE
  # ),
  # 
  # # 2. parameter search ~ modularity
  # 
  # ip_module_stats = target(
  #   wgcna_parameter_search(similarity=similarity,
  #                          corr_like=corr_like,
  #                          method=method,
  #                          powers=softpower,  #don't use looping over powers, to save total memory requirement per SLURM job
  #                          # powers=powers,   # looping over powers for same similarity, saves on number of large object imports (similarity and corr_like) required
  #                          cut_pars=expand.grid(deep =c(4, 1),
  #                                               minModuleSize=c(3,5,10,20)),
  #                          create_plot=FALSE,
  #   ),
  #   # transform = map(similarity=c(similarity_signed_hybrid_ips_pearson, similarity_signed_hybrid_ips_bicor))
  #   transform = map(.data=!!ip_network_pars,
  #                   .id = c(nwtype, corr_like, softpower)),
  #   resources = list(ncpus=1, memory=1024*220, walltime=3600*48)
  # ),
  # 
  # ip_module_stats_combined = target(
  #   bind_rows_plan(ip_module_stats, drakeplan=ignore(plan_slurm), vars=c("corfnc", "nwtype", "method", "softpower")),
  #   transform=combine(ip_module_stats),
  #   resources= list(ncpus=1),
  #   hpc=FALSE
  # ),
  # combined_bare_rds = target(
  #   saveRDS(ip_module_stats_combined, file = file_out("results/ip_module_stats_combined.RDS")),
  #   hpc=FALSE
  # ),

  #
  # network_module_heatmaps = make_module_heatmap_pdf(plots = glycan_module_stats$plots,
  #                                                   file_path = file_out("results/figures/WGCNA_pars/glycans/overview.pdf"),
  #                                                   pointsize = 2,
  #                                                   width = 10,
  #                                                   height=10,
  #                                                   network_pars = glycan_network_pars),
  #
  # glycan_modularity_overview = write_csv(ip_module_stats$modularity_overview, path = file_out("results/ip_modularity_overview.csv")),
  trace = TRUE
)

plan_slurm$resources <- lapply(plan_slurm$resources, eval)