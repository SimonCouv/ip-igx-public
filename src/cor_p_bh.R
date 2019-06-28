cor_p_bh <- function(x, y, xname="glycan_module", yname="IP", ...){
  
  res <- WGCNA::corAndPvalue(x, y, ...)
  
  cor_dt_m <- res$cor %>% 
    data.table(., keep.rownames = TRUE) %>% 
    setnames(., c(xname, names(.)[-1])) %>% 
    .[glycan_module!="grey",] %>% 
    melt(id.vars=xname,
         variable.name= yname,
         value.name="correlation"
    ) %>% 
    setkeyv(c(xname,yname))
  
  p_dt_m <- res$p %>% 
    data.table(., keep.rownames = TRUE) %>% 
    setnames(., c(xname, names(.)[-1])) %>% 
    .[glycan_module!="grey",] %>% 
    melt(id.vars=xname,
         variable.name= yname,
         value.name="p.val"
    ) %>% 
    .[, p.val_bh := p.adjust(.[,p.val], method = "BH")] %>%
    .[, p.val_storey := qvalue(p.val, pi0.method="smoother")$qvalues] %>% 
    setkeyv(c(xname,yname))
  
  # join and order
  cor_dt_m[p_dt_m] %>% 
    .[order(p.val_bh)]
}


