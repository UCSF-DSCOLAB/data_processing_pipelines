
isotype_ctl_plot = function(ADT_counts, isotype_ctl_data, params){
  iso_max = apply(isotype_ctl_data, 2, max)
  adt_tab = tibble("iso_max"=iso_max+0.1, "nadt"=colSums(ADT_counts))

  nabove = adt_tab %>% filter(iso_max > params["ADT_isotype_ctl.upper",]) %>% nrow() # 41

  p= ggplot(adt_tab, aes(x=nadt, y=iso_max))+
    geom_point(size=0.1,alpha=0.1) +
    geom_hex(bins=100) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_classic() +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    geom_hline(yintercept=params["ADT_isotype_ctl.upper",], alpha=0.5)+
    geom_vline(xintercept=params["nCount_ADT.upper",], alpha=0.5)+
    ylab("isotype control max")+
    xlab("nCount_ADT")+
    geom_text(data = data.frame(x=1, y=params["nCount_ADT.upper",],text=sprintf("%s droplets > isotype ctl max", nabove)), 
              aes(x=x,y=y,label=text), hjust=0, vjust=1, color="darkred", size=4)
  return(ggExtra::ggMarginal(p, type = "histogram", bins=50))
}

background_plot <- function(gex_data, ab_data, params){
  tab = tibble(
  "nCount_RNA" = colSums(gex_data),
  "nCount_ADT" = colSums(ab_data)
  ) %>% filter(nCount_RNA > 0, nCount_ADT > 0)
  if (!any(str_detect(rownames(params), "background"))){
    adt_lower=10^1.5
    adt_upper=params['nCount_ADT.upper',]
    rna_upper=max(300, params['nFeature_RNA.lower',])
  } else {
    adt_lower=params['background_ADT.lower',]
    rna_upper=params['background_RNA.upper',]
    adt_upper=params['background_ADT.upper',]
  }

  nbackground = tab %>% filter(nCount_RNA < rna_upper,
                nCount_ADT < adt_upper,
                nCount_ADT > adt_lower) %>% nrow() 

  p = ggplot(tab, aes(x=nCount_RNA, y=nCount_ADT))+
    geom_point(size=0.1,alpha=0.1) +
    geom_hex(bins=100) +
    geom_vline(xintercept=rna_upper, alpha=0.5)+
    geom_hline(yintercept=adt_lower, alpha=0.5)+
    geom_hline(yintercept=adt_upper, alpha=0.5)+
    scale_fill_distiller(palette = "RdYlBu") +
    theme_classic() +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    geom_text(data = data.frame(x=1, y=adt_upper,text=sprintf("%s droplets (%s%%)", nbackground, round(nbackground/nrow(tab)*100,1))), 
                                aes(x=x,y=y,label=text), hjust=0, vjust=1, color="darkred", size=4)

  return(ggExtra::ggMarginal(p, type = "histogram", bins=50))
}

adt_plot <- function(sobj){
  adt_dat = sobj@assays$ADT@data %>% as_tibble(rownames="adt")
  adt_dat %>% 
    pivot_longer(-adt, names_to="cell", values_to="val") %>%
    ggplot(aes(x=val))+
    geom_histogram()+
    theme_classic()+
    theme(axis.title.x=element_blank())+
    ylab("number of cells")+
    facet_wrap(~adt, scales="free")
}

scatterhist <- function(x, y, sobj, params, add_stats=T, log.scale.x=F, log.scale.y=F, group_col=NULL, legend_only=F) {

  if (!is.null(group_col)){
    p = ggplot(sobj@meta.data, aes(x=!!sym(x)+0.1, y=!!sym(y)+0.1, col=!!sym(group_col))) +
      scale_color_manual()+
      geom_point(size=0.01)+
      scale_color_manual(values=dittoColors())

    if (str_detect(group_col, "DROPLET.TYPE")){
      doublet_colors = c("gray", dittoColors()[c(3,2,1,4)])
      names(doublet_colors)= c("AMB", "Intra.DBL", "Inter.DBL", "SNG", "DBL")
      p = p+scale_color_manual(values=doublet_colors)
    }

  } else {
    p = ggplot(sobj@meta.data, aes(x=!!sym(x)+0.1, y=!!sym(y)+0.1)) +
      geom_point(size=0.1,alpha=0.1) +
      geom_hex(bins=100)+
      scale_fill_distiller(palette = "RdYlBu") 
  }

  p = p + 
    theme_classic() +
    xlab(x)+
    ylab(y)
  if (legend_only & !is.null(group_col)){
    return(get_legend(p))
  }
  if (add_stats){
    x_upper = as.numeric(params[paste0(x,".upper"),])
    x_lower = as.numeric(params[paste0(x,".lower"),])
    y_upper = as.numeric(params[paste0(y,".upper"),])
    y_lower = as.numeric(params[paste0(y,".lower"),])
    xy_data = FetchData(sobj, vars=c(x, y))
    filter_cells = xy_data[,1] <= x_upper & 
      xy_data[,1] >= x_lower & 
      xy_data[,2] <= y_upper & 
      xy_data[,2] >= y_lower
    filter_cell_percent = round( sum(filter_cells) / ncol(sobj) ,3) * 100
    p = p+
      geom_vline(xintercept = x_upper) +
      geom_vline(xintercept = x_lower) +
      geom_hline(yintercept = y_upper) +
      geom_hline(yintercept = y_lower) + 
      geom_rect(aes(xmin=x_lower, xmax=x_upper, ymin=y_lower, ymax=y_upper), color="red", alpha=0) +
      geom_text(data = data.frame(x=x_lower, y=y_upper,text=paste0(filter_cell_percent, "%")), aes(x=x,y=y,label=text), hjust=0, vjust=1, color="darkred", size=8)
  }
   if (log.scale.x){
    p = p+ scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
          labels = scales::trans_format("log10", scales::math_format(10^.x))) 
  } 
  if (log.scale.y){
    p = p+ scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
          labels = scales::trans_format("log10", scales::math_format(10^.x))) 
  } 
  if (!is.null(group_col)){
    p = p+theme(legend.position="none")
    return( ggExtra::ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE))
  } else {
    return( ggExtra::ggMarginal(p, type = "histogram", bins=50) )
  }
}

load_params <- function(params_df){
  params_df %>% 
  filter(parameter!="reviewed") %>%  
  mutate(value=case_when(
    !is.na(value) ~ as.numeric(value),
    str_detect(parameter, "upper") ~ Inf,
    str_detect(parameter, "lower") ~ -Inf
  )) %>%
  column_to_rownames("parameter") 
}

quantile_frac_table <- function(sobj, adt.present){
  df <- data.frame(
          cell_counts=seq(0, 1.01, 0.1)*dim(sobj@meta.data)[1],
          percent.mt=quantile(sobj@meta.data[["percent.mt"]], seq(0, 1.01, 0.1)),
          percent.ribo=quantile(sobj@meta.data[["percent.ribo"]], seq(0, 1.01, 0.1)),
          nFeature_RNA=quantile(sobj@meta.data[["nFeature_RNA"]], seq(0, 1.01, 0.1)),
          nCount_RNA=quantile(sobj@meta.data[["nCount_RNA"]], seq(0, 1.01, 0.1)),
          row.names=seq(0, 1.01, 0.1)
          )
  if (adt.present){
    df["nCount_ADT"] <- quantile(sobj@meta.data[["nCount_ADT"]], seq(0, 1.01, 0.1))
    df["nFeature_ADT"] <- quantile(sobj@meta.data[["nFeature_ADT"]], seq(0, 1.01, 0.1))
    df["isotype_ctl_max"] <- quantile(sobj@meta.data[["isotype_ctl_max"]], seq(0, 1.01, 0.1))
  }
  return(df %>% 
    mutate(cell_counts=floor(cell_counts)) %>%
    mutate(across(contains("percent"), ~round(.x, 2))) %>%
    mutate(across(contains("nCount"), ~floor(.x))) %>% 
    mutate(across(contains("nFeature"), ~floor(.x))))
}

filter_cells = function(sobj, params, adt.present){
  sobj = subset(sobj, 
              nFeature_RNA <= params["nFeature_RNA.upper",] & 
                nFeature_RNA >= params["nFeature_RNA.lower",] &
                nCount_RNA <= params["nCount_RNA.upper",] &
                nCount_RNA >= params["nCount_RNA.lower",] &
                percent.mt <= params["percent.mt.upper",] &
                percent.mt >= params["percent.mt.lower",] &
                percent.ribo <= params["percent.ribo.upper",] &
                percent.ribo >= params["percent.ribo.lower",] )
  if (adt.present){
    sobj = subset(sobj, 
                  nFeature_ADT <= params["nFeature_ADT.upper",] &
                    nFeature_ADT >= params["nFeature_ADT.lower",] & 
                    nCount_ADT <= params["nCount_ADT.upper",] &
                    nCount_ADT >= params["nCount_ADT.lower",] &
                    isotype_ctl_max <= params["ADT_isotype_ctl.upper",] )

  }
  return(sobj)
}

load_clonotypes <- function(library, data_type, contig_path) {
  vdj_library=str_replace(library, "SCG", sprintf("SC%s", substr(data_type, 1,1) ))

  annot = read.csv(contig_path) %>%
      filter(
          high_confidence=="true",
          full_length=="true",
          productive=="true",
          raw_clonotype_id!="" & !is.na(raw_clonotype_id) # target=NA, but may be turned into ""
      )
  
  # The TCR/BCR sequencies for each barcode are split into different chains across multiple lines, but all those rows have same raw_clonotype_id values. Extracting the barcode and raw_clonotype_id and removing the redundancies.
  annot_set = annot %>% select(barcode, raw_clonotype_id) %>% unique()
  
  # Get AA sequences and other info per contig and join them together per clonotype_id
  columns_grab <- c("v_gene","d_gene","j_gene","c_gene","fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4","reads","umis")
  clonotype_data <- do.call(rbind, lapply(
      1:nrow(annot_set),
      function(i) {
          this_annot <- annot[annot$barcode==annot_set$barcode[i] & annot$raw_clonotype_id==annot_set$raw_clonotype_id[i],]
          this_clonotype <- data.frame(barcode = annot_set$barcode[i])
          for (col in columns_grab) {
              this_clonotype[[col]] <- paste0(this_annot[,col], collapse = ";")
          }
          this_clonotype
      }
  ))
  
  # Make barcode the rownames and adjust the column names.
  clonotype_data = clonotype_data %>% 
      column_to_rownames("barcode") %>% 
      setNames(paste0(data_type, ".", columns_grab))
  
  return(clonotype_data)
}

process_adt <- function(ADT_counts, ADT_background){
  
  ADT_background <- ADT_background[rownames(ADT_counts), !colnames(ADT_background) %in% 
                                     colnames(ADT_counts)]
  isotype_ctls = rownames(ADT_counts)[str_detect(tolower(rownames(ADT_counts)), "iso")]
  adt_norm = DSBNormalizeProtein(
    # cell-containing droplet raw protein count matrix
    cell_protein_matrix = ADT_counts, 
    # empty/background droplet raw protein counts
    empty_drop_matrix = ADT_background,
    # recommended step: model + remove the technical component of each cell's protein library
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    quantile.clipping = TRUE, # adding this in to clip extreme outliers
    quantile.clip = c(0.01, 0.99), # slightly more stringent than the default to clip large outliers
    isotype.control.name.vec = isotype_ctls
  )
  return(adt_norm)
}

loadFreemuxletData <- function(sObj, freemuxletSampleF) {
  fml0 = ""
  if(grepl("gz$", freemuxletSampleF)) {
    zz = gzfile(freemuxletSampleF, 'rt')
    fml0 = read.table(zz, header=T, row.names=1)
    close(zz)
  } else {
    fml0 = read.table(freemuxletSampleF, header=T, row.names=1)
  }
  fml = sObj@meta.data %>% 
    as_tibble(rownames="cell_id") %>% 
    left_join(fml0 %>% as_tibble(rownames="cell_id"), by=c("cell_id")) 
  tbl = table(fml[, "DROPLET.TYPE"])
  sObj@misc$scStat$fmlDropletTypeComp = as.matrix(tbl)
  sObj@misc$scStat$fmlDropletTypeProp = as.matrix(prop.table(tbl))
  
  sObj = AddMetaData(sObj, fml$`DROPLET.TYPE`, "DROPLET.TYPE")  
  sObj = AddMetaData(sObj, fml$`BEST.GUESS`, "BEST.GUESS")
  
  return(sObj)
}

make_doublet_plot = function(sobj, adt.present=F){
  doublet_colors = c("gray", dittoColors()[c(3,2,1)])
  sobj = subset(sobj, DROPLET.TYPE.FINAL %in% c("AMB", "Intra.DBL", "Inter.DBL", "SNG"))
  names(doublet_colors)= c("AMB", "Intra.DBL", "Inter.DBL", "SNG")

  if (adt.present){
    p = ggplot(sobj@meta.data, aes(x=nCount_RNA, y=nCount_ADT, col=DROPLET.TYPE.FINAL))
  } else {
    p = ggplot(sobj@meta.data, aes(x=nCount_RNA, y=percent.mt, col=DROPLET.TYPE.FINAL))
  }
  p = p +
        geom_point(size=0.01)+
        scale_color_manual(values=doublet_colors)+
        theme_bw()+
        labs(col="DROPLET.TYPE")+
        theme(panel.grid=element_blank(),
              axis.title.y=element_blank())
  p = p+theme(legend.position="none") 
  p = p %>% ggMarginal(type = "density", groupColour = TRUE, groupFill = TRUE)
  if (adt.present){
    p2 = p+scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))
  } else {
    p2 = p
  }
  p2 = p2 +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^(x)),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) 
  legend = get_legend(p2)
  p2 = p2+theme(legend.position="none")
  p2 = p2 %>% ggMarginal(type = "density", groupColour = TRUE, groupFill = TRUE)

  return(plot_grid(p, p2, legend, rel_widths = c(1, 0.9, 0.35), ncol=3))
}


print_message <- function(...) {
  cat("[", format(Sys.time()), "] ", ..., "\n", sep="")
}

## from aarao_scripts
TriPlot <- function(sobj, features, reduction.use="umap", group.by="seurat_clusters", jitter=TRUE) {
  if (jitter){
    pt.size = 1
  } else {
    pt.size = 0
  }
  plots <- list()
  layout = rbind(c(1,1,2,2),
                 c(1,1,2,2),
                 c(3,3,3,3),
                 c(3,3,3,3))
  for (f in features){
    tmp_plots <- list(p1=DimPlot(sobj,
                                 group.by=group.by, 
                                 reduction=reduction.use, 
                                 label=TRUE) + NoLegend(),
                      p2=FeaturePlot(sobj,
                                     features=f,
                                     reduction=reduction.use),
                      p3=VlnPlot(sobj, 
                                 features=f,
                                 group.by=group.by,
                                 pt.size = pt.size) + NoLegend()
    )
    plots[[f]] <- grid.arrange(grobs=tmp_plots, layout_matrix=layout)
  }
  plots
}


make_plots = function(sobj, params, adt.present=F, add_stats=T, group=NULL) {
  plot_list = list()
  plot_list[[length(plot_list)+1]] = scatterhist("percent.ribo","percent.mt",sobj,params, add_stats, group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.mt",sobj,params, add_stats, log.scale.x=F,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.mt",sobj,params, add_stats, log.scale.x=T, group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.mt",sobj,params, add_stats, log.scale.x=T,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.mt",sobj,params, add_stats, log.scale.x=F,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.ribo",sobj,params, add_stats, log.scale.x=T,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.ribo",sobj,params, add_stats, log.scale.x=F,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.ribo",sobj,params, add_stats, log.scale.x=T,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.ribo",sobj,params, add_stats, log.scale.x=F,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","nFeature_RNA",sobj,params, add_stats, log.scale.x=T, log.scale.y=T,  group_col=group)
  plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","nFeature_RNA",sobj,params, add_stats, log.scale.x=F, log.scale.y=F,  group_col=group)
  if (adt.present){
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_ADT","percent.mt",sobj,params, add_stats, log.scale.x=T,  group_col=group)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_ADT","percent.mt",sobj,params, add_stats, log.scale.x=F,  group_col=group)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_ADT","percent.ribo",sobj,params, add_stats, log.scale.x=T,  group_col=group)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_ADT","percent.ribo",sobj,params, add_stats, log.scale.x=F,  group_col=group)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","nCount_ADT",sobj,params, add_stats, log.scale.x=T,log.scale.y=T, group_col=group)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","nCount_ADT",sobj,params, add_stats, log.scale.x=F,log.scale.y=F,  group_col=group)
  }
  if (!is.null(group)){
    plot_list[[length(plot_list)+1]] = scatterhist("percent.ribo","percent.mt",sobj,params, group_col=group, legend_only=T)
  }
  # add a legend
  return(plot_list)
}