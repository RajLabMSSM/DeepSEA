## ~~~~~~~~~~~~~~~~~~        ~~~~~~~~~~~~~~~~~~ 
# ~~~~~~~~~~~~~~~~~~ DeepSEA ~~~~~~~~~~~~~~~~~~ 
## ~~~~~~~~~~~~~~~~~~        ~~~~~~~~~~~~~~~~~~ 

# Results generated using DeepSEA online server(http://deepsea.princeton.edu/help/)
## or the selene python module (https://github.com/FunctionLab/selene/blob/master/tutorials/analyzing_mutations_with_trained_models/analyzing_mutations_with_trained_models.ipynb).

library(dplyr)
library(ggplot2)
library(ggrepel)

printer <- function(..., v=T){if(v){print(paste(...))}}
root <- "~/Desktop/Fine_Mapping"

DeepSEA.subset_kunkle <- function(locally=T){
  # Summary stats
  if(locally){
    kunk <- data.table::fread(filepath(root,"Data/GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt.gz", nThread = 4))
  } else {
    kunk <- data.table::fread("/sc/orga/projects/ad-omics/data/AD_GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt", 
                              nThread = 4) 
  }  
  kunk$Pvalue <- as.numeric(kunk$Pvalue) 
  kunk.sig <- subset(kunk, Pvalue<=0.05)
  # Locus coordinates 
  kunkle_loci <- readxl::read_excel(file.path(root, "Data/GWAS/Kunkle_2019/Kunkle2019_supplementary_tables.xlsx"), 
                                    sheet = "Supplementary Table 8") %>%
    dplyr::select(Locus=`Lead SNV Gene`,SNP=`Top Associated SNV`, LD_block=`LD Block (GRCh37)`) %>% 
    data.table::data.table() %>% 
    tidyr::separate(LD_block, c("CHR","Start","End"), sep=":|-")
  
  kunkle.merge <- lapply(kunkle_loci$Locus, function(locus){
    printer("Merging data for locus:",locus)
    locus.sub <- subset(kunkle_loci, Locus==locus)
    kunk.sub <- subset(kunk.sig, (Chromosome==locus.sub$CHR) & (Position >= locus.sub$Start) & (Position <= locus.sub$End))
    if(nrow(kunk.sub)>0){
      kunk.sub <- cbind(Locus=locus.sub$Locus, kunk.sub)
      return(kunk.sub)
    } else {printer("No SNPs identified"); return(NULL)}
  }) %>% data.table::rbindlist(fill=T)
  
  # Tally before sig subset
  kunkle.merge %>% dplyr::group_by(Locus) %>% tally()
  kunkle.merge <- kunkle.merge %>% dplyr::mutate(FDR=p.adjust(Pvalue, method = "fdr", n = nrow(kunk)))
  # Tally after sig subset 
  subset(kunkle.merge, FDR<=0.05)  %>% dplyr::group_by(Locus) %>% tally() %>% data.frame()
  kunkle.merge.sig <- subset(kunkle.merge, FDR<=0.05) 
  
  data.table::fwrite(kunkle.merge.sig,file.path(root, "Data/GWAS/Kunkle_2019/Kunkle_Stage1_sig.txt"), sep="\t")
}




DeepSEA.vcf_subset <- function(locally=T){
  locus_list <- c("BIN1",
                  "TREM2",
                  "NYAP1",
                  # "LRRK2",
                  "PTK2B",
                  "SPI1",
                  "MS4A2",
                  "ABCA7")
  kunkle.sig <- data.table::fread(file.path(root,"./Data/GWAS/Kunkle_2019/Kunkle_Stage1_sig.txt"))
  kunkle_loci <- readxl::read_excel(file.path(root,"./Data/GWAS/Kunkle_2019/Kunkle2019_supplementary_tables.xlsx"), 
                     sheet = "Supplementary Table 8")
  if(locally){
    vcf_folder <- "./ROSMAP"
  } else{
    vcf_folder <- "/sc/orga/projects/ad-omics/brian/ROSMAP"
  }
  
  for(locus in locus_list){  
    printer("DeepSEA:: Preparing locus",locus,"...")
    regions.file <- file.path(vcf_folder, paste0(locus,".regions.txt"))
    if(locus=="LRRK2"){
      finemap_DT <- data.table::fread(file.path(root,"Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt"))
      coords <- data.frame(CHROM=unique(finemap_DT$CHR),POS=sort(subset(finemap_DT, P<=5e-8)$POS) )
    }else { 
      kunkle.sub <- subset(kunkle.merge.sig, Locus==locus) 
      coords <- data.frame(CHROM=unique(kunkle.sub$Chromosome)[1], POS=kunkle.sub$Position) 
    }
    data.table::fwrite(coords, regions.file, sep="\t", col.names = F)
    
    chr <- coords$CHROM[1]
    vcf.gz <- file.path(vcf_folder, paste0("DEJ_11898_B01_GRM_WGS_2017-05-15_",chr,".recalibrated_variants.vcf.gz"))
    locus.vcf <- file.path(vcf_folder,locus, paste0(locus,".vcf"))
    dir.create(dirname(locus.vcf),showWarnings = F, recursive = T)
    # Subset vcf with tabix ==> 
    ## ==> Force multiallelic sites to become biallelic ==> 
    ## filter snps with alt alleles longer than 100bp (DeepSEA's upper limit) ==>
    ## ==> Remove all but the first 5 cols  ==>
    ## ==> replace * with "'==>

    ## Important: need to include '-h' flag in Tabix to include the header so that bcftools can understand the file format.  
    cmd <- capture.output( 
      cat("tabix",vcf.gz,"-h -R",regions.file,
                  "| bcftools norm -m-", 
                  # "| bcftools filter -e '(STRLEN(ALT[0])>=90)'",
                    # "| bcftools query -H -f '%CHROM\\t%POS\\t%TYPE\\t%REF\\t%ALT\\n'", # Needs to be last bcftool command bc erases header
                  "| bcftools view -e '(STRLEN(REF)>100) | (STRLEN(ALT)>100)'",
                  "| cut -f1-5",
                  "| sed 's/\\*//g'",
                  "> ",locus.vcf)  
    )
    system(cmd) 
    cat("\n\n")
  } 
   
  results_links <- list(LRRK2="http://deepsea.princeton.edu/job/analysis/results/ea111803-8f88-4299-8a0e-872d8594b854",
                        ABCA7="http://deepsea.princeton.edu/job/analysis/results/5fa99587-1051-4a29-b961-5182a12ec81c",
                        BIN1="http://deepsea.princeton.edu/job/analysis/results/8753189d-75f7-47af-afd3-6ee7a21b6b66",
                        MS4A2="http://deepsea.princeton.edu/job/analysis/results/1732d2a7-0e0b-4c0b-aa71-cdc8e677aa79",
                        NYAP1="http://deepsea.princeton.edu/job/analysis/results/8682bd9e-6431-4b12-9c57-de9ab812d51d",
                        PTK2B="http://deepsea.princeton.edu/job/analysis/results/e2c72379-f163-4571-9e70-a19b828943b5",
                        SPI1="http://deepsea.princeton.edu/job/analysis/results/67c215c2-7474-46b9-b732-766f99391f35",
                        TREM1="http://deepsea.princeton.edu/job/analysis/results/19f228f9-ca64-4a62-bda2-f3e23d7b374f") 
}


DeepSEA.gather_predictions <- function(deepsea_path){ 
  printer("DeepSEA:: Gathering prediction results...")
  loci <- dirname(list.files(deepsea_path,"infile.vcf.out.funsig", recursive = T, full.names = F))
  DS.predict <- lapply(loci, function(locus){ 
    # SNP-wise 'Functional Significant Scores'
    funsig <- data.table::fread(file.path(deepsea_path,locus,"infile.vcf.out.funsig")) 
    # SNP-wise probabilites of being eQTL, GWAS, or HGMD (The Human Gene Mutation Database) hit
    snpclass <- data.table::fread(file.path(deepsea_path,locus,"infile.vcf.out.snpclass"))  
    deepsea.dat <- cbind(snpclass, `Functional significance score`=funsig$`Functional significance score`) 
    return(deepsea.dat)
  }) %>% data.table::rbindlist(fill=T)
  return(DS.predict)
}
 
DeepSEA.gather_enrichment <- function(deepsea_path, top_annots=F){
  printer("DeepSEA:: Gathering enrichment results...")
  enrich.files <- list.files(deepsea_path, "infile.vcf.out.logfoldchange", recursive = T, full.names = T)
  enrich.DF <- lapply(enrich.files, function(f){
    locus <- basename(dirname(f))
    enrich <- data.table::fread(f) 
    enrich <- cbind(Locus=locus, enrich)
    return(enrich)
  }) %>% data.table::rbindlist(fill=T) 
  DS.enrich <- data.table:::melt.data.table(enrich.DF, 
                                              id.vars = c("Locus","V1","chr","pos","name","ref","alt"),
                                              variable.name = "Annotation",
                                              value.name = "logFC") 
  if(top_annots!=F){
    DS.enrich <- DS.enrich %>% dplyr::group_by(Locus, chr, pos, Annotation) %>%  
      dplyr::group_by(Locus, chr, pos) %>%
      top_n(n=top_annots, wt=logFC)
  }
  return(DS.enrich)
}


DeepSEA.prepare_data <- function(root="~/Desktop/Fine_Mapping",
                                 deepsea_path="./ROSMAP",
                                 GWAS_path=file.path(root,"Data/GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt.gz") ){
  # merge DeepSEA results
  DS.predict <- DeepSEA.gather_predictions(deepsea_path)
  DS.enrich <- DeepSEA.gather_enrichment(deepsea_path, top_annots=1)
  deepsea.DAT <- data.table:::merge.data.table(DS.predict,
                                               DS.enrich,
                                               by= c("chr","pos","ref","alt"), all.x = T) %>% 
    dplyr::mutate(Chromosome = as.numeric(gsub("chr","",chr)), Position=pos)
 
  # merge with Kunkle GWAS
  printer("DeepSEA:: Merging with GWAS data...")
  kunk <- data.table::fread(GWAS_path, nThread = 4)
  deepsea.DAT <- data.table:::merge.data.table(deepsea.DAT,
                                                kunk,
                                                by= c("Chromosome","Position"), all.x = T) %>% 
    dplyr::rename(SNP=MarkerName) 
  remove(kunk, DS.predict, DS.enrich)
  
  # Fill NA SNPs
  deepsea.DAT[is.na(deepsea.DAT$SNP), ]$SNP <- paste0(deepsea.DAT[is.na(deepsea.DAT$SNP), ]$chr,":", 
                                                      deepsea.DAT[is.na(deepsea.DAT$SNP), ]$pos)
  return(deepsea.DAT)
} 




DeepSEA.corrplot <- function(deepsea_path="./ROSMAP", deepsea.DAT, save_path=F){ 
  # PP.cols <- c(grep(".Probability",colnames(deepsea.DAT), value = T), "mean.PP")
  cor.dat <- deepsea.DAT[,c("eQTL-probability","GWAS-probability","HGMD-probability","Functional significance score")]
  cor.df <- cor(cor.dat)
  cor.df[is.na(cor(cor.dat))] <- 0 
  if(save_path!=F){png(file = file.path(deepsea_path,"_plots","DeepSEA.corrplot.png") )}
  corrplot::corrplot(cor.df,
                     # title = "Correlations Between Fine-mapping PPs and DeepSEA Predictions",
                     bg = "black",
                     order = "hclust",
                     # addrect = 2,
                     type = "upper",
                     tl.col = "black", 
                     method = "ellipse",
                     addCoef.col = "cyan",
                     number.cex = 1,
                     tl.srt = 45,
                     tl.cex = .7
                     # col = colorRampPalette(c("red","white","blue"))(200)
                     # col = RColorBrewer::brewer.pal(5,"RdBu")
                     ) 
  dev.off()
}



 
DeepSEA.plot_predictions <- function(deepsea.DAT, 
                                     label.snps,
                                      keep_cols, 
                                      facet_xlabs = T, 
                                      custom_ylab=F){
  deepsea.melt <- data.table::melt.data.table(data.table::data.table(deepsea.DAT),
                                              id.vars = c("Locus","SNP","chr","pos"),
                                              measure.vars = c("Kunkle(2019) GWAS",
                                                               "eQTL-probability",
                                                               "GWAS-probability",
                                                               "HGMD-probability",
                                                               "Functional significance score"),
                                              variable.name = "Prediction",
                                              value.name = "Probability") 
  deepsea.melt <- subset(deepsea.melt, Prediction %in% keep_cols, drop=F)
  deepsea.melt$Prediction <- gsub("-| ","\n", deepsea.melt$Prediction)
  deepsea.melt <- deepsea.melt %>% arrange(pos)
  deepsea.melt$Prediction <- factor(deepsea.melt$Prediction, levels = unique(deepsea.melt$Prediction) )
  deepsea.melt <- deepsea.melt %>% arrange(pos)
  deepsea.melt$Mb <- deepsea.melt$pos/1000000
  
  # Plot
  gp <- ggplot(data=deepsea.melt, aes(x=Mb, y=Probability, color=Probability)) +
    geom_point() + 
    facet_grid(Prediction~Locus, scales = "free", drop = T) + 
    theme_dark() +
    scale_color_viridis_c() +  
    geom_hline(yintercept = 0) +  
    geom_point(data=subset(deepsea.melt, SNP %in% label.snps ), 
               pch=21, fill=NA, size=4, color="blue", stroke=1.5, alpha=0.8,  label.padding = .25,) +
    geom_label_repel(data=subset(deepsea.melt, SNP %in% label.snps ), aes(label=SNP), 
                     alpha=.7,
                     # aes(color=Probability),  
                     nudge_x = .5,  
                     label.padding = .25,
                     label.size=NA, 
                     seed = 1, 
                     size = 3,
                     min.segment.length = 1) + 
    theme(#axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          strip.text.y = element_text(angle = 0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())  +
    scale_x_continuous(labels=function(x) sprintf("%.2f", x))
  if(facet_xlabs==F){
    gp <- gp + 
      theme( strip.background.x = element_blank(), 
             strip.text.x = element_blank())
  }
  if(custom_ylab!=F){
    gp <- gp + labs(y=custom_ylab, color=custom_ylab) 
  } else {gp <- gp +scale_y_continuous(limits = c(0,1.1), breaks = c(0,.5, 1))  }
  return(gp)
}



DeepSEA.plot_enrichment <- function(deepsea.DAT, label.snps){
  dat <- subset(deepsea.DAT, SNP %in% label.snps)
  # infile.vcf.out.logfoldchange: Chromatin feature probability log fold changes log(p_alt/(1-p_alt))-log(p_ref/(1-p_ref)) for each variant. Computed by comparing the 'ref' and 'alt' files. Note that log fold change can be noisy when probabilities are small, therefore this score should be considered together with probability differences and E-values when evaluating the impact of a variant. 
  enrich.plots <- lapply(sort(unique(dat$Locus)), function(locus){
    printer("DeepSEA:: Enrichment plot for",locus) 
    # Average by Annotation
    locus.annot <- subset(dat, Locus==locus, .drop=T) %>% 
      dplyr::group_by(SNP) %>%
      top_n(n=1, wt=logFC)
    # ------- plot --------
    enrich.plot <-  ggplot(locus.annot, aes(y=logFC, x=SNP, fill=logFC)) +  
      # geom_col(position = "dodge", show.legend = F) +
      geom_bar(stat="identity", show.legend = F ) +
      coord_flip() +
      facet_grid(~Locus, drop=T, scales="free", space = "free") +
      geom_text(aes(y=0, label=paste0(SNP,"\n",Annotation), hjust=.5, srt=0), 
                size=3, color="magenta", nudge_y = -.005) +
      theme_dark() +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            rect = element_rect(fill = "transparent"),
            panel.background = element_rect(fill = "transparent"),
            # strip.text.y = element_text(angle = 0),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +  
      # labs(subtitle=paste0(locus,"\nEnrichment")) +
      scale_fill_viridis_c() + 
      ylim(c(min(dat$logFC)-2, max(dat$logFC)))
    # print(enrich.plot)
    return(enrich.plot)
  })
}



DeepSEA.stacked_plots <- function(root="~/Desktop/Fine_Mapping/",
                                  deepsea_path="./ROSMAP",
                                  deepsea.DAT=DeepSEA.prepare_data(deepsea_path=deepsea_path,
                                                                   GWAS_path=file.path(root,"Data/GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt.gz") ),
                                  label_topn=3){   
  # Get the top n SNPs to label
  label.snps <- (deepsea.DAT %>% dplyr::select(Locus, SNP, `Functional significance score`) %>%
    dplyr::group_by(Locus) %>%
    top_n(n=label_topn, wt=`Functional significance score`))$SNP
  deepsea.DAT <- subset(deepsea.DAT, !(Locus %in% c("LRRK2",NA)) ) %>% 
    dplyr::mutate("Kunkle(2019) GWAS"= -log10(as.numeric(Pvalue)))
  

  # GWAS row  
  gp1 <- DeepSEA.plot_predictions(deepsea.DAT, 
                                  label.snps = label.snps,
                                  keep_cols = grep("Kunkle",colnames(deepsea.DAT), value=T ), 
                                  custom_ylab = "-log10(P-value)") + 
    scale_color_gradient(low="blue", high="red")  
  # DeepSEA rows 
  gp2 <- DeepSEA.plot_predictions(deepsea.DAT, 
                                  label.snps = label.snps,
                                   keep_cols = grep("probability|Functional ",colnames(deepsea.DAT), value=T ), 
                                   facet_xlabs = T) 
  # Enrichment plots
  enrich.plots <- DeepSEA.plot_enrichment(deepsea.DAT, label.snps)
  ep <- cowplot::plot_grid(plotlist = enrich.plots, nrow=1)
  
  # merge all plots
  # cowplot::plot_grid(gp1, gp2, ep, ncol = 1)
  gg <- ggpubr::ggarrange(gp1, gp2, ep, nrow = 3, align = c("hv"), 
                          widths = c(1,.8,.5), 
                          heights = c(.5,1.5,.75))
  # gg
  if(save_path!=F){
    ggsave(filename = file.path(deepsea_path,"_plots","DeepSEA.predict.enrich.png"),
           plot = gg, dpi = 600, height=10, width=25)
  } 
  return(gp)
}



