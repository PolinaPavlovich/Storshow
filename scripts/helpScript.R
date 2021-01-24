#library(rjson)
library(yaml)

#config_json <- fromJSON(file = file.path(inputDir, "config.json"))
config_file <- yaml.load_file(file.path(inputDir, "config.yaml"))

ShinyInternal_script <- config_file[["ShinyInternal_script"]]
storshow_functions <- config_file[["storshow_functions"]]
count_matrix_fn <- config_file[["count_matrix_fn"]]
extra_data_dds <- config_file[["extra_data_dds"]]
dataDirectory_orig <- config_file[["dataDirectory_orig"]]
samp <- config_file[["samp"]]
avg_cond_thr <- config_file[["avg_cond_thr"]]
treat_cond <- config_file[["treat_cond"]]
ctr_cond <- config_file[["ctr_cond"]]
DESeq2_padj_thr <- config_file[["DESeq2_padj_thr"]]
DESeq2_absLog2FC_thr <- config_file[["DESeq2_absLog2FC_thr"]]
MA_plots_ylim <- config_file[["MA_plots_ylim"]]
ORA_padj_thr <- config_file[["ORA_padj_thr"]]
maxCores <- config_file[["maxCores"]]
GSEA_padj_thr <- config_file[["GSEA_padj_thr"]]
gmt.file <- config_file[["gmt.file"]]
avg_pair_overl_thr <- config_file[["avg_pair_overl_thr"]]
levels_ordered <- config_file[["levels_ordered"]]
colors_ordered <- config_file[["colors_ordered"]]
titlePage <- config_file[["titlePage"]]
first_gene_to_show <- config_file[["first_gene_to_show"]]
first_custom_GSEA_to_show <- config_file[["first_custom_GSEA_to_show"]]
first_allPath_GSEA_to_show <- config_file[["first_allPath_GSEA_to_show"]]
extra_data_comment <- config_file[["extra_data_comment"]]
front_target_barplot_height <- config_file[["front_target_barplot_height"]]
front_target_barplot_width <- config_file[["front_target_barplot_width"]]
front_extra_barplot_height <- config_file[["front_extra_barplot_height"]]
front_extra_barplot_width <- config_file[["front_extra_barplot_width"]]
showCateg <- config_file[["showCateg"]]
showCategCnet <- config_file[["showCategCnet"]]
maxNgenes <- config_file[["maxNgenes"]]
DE_df_height <- config_file[["DE_df_height"]]
MA_plots_width <- config_file[["MA_plots_width"]]
Volcano_plots_width <- config_file[["Volcano_plots_width"]]
Volcano_nShow_geneNames <- config_file[["Volcano_nShow_geneNames"]]
betaPrior <- as.logical(config_file[["betaPrior"]])
organism <- config_file[["organism"]] # human, mouse, dmelan (fruit fly)
rm_samples <- config_file[["rm_samples"]] # samples to remove before DESeq2

if(organism == "mouse"){
    mart_dataset <- "mmusculus_gene_ensembl"
    OrgDb <- "org.Mm.eg.db"
    KEGG_organism <- 'mmu'
    
}else if(organism == "human"){
    mart_dataset <- "hsapiens_gene_ensembl"
    OrgDb <- "org.Hs.eg.db"
    KEGG_organism <- 'hsa'
        
}else if(organism == "dmelan"){
    mart_dataset <- "dmelanogaster_gene_ensembl"
    OrgDb <- "org.Dm.eg.db"
    KEGG_organism <- 'dme'
        
}else{
    mart_dataset <- config_file[["mart_dataset"]]
    OrgDb <- config_file[["OrgDb"]]
    KEGG_organism <- config_file[["KEGG_organism"]]
    
}



