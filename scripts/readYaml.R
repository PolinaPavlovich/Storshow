require(yaml)

config_file <- yaml.load_file(file.path(workdir, "defaults.yaml"))

# Read params from yaml
filePrefix <- config_file[["fPrefix"]]
geneFilt <- config_file[["geneFilt"]]
#compar_file <- config_file[['compar']]
DESeq2_padj_thr <- config_file[["DESeq2_padj"]]
DESeq2_absLog2FC_thr <- config_file[["DESeq2_absLog2FC"]]

ORA_padj_thr <- config_file[["ORA_padj"]]
GSEA_padj_thr <- config_file[["GSEA_padj"]]
avg_pair_overl_thr <- config_file[["avg_pair_overl"]]
#style_file <- config_file[["style"]]

#betaPrior <- as.logical(config_file[["betaPrior"]])
betaPrior <- TRUE
organism <- config_file[["organism"]] # human, mouse, dmelan (fruit fly)
rm_samples <- config_file[["rm_samples"]] # samples to remove before DESeq2

# Set params that can be changed inside the app
MA_plots_ylim <- c(-3,3)

first_gene_to_show <- "Cyp26b1"
first_custom_GSEA_to_show <- "clus31_MPPs"
first_allPath_GSEA_to_show <- "5991454_M_Phase"

extra_data_comment <- ""

front_target_barplot_height <- "auto"
front_target_barplot_width <- "auto"
front_extra_barplot_height <- "280px"
front_extra_barplot_width <- "600px"
showCateg <- 25
showCategCnet <- 6
maxNgenes <- 7
DE_df_height <- "400px"
MA_plots_width <- "600px"
Volcano_plots_width <- "600px"
Volcano_nShow_geneNames <- 25

# Parse some parameters
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

if(filePrefix != ""){
    filePrefix <- paste0(filePrefix, "_")
}

ctr_treat <- read.delim(as.character(compar_file), stringsAsFactors = F)
ctr_cond <- ctr_treat[, 1]
treat_cond <- ctr_treat[, 2]