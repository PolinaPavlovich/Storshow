#!/usr/bin/env Rscript
.libPaths(R.home("library"))

workdir <- snakemake@params[['workdir']]
compar_file <- snakemake@params[['compar']]
outdir <- snakemake@params[['outdir']]
log_fn <- snakemake@log[['out']]
storshow_dir <- snakemake@params[['storshow_dir']]
deseq2_dir <- snakemake@params[['deseq2_dir']]
filePrefix <- snakemake@params[['filePrefix']]

# Set up sink for logging
logfile <- file(log_fn, open = "w+")
sink(logfile, type = c("output", "message"))

if(filePrefix != ""){
    filePrefix <- paste0(filePrefix, "_")
}

source(file.path(storshow_dir, "scripts/storshow_functions.R"))

# Read file with shared output file names
fn_df <- read.delim(file.path(workdir, "DESeq/fn_geneCl.tsv"),
                    stringsAsFactors = FALSE)

# Load DESeq2 object
dds <- readRDS(file.path(deseq2_dir, fn_df$Name[fn_df$Type == "dds"]))

# Read the file with comparisons
ctr_treat <- read.delim(as.character(compar_file), stringsAsFactors = F)
ctr_cond <- ctr_treat[, 1]
treat_cond <- ctr_treat[, 2]

print("Checkpoint 1")
# Save DE genes showing the same trend
padj_df_dir <- workdir
allG_dir <- file.path(deseq2_dir, "DESeq2_results_allG_rds")
categ_df_fn <- save3cl(padj_df_dir = padj_df_dir, 
                       fn_df = fn_df, 
                       allG_dir = allG_dir, 
                       outputDir = outdir, 
                       filePrefix = filePrefix,
                       pairsToInclude = "all",
                       ctr_cond = ctr_cond, 
                       treat_cond = treat_cond)

if(all(!grepl("categ_df", fn_df$Type))){ # if categ_df is absent -> add
    fn_df <- rbind(fn_df, c("categ_df", categ_df_fn))
}else{
    fn_df$Name[fn_df$Type == "categ_df"] <- categ_df_fn
}

print("Checkpoint 2")
# Advanced gene clusters with trends "allUp" in treated, "allDown" and "so-so"
cl3_df <- read.delim(file.path(workdir, 
                               fn_df$Name[fn_df$Type == "categ_df"]),
                     stringsAsFactors = FALSE)
geneClusters_df <- read.delim(file.path(workdir, 
                                        fn_df$Name[fn_df$Type == "Gene_clusters"]),
                              stringsAsFactors = FALSE)
advanced_geneCl_fn <- advancedGeneClusters(cl3_df, geneClusters_df,
                                           outputDir = workdir, filePrefix, fn_affix = "")

if(all(!grepl("Gene_cl_trend", fn_df$Type))){
    fn_df <- rbind(fn_df, c("Gene_cl_trend", advanced_geneCl_fn))
}else{
    fn_df$Name[fn_df$Type == "Gene_cl_trend"] <- advanced_geneCl_fn
}
write.table(fn_df, file=file.path(workdir, "DESeq/fn_geneClTrend.tsv"),
            quote = FALSE, sep = "\t", row.names = F)

print(warnings())

print(sessionInfo())
# Close sink
close(logfile)
