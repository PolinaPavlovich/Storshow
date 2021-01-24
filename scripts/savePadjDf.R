#!/usr/bin/env Rscript
.libPaths(R.home("library"))

suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel, warn.conflicts = F))
suppressMessages(require(doParallel))

storshow_dir <- snakemake@params[['storshow_dir']]
compar_file <- snakemake@params[['compar']]
deseq2_dir <- snakemake@params[['deseq2_dir']]
#filePrefix <- snakemake@params[['filePrefix']]
outdir <- snakemake@params[['outdir']]

logfile <- file(snakemake@log[['out']], open="w+")
sink(logfile, type = c("output", "message"))

print(paste0("storshow_dir is ", storshow_dir))
source(file.path(storshow_dir, "scripts/storshow_functions.R"))

# Load DESeq2 object
dds <- readRDS(file.path(deseq2_dir, "DDSobject.rds"))

# Read the file with comparisons
ctr_treat <- read.delim(as.character(compar_file), stringsAsFactors = F)

# Save padj_df
inputDir <- file.path(deseq2_dir, "DESeq2_results_sign_DE")

# the 1st column is ctr, the 2nd column is treated condition
fns <- reorderFns(list.files(inputDir), ctr_cond = ctr_treat[, 1], treat_cond = ctr_treat[, 2])
results_list <- vector("list", length = length(fns))
for(i in 1:length(fns)){
    results_list[[i]] <- read.delim(file.path(inputDir, fns[i]), stringsAsFactors = FALSE)
}

all_deG <- c()
for(i in 1:length(fns)){
    all_deG <- c(all_deG, results_list[[i]]$Gene)
}
all_deG <- unique(all_deG)
padj_df <- data.frame(Gene = all_deG)

for(i in 1:length(fns)){
    padj_df <- merge(padj_df, results_list[[i]][, c("padj", "Gene")],
                     by.x="Gene", by.y="Gene",
                     all.x=TRUE, all.y=TRUE, sort=FALSE)
    colnames(padj_df)[ncol(padj_df)] <- gsub("_sign_.*tsv", "", fns[i], perl = TRUE)
    colnames(padj_df)[ncol(padj_df)] <- gsub(".*_DESeq2_res_", "", 
                                             colnames(padj_df)[ncol(padj_df)], perl = TRUE)
}

#f(filePrefix != ""){
#    padj_df_fn <- paste0(samp, "_padj_df.tsv")
#}else{
#    padj_df_fn <- paste0("padj_df.tsv")
#}

# Would be cool to add filePrefix
padj_df_fn <- paste0("padj_df.tsv")

dir.create(outdir, showWarnings = F)
write.table(padj_df, file.path(outdir, padj_df_fn),
            quote = FALSE, row.names = FALSE, sep = "\t")

print(sessionInfo())

# Close sink
close(logfile)


