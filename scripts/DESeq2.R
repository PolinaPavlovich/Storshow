#!/usr/bin/env Rscript
.libPaths(R.home("library"))

# Read config jsons to parse vars, for now hardcoded:
betaPrior <- TRUE

# Read parameters/set up sink
indir <- snakemake@params[['indir']]
compar_file <- snakemake@params[['compar']]
outdir <- snakemake@params[['outdir']]
#batch <- snakemake@params[['batch']]
storshow_dir <- snakemake@params[['storshow_dir']]
filePrefix <- snakemake@params[['filePrefix']]
DESeq2_padj_thr <- as.numeric(snakemake@params[['DESeq2_padj']])
DESeq2_absLog2FC_thr <- as.numeric(snakemake@params[['DESeq2_absLog2FC']])
workdir <- snakemake@params[['workdir']]

logfile <- file(snakemake@log[['out']], open="w+")
sink(logfile, type = c("output", "message"))

print(paste0("storshow_dir is ", storshow_dir))
source(file.path(storshow_dir, "scripts/storshow_functions.R"))

#print(batch)
fn_df <- read.delim(file.path(workdir, "preProc/fn_preProc.tsv"),
                    stringsAsFactors = FALSE)

counts_cells <- read.delim(file.path(indir, fn_df$Name[fn_df$Type == "Raw_count_matrix"]), 
                           stringsAsFactors = F)
rowdata <- read.delim(file.path(indir, fn_df$Name[fn_df$Type == "Rowdata"]), 
                      stringsAsFactors = F)
coldata <- read.delim(file.path(indir, fn_df$Name[fn_df$Type == "Coldata"]), 
                      stringsAsFactors = F)

# Create a DESeq2 object
dds <- DESeq2::DESeqDataSetFromMatrix(counts_cells, coldata, 
                                     design = ~ state, rowData = rowdata)
dds <- DESeq2::DESeq(dds, betaPrior = betaPrior)
vsd <- DESeq2::varianceStabilizingTransformation(dds)

if(filePrefix != ""){
    filePrefix <- paste0(filePrefix, "_")
}

dds_fn <- paste0(filePrefix, "dds_", ncol(counts_cells), "s_", 
                 nrow(counts_cells),"g.rds")

saveRDS(dds, file = file.path(outdir, dds_fn))

fn_df$Name[fn_df$Type == "dds"] <- dds_fn

# Read the file with comparisons
ctr_treat <- read.delim(as.character(compar_file), stringsAsFactors = F)
ctr_cond <- ctr_treat[, 1]
treat_cond <- ctr_treat[, 2]

outputDirDE <- file.path(outdir, "DESeq2_results_sign_DE")
outputDirRDS <- file.path(outdir, "DESeq2_results_allG_rds")

res <- list()

for(i in 1:nrow(ctr_treat)){
    # the 1st column is ctr, the 2nd column is treated condition
    cat(paste0("treat", treat_cond[i], " vs ctr ", ctr_cond[i], "\n"))
    res[[i]] <- compare2cond(dds, treat_cond = treat_cond[i], ctr_cond = ctr_cond[i], 
                             samp = filePrefix,
                             outputDirDE = outputDirDE, outputDirRDS = outputDirRDS,
                             DESeq2_padj_thr = DESeq2_padj_thr, 
                             DESeq2_absLog2FC_thr = DESeq2_absLog2FC_thr)
}

# Save list of DE genes for DAVID
dir.create(file.path(workdir, "DESeq/for_DAVID"), showWarnings = FALSE)
inputDir <- file.path(workdir, "DESeq/DESeq2_results_sign_DE")

fns <- reorderFns(list.files(inputDir), 
                  ctr_cond = ctr_cond, 
                  treat_cond = treat_cond)
fns_prefix <- gsub(".tsv", "", fns)

for(i in 1:length(fns)){
    de_genes <- read.delim(file.path(inputDir, fns[i]), stringsAsFactors = FALSE)
    # UP
    de_genes_up <- de_genes[de_genes$log2FoldChange > 0, ]
    write.table(as.data.frame(de_genes_up$Gene), 
                file=file.path(file.path(outdir, "for_DAVID"), paste0(fns_prefix[i], "_", 
                                              nrow(de_genes_up), "g_UP.tsv")),
                quote = FALSE, sep = "\t", row.names = F, col.names = F)
    # DOWN
    de_genes_down <- de_genes[de_genes$log2FoldChange < 0, ]
    write.table(as.data.frame(de_genes_down$Gene), 
                file=file.path(file.path(outdir, "for_DAVID"), paste0(fns_prefix[i], "_", 
                                              nrow(de_genes_down), "g_DOWN.tsv")),
                quote = FALSE, sep = "\t", row.names = F, col.names = F)
}

# Save ranked file for GSEA
dir.create(file.path(workdir, "DESeq/for_GSEA_preranked"), showWarnings = F)
for(i in 1:nrow(ctr_treat)){
    DESeq_result <- res[[i]][[1]]
    deseq_res_df = as.data.frame(DESeq_result)
    
    for_GSEA = data.frame(Gene = rownames(deseq_res_df), log2FC = deseq_res_df$log2FoldChange)
    for_GSEA = for_GSEA[order(for_GSEA$log2FC, decreasing = TRUE),]
    for_GSEA = for_GSEA[!is.na(for_GSEA$log2FC),]
    for_GSEA = for_GSEA[for_GSEA$Gene != "",]
    for_GSEA_fn <- paste0(filePrefix, "contrast_", treat_cond[i], "_", ctr_cond[i],
                          "_posIsUpInTreated-firstInRank.rnk")
    
    write.table(for_GSEA, file.path(workdir, "DESeq/for_GSEA_preranked", for_GSEA_fn),
                quote = FALSE, sep = "\t", row.names = FALSE)
}

# Save normalized counts
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
norm_counts <- round(norm_counts, 1)

norm_counts_fn <- paste0(filePrefix, "norm_counts_", nrow(norm_counts), "g_", 
                         ncol(norm_counts), "s_avg_cond_thr100_SYMBOL.tsv")

write.table(norm_counts, file.path(outdir, norm_counts_fn),
            quote = FALSE, sep = "\t")

fn_df$Name[fn_df$Type == "Normalized_counts"] <- norm_counts_fn

# If Normalized_counts output filename is present in this file,
# then the rule is completed successfully
write.table(fn_df, file=file.path(workdir, "DESeq/fn_deSeq2.tsv"),
            quote = FALSE, sep = "\t", row.names = F)

print(warnings())
print(sessionInfo())

# Close sink
close(logfile)
