#!/usr/bin/env Rscript
.libPaths(R.home("library"))

workdir <- snakemake@params[['workdir']]
compar_file <- snakemake@params[['compar']]
outdir <- snakemake@params[['outdir']]
log_fn <- snakemake@log[['out']]
storshow_dir <- snakemake@params[['storshow_dir']]
deseq2_dir <- snakemake@params[['deseq2_dir']]
GSEA_padj_thr <- as.numeric(snakemake@params[['GSEA_padj']])

# Set up sink for logging
logfile <- file(log_fn, open = "w+")
sink(logfile, type = c("output", "message"))

source(file.path(storshow_dir, "scripts/storshow_functions.R"))

fn_df <- read.delim(file.path(workdir, "DESeq/fn_deSeq2.tsv"),
                    stringsAsFactors = FALSE)

# Load DESeq2 object
dds <- readRDS(file.path(deseq2_dir, fn_df$Name[fn_df$Type == "dds"]))

# Read the file with comparisons
ctr_treat <- read.delim(as.character(compar_file), stringsAsFactors = F)

resRDSdir <- file.path(deseq2_dir, "DESeq2_results_allG_rds")

suppressMessages(require(fgsea))
gmt.file <- file.path(storshow_dir, "scripts/signatures",
                      "Signatures_for_GSEA_updated_20200714.gmt")
pathways <- gmtPathways(gmt.file)

saveFgsea(dds, pathways, "Symbol", 
          resRDSdir,
          outdir, 
          fnPart = "custom",
          GSEA_padj_thr = GSEA_padj_thr,
          ctr_cond = ctr_treat[, 1], 
          treat_cond = ctr_treat[, 2],
          organism = organism)

outTable <- data.frame(successMessage = 
                           paste0("saveFgseaCustom is completed successfully.\n",
                                  "The output is in ", outdir))
write.table(outTable, file=file.path(workdir, "FGSEA/saveFgseaCustom.tsv"),
            quote = FALSE, sep = "\t", row.names = F, col.names = F)

print(warnings())

print(sessionInfo())
# Close sink
close(logfile)
