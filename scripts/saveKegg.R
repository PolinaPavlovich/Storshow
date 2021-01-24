#!/usr/bin/env Rscript
.libPaths(R.home("library"))

workdir <- snakemake@params[['workdir']]
compar_file <- snakemake@params[['compar']]
outdir <- snakemake@params[['outdir']]
log_fn <- snakemake@log[['out']]
storshow_dir <- snakemake@params[['storshow_dir']]
organism <- snakemake@params[['organism']]
deseq2_dir <- snakemake@params[['deseq2_dir']]
ORA_padj_thr <- as.numeric(snakemake@params[['ORA_padj']])

# Set up sink for logging
logfile <- file(log_fn, open = "w+")
sink(logfile, type = c("output", "message"))

if(organism == "mouse"){
    OrgDb <- "org.Mm.eg.db"
    KEGG_organism <- 'mmu'
    
}else if(organism == "human"){
    OrgDb <- "org.Hs.eg.db"
    KEGG_organism <- 'hsa'
    
}else if(organism == "dmelan"){
    OrgDb <- "org.Dm.eg.db"
    KEGG_organism <- 'dme'
    
}else{
    stop("Unknown organism.")
}

source(file.path(storshow_dir, "scripts/storshow_functions.R"))

fn_df <- read.delim(file.path(workdir, "DESeq/fn_deSeq2.tsv"),
                    stringsAsFactors = FALSE)

# Load DESeq2 object
dds <- readRDS(file.path(deseq2_dir, fn_df$Name[fn_df$Type == "dds"]))

# Read the file with comparisons
ctr_treat <- read.delim(as.character(compar_file), stringsAsFactors = F)

resDEdir <- file.path(deseq2_dir, "DESeq2_results_sign_DE")

saveKegg(dds, resDEdir, outdir,
         KEGG_organism = KEGG_organism, OrgDb = OrgDb, 
         padjCutoff = ORA_padj_thr,
         ctr_cond = ctr_treat[, 1], treat_cond = ctr_treat[, 2])

outTable <- data.frame(successMessage = 
                           paste0("saveKegg is completed successfully.\n",
                                  "The output is in ", outdir))
write.table(outTable, file=file.path(workdir, "KEGG/saveKegg.tsv"),
            quote = FALSE, sep = "\t", row.names = F, col.names = F)

print(sessionInfo())

# Close sink
close(logfile)
