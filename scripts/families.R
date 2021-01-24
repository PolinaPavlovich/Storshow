#!/usr/bin/env Rscript
.libPaths(R.home("library"))

workdir <- snakemake@params[['workdir']]
compar_file <- snakemake@params[['compar']]
outdir <- snakemake@params[['outdir']]
log_fn <- snakemake@log[['out']]
storshow_dir <- snakemake@params[['storshow_dir']]
deseq2_dir <- snakemake@params[['deseq2_dir']]
avg_pair_overl_thr <- as.numeric(snakemake@params[['avg_pair_overl']])
ORA_padj_thr <- as.numeric(snakemake@params[['ORA_padj']])

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
ctr_cond <- ctr_treat[, 1]
treat_cond <- ctr_treat[, 2]

# Save families
## Create output folders
createFamilyDirectories(workdir)

# Upload FGSEA custom
inputDir <- file.path(workdir, "FGSEA/custom")
fns <- reorderFns(list.files(inputDir, pattern = "fgsea_custom.rds"), 
                  ctr_cond = ctr_cond, treat_cond = treat_cond)

fgsea_list <- lapply(fns, function(x) readRDS(file.path(inputDir, x)))

# Save families FGSEA Custom
outputDir_up <- file.path(workdir, "families/FGSEA/custom/up")
outputDir_down <- file.path(workdir, "families/FGSEA/custom/down")
outfile_prefix <- gsub(".rds", "", fns, perl = FALSE)

suppressMessages(library(dbscan))

## Save families FGSEA Custom UP
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_up)
for(i in 1:length(fgsea_list)){
    fgsea_res <- fgsea_list[[i]]
    fgsea_pos <- fgsea_res[fgsea_res$NES > 0, ]
    
    saveFGSEAFamilies(fgsea_input = fgsea_pos, 
                      outputDir = outputDir_up, 
                      outfile_prefix = outfile_prefix[i],
                      res_fn = "families_df_up.rds",
                      avg_pair_overl_thr = avg_pair_overl_thr)
}

## Save families FGSEA Custom DOWN
setFamDir(outputDir_down)
for(i in 1:length(fgsea_list)){
    fgsea_res <- fgsea_list[[i]]
    fgsea_neg <- fgsea_res[fgsea_res$NES < 0, ]
    
    saveFGSEAFamilies(fgsea_neg, outputDir_down, 
                      outfile_prefix[i],
                      res_fn = "families_df_down.rds",
                      avg_pair_overl_thr = avg_pair_overl_thr)
}

# Upload FGSEA allPath
inputDir <- file.path(workdir, "FGSEA/allPath")
fns <- reorderFns(list.files(inputDir, pattern = "fgsea_allPath.rds"), 
                  ctr_cond = ctr_cond, treat_cond = treat_cond)

fgsea_list <- lapply(fns, function(x) readRDS(file.path(inputDir, x)))

# Save families FGSEA allPath
outputDir_up <- file.path(workdir, "families/FGSEA/allPath/up")
outputDir_down <- file.path(workdir, "families/FGSEA/allPath/down")
outfile_prefix <- gsub(".rds", "", fns, perl = FALSE)

## Save families FGSEA allPath UP
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_up)
for(i in 1:length(fgsea_list)){
    fgsea_res <- fgsea_list[[i]]
    fgsea_pos <- fgsea_res[fgsea_res$NES > 0, ]
    
    saveFGSEAFamilies(fgsea_input = fgsea_pos, 
                      outputDir = outputDir_up, 
                      outfile_prefix = outfile_prefix[i],
                      res_fn = "families_df_up.rds",
                      avg_pair_overl_thr = avg_pair_overl_thr)
}

## Save families FGSEA allPath DOWN
setFamDir(outputDir_down)
for(i in 1:length(fgsea_list)){
    fgsea_res <- fgsea_list[[i]]
    fgsea_neg <- fgsea_res[fgsea_res$NES < 0, ]
    
    saveFGSEAFamilies(fgsea_neg, outputDir_down, 
                      outfile_prefix[i],
                      res_fn = "families_df_down.rds",
                      avg_pair_overl_thr = avg_pair_overl_thr)
}

# Upload BP
inputDir <- file.path(workdir, "enrichGO_simpl/BP")

fns_pos <- reorderFns(list.files(inputDir, pattern = "pos.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_pos <- lapply(fns_pos, function(x) readRDS(file.path(inputDir, x)))

fns_neg <- reorderFns(list.files(inputDir, pattern = "neg.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_neg <- lapply(fns_neg, function(x) readRDS(file.path(inputDir, x)))

# Save families BP
outputDir_up <- file.path(workdir, "families/BP/up")
outputDir_down <- file.path(workdir, "families/BP/down")
outfile_prefix_up <- gsub(".rds", "", fns_pos, perl = FALSE)
outfile_prefix_down <- gsub(".rds", "", fns_neg, perl = FALSE)

## Save families BP UP
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_up)
for(i in 1:length(ora_pos)){
    cat(i, " ")
    saveORAfamilies(ora_input = ora_pos[[i]], 
                    outputDir = outputDir_up, 
                    outfile_prefix = outfile_prefix_up[i],
                    res_fn = "families_df_up.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

## Save families BP DOWN
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_down)
for(i in 1:length(ora_neg)){
    cat(i, " ")
    saveORAfamilies(ora_neg[[i]], outputDir_down, 
                    outfile_prefix_down[i],
                    res_fn = "families_df_down.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

# Upload MF
inputDir <- file.path(workdir, "enrichGO_simpl/MF")

fns_pos <- reorderFns(list.files(inputDir, pattern = "pos.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_pos <- lapply(fns_pos, function(x) readRDS(file.path(inputDir, x)))

fns_neg <- reorderFns(list.files(inputDir, pattern = "neg.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_neg <- lapply(fns_neg, function(x) readRDS(file.path(inputDir, x)))

# Save families MF
outputDir_up <- file.path(workdir, "families/MF/up")
outputDir_down <- file.path(workdir, "families/MF/down")
outfile_prefix_up <- gsub(".rds", "", fns_pos, perl = FALSE)
outfile_prefix_down <- gsub(".rds", "", fns_neg, perl = FALSE)

## Save families MF UP
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_up)
for(i in 1:length(ora_pos)){
    cat(i, " ")
    saveORAfamilies(ora_pos[[i]], outputDir_up, 
                    outfile_prefix_up[i],
                    res_fn = "families_df_up.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

## Save families MF DOWN
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_down)
for(i in 1:length(ora_neg)){
    cat(i, " ")
    saveORAfamilies(ora_neg[[i]], outputDir_down, 
                    outfile_prefix_down[i],
                    res_fn = "families_df_down.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

# Upload CC
inputDir <- file.path(workdir, "enrichGO_simpl/CC")

fns_pos <- reorderFns(list.files(inputDir, pattern = "pos.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_pos <- lapply(fns_pos, function(x) readRDS(file.path(inputDir, x)))

fns_neg <- reorderFns(list.files(inputDir, pattern = "neg.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_neg <- lapply(fns_neg, function(x) readRDS(file.path(inputDir, x)))

# Save families CC
outputDir_up <- file.path(workdir, "families/CC/up")
outputDir_down <- file.path(workdir, "families/CC/down")
outfile_prefix_up <- gsub(".rds", "", fns_pos, perl = FALSE)
outfile_prefix_down <- gsub(".rds", "", fns_neg, perl = FALSE)

## Save families CC UP
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_up)
for(i in 1:length(ora_pos)){
    cat(i, " ")
    saveORAfamilies(ora_pos[[i]], outputDir_up, 
                    outfile_prefix_up[i],
                    res_fn = "families_df_up.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

## Save families CC DOWN
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_down)
for(i in 1:length(ora_neg)){
    cat(i, " ")
    saveORAfamilies(ora_neg[[i]], outputDir_down, 
                    outfile_prefix_down[i],
                    res_fn = "families_df_down.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

# Upload KEGG
inputDir <- file.path(workdir, "KEGG")

fns_pos <- reorderFns(list.files(inputDir, pattern = "pos.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_pos <- lapply(fns_pos, function(x) readRDS(file.path(inputDir, x)))

fns_neg <- reorderFns(list.files(inputDir, pattern = "neg.rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
ora_neg <- lapply(fns_neg, function(x) readRDS(file.path(inputDir, x)))

# Save families KEGG
outputDir_up <- file.path(workdir, "families/KEGG/up")
outputDir_down <- file.path(workdir, "families/KEGG/down")
outfile_prefix_up <- gsub(".rds", "", fns_pos, perl = FALSE)
outfile_prefix_down <- gsub(".rds", "", fns_neg, perl = FALSE)

## Save families KEGG UP
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_up)
for(i in 1:length(ora_pos)){
    cat(i, " ")
    saveORAfamilies(ora_pos[[i]], outputDir_up, 
                    outfile_prefix_up[i],
                    res_fn = "families_df_up.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

## Save families KEGG DOWN
# Create output folders if did not exist
# and delete old files there there were any
setFamDir(outputDir_down)
for(i in 1:length(ora_neg)){
    cat(i, " ")
    saveORAfamilies(ora_neg[[i]], outputDir_down, 
                    outfile_prefix_down[i],
                    res_fn = "families_df_down.rds",
                    avg_pair_overl_thr = avg_pair_overl_thr,
                    ORA_padj_thr = ORA_padj_thr)
}

outTable <- data.frame(successMessage = 
                           paste0("saveFamilies is completed successfully.\n",
                                  "The output is in ", outdir))
write.table(outTable, file=file.path(workdir, "FGSEA/families.tsv"),
            quote = FALSE, sep = "\t", row.names = F, col.names = F)

print(warnings())

print(sessionInfo())
# Close sink
close(logfile)
