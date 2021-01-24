#!/usr/bin/env Rscript
.libPaths(R.home("library"))

# Read up snake params
cmatloc <- snakemake@params[['inc']]
col_fn <- snakemake@params[['coldata']]
outdir <- snakemake@params[['outdir']]
organism <- snakemake@params[['organism']]
rm_samples <- snakemake@params[['rm_samples']]
geneFilt <- snakemake@params[['geneFilt']]
log_fn <- snakemake@log[['out']]
filePrefix <- snakemake@params[['filePrefix']]
workdir <- snakemake@params[['workdir']]

# Set up sink for logging
logfile <- file(log_fn, open="w+")
sink(logfile, type = c("output", "message"))

rm_samples <- unlist(strsplit(rm_samples, ","))

# Sanity check for geneFilt format
## $ means that it should end with a digit
if(!grepl("avgCond_\\d+$", geneFilt) & 
   !grepl("total_\\d+$", geneFilt)){
    
    messagePrint <- paste0("Wrong geneFilt format. Check help info.")
    cat(paste0("Error: ", messagePrint), file = logfile)
    stop(messagePrint)
}

# Parge geneFilt
geneFilt <- unlist(strsplit(geneFilt, "_"))

# Parse organism
if(organism == "mouse"){
    mart_dataset <- "mmusculus_gene_ensembl"
    
}else if(organism == "human"){
    mart_dataset <- "hsapiens_gene_ensembl"
    
}else if(organism == "dmelan"){
    mart_dataset <- "dmelanogaster_gene_ensembl"
    
}else{
    stop("Unknown organism.")
    
}

# Debug info
print(rm_samples)
print(cmatloc)
print(col_fn)
print(organism)

# Read Data
counts_cells <- read.delim(as.character(cmatloc), stringsAsFactors = FALSE)

# Sanity check for a correct organism (only for mouse so far)
n_mouseEns <- length(counts_cells[, "X"][grepl("ENSMUSG", counts_cells[, "X"])])
if(organism == "mouse"){
    if(n_mouseEns < nrow(counts_cells)){
        messagePrint <- paste0("Only ", n_mouseEns, " out of ", nrow(counts_cells), 
                               " feature names start with mouse ENSMUSG Ensembl ID.")
        cat(paste0("Error: ", messagePrint), file = logfile)
        stop(messagePrint)
    }
}

if((organism != "mouse") & n_mouseEns > 0){
    messagePrint <- paste0("The specified organism is ", organism, ", but ", 
                           n_mouseEns, " out of ", nrow(counts_cells), 
                           " feature names start with mouse ENSMUSG Ensembl ID.")
    cat(paste0("Error: ", messagePrint), file = logfile)
    stop(messagePrint)
}

# Create gene database
ensembl <- biomaRt::useMart(biomart = "ensembl", dataset = mart_dataset)
database <- biomaRt::getBM(attributes=c("ensembl_gene_id", 
                               "external_gene_name",
                               "description",
                               "gene_biotype",
                               "entrezgene_id",
                               "chromosome_name"),  mart=ensembl)

colnames(database)[colnames(database) == "ensembl_gene_id"] <- "Ensembl"
colnames(database)[colnames(database) == "external_gene_name"] <- "Symbol"
colnames(database)[colnames(database) == "entrezgene_id"] <- "Entrez.Gene.ID"
colnames(database)[colnames(database) == "description"] <- "Name"
colnames(database)[colnames(database) == "gene_biotype"] <- "Feature.Type"

database <- database[!duplicated(database$Ensembl), ]

rownames(counts_cells) <- counts_cells[, "X"]
counts_cells <- counts_cells[, colnames(counts_cells) != "X"]
counts_cells <- as.matrix(counts_cells)

# Coldata
coldata <- read.delim(as.character(col_fn), stringsAsFactors = FALSE)

# Rowdata
rowdata <- data.frame(fullGeneID = rownames(counts_cells))
rowdata$fullGeneID <- gsub("[.]", "_", rowdata$fullGeneID, perl = TRUE)
rowdata$ENSEMBL <- do.call(rbind, strsplit(as.character(rowdata$fullGeneID), "_"))[,1]
rowdata <- merge(rowdata, database,
                by.x = "ENSEMBL", by.y = "Ensembl",
                all.x = TRUE, all.y = FALSE, sort = FALSE)
# Optional remove, needs arguments.

rownames(counts_cells) <- gsub("[.]", "_", rownames(counts_cells), perl = TRUE)

# Removing samples is needed
counts_cells <- counts_cells[ , !(colnames(counts_cells) %in% rm_samples)]
coldata <- coldata[!(coldata$sampleName %in% rm_samples), ]

# Select which genes to keep
if(geneFilt[1] == "avgCond"){
    conditions <- unique(coldata$state)
    
    counts_split <- lapply(conditions, function(cond){
        counts_cells[, grepl(cond, colnames(counts_cells))]
    })
    
    names(counts_split) <- conditions
    
    rowSums_conds <- lapply(counts_split, function(df){
        rowSums(df) / ncol(df)
    })
    
    rowSums_conds_df <- do.call(cbind, rowSums_conds)
    stopifnot(ncol(rowSums_conds_df) == length(conditions))
    max_avg_cond <- apply(rowSums_conds_df, 1, max)
    sel_genes <- names(max_avg_cond[max_avg_cond > as.numeric(geneFilt[2])])
    
    max_avg_cond <- as.data.frame(max_avg_cond)
    max_avg_cond$fullGeneID = rownames(max_avg_cond)
    
    rowdata <- merge(rowdata, max_avg_cond,
                     by.x = "fullGeneID", by.y = "fullGeneID",
                     all.x = TRUE, all.y = TRUE, sort = FALSE)
    
}else if(geneFilt[1] == "total"){
    readsSums <- data.frame(fullGeneID = rownames(counts_cells),
                           readsSum = rowSums(counts_cells))
    
    sel_genes <- readsSums$fullGeneID[readsSums$readsSum > as.numeric(geneFilt[2])]
    
    rowdata <- merge(rowdata, readsSums,
                     by.x = "fullGeneID", by.y = "fullGeneID",
                     all.x = TRUE, all.y = TRUE, sort = FALSE)
}

# Filter genes in counts_cells and rowdata
counts_cells <- counts_cells[order(rownames(counts_cells)),]
rowdata <- rowdata[order(rowdata$fullGeneID),]
stopifnot(all(rowdata$fullGeneID == rownames(counts_cells)))
counts_cells <- counts_cells[rownames(counts_cells) %in% sel_genes, ]
rowdata <- rowdata[rowdata$fullGeneID %in% sel_genes, ]

# Substitute Ensembl IDs with Symbols where possible
rownames(counts_cells) = do.call(rbind, strsplit(as.character(rownames(counts_cells)), "_"))[,1]
counts_cells <- counts_cells[order(rownames(counts_cells)),]
rowdata <- rowdata[order(rowdata$ENSEMBL),]
stopifnot(all(rownames(counts_cells) == rowdata$ENSEMBL))

rownames(counts_cells)[!is.na(rowdata$Symbol) & 
                           !(duplicated(rowdata$Symbol)) & 
                           (rowdata$Symbol != "")] <- 
    rowdata$Symbol[!is.na(rowdata$Symbol) & 
                           !(duplicated(rowdata$Symbol)) & 
                           (rowdata$Symbol != "")]

# Write out
if(filePrefix != ""){
    filePrefix <- paste0(filePrefix, "_")
}

rowdata_fn <- paste0(filePrefix, "rowdata_", nrow(rowdata), "g.tsv")
coldata_fn <- paste0(filePrefix, "coldata_", nrow(coldata), "s.tsv")
counts_cells_fn <- paste0(filePrefix, "counts_cells_", 
                          nrow(counts_cells), "g_", ncol(counts_cells), 
                          "s_", geneFilt[1], "_", as.numeric(geneFilt[2]), "_SYMBOL.tsv")

#if(organism == "human"){
#    rownames(counts_cells) <- stringr::str_to_title(rownames(counts_cells))
#    rowdata$Symbol <- stringr::str_to_title(rowdata$Symbol)
#}

dir.create(outdir, showWarnings = F)
write.table(rowdata, file=file.path(outdir, rowdata_fn),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(coldata, file=file.path(outdir, coldata_fn),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(counts_cells, file=file.path(outdir, counts_cells_fn),
            quote = FALSE, sep = "\t")

# Save output file names
fn_types <- c("Raw_count_matrix",
              "Coldata",
              "Rowdata",
              "dds",
              "Normalized_counts",
              "padj_df",
              "Gene_clusters")
fn_df <- data.frame(Type = fn_types,
                    Name = c(counts_cells_fn,
                             coldata_fn,
                             rowdata_fn,
                             rep(NA, length(fn_types) - 3)))

write.table(fn_df, file=file.path(workdir, "preProc/fn_preProc.tsv"),
            quote = FALSE, sep = "\t", row.names = F)

print(warnings())
print(sessionInfo())
# Close sink
close(logfile)

