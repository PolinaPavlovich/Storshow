# deseq2 rule
compare2cond <- function(dds, treat_cond, ctr_cond, samp, 
                         outputDirDE, outputDirRDS, 
                         DESeq2_padj_thr = 0.1, DESeq2_absLog2FC_thr = 0.5,
                         decreasingOrder = FALSE,
                         nameAddColumn = "Gene",
                         baseMean_thr = -0.1,
                         featureFn = "g"){
    (DESeq_result <- DESeq2::results(dds, contrast = c("state", treat_cond, ctr_cond)))
    deseq_res_df <- as.data.frame(DESeq_result)
    deseq_res_df[nameAddColumn] <- rownames(deseq_res_df)
    
    if(is.na(DESeq2_padj_thr) | is.null(DESeq2_padj_thr) | (DESeq2_padj_thr > 1)){
        deseq_res_sign <- deseq_res_df[!is.na(deseq_res_df$log2FoldChange) & 
                                           (abs(deseq_res_df$log2FoldChange) > 
                                                DESeq2_absLog2FC_thr) &
                                           (deseq_res_df$baseMean > baseMean_thr),]
    }else{
        deseq_res_sign <- deseq_res_df[!is.na(deseq_res_df$padj) & 
                                           (deseq_res_df$padj < DESeq2_padj_thr) &
                                           !is.na(deseq_res_df$log2FoldChange) &
                                           (abs(deseq_res_df$log2FoldChange) > 
                                                DESeq2_absLog2FC_thr) &
                                           (deseq_res_df$baseMean > baseMean_thr),]
    }
    
    deseq_res_sign <- deseq_res_sign[order(deseq_res_sign$log2FoldChange,
                                           decreasing = decreasingOrder),]
    
    dir.create(outputDirDE, showWarnings = FALSE)
    dir.create(outputDirRDS, showWarnings = FALSE)
    
    if(samp != ""){
        fAffix <- "_DESeq2_res_"
    }else{
        fAffix <- "DESeq2_res_"
    }
    
    if(baseMean_thr == -0.1){
        deseq_res_sign_fn <- paste0(samp, fAffix, ctr_cond, "_vs_", 
                                    treat_cond,"_sign_", 
                                    nrow(deseq_res_sign), 
                                    featureFn, "_padj", DESeq2_padj_thr, "_logFC", 
                                    DESeq2_absLog2FC_thr, ".tsv")
    }else{
        deseq_res_sign_fn <- paste0(samp, fAffix, ctr_cond, "_vs_", 
                                    treat_cond,"_sign_", 
                                    nrow(deseq_res_sign), 
                                    featureFn, "_padj", DESeq2_padj_thr, "_logFC", 
                                    DESeq2_absLog2FC_thr, "_baseMean", 
                                    baseMean_thr, ".tsv")
    }
    
    write.table(deseq_res_sign, 
                file = file.path(outputDirDE, deseq_res_sign_fn),
                quote = FALSE, sep = "\t")
    write.table(as.data.frame(DESeq_result), 
                file = file.path(outputDirRDS, 
                                 paste0(samp, fAffix, ctr_cond, "_vs_",
                                        treat_cond, "_all_", 
                                        nrow(deseq_res_df), featureFn, ".tsv")),
                quote = FALSE, sep = "\t")
    
    saveRDS(DESeq_result, file.path(outputDirRDS, 
                                    paste0(samp, fAffix, ctr_cond, "_vs_",
                                           treat_cond, "_all_", 
                                           nrow(deseq_res_df), featureFn, ".rds"))) 
    return(list(DESeq_result, deseq_res_sign))
}

# everywhere starting from the padj_df rule
reorderFns <- function(fns, ctr_cond, treat_cond){
    ordered_patterns <- paste0(ctr_cond, "_vs_", treat_cond)
    reference_order <- data.frame(ref_cond = ordered_patterns,
                                  desired = c(1:length(ordered_patterns)), 
                                  real = NA)
    
    for(i in 1:nrow(reference_order)){
        #cat(i, "\n")
        reference_order$real[i] <- 
            which(grepl(paste0(reference_order$ref_cond[i]), fns))
    }
    
    fns <- fns[reference_order$real]
    return(fns)
}

# saving geneClusters
saveGeneClusters <- function(padj_df, filePrefix, outputDir){
    grid_df <- padj_df
    rownames(grid_df) <- grid_df$Gene
    grid_df <- grid_df[-1]
    grid_df[!is.na(grid_df)] <- 1
    grid_df[is.na(grid_df)] <- 0
    
    geneClusters <- plyr::count(grid_df, vars = colnames(grid_df))
    
    geneClusters$mask <- apply(geneClusters[1:ncol(grid_df)],1,paste,collapse="")
    
    geneClusters$genes <- ""
    
    stopifnot(all(colnames(grid_df) == colnames(geneClusters)[1:ncol(grid_df)]))
    for(i in 1:nrow(grid_df)){
        geneClusters$genes[geneClusters$mask == paste0(grid_df[i, ], collapse = "")] <- 
            paste0(geneClusters$genes[geneClusters$mask == paste0(grid_df[i, ], collapse = "")], 
                   rownames(grid_df)[i], sep = ",")
    }
    
    geneClusters <- geneClusters[order(geneClusters$freq, decreasing = TRUE),]
    geneClusters$genes <- gsub(",[^,]*$", "", geneClusters$genes) # delete the last comma
    geneClusters$mask <- paste0("M", geneClusters$mask)
    geneClusters <- geneClusters[order(rowSums(geneClusters[1:ncol(grid_df)]), 
                                       decreasing = TRUE), ]
    
    geneClusters_fn <- paste0(filePrefix, nrow(geneClusters), "_nonZero_geneClusters.tsv")
    write.table(geneClusters, file.path(outputDir, geneClusters_fn),
                sep = "\t", quote = FALSE, row.names = FALSE)
    return(geneClusters_fn)
}

# saving simplified enrichGO results for all comparisons
saveEgoS <- function(dds, resDEdir, 
                     outputDir, 
                     OrgDb,
                     pAdjustMethod = "BH", 
                     padjCutoff = 0.1,
                     maxCores = 10,
                     simplify_cutoff = 0.7,
                     simplify_by = "p.adjust",
                     ctr_cond, 
                     treat_cond){
    
    suppressMessages(library(clusterProfiler)) 
    suppressMessages(library(foreach))
    #require(OrgDb)
    
    giveEgoS <- function(gene, ontol, pAdjustMethod, padjCutoff){
        if(length(gene) < 3){
            library(clusterProfiler)
            result <- data.frame(ID = character(), Description = character(), 
                                 GeneRatio = character(), BgRatio = character(), 
                                 pvalue = numeric(), p.adjust = numeric(), 
                                 qvalue = numeric(), geneID = character(), 
                                 Count = numeric())
            ego_s <- new("enrichResult",
                         result         = result,
                         pvalueCutoff   = padjCutoff,
                         pAdjustMethod  = pAdjustMethod,
                         gene           = as.character(gene),
                         ontology       = ontol,
                         readable       = TRUE)
        }else{
            ego <- clusterProfiler::enrichGO(gene          = gene,
                                             OrgDb         = OrgDb,
                                             ont           = ontol,
                                             pAdjustMethod = pAdjustMethod,
                                             readable      = TRUE,
                                             pvalueCutoff  = padjCutoff)
            ego@result = ego@result[ego@result$p.adjust < padjCutoff,]
            ego@geneSets = data.frame()
            ego_s <- clusterProfiler::simplify(ego, cutoff=simplify_cutoff, 
                                               by=simplify_by, select_fun=min)
        }
        
        return(ego_s)
    }
    
    dir.create(outputDir, showWarnings = FALSE)
    
    rowdata <- as.data.frame(SummarizedExperiment::rowData(dds))
    rowdata <- rowdata[!duplicated(rowdata$Entrez.Gene.ID),]
    
    fns_DE <- reorderFns(list.files(resDEdir), 
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    
    # * 2 is because there are UP and DOWN directions for each file
    cores <- min(length(fns_DE) * 2, maxCores, 
                 parallel::detectCores(all.tests = FALSE, logical = TRUE))
    
    results_list <- vector("list", length = length(fns_DE))
    for(i in 1:length(fns_DE)){
        results_list[[i]] <- read.delim(file.path(resDEdir, fns_DE[i]), stringsAsFactors = FALSE)
    }
    filePrefix <- do.call(rbind, strsplit(fns_DE, "\\.[^\\.]*$"))[,1]
    ont_vec <- c("CC", "MF", "BP")
    
    myCluster <- parallel::makeCluster(cores, # number of cores to use
                                       type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)
    
    for(ont_i in 1:length(ont_vec)){
        #foreach::foreach(ont_i = 1:length(ont_vec)) %dopar% {
        dir.create(file.path(outputDir, ont_vec[ont_i]), showWarnings = FALSE)
        #for(res_i in 1:length(fns_DE)){
        foreach::foreach(res_i = 1:length(fns_DE)) %dopar% {
            cat(ont_vec[ont_i], " ", res_i, "up\n")
            gene <- results_list[[res_i]]$Gene[results_list[[res_i]]$log2FoldChange > 0]
            cat(length(gene), "\n")
            gene <- rowdata$Entrez.Gene.ID[rowdata$Symbol %in% gene]
            cat(length(gene), "\n")
            
            ego_s <- giveEgoS(gene, ont_vec[ont_i], pAdjustMethod, 
                              padjCutoff)
            saveRDS(ego_s, file.path(outputDir, paste0(ont_vec[ont_i], "/", 
                                                       filePrefix[res_i], "_ego_s_ont", 
                                                       ont_vec[ont_i], "_logFC_pos.rds")))
            
            cat(ont_vec[ont_i], " ", res_i, "down\n")
            gene <- results_list[[res_i]]$Gene[results_list[[res_i]]$log2FoldChange < 0]
            cat(length(gene), "\n")
            gene <- rowdata$Entrez.Gene.ID[rowdata$Symbol %in% gene]
            cat(length(gene), "\n")
            
            ego_s <- giveEgoS(gene, ont_vec[ont_i], pAdjustMethod, 
                              padjCutoff)
            saveRDS(ego_s, file.path(outputDir, paste0(ont_vec[ont_i], "/", 
                                                       filePrefix[res_i], "_ego_s_ont", 
                                                       ont_vec[ont_i], "_logFC_neg.rds")))
        }
    }
    parallel::stopCluster(myCluster)
}

# saving KEGG results for all comparisons
saveKegg <- function(dds, resDEdir, outputDir, 
                     KEGG_organism, OrgDb,
                     padjCutoff = 0.1,
                     ctr_cond, 
                     treat_cond,
                     pAdjustMethod = "BH"){
    
    suppressMessages(library(DOSE))
    suppressMessages(library(clusterProfiler))
    suppressMessages(library(SummarizedExperiment))
    
    giveKegg <- function(gene, padjCutoff, pAdjustMethod){
        if(length(gene) < 3){
            library(clusterProfiler)
            result <- data.frame(ID = character(), Description = character(), 
                                 GeneRatio = character(), BgRatio = character(), 
                                 pvalue = numeric(), p.adjust = numeric(), 
                                 qvalue = numeric(), geneID = character(), 
                                 Count = numeric())
            kk_r <- new("enrichResult",
                         result         = result,
                         pvalueCutoff   = padjCutoff,
                         pAdjustMethod  = pAdjustMethod,
                         gene           = as.character(gene),
                         readable       = TRUE)
        }else{
            kk <- clusterProfiler::enrichKEGG(gene         = gene,
                                              organism     = KEGG_organism,
                                              pvalueCutoff = padjCutoff)
            if(!is.null(kk)){
                kk_r <- setReadable(kk, OrgDb = OrgDb, keyType = "ENTREZID")
                kk_r@result = kk_r@result[!(is.na(kk_r@result$qvalue)) &
                                              kk_r@result$p.adjust < padjCutoff,]
                kk_r@geneSets = data.frame()
            }else{
                kk_r <- kk
            }
        }
        return(kk_r)
    }
    
    dir.create(outputDir, showWarnings = FALSE)
    
    rowdata <- as.data.frame(SummarizedExperiment::rowData(dds))
    rowdata <- rowdata[!duplicated(rowdata$Entrez.Gene.ID),]
    
    fns_DE <- reorderFns(list.files(resDEdir), 
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    results_list <- vector("list", length = length(fns_DE))
    for(i in 1:length(fns_DE)){
        results_list[[i]] <- read.delim(file.path(resDEdir, fns_DE[i]), stringsAsFactors = FALSE)
    }
    filePrefix <- do.call(rbind, strsplit(fns_DE, "\\.[^\\.]*$"))[,1]
    
    for(res_i in 1:length(fns_DE)){
        cat(res_i, "up\n")
        gene <- results_list[[res_i]]$Gene[results_list[[res_i]]$log2FoldChange > 0]
        cat(length(gene), "\n")
        gene <- rowdata$Entrez.Gene.ID[rowdata$Symbol %in% gene]
        cat(length(gene), "\n")
        
        kegg_res <- giveKegg(gene, padjCutoff, pAdjustMethod)
        saveRDS(kegg_res, file.path(outputDir, paste0(filePrefix[res_i], "_kegg_logFC_pos.rds")))
        
        cat(res_i, "down\n")
        gene <- results_list[[res_i]]$Gene[results_list[[res_i]]$log2FoldChange < 0]
        cat(length(gene), "\n")
        gene <- rowdata$Entrez.Gene.ID[rowdata$Symbol %in% gene]
        cat(length(gene), "\n")
        
        kegg_res <- giveKegg(gene, padjCutoff, pAdjustMethod)
        saveRDS(kegg_res, file.path(outputDir, paste0(filePrefix[res_i], "_kegg_logFC_neg.rds")))
    }
}

# saving FGSEA results for all comparisons
saveFgsea <- function(dds, pathways, geneMode, #geneMode "Symbol" or "EntrezID"
                      resRDSdir,
                      outputDir, 
                      fnPart = "custom",
                      GSEA_padj_thr = 0.1,
                      ctr_cond, 
                      treat_cond,
                      organism = "mouse"){
    
    suppressMessages(require(fgsea))
    
    dir.create(outputDir, showWarnings = FALSE)
    outputDir_sign <- file.path(outputDir, fnPart)
    
    # create output directories
    dir.create(outputDir_sign, showWarnings = FALSE)
    
    # delete old files if any
    oldFiles <- list.files(outputDir_sign)
    deletedFiles <- sapply(oldFiles, function(fileName) 
        file.remove(file.path(outputDir_sign, fileName)))
    
    rowdata <- as.data.frame(SummarizedExperiment::rowData(dds))
    rowdata <- rowdata[!duplicated(rowdata$Entrez.Gene.ID),]
    
    fns_allG <- reorderFns(list.files(resRDSdir, pattern = ".rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    filePrefix <- gsub("_all.*", "", fns_allG)
    results_allG <- vector("list", length = length(fns_allG))
    geneList <- vector("list", length = length(fns_allG))
    symbList <- vector("list", length = length(fns_allG))
    for(i in 1:length(fns_allG)){
        results_allG[[i]] <- as.data.frame(readRDS(file.path(resRDSdir, fns_allG[i])))
        results_allG[[i]]$Gene <- rownames(results_allG[[i]])
        results_allG[[i]] <- merge(results_allG[[i]], rowdata[,c("Symbol", "Entrez.Gene.ID")],
                                   by.x="Gene", by.y="Symbol",
                                   all.x=FALSE, all.y=FALSE, sort=FALSE)
        
        geneList[[i]] <- results_allG[[i]]$log2FoldChange
        names(geneList[[i]]) <- results_allG[[i]]$Entrez.Gene.ID
        geneList[[i]] <- sort(geneList[[i]], decreasing=TRUE)
        
        symbList[[i]] <- results_allG[[i]]$log2FoldChange
        names(symbList[[i]]) <- results_allG[[i]]$Gene
        symbList[[i]] <- sort(symbList[[i]], decreasing=TRUE)
        
    }
    
    giveFgsea <- function(pathways, stats, GSEA_padj_thr){
        fgseaRes <- fgseaMultilevel(pathways = pathways, 
                          stats = stats,
                          minSize=15,
                          maxSize=500,
                          eps = 0.0)
        #fgseaRes$padj2x <- pmin(1, fgseaRes$padj * nPairs)
        #fgseaRes <- fgseaRes[fgseaRes$padj2x < GSEA_padj_thr,]
        fgseaRes <- fgseaRes[fgseaRes$padj < GSEA_padj_thr,]
        fgseaRes <- fgseaRes[order(fgseaRes$NES),]
        #fgseaRes <- fgseaRes[, c(1:3,9,4:8)]
        
        #return(list(fgseaRes[fgseaRes$padj2x < GSEA_padj_thr, ],
        #            fgseaRes))
        return(fgseaRes)
    }
    
    entrezToSymbols <- function(fgseaRes){
        if(nrow(fgseaRes) > 0){
            for(i in 1:nrow(fgseaRes)){
                tmp <- fgseaRes[i,"leadingEdge"]
                tmp_df <- data.frame(EntrezID = tmp[[1]][[1]],
                                     order = 1:length(tmp[[1]][[1]]))
                tmp_df <- merge(tmp_df, rowdata[,c("Symbol", "Entrez.Gene.ID")],
                                by.x="EntrezID", by.y="Entrez.Gene.ID",
                                all.x=TRUE, all.y=FALSE, sort=FALSE)
                tmp_df <- tmp_df[order(tmp_df$order),]
                tmp[[1]][[1]] <- tmp_df$Symbol
                fgseaRes[i, "leadingEdge"] <- tmp
            }
            return(fgseaRes)
        }else{
            return(fgseaRes)
        }
    }
    
    for(res_i in 1:length(fns_allG)){
        cat(res_i, " ")
        if(geneMode == "Symbol"){
            stats <- symbList[[res_i]]
        }else if(geneMode == "EntrezID"){
            stats <- geneList[[res_i]]
        }else{
            cat("Wrong geneMode. Please choose between 'Symbol' and 'EntrezID'.")
        }
        
        
        fgseaRes <- giveFgsea(pathways, stats, GSEA_padj_thr)
        
        if(geneMode == "EntrezID"){
            fgseaRes <- entrezToSymbols(fgseaRes)
        }else{
            fgseaRes <- fgseaRes
        }
        
        fgseaRes <- fgseaRes[order(fgseaRes$NES), ]
        
        saveRDS(fgseaRes, file.path(outputDir, fnPart, paste0(filePrefix[res_i], 
                                                              "_fgsea_", fnPart, ".rds")))
    }
}

createFamilyDirectories <- function(working_dir){
    dir.create(file.path(working_dir, "families"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/FGSEA"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/BP"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/MF"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/CC"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/KEGG"), showWarnings = FALSE)
    
    dir.create(file.path(working_dir, "families/FGSEA/custom"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/FGSEA/allPath"), showWarnings = FALSE)
    
    dir.create(file.path(working_dir, "families/FGSEA/custom/up"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/FGSEA/custom/down"), showWarnings = FALSE)
    
    dir.create(file.path(working_dir, "families/FGSEA/allPath/up"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/FGSEA/allPath/down"), showWarnings = FALSE)
    
    dir.create(file.path(working_dir, "families/BP/up"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/BP/down"), showWarnings = FALSE)
    
    dir.create(file.path(working_dir, "families/MF/up"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/MF/down"), showWarnings = FALSE)
    
    dir.create(file.path(working_dir, "families/CC/up"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/CC/down"), showWarnings = FALSE)
    
    dir.create(file.path(working_dir, "families/KEGG/up"), showWarnings = FALSE)
    dir.create(file.path(working_dir, "families/KEGG/down"), showWarnings = FALSE)
}

saveFGSEAFamilies <- function(fgsea_input, outputDir, 
                              outfile_prefix,
                              res_fn = "families_df",
                              avg_pair_overl_thr = 70){ #70 recommended
    
    fgsea_input <- fgsea_input[order(abs(fgsea_input$NES), decreasing = TRUE), ]
    
    if(nrow(fgsea_input) == 0){
        colnames_result <- c("Pathways", "Avg_pval", "Avg_padj", #"Avg_padj2x",
                             "Avg_ES", "Avg_NES", "Family_size", "Overlap_size",
                             "Overlap_pct",
                             "Overlap_leadingEdge", 
                             "Add_leadingEdge")
        result <- as.data.frame(matrix(NA, ncol = length(colnames_result), 
                                       nrow = 1))
        colnames(result) <- colnames_result
    }else if(nrow(fgsea_input) == 1){
        colnames_result <- c("Pathways", "Avg_pval", "Avg_padj", #"Avg_padj2x",
                             "Avg_ES", "Avg_NES", "Family_size", "Overlap_size",
                             "Overlap_pct",
                             "Overlap_leadingEdge", 
                             "Add_leadingEdge")
        result <- as.data.frame(matrix(NA, ncol = length(colnames_result), 
                                       nrow = nrow(fgsea_input)))
        colnames(result) <- colnames_result
        result[, c(1:6)] <- fgsea_input[, c(1:6)]
        result$Family_size <- 1
        result$Overlap_size <- length(fgsea_input$leadingEdge[[1]])
        result$Overlap_pct <- 100
        result$Overlap_leadingEdge <- fgsea_input$leadingEdge
        result$Add_leadingEdge <- list(character(0))
    }else{
        input_list <- fgsea_input$leadingEdge
        names(input_list) <- paste0("sign_", 1:length(input_list))
        fgsea_input$sign_name <- names(input_list)
        
        icpt_df <- matrix(data = NA, nrow = length(input_list), ncol = length(input_list))
        rownames(icpt_df) <- names(input_list)
        colnames(icpt_df) <- names(input_list)
        
        for(row in 1:length(input_list)){
            for(col in 1:length(input_list)){
                #cat("row ", row, " col ", col, "\n")
                icpt_df[row, col] <- round(100 * 2 * length(intersect(input_list[[row]], 
                                                                      input_list[[col]])) / 
                                               sum(length(input_list[[row]]), length(input_list[[col]])), 1)
            }
        }
        
        # Find optimal epsilon automatically
        kNNdist_vec <- sort(kNNdist(icpt_df, k=1))
        #eps_opt <- (max(kNNdist_vec) + min(kNNdist_vec)) / 2
        eps_opt <- max(kNNdist_vec) / 2
        
        clusters <- dbscan(icpt_df, eps = eps_opt, minPts = 1)$cluster
        length(unique(clusters))
        clusters
        
        clus_df <- data.frame(sign = names(input_list), 
                              clusters = clusters,
                              stringsAsFactors = FALSE)
        
        icpt_df2 <- icpt_df[order(clus_df$clusters),]
        icpt_df2 <- icpt_df2[,order(clus_df$clusters)]
        
        saveSignSimilarity(icpt_df, clus_df, outputDir, outfile_prefix)
        
        simMed <- mkSimMed(as.matrix(icpt_df), clusters)
        pheatmapObject <- pheatmap::pheatmap(simMed)
        
        # Save pictures
        onefile = FALSE
        clustersNumber <- length(unique(clusters))
        
        ## Save simMed heatmap
        pdf(file=file.path(outputDir, "pictures",
                           paste(outfile_prefix, "families_intersection", 
                                 clustersNumber,"clusters.pdf", sep="_")),
            width=max(4, ncol(simMed) / 4), height=max(3.5, ncol(simMed) / 4.5), onefile=onefile)
        grid::grid.newpage()
        grid::grid.draw(pheatmapObject$gtable)
        dev.off()
        
        ## Save simMed scatterplot
        pdf(file=file.path(outputDir, "pictures", paste(outfile_prefix,
                                                        "families_intersection", clustersNumber,
                                                        "clusters_scatterplot.pdf", sep="_")),
            width=5, height=4, onefile=onefile # not changed by default
        )
        plot(diag(simMed), ylim = c(min(70, min(diag(simMed))), 100), xlab = "Signature family name", 
             ylab = "Avg intersection within a family", yaxt = "n", xaxt = "n")
        # avg_pair_overl_thr is the threshold for avg intercept
        abline(h=avg_pair_overl_thr, lty=2, col = "green") 
        axis(2, las = 2)
        axis(1, at = seq(1, length(diag(simMed)), by = 1), las = 0)
        dev.off()
        
        ## Save kNNdist plot
        pdf(file=file.path(outputDir, "pictures", paste(outfile_prefix,
                                                        "kNNdistancePlot_k1_optimalEps.pdf", sep="_")),
            width=5, height=4, onefile=onefile # not changed by default
        )
        plot(sort(kNNdist_vec), type="l", ylab=paste(1, "-NN distance", sep=""),
             xlab = "Points (sample) sorted by distance", yaxt = "n", xaxt = "n")
        abline(h = eps_opt, lty=2)
        axis(2, las = 2)
        axis(1, las = 0)
        text(1/3*length(kNNdist_vec), 1.1 * eps_opt,
             paste0("Optimal epsilon = ", round(eps_opt, 0)))
        dev.off()
        
        # Split outliers if they cluster together
        if(any(diag(simMed) < avg_pair_overl_thr)){
            outlier_clusters <- rownames(simMed)[diag(simMed) < avg_pair_overl_thr]
            clus_df$clusters_approved <- clus_df$clusters
            for(k in 1:length(outlier_clusters)){
                length_outl <- length(clus_df$clusters[clus_df$clusters == 
                                                           outlier_clusters[k]])
                clus_df$clusters_approved[clus_df$clusters == outlier_clusters[k]] <-
                    c((max(clus_df$clusters_approved) + 1) :
                          (max(clus_df$clusters_approved) + length_outl))
            }
            clus_df$clusters <- clus_df$clusters_approved
            clus_df <- clus_df[, colnames(clus_df) != "clusters_approved"]
        }
        
        # Calculate the result table
        ## Order families (clusters) by size
        clus_freq <- as.data.frame(table(clus_df$clusters))
        clus_freq <- clus_freq[order(clus_freq$Freq, decreasing = TRUE),]
        colnames(clus_freq)[1] <- "clusters"
        
        fgsea_input <- merge(fgsea_input, clus_df,
                             by.x="sign_name", by.y="sign",
                             all.x=TRUE, all.y=TRUE, sort=FALSE)
        
        colnames_result <- c("Pathways", "Avg_pval", "Avg_padj", #"Avg_padj2x",
                             "Avg_ES", "Avg_NES", "Family_size", "Overlap_size",
                             "Overlap_pct",
                             "Overlap_leadingEdge", 
                             "Add_leadingEdge")
        result <- as.data.frame(matrix(NA, ncol = length(colnames_result), 
                                       nrow = nrow(clus_freq)))
        colnames(result) <- colnames_result
        for(i in 1:nrow(clus_freq)){
            pathways_list <- as.list(fgsea_input$pathway[fgsea_input$clusters == 
                                                             clus_freq$clusters[i]])
            names(pathways_list) <- fgsea_input$sign_name[fgsea_input$clusters == 
                                                              clus_freq$clusters[i]]
            result$Pathways[i] <- list(pathways_list)
            
            result$Avg_pval[i] <- mean(fgsea_input$pval[fgsea_input$clusters == 
                                                            clus_freq$clusters[i]])
            
            result$Avg_padj[i] <- mean(fgsea_input$padj[fgsea_input$clusters ==
                                                            clus_freq$clusters[i]])
            
            #result$Avg_padj2x[i] <- mean(fgsea_input$padj2x[fgsea_input$clusters ==
            #                                                    clus_freq$clusters[i]])
            
            result$Avg_ES[i] <- mean(fgsea_input$ES[fgsea_input$clusters ==
                                                        clus_freq$clusters[i]])
            
            result$Avg_NES[i] <- mean(fgsea_input$NES[fgsea_input$clusters ==
                                                          clus_freq$clusters[i]])
            
            result$Family_size[i] <- length(pathways_list)
            
            val_list <- fgsea_input$leadingEdge[fgsea_input$clusters ==
                                                    clus_freq$clusters[i]]
            overlapGenes <- list(Reduce(intersect, val_list))
            
            result$Overlap_size[i] <- length(overlapGenes[[1]])
            
            result$Overlap_pct[i] <- round(100 * length(val_list) * result$Overlap_size[i] / 
                                               sum(lengths(val_list)), 0)
            
            result$Overlap_leadingEdge[i] <- overlapGenes
            
            add_list <- fgsea_input$leadingEdge[fgsea_input$clusters == 
                                                    clus_freq$clusters[i]]
            add_list <- lapply(add_list, function(x) x[!(x %in% result$Overlap_leadingEdge[i][[1]])])
            names(add_list) <- fgsea_input$sign_name[fgsea_input$clusters == 
                                                         clus_freq$clusters[i]]
            result$Add_leadingEdge[i] <- list(add_list)
        }
    }
    result <- result[order(abs(result$Avg_NES), decreasing = TRUE), ]
    
    ## Save result table
    saveRDS(result, file.path(outputDir, "tables", paste0(outfile_prefix, "_", res_fn)))
}

saveORAfamilies <- function(ora_input, outputDir, 
                            outfile_prefix,
                            res_fn = "families_df",
                            avg_pair_overl_thr = 70, #70 recommended
                            ORA_padj_thr = 0.01){ 
    
    if(is.null(ora_input)){
        colnames_result <- c("Descriptions", "Avg_pval", 
                             "Avg_padj", 
                             "Avg_qvalue", "Family_size", "Overlap_size",
                             "Overlap_pct",
                             "Overlap_geneID", 
                             "Add_geneID", "IDs")
        result <- as.data.frame(matrix(NA, ncol = length(colnames_result), nrow = 1))
        colnames(result) <- colnames_result
    }else{
        ora_input <- ora_input@result
        ora_input <- ora_input[!is.na(ora_input$Description), ]
        ora_input <- ora_input[!is.na(ora_input$p.adjust) & 
                                   ora_input$p.adjust < ORA_padj_thr, ]
        ora_input <- ora_input[order(ora_input$p.adjust), ]
        
        if(nrow(ora_input) == 0){
            colnames_result <- c("Descriptions", "Avg_pval", 
                                 "Avg_padj", 
                                 "Avg_qvalue", "Family_size", "Overlap_size",
                                 "Overlap_pct",
                                 "Overlap_geneID", 
                                 "Add_geneID", "IDs")
            result <- as.data.frame(matrix(NA, ncol = length(colnames_result), 
                                           nrow = 1))
            colnames(result) <- colnames_result
        }else if(nrow(ora_input) == 1){
            colnames_result <- c("Descriptions", "Avg_pval", 
                                 "Avg_padj", 
                                 "Avg_qvalue", "Family_size", "Overlap_size",
                                 "Overlap_pct",
                                 "Overlap_geneID", 
                                 "Add_geneID", "IDs")
            result <- as.data.frame(matrix(NA, ncol = length(colnames_result), 
                                           nrow = nrow(ora_input)))
            colnames(result) <- colnames_result
            result[, c("Descriptions", "Avg_pval",
                       "Avg_padj", "Avg_qvalue", "IDs")] <- 
                ora_input[, c("Description", "pvalue", "p.adjust", "qvalue", "ID")]
            #result[, c(1:4, 10)] <- ora_input[, c(2, 5, 6, 7, 1)]
            result$Family_size <- 1
            result$Overlap_size <- length(strsplit(ora_input$geneID[[1]], "/"))
            result$Overlap_pct <- 100
            result$Overlap_geneID <- strsplit(ora_input$geneID, "/")
            result$Add_geneID <- list(character(0))
        }else{
            input_list <- strsplit(ora_input$geneID, "/")
            names(input_list) <- paste0("sign_", 1:length(input_list))
            ora_input$sign_name <- names(input_list)
            
            icpt_df <- matrix(data = NA, nrow = length(input_list), ncol = length(input_list))
            rownames(icpt_df) <- names(input_list)
            colnames(icpt_df) <- names(input_list)
            
            for(row in 1:length(input_list)){
                for(col in 1:length(input_list)){
                    #cat("row ", row, " col ", col, "\n")
                    icpt_df[row, col] <- round(100 * 2 * length(intersect(input_list[[row]], 
                                                                          input_list[[col]])) / 
                                                   sum(length(input_list[[row]]), length(input_list[[col]])), 1)
                }
            }
            
            # Find optimal epsilon automatically
            kNNdist_vec <- sort(kNNdist(icpt_df, k=1))
            #eps_opt <- (max(kNNdist_vec) + min(kNNdist_vec)) / 2
            eps_opt <- max(kNNdist_vec) / 2
            
            clusters <- dbscan(icpt_df, eps = eps_opt, minPts = 1)$cluster
            length(unique(clusters))
            clusters
            
            clus_df <- data.frame(sign = names(input_list), 
                                  clusters = clusters,
                                  stringsAsFactors = FALSE)
            
            icpt_df2 <- icpt_df[order(clus_df$clusters),]
            icpt_df2 <- icpt_df2[,order(clus_df$clusters)]
            
            # Save pictures
            onefile = FALSE
            clustersNumber <- length(unique(clusters))
            
            if(length(unique(clusters)) > 1){
                saveSignSimilarity(icpt_df, clus_df, outputDir, outfile_prefix)
                
                simMed <- mkSimMed(as.matrix(icpt_df), clusters)
                pheatmapObject <- pheatmap::pheatmap(simMed)
                
                ## Save simMed heatmap
                pdf(file=file.path(outputDir, "pictures",
                                   paste(outfile_prefix, "families_intersection", 
                                         clustersNumber,"clusters.pdf", sep="_")),
                    width=max(4, ncol(simMed) / 4), height=max(3.5, ncol(simMed) / 4.5), onefile=onefile)
                grid::grid.newpage()
                grid::grid.draw(pheatmapObject$gtable)
                dev.off()
                
                ## Save simMed scatterplot
                pdf(file=file.path(outputDir, "pictures", paste(outfile_prefix,
                                                                "families_intersection", clustersNumber,
                                                                "clusters_scatterplot.pdf", sep="_")),
                    width=5, height=4, onefile=onefile # not changed by default
                )
                plot(diag(simMed), ylim = c(min(70, min(diag(simMed))), 100), xlab = "Signature family name", 
                     ylab = "Avg intersection within a family", yaxt = "n", xaxt = "n")
                # avg_pair_overl_thr is the threshold for avg intercept
                abline(h=avg_pair_overl_thr, lty=2, col = "green") 
                axis(2, las = 2)
                axis(1, at = seq(1, length(diag(simMed)), by = 1), las = 0)
                dev.off()
            }else{
                message("Signatures form only one cluster: simMat and simMed won't be plotted,\n", 
                        "but the result table will be saved. Everything is ok.")
                val_list <- strsplit(ora_input$geneID, "/")
                overlapSize <- length(Reduce(intersect, val_list))
                simMed <- matrix(round(100 * length(val_list) * overlapSize / 
                                           sum(lengths(val_list)), 0))
                colnames(simMed) <- c("1")
                rownames(simMed) <- c("1")
            }
            
            ## Save kNNdist plot
            pdf(file=file.path(outputDir, "pictures", paste(outfile_prefix,
                                                            "kNNdistancePlot_k1_optimalEps.pdf", sep="_")),
                width=5, height=4, onefile=onefile # not changed by default
            )
            plot(sort(kNNdist_vec), type="l", ylab=paste(1, "-NN distance", sep=""),
                 xlab = "Points (sample) sorted by distance", yaxt = "n", xaxt = "n")
            abline(h = eps_opt, lty=2)
            axis(2, las = 2)
            axis(1, las = 0)
            text(1/3*length(kNNdist_vec), 1.1 * eps_opt,
                 paste0("Optimal epsilon = ", round(eps_opt, 0)))
            dev.off()
            
            # Split outliers if they cluster together
            if(any(diag(simMed) < avg_pair_overl_thr)){
                outlier_clusters <- rownames(simMed)[diag(simMed) < avg_pair_overl_thr]
                clus_df$clusters_approved <- clus_df$clusters
                for(k in 1:length(outlier_clusters)){
                    length_outl <- length(clus_df$clusters[clus_df$clusters == 
                                                               outlier_clusters[k]])
                    clus_df$clusters_approved[clus_df$clusters == outlier_clusters[k]] <-
                        c((max(clus_df$clusters_approved) + 1) :
                              (max(clus_df$clusters_approved) + length_outl))
                }
                clus_df$clusters <- clus_df$clusters_approved
                clus_df <- clus_df[, colnames(clus_df) != "clusters_approved"]
            }
            
            # Calculate the result table
            ## Order families (clusters) by size
            clus_freq <- as.data.frame(table(clus_df$clusters))
            clus_freq <- clus_freq[order(clus_freq$Freq, decreasing = TRUE),]
            colnames(clus_freq)[1] <- "clusters"
            
            ora_input <- merge(ora_input, clus_df,
                               by.x="sign_name", by.y="sign",
                               all.x=TRUE, all.y=TRUE, sort=FALSE)
            
            colnames_result <- c("Descriptions", "Avg_pval", 
                                 "Avg_padj", 
                                 "Avg_qvalue", "Family_size", "Overlap_size",
                                 "Overlap_pct",
                                 "Overlap_geneID", 
                                 "Add_geneID", "IDs")
            result <- as.data.frame(matrix(NA, ncol = length(colnames_result), 
                                           nrow = nrow(clus_freq)))
            colnames(result) <- colnames_result
            for(i in 1:nrow(clus_freq)){
                pathways_list <- as.list(ora_input$Description[ora_input$clusters == 
                                                                   clus_freq$clusters[i]])
                names(pathways_list) <- ora_input$sign_name[ora_input$clusters == 
                                                                clus_freq$clusters[i]]
                result$Descriptions[i] <- list(pathways_list)
                
                result$Avg_pval[i] <- mean(ora_input$pvalue[ora_input$clusters == 
                                                                clus_freq$clusters[i]])
                
                result$Avg_padj[i] <- mean(ora_input$p.adjust[ora_input$clusters ==
                                                                  clus_freq$clusters[i]])
                
                result$Avg_qvalue[i] <- mean(ora_input$qvalue[ora_input$clusters ==
                                                                  clus_freq$clusters[i]])
                
                result$Family_size[i] <- length(pathways_list)
                
                val_list <- strsplit(ora_input$geneID[ora_input$clusters ==
                                                          clus_freq$clusters[i]], "/")
                overlapGenes <- list(Reduce(intersect, val_list))
                
                result$Overlap_size[i] <- length(overlapGenes[[1]])
                
                result$Overlap_pct[i] <- round(100 * length(val_list) * result$Overlap_size[i] / 
                                                   sum(lengths(val_list)), 0)
                
                result$Overlap_geneID[i] <- overlapGenes
                
                add_list <- strsplit(ora_input$geneID[ora_input$clusters == 
                                                          clus_freq$clusters[i]], "/")
                add_list <- lapply(add_list, function(x) x[!(x %in% result$Overlap_geneID[i][[1]])])
                names(add_list) <- ora_input$sign_name[ora_input$clusters == 
                                                           clus_freq$clusters[i]]
                result$Add_geneID[i] <- list(add_list)
                
                ids_list <- as.list(ora_input$ID[ora_input$clusters == 
                                                     clus_freq$clusters[i]])
                
                names(ids_list) <- ora_input$sign_name[ora_input$clusters == 
                                                           clus_freq$clusters[i]]
                result$IDs[i] <- list(ids_list)
            }
        }
    }
    
    result <- result[order(result$Avg_qvalue), ]
    
    ## Save result table
    saveRDS(result, file.path(outputDir, "tables", paste0(outfile_prefix, "_", res_fn)))
}

setFamDir <- function(outputDir){
    # Create output directory for pictures if does not exist
    dir.create(file.path(outputDir, "pictures"), showWarnings = FALSE)
    
    # Delete old files in output directory if any
    oldFiles <- list.files(file.path(outputDir, "pictures"))
    deletedFiles <- sapply(oldFiles, function(fileName) 
        file.remove(file.path(file.path(outputDir, "pictures"), fileName)))
    
    ## Create output directory if does not exist
    dir.create(file.path(outputDir, "tables"), showWarnings = FALSE)
    
    ## Delete old files in output directory if any
    oldFiles <- list.files(file.path(outputDir, "tables"))
    deletedFiles <- sapply(oldFiles, function(fileName) 
        file.remove(file.path(file.path(outputDir, "tables"), fileName)))
}

saveSignSimilarity <- function(signSimilarityMatrix, clus_df, dataDirectory,
                               samp, colorPalette="default",
                               statePalette="default", clusteringMethod="ward.D2",
                               plot = TRUE,
                               returnPlot = FALSE,
                               width = max(5, ncol(signSimilarityMatrix) / 8), 
                               height = max(4.5, ncol(signSimilarityMatrix) / 9), 
                               onefile = FALSE, #pdf
                               color = colorRampPalette(rev( #pheatmap
                                   RColorBrewer::brewer.pal(n = 7,
                                                            name = "RdYlBu")))(100),
                               breaks = NA,
                               border_color = "grey60",
                               cellwidth = NA, cellheight = NA,
                               scale = "none",
                               cluster_rows = FALSE, #not default
                               cluster_cols = FALSE, #not default
                               clustering_callback = identity2,
                               cutree_rows = NA, cutree_cols = NA,
                               treeheight_row = ifelse((
                                   class(cluster_rows) == "hclust") ||
                                       cluster_rows, 50, 0),
                               treeheight_col = ifelse((
                                   class(cluster_cols) == "hclust") ||
                                       cluster_cols, 50, 0),
                               annotation_row = NA,
                               annotation_col = NA,
                               annotation_colors = NA,
                               show_rownames = TRUE, #not default
                               show_colnames = TRUE, #not default
                               fontsize = 7.5, #not default (10)
                               fontsize_row = fontsize, #not default
                               fontsize_number = 0.8 * fontsize,
                               main = paste0("Signatures similarity matrix ",
                                             ncol(signSimilarityMatrix),
                                             " columns, ",
                                             nrow(signSimilarityMatrix),
                                             " rows.")
){
    # plots cells correlation matrix gained form
    # clusterCellsInternal() function as the result of DBSCAN
    
    colData <- clus_df
    clustersNumber <- length(unique(clus_df$clusters))
    rownames(colData) <- colData$sign
    colData$clusters <- as.factor(colData$clusters)
    
    #distanceMatrix <- as.dist(sqrt((100-signSimilarityMatrix)/2))
    distanceMatrix <- dist(signSimilarityMatrix)
    clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
    cluster_cols <- clusteringTree
    cluster_rows <- clusteringTree
    
    annotationColors <- generateAnnotationColors(colData, colorPalette,
                                                 statePalette)
    columnsToPlot <- switch(is.null(colData$state) + 1, c("clusters", "state"),
                            c("clusters"))
    
    pheatmapObject <- pheatmap::pheatmap(signSimilarityMatrix,
                                         show_colnames=show_colnames,
                                         show_rownames=show_rownames,
                                         annotation_col=as.data.frame(colData[columnsToPlot]),
                                         annotation_colors=annotationColors,
                                         fontsize_row=fontsize_row,
                                         cluster_cols=clusteringTree,
                                         cluster_rows=clusteringTree,
                                         fontsize=fontsize,
                                         main = main,
                                         color = color # not changed by default
    )
    
    pdf(file=file.path(dataDirectory, "pictures", paste(samp,
                                                        "signatures_intersection", clustersNumber,
                                                        "clusters.pdf", sep="_")),
        width=width, height=height, onefile=onefile # not changed by default
    )
    grid::grid.newpage()
    grid::grid.draw(pheatmapObject$gtable)
    dev.off()
    
    if(returnPlot){
        return(pheatmapObject)
    }
}

generateAnnotationColors <- function(colData, colorPaletteParameter,
                                     statePalette){
    
    clusters <- levels(colData$clusters)
    states <- unique(colData$state)
    clusterNumber <- length(unique(colData$clusters))
    
    colorAnnotationClusters <- choosePalette(colorPaletteParameter, clusterNumber)
    #colorAnnotationState <- chooseStatePalette(length(states))
    colorAnnotationState <- choosePalette(statePalette, length(states))
    names(colorAnnotationState) <- states
    names(colorAnnotationClusters) <- clusters
    
    return(list(state=colorAnnotationState, clusters=colorAnnotationClusters))
}

# A function from CONCLUS r package
#' Choose palette for a plot.
#'
#' It is an internal function usually applied for choosing the palette for clusters.
#' Depending if the number of clusters is more than 12 or not, one of two built-in palettes will be applied.
#' If you give your vector of colors, the function will not change them.
#' If the number of clusters is more than 26, it will copy colors to get the needed length of the palette.
#'
#' @param colorPalette Either "default" or a vector of colors, for example c("yellow", "#CC79A7"). 
#' @param clustersNumber number of clusters in the output palette.
#'
#' @return Color palette with the number of colors equal to the clusterNumber parameter.
#' @export
choosePalette <- function(colorPalette, clustersNumber){
    
    colorPalette26 <- c( "yellow", "darkgoldenrod1", "coral1", "deeppink",
                         "indianred", "coral4", "darkblue", "darkmagenta",
                         "darkcyan", "mediumorchid", "plum2", "gray73",
                         "cadetblue3", "khaki",
                         "darkolivegreen1", "lightgreen", "limegreen",
                         "darkolivegreen4", "green", "#CC79A7", "violetred3",
                         "brown3", "darkslategray1", "gray51", "slateblue2",
                         "blue")
    
    pickDefaultPalette <- function(clustersNumber, colorPalette26){
        if(clustersNumber < 13) return(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                                         "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                         "#CAB2D6", "#6A3D9A", "#FFFF99",
                                         "#B15928")[1:clustersNumber])
        return(rep(colorPalette26,
                   round(clustersNumber/length(colorPalette26))+1)[1:clustersNumber])
    }
    
    if(length(colorPalette) == 1){
        if(colorPalette == "default"){
            return(pickDefaultPalette(clustersNumber, colorPalette26))
        }
    }
    
    if(clustersNumber > length(colorPalette)){
        message("The number of clusters is greater than the number of colors.
                Using default CONCLUS palette instead.")
        return(pickDefaultPalette(clustersNumber, colorPalette26))
    }
    
    return(colorPalette[1:clustersNumber])
}

mkSimMed <- function(simMat, clusters){
    
    clusMed <- matrix(ncol=length(unique(clusters)), nrow=nrow(simMat))
    clusterNames <- unique(clusters)
    
    for(i in 1:ncol(clusMed)){
        if(length(clusters[clusters == clusterNames[i]]) == 1){
            clusMed[,i] <- simMat[,clusters == clusterNames[i]]
        }else{
            clusMed[,i] <- round(rowMeans(simMat[,clusters == clusterNames[i]]), 1)
        }
    }
    
    clusMed <- t(clusMed)
    
    simMed <- matrix(ncol=length(unique(clusters)),
                     nrow=length(unique(clusters)))
    
    for(i in 1:ncol(simMed)){
        if(length(clusters[clusters == clusterNames[i]]) == 1){
            simMed[,i] <- clusMed[,clusters == clusterNames[i]]
        }else{
            simMed[,i] <- round(rowMeans(clusMed[,clusters == clusterNames[i]]), 1)
        }
    }
    
    # colnames(simMed) = 1:length(unique(clusters))
    # rownames(simMed) = 1:length(unique(clusters))
    
    colnames(simMed) <- clusterNames
    rownames(simMed) <- clusterNames
    
    return(simMed)
}

save3cl <- function(padj_df_dir, fn_df, 
                    allG_dir, outputDir, filePrefix,
                    pairsToInclude = "all",
                    ctr_cond, 
                    treat_cond){
    
    padj_df <- read.delim(file.path(padj_df_dir, fn_df$Name[fn_df$Type == "padj_df"]),
                          stringsAsFactors = FALSE)
    if(pairsToInclude == "all"){
        pairsToInclude <- 1:(ncol(padj_df)-1)
    }
    
    fns_allG <- reorderFns(list.files(allG_dir, pattern = ".rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    (fns_allG <- fns_allG[pairsToInclude])
    #deseq2_filePrefix <- gsub("_all.*", "", fns_allG)
    results_allG <- vector("list", length = length(fns_allG))
    
    for(i in 1:length(fns_allG)){
        results_allG[[i]] <- as.data.frame(readRDS(file.path(allG_dir, fns_allG[i])))
        results_allG[[i]]$Gene <- rownames(results_allG[[i]])
        results_allG[[i]]$Index <- i
    }
    
    all_res_df <- results_allG[[1]]
    if(length(pairsToInclude) > 1){
        for(i in 2:length(fns_allG)){
            all_res_df <- rbind(all_res_df, results_allG[[i]])
        } 
    }
    
    categ_df <- data.frame(Gene = padj_df$Gene, categ = NA)
    for(i in 1:nrow(categ_df)){
        if(all(all_res_df$log2FoldChange[all_res_df$Gene == categ_df$Gene[i]] > 0)){
            categ_df$categ[i] <- "allDown"
        }else if(all(all_res_df$log2FoldChange[all_res_df$Gene == categ_df$Gene[i]] < 0)){
            categ_df$categ[i] <- "allUp"
        }else{
            categ_df$categ[i] <- "so-so"
        }
    }
    
    categ_df <- categ_df[order(categ_df$categ), ]
    categ_df_fn <- paste0(filePrefix, "3categ_", 
                          nrow(categ_df), 
                          "_DE_genes_pairs", 
                          paste0(pairsToInclude, collapse = ""), ".tsv")
    write.table(categ_df, file.path(outputDir, categ_df_fn),
                quote = FALSE, sep = "\t", row.names = FALSE)
    return(categ_df_fn)
}

advancedGeneClusters <- function(cl3_df, geneClusters_df,
                                 outputDir, filePrefix, fn_affix = ""){
    result <- geneClusters_df[, 1:(ncol(geneClusters_df) - 1)]
    unique3cl <- unique(cl3_df$categ)
    result[, unique3cl[2]] <- NA
    result[, unique3cl[1]] <- NA
    if(length(unique3cl) == 3){
        result[, unique3cl[3]] <- NA 
    }
    
    for(i in 1:nrow(geneClusters_df)){
        #cat(i , " ")
        tmp <- unlist(strsplit(geneClusters_df$genes[i], ","))
        result[, unique3cl[2]][i] <- 
            paste0(tmp[tmp %in% cl3_df$Gene[cl3_df$categ == unique3cl[2]]], collapse = ",")
        
        result[, unique3cl[1]][i] <- 
            paste0(tmp[tmp %in% cl3_df$Gene[cl3_df$categ == unique3cl[1]]], collapse = ",")
        
        if(length(unique3cl) == 3){
            result[, unique3cl[3]][i] <- 
                paste0(tmp[tmp %in% cl3_df$Gene[cl3_df$categ == unique3cl[3]]], collapse = ",")
        }
    }
    advanced_geneCl_fn <- paste0(filePrefix, nrow(geneClusters_df), "_nonZero_geneClusters_",
                                 "trend", fn_affix, ".tsv")
    write.table(result, file.path(outputDir, 
                                  advanced_geneCl_fn),
                quote = FALSE, sep = "\t")
    return(advanced_geneCl_fn)
}

# For debugging parseStyle
#inputdir <- "/data/cabezas/group/pavlovich/data/shared/storshow/runs/input"
#dds <- readRDS("/data/cabezas/group/pavlovich/data/shared/storshow/runs/outputdir137/DESeq/dds_16s_12351g.rds")
#coldata <- as.data.frame(SummarizedExperiment::colData(dds))
#source("/data/cabezas/group/pavlovich/data/shared/storshow/R/dev_Polina/Storshow/scripts/storshow_functions.R")
#style_file <- file.path(inputdir, "CypKO_style_noColNames.tsv")

parseStyle <- function(style_file, coldata){
    if(style_file == ""){
        levels_ordered <- unique(coldata$state)
        colors_ordered <- choosePalette("default", coldata_nstates)
    }else{
        style_df <- read.delim(style_file, stringsAsFactors = F, header = TRUE)
        
        if(colnames(style_df)[1] %in% unique(coldata$state)){
            style_df <- read.delim(style_file, stringsAsFactors = F, header = FALSE)
        }
        
        coldata_nstates <- length(unique(coldata$state))
        if(nrow(style_df) != coldata_nstates){
            warning("The number of conditions in the style.tsv file does not match the dds object.",
                    "Levels and colors will be assigned automatically.")
            levels_ordered <- unique(coldata$state)
            colors_ordered <- choosePalette("default", coldata_nstates)
        }
        if(ncol(style_df) >= 2){
            levels_ordered <- style_df[, 1]
            colors_ordered <- style_df[, 2]
        }else if(ncol(style_df) == 1){
            levels_ordered <- style_df[, 1]
            colors_ordered <- choosePalette("default", length(levels_ordered))
        }else{
            warning("Style.tsv is empty. Levels and colors will be assigned automatically.")
            levels_ordered <- unique(coldata$state)
            colors_ordered <- choosePalette("default", length(levels_ordered))
        }
    }
    
    return(list(levels_ordered, colors_ordered))
}

# for ShinyInternal
# giveMult counts if a pathway was significant in more than one comparison
giveMult <- function(GO_resList, pairs_df){
    tmp_list <- vector("list", length = nrow(pairs_df))
    for(i in 1:nrow(pairs_df)){
        tmp_list[[i]] <- GO_resList[[i]]@result
        if(nrow(tmp_list[[i]]) > 0){
            tmp_list[[i]]$pairNum <- pairs_df$pairNum[i] 
            tmp_list[[i]]$cond1 <- pairs_df$cond1[i]
            tmp_list[[i]]$cond2 <- pairs_df$cond2[i]
        }
    }
    tmp <- tmp_list[[1]]
    if(length(tmp_list) > 1){
        for(i in 2:nrow(pairs_df)){
            if(nrow(tmp_list[[i]]) > 0){
                tmp <- rbind(tmp, tmp_list[[i]])
            }
        }
    }
    
    tmp <- tmp[!is.na(tmp$Description),]
    path_counts <- as.data.frame(table(tmp$Description))
    if(nrow(path_counts) > 0){
        path_counts <- path_counts[order(path_counts$Freq, decreasing = TRUE),]
        path_counts <- path_counts[path_counts$Freq > 1,]
        if(nrow(path_counts) > 0){
            tmp <- tmp[tmp$Description %in% path_counts$Var1,]
            tmp <- merge(tmp, path_counts,
                         by.x="Description", by.y="Var1",
                         all.x=TRUE, all.y=TRUE, sort=FALSE)
            tmp <- tmp[,c(ncol(tmp), 1:(ncol(tmp)-1))]
            tmp <- tmp[order(tmp$Freq, tmp$Description, decreasing = TRUE),]
            rownames(tmp) <- 1:nrow(tmp)
            return(tmp)
        }else{
            return(data.frame(Data = "Zero repeating signatures"))
        }
    }else{
        return(data.frame(Data = "Zero repeating signatures"))
    }
}

giveMultFGSEA <- function(FGSEA_resList, pairs_df){
    tmp_list <- vector("list", length = nrow(pairs_df))
    for(i in 1:nrow(pairs_df)){
        tmp_list[[i]] <- as.data.frame(FGSEA_resList[[i]])
        if(nrow(tmp_list[[i]]) > 0){
            tmp_list[[i]]$pairNum <- pairs_df$pairNum[i] 
            tmp_list[[i]]$cond1 <- pairs_df$cond1[i]
            tmp_list[[i]]$cond2 <- pairs_df$cond2[i]
        }
    }
    tmp <- tmp_list[[1]]
    if(length(tmp_list) > 1){
        for(i in 2:nrow(pairs_df)){
            if(nrow(tmp_list[[i]]) > 0){
                tmp <- rbind(tmp, tmp_list[[i]])
            }
        }
    }
    path_counts <- table(tmp$pathway)
    path_counts <- path_counts[order(path_counts, decreasing = TRUE)]
    path_counts <- path_counts[path_counts > 1]
    tmp <- tmp[tmp$pathway %in% names(path_counts),]
    
    if(nrow(tmp) == 0){
        tmp_down <- data.frame(Data = "Zero repeating signatures")
        tmp_up <- data.frame(Data = "Zero repeating signatures")
    }else{
        path_counts <- as.data.frame(path_counts)
        tmp <- merge(tmp, path_counts,
                     by.x="pathway", by.y="Var1",
                     all.x=TRUE, all.y=TRUE, sort=FALSE)
        tmp <- tmp[,c(ncol(tmp), 1:(ncol(tmp)-1))]
        tmp <- tmp[order(tmp$Freq, tmp$pathway, decreasing = TRUE),]
        keepPaths_up <- c()
        keepPaths_down <- c()
        for(k in 1:length(unique(tmp$pathway))){
            onePath_df <- tmp[(tmp$pathway == unique(tmp$pathway)[k]) &
                                  !(tmp$pairNum %in% c(7)),] # exclude Rescue pairs
            if(all(onePath_df$ES > 0)){
                keepPaths_up <- c(keepPaths_up, unique(tmp$pathway)[k])
            }else if(all(onePath_df$ES < 0)){
                keepPaths_down <- c(keepPaths_down, unique(tmp$pathway)[k])
            }
        }
        tmp_down <- tmp[tmp$pathway %in% keepPaths_down,]
        tmp_up <- tmp[tmp$pathway %in% keepPaths_up,]
        
        if(nrow(tmp_down) > 0){
            rownames(tmp_down) <- 1:nrow(tmp_down)
        }
        if(nrow(tmp_up) > 0){
            rownames(tmp_up) <- 1:nrow(tmp_up)
        }
    }
    return(list(tmp_down, tmp_up))
}

# for ShinyInternal
limitNgenes <- function(enrichObj, maxNgenes = 7){
    tmp <- enrichObj@result$geneID
    tmp <- unname(sapply(tmp, function(x) strsplit(x, split = "/"))) 
    tmp <- unlist(lapply(tmp, function(x) paste0(x[1:min(maxNgenes, length(x))], collapse = "/")))
    enrichObj@result$geneID <- tmp
    return(enrichObj)
}

# for ShinyInternal
# showStates is a list where each element is a vector of length 2
# levels_ordered is a vector of levels in a desired order
plotBarPlots <- function(geneName, dds, padj_df, 
                         levels_ordered, 
                         colorSwitch, # a data frame
                         showStates, # a list
                         cex.empty = 0.7){
    if(geneName == "coming"){
        return(plotOnlyText(textToPlot = "If you selected a signature, please wait."),
               cex = cex.empty)
    }else{
        if(!is.vector(levels_ordered)){
            stop("levels_ordered is not a vector.")
        }
        if(!is.list(showStates)){
            stop("showStates is not a list.")
        }
        #if(!is.function(colorSwitch)){
        #    stop("colorSwitch is not a function.")
        #}
        if(ncol(colorSwitch) != 2 | 
           nrow(colorSwitch) != length(unique(colData(dds)$state))){
            stop(paste0("colorSwitch has incorrect dimenstions. ",
                        "It must have two columns and ",
                        length(unique(colData(dds)$state)), " rows."))
        }
        p <- list()
        for(i in 1:length(showStates)){
            p[[i]] <- giveBarPlot(geneName, dds, showStates[[i]], padj_df,
                                  levels_ordered, colorSwitch,
                                  titleAffix = "")
        }
        
        return(plot_grid(plotlist = p, align = "hv"))
    }
}

# for ShinyInternal
# returns a barplot for VitA data
giveBarPlot <- function(geneName, dds, showStates_i, padj_df, 
                        levels_ordered, colorSwitch, #colorSwitch is a data frame
                        titleAffix = ""){
    
    df <- data.frame(state = colData(dds)$state[colData(dds)$state %in% showStates_i],
                     expression = counts(dds, 
                                         normalized = TRUE)[rownames(dds) == geneName,
                                                            colData(dds)$state %in% showStates_i])
    
    
    if(length(levels_ordered) == length(unique(colData(dds)$state))){
        df <- within(df, state <- factor(state, levels=levels_ordered))
    }else{
        warning(paste0("The number of conditions in levels_ordered does not match ",
                       "the number of states in the DESeq2 object. The order of ",
                       "levels will be defined automatically."))
    }
    
    #df$state <- as.factor(df$state)
    
    colorPalette <- c()
    for(col_i in 1:length(showStates_i)){
        #colorPalette <- c(colorPalette, colorSwitch(showStates_i[col_i]))
        colorPalette <- c(colorPalette, 
                          colorSwitch$colors_ordered[colorSwitch$levels_ordered == 
                                                         showStates_i[col_i]])
    }
    
    df.summary <- df %>%
        group_by(state) %>%
        summarise(
            sd = sd(expression, na.rm = TRUE),
            expression = mean(expression)
        )
    df.summary
    
    y <- max(df$expression)
    
    giveLabelSets <- function(geneName, padj_df, df, index1, index2){
        sizeL <- 4
        posYL <- y + 0.14*y
        label <- padj_df[padj_df$Gene == geneName, 
                         grepl(unique(df$state)[index1], colnames(padj_df)) & 
                             grepl(unique(df$state)[index2], colnames(padj_df))]
        if(length(label) == 0){
            label <- "ns"
            sizeL <- 3
            posYL <- y + 0.18*y
        }else if(is.na(label)){ 
            label <- "ns"
            sizeL <- 3
            posYL <- y + 0.18*y
        }else if(label < 0.0001){label <- "****"
        }else if(label < 0.001){label <- "***"
        }else if(label < 0.01){label <- "**"
        }else if(label < 0.1){label <- "*"}
        return(list(posYL, label, sizeL))
    }
    labelSets1 <- giveLabelSets(geneName, padj_df, df, 1, 2)
    
    # (2) Bar plots of means + individual jitter points + errors
    p <- ggplot(df, aes(state, expression)) +
        geom_bar(stat = "identity", data = df.summary, fill = colorPalette,
                 color = "black") +
        geom_jitter(position = position_jitter(0.15), color = "black") + 
        geom_errorbar(aes(ymin = expression-sd, ymax = expression+sd),
                      data = df.summary, width = 0.2) +
        theme_classic() + 
        ggtitle(paste0(geneName, titleAffix)) +
        labs(y="norm. expr") +
        theme(text = element_text(size=10),
              axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 11),
              axis.text.y = element_text(colour = "black", size = 9),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = 11),
              plot.title = element_text(size = 11))  +
        geom_text(aes(x = c(1.5), y = labelSets1[[1]], label = labelSets1[[2]]), size = labelSets1[[3]]) +
        geom_segment(aes(x = 1.25, xend = 1.75, y = y + 0.1*y, yend = y + 0.1*y))
    
    #annotate("text", label = "**", x = 2.5, y = 15, size = 5, colour = "red") + 
    if(length(unique(df$state)) == 3){
        labelSets2 <- giveLabelSets(geneName, padj_df, df, 2, 3)
        labelSets3 <- giveLabelSets(geneName, padj_df, df, 1, 3)
        p <- p + 
            coord_cartesian(ylim = c(0, y + 0.35*y)) + 
            geom_segment(aes(x = 2.25, xend = 2.75, y = y + 0.05*y, yend = y + 0.05*y)) +
            geom_segment(aes(x = 1.25, xend = 2.75, y = y + 0.22*y, yend = y + 0.22*y)) +
            geom_text(aes(x = 2.5, y = labelSets2[[1]], label = labelSets2[[2]]), 
                      size = labelSets2[[3]]) +
            geom_text(aes(x = 2, y = labelSets3[[1]] + 0.18*y, label = labelSets3[[2]]), 
                      size = labelSets3[[3]])
    }else if(length(unique(df$state)) == 2){
        p <- p + coord_cartesian(ylim = c(0, y + 0.15*y))
    }
    
    return(p)
}

# for ShinyInternal
plotOnlyText <- function(textToPlot = "Zero significant", cex = 1.6){
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste(textToPlot), 
         cex = cex, col = "gray40")
}

# for ShinyInternal
giveVolcanoPlot <- function(results_df, cond1, cond2, 
                            nShow = 25, 
                            DESeq2_padj_thr = 0.1, 
                            DESeq2_absLog2FC_thr = 0.5){
    #library(ggrepel)
    results_df <- results_df[order(results_df$padj), ]
    results_df <- results_df[!is.na(results_df$padj), ]
    results_df$Significant <- "ns"
    results_df$Significant[(results_df$padj < DESeq2_padj_thr) & 
                               (results_df$log2FoldChange > DESeq2_absLog2FC_thr)] <- "Up"
    results_df$Significant[(results_df$padj < DESeq2_padj_thr) & 
                               (results_df$log2FoldChange < -DESeq2_absLog2FC_thr)] <- "Down"
    
    colorSwitch <- function(inputVec){
        res <- c()
        for(cond in inputVec){
            res <- c(res, switch(cond, 
                                 "Down" = "blue", 
                                 "Up" = "red", 
                                 "ns" = "gray32"))
        }
        return(res)
    }
    
    results_df$Significant <- factor(results_df$Significant)
    levels_sign <- levels(results_df$Significant)
    
    subset_up <- subset(results_df, (padj < DESeq2_padj_thr) & log2FoldChange > DESeq2_absLog2FC_thr)
    subset_down <- subset(results_df, (padj < DESeq2_padj_thr) & log2FoldChange < -DESeq2_absLog2FC_thr)
    
    p <- ggplot(results_df, aes(x = -log10(padj), y = log2FoldChange)) +
        geom_point(aes(color = Significant), size = 1) +
        scale_color_manual(values = colorSwitch(levels_sign)) + # c("blue", "gray32", "red")
        theme_bw(base_size = 12) + theme(legend.position = "bottom") +
        geom_text_repel(
            data = subset_up[1:min(nShow, nrow(subset_up)), ],
            aes(label = Gene),
            size = 3,
            colour = "#6e0707", #"red",
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")) +
        geom_text_repel(
            data = subset_down[1:min(nShow, nrow(subset_down)), ],
            aes(label = Gene),
            size = 3,
            colour = "#050563", #"blue",
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        )+
        ggtitle(paste0("Pos log2FC is up in ", 
                       cond2, ", neg is up in ", cond1))
    return(p)
}

plotVolcanoPlots <- function(results_allG_list, pairs_df,
                             nShow = 25, DESeq2_padj_thr = 0.1, DESeq2_absLog2FC_thr = 0.5){
    #library(ggpubr)
    plots <- vector("list", length = length(results_allG_list))
    for(i in 1:length(results_allG_list)){
        plots[[i]] <- giveVolcanoPlot(results_allG_list[[i]], 
                                      cond1 = pairs_df$cond1[i],
                                      cond2 = pairs_df$cond2[i],
                                      nShow = nShow, 
                                      DESeq2_padj_thr = DESeq2_padj_thr, 
                                      DESeq2_absLog2FC_thr = DESeq2_absLog2FC_thr)
    }
    return(ggarrange(plotlist = plots, ncol = 1, nrow = length(plots)))
}

# for ShinyInternal
# MA plots
giveMAplots <- function(inputDir, samp,
                        ylim = c(-2,2), DESeq2_padj_thr = 0.1, 
                        ctr_cond, treat_cond){
    fns <- reorderFns(list.files(inputDir, pattern = ".rds"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    results_list <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        results_list[[i]] <- readRDS(file.path(inputDir, fns[i]))
    }
    #ctr_cond <- gsub(paste0(samp, "_DESeq2_res_"), "", fns)
    #ctr_cond <- gsub("_vs_.*", "", ctr_cond, perl = TRUE)
    #treat_cond <- gsub(paste0(samp, "_DESeq2_res_.*_vs_"), "", fns, perl = T)
    #treat_cond <- gsub("_all.*", "", treat_cond)
    par(mfrow=c(length(fns),1), mar = c(6, 6, 5, 2) + 0.1)
    for(i in 1:length(fns)){
        DESeq2::plotMA(results_list[[i]], ylim=ylim, alpha = DESeq2_padj_thr,
                       cex = 0.5, cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2,
                       main = paste0("Pos log2FC is up in ", 
                                     treat_cond[i], ", neg is up in ", ctr_cond[i]),) 
    }
}

#for ShinyInternal
# DE counts
giveDEcounts <- function(inputDir, ctr_cond, treat_cond){
    fns <- reorderFns(list.files(inputDir), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    results_list <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        results_list[[i]] <- read.delim(file.path(inputDir, fns[i]), stringsAsFactors = FALSE)
    }
    DE_counts_df <- as.data.frame(matrix(ncol = 2, nrow = length(fns)))
    colnames(DE_counts_df) <- c("up_cond1", "up_cond2")
    for(i in 1:length(fns)){
        DE_counts_df$up_cond1[i] <- nrow(results_list[[i]][results_list[[i]]$log2FoldChange < 0,])
        DE_counts_df$up_cond2[i] <- nrow(results_list[[i]][results_list[[i]]$log2FoldChange > 0,])
    }
    #DE_counts_df <- t(DE_counts_df)
    #colnames(DE_counts_df) <- c(1:length(fns))
    DE_counts_df$pairNum <- c(1:length(fns))
    return(DE_counts_df)
}

# for ShinyInternal
fgseaToEnrichResult <- function(fgseaInput, 
                                enrichExample,
                                dataDirectory,
                                rowdata,
                                forPlot = "cnetplot", # or "barplot"#
                                maxNgenes = 7){
    #enrichExample <- bp_ego_s_up[[1]]
    
    enrichRes_custom <- enrichExample
    
    #rowdata <- read.delim(file.path(dataDirectory, "VitA/FourExp_rowdata_22970g.tsv"),
    #                      stringsAsFactors = FALSE)
    
    fgseaConvertedToEnRes <- as.data.frame(matrix(NA, ncol = ncol(enrichRes_custom@result),
                                                  nrow = nrow(fgseaInput)))
    colnames(fgseaConvertedToEnRes) <- colnames(enrichRes_custom@result)
    
    rownames(fgseaConvertedToEnRes) <- 1:nrow(fgseaConvertedToEnRes) 
    fgseaConvertedToEnRes$ID <- 1:nrow(fgseaConvertedToEnRes)
    fgseaConvertedToEnRes$Description <- substr(fgseaInput$pathway, 1, 50)
    
    fgseaConvertedToEnRes$GeneRatio <- rep("", nrow(fgseaConvertedToEnRes))
    fgseaConvertedToEnRes$BgRatio <- rep("", nrow(fgseaConvertedToEnRes))
    
    fgseaConvertedToEnRes$pvalue <- fgseaInput$pval
    fgseaConvertedToEnRes$p.adjust <- fgseaInput$padj
    fgseaConvertedToEnRes$qvalue <- fgseaInput$padj
    
    if(forPlot == "barplot"){
        fgseaConvertedToEnRes$geneID <- unlist(lapply(fgseaInput$leadingEdge, 
                                                      function(x) 
                                                          paste(x, 
                                                                collapse = "/")))
    }else if(forPlot == "cnetplot"){
        fgseaConvertedToEnRes$geneID <- unlist(lapply(fgseaInput$leadingEdge, 
                                                      function(x) 
                                                          paste(x[1:min(length(x), maxNgenes)], 
                                                                collapse = "/")))
    }else{
        message("Only 'cnetplot' and 'barplot' are options for the 'forPlot' argument.")
    }
    
    fgseaConvertedToEnRes$Count <- lengths(strsplit(fgseaConvertedToEnRes$geneID, "/"))
    
    enrichRes_custom@result <- fgseaConvertedToEnRes
    
    genesToUse <- unique(unlist(strsplit(enrichRes_custom@result$geneID, "/")))
    names(genesToUse) <- rowdata$Entrez.Gene.ID[rowdata$Symbol %in% genesToUse]
    enrichRes_custom@gene2Symbol <- genesToUse
    
    enrichRes_custom@gene <- names(genesToUse)
    return(enrichRes_custom)
}