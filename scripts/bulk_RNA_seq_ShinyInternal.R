#detach("package:rlang", unload=TRUE)
#library(rlang)
require(DESeq2)
#library(Rcpp)
require(dplyr)

#install.packages("shiny")
require(shiny) # version >= 1.2.0 required
require(rstudioapi)
require(ggplot2)
require(knitr)
require(kableExtra)
#library(xtable)
require(DT)
require(cowplot)
require(enrichplot)
require(DOSE)
require(fgsea)
require(ggrepel)
require(ggpubr)

# Shiny Modules
cnetplotsUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            plotOutput(NS(id, sprintf("Cnet%01d",i)))
        )
    })
}

cnetplotsServer <- function(id, pairs_df, data, 
                            maxNgenes, showCategCnet, geneList) {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            if(is.null(data[[i]])){
                output[[sprintf("Cnet%01d",i)]] <- renderPlot(plotOnlyText())
            }else if(nrow(data[[i]]@result[!is.na(data[[i]]@result$Description),]) == 0){
                output[[sprintf("Cnet%01d",i)]] <- renderPlot(plotOnlyText())
            }else{
                dataToShow <- data[[i]]
                dataToShow <- limitNgenes(dataToShow, maxNgenes = maxNgenes)
                output[[sprintf("Cnet%01d",i)]] <- 
                    renderPlot(cnetplot(dataToShow, showCategory = showCategCnet, 
                                        foldChange = geneList[[i]]))
            }
        })
    })
}

tablesUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            tableOutput(NS(id, sprintf("Table%01d",i)))
        )
    })
}

tablesServer <- function(id, pairs_df, data, 
                            styling_opt, DE_df_height) {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            if(nrow(data[[i]]) == 0 | all(is.na(data[[i]]@result))){
                dataToShow <- data.frame(Data = "Zero significant")
            }else{
                dataToShow <- cbind(data.frame(N = 1:nrow(data[[i]]@result)), 
                                    data[[i]]@result)
            }
            output[[sprintf("Table%01d",i)]] <- function() {
                dataToShow %>%
                    knitr::kable(format = "html", row.names = TRUE) %>%
                    kable_styling(styling_opt, full_width = F) %>%
                    scroll_box(width = "100%", height = DE_df_height)
            }
        })
    })
}

barplotsUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            plotOutput(NS(id, sprintf("Barplot%01d",i)))
        )
    })
}

barplotsServer <- function(id, pairs_df, data, showCateg) {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            if(is.null(data[[i]])){
                output[[sprintf("Barplot%01d",i)]] <- renderPlot(plotOnlyText())
            }else if(nrow(data[[i]]@result[!is.na(data[[i]]@result$Description),]) == 0){
                output[[sprintf("Barplot%01d",i)]] <- renderPlot(plotOnlyText())
            }else{
                output[[sprintf("Barplot%01d",i)]] <- 
                    renderPlot(suppressMessages(barplot(data[[i]], showCategory=showCateg) +
                                                    scale_fill_continuous(low="red", high="gray", 
                                                                          name="p.adjust", 
                                                                          guide=guide_colorbar(reverse=TRUE))))
            }
        })
    })
}

famTablesUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            tableOutput(NS(id, sprintf("famTable%01d",i)))
        )
    })
}

famTablesServer <- function(id, pairs_df, data, 
                         styling_opt, DE_df_height) {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            if(nrow(data[[i]]) == 0 | all(is.na(data[[i]]))){
                dataToShow <- data.frame(Data = "Zero significant")
            }else{
                dataToShow <- cbind(data.frame(N = 1:nrow(data[[i]])), 
                                    as.data.frame(data[[i]]))
            }
            output[[sprintf("famTable%01d",i)]] <- function() {
                dataToShow %>%
                    knitr::kable(format = "html", row.names = FALSE) %>%
                    kable_styling(styling_opt, full_width = T) %>%
                    scroll_box(width = "100%", height = DE_df_height)
            }
        })
    })
}

famCnetplotsUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            plotOutput(NS(id, sprintf("famCnet%01d",i)))
        )
    })
}

famCnetplotsServer <- function(id, pairs_df, data, famData,
                            maxNgenes, showCategCnet, geneList) {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            if(is.null(data[[i]])){
                output[[sprintf("famCnet%01d",i)]] <- renderPlot(plotOnlyText())
            }else if(nrow(data[[i]]@result[!is.na(data[[i]]@result$Description),]) == 0){
                output[[sprintf("famCnet%01d",i)]] <- renderPlot(plotOnlyText())
            }else{
                dataToShow <- data[[i]]
                signToShow <- unname(unlist(lapply(famData[[i]]$Descriptions, 
                                                   function(x) x[1])))
                dataToShow@result <- dataToShow@result[dataToShow@result$Description %in% 
                                                           signToShow, ]
                dataToShow <- limitNgenes(dataToShow, maxNgenes = maxNgenes)
                
                output[[sprintf("famCnet%01d",i)]] <- 
                    renderPlot(cnetplot(dataToShow, showCategory = showCategCnet, 
                                        foldChange = geneList[[i]]))
            }
        })
    })
}

famBarplotsUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            plotOutput(NS(id, sprintf("famBarplot%01d",i)))
        )
    })
}

famBarplotsServer <- function(id, pairs_df, data, famData, showCateg) {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            if(is.null(data[[i]])){
                output[[sprintf("famBarplot%01d",i)]] <- renderPlot(plotOnlyText())
            }else if(nrow(data[[i]]@result[!is.na(data[[i]]@result$Description),]) == 0){
                output[[sprintf("famBarplot%01d",i)]] <- renderPlot(plotOnlyText())
            }else{
                dataToShow <- data[[i]]
                signToShow <- unname(unlist(lapply(famData[[i]]$Descriptions, 
                                                   function(x) x[1])))
                dataToShow@result <- dataToShow@result[dataToShow@result$Description %in% 
                                                           signToShow, ]
                output[[sprintf("famBarplot%01d",i)]] <- 
                    renderPlot(suppressMessages(barplot(dataToShow, showCategory=showCateg) +
                                                    scale_fill_continuous(low="red", high="gray", 
                                                                          name="p.adjust", 
                                                                          guide=guide_colorbar(reverse=TRUE))))
            }
        })
    })
}

famTablesUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            tableOutput(NS(id, sprintf("famTable%01d",i)))
        )
    })
}

famTablesServer <- function(id, pairs_df, data, 
                            styling_opt, DE_df_height) {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            if(nrow(data[[i]]) == 0 | all(is.na(data[[i]]))){
                dataToShow <- data.frame(Data = "Zero significant")
            }else{
                dataToShow <- cbind(data.frame(N = 1:nrow(data[[i]])), 
                                    as.data.frame(data[[i]]))
            }
            output[[sprintf("famTable%01d",i)]] <- function() {
                dataToShow %>%
                    knitr::kable(format = "html", row.names = FALSE) %>%
                    kable_styling(styling_opt, full_width = T) %>%
                    scroll_box(width = "100%", height = DE_df_height)
            }
        })
    })
}

famFGSEACnetplotsUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            plotOutput(NS(id, sprintf("famFGSEACnet%01d",i)))
        )
    })
}

famFGSEACnetplotsServer <- function(id, pairs_df, fgsea_obj, fgsea_table, 
                               enrich_example, dds,
                               maxNgenes, showCategCnet, geneList, upDown = "Down") {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            initialData <- fgsea_obj[[i]]
            if(upDown == "Down"){
                initialData <- initialData[initialData$NES < 0, ]
            }else if(upDown == "Up"){
                initialData <- initialData[initialData$NES > 0, ]
            }else{
                stop("Wrong direction. Only Up and Down are valid.")
            }
            
            if(is.null(initialData)){
                output[[sprintf("famFGSEACnet%01d",i)]] <- renderPlot(plotOnlyText())
            }else if(nrow(initialData[!is.na(initialData$pathway),]) == 0){
                output[[sprintf("famFGSEACnet%01d",i)]] <- renderPlot(plotOnlyText())
            }else{
                signToShow <- unname(unlist(lapply(fgsea_table[[i]]$Pathways, 
                                                   function(x) x[1])))
                initialData <- initialData[initialData$pathway %in% 
                                               signToShow, ]
                
                output[[sprintf("famFGSEACnet%01d",i)]] <- 
                    renderPlot(cnetplot(fgseaToEnrichResult(initialData, 
                                                            enrichExample = enrich_example[[1]],
                                                            workdir,
                                                            as.data.frame(rowData(dds)),
                                                            forPlot = "cnetplot",
                                                            maxNgenes = maxNgenes), 
                                        showCategory = showCategCnet, 
                                        foldChange = geneList[[i]]))
            }
        })
    })
}

famFGSEAbarplotsUI <- function(id, pairs_df, upDown = "Down") {
    lapply(1:nrow(pairs_df), function(i){
        tagList(
            h5(paste0(upDown, " in ", pairs_df$cond2[i], 
                      " vs ", pairs_df$cond1[i], ".")),
            plotOutput(NS(id, sprintf("famFGSEAbarplot%01d",i)))
        )
    })
}

famFGSEAbarplotsServer <- function(id, pairs_df, fgsea_obj, fgsea_table, 
                                    enrich_example, dds,
                                    maxNgenes, showCategCnet, geneList, upDown = "Down") {
    moduleServer(id, function(input, output, session) {
        lapply(1:nrow(pairs_df), function(i){
            initialData <- fgsea_obj[[i]]
            if(upDown == "Down"){
                initialData <- initialData[initialData$NES < 0, ]
            }else if(upDown == "Up"){
                initialData <- initialData[initialData$NES > 0, ]
            }else{
                stop("Wrong direction. Only Up and Down are valid.")
            }
            
            initialData <- initialData[order(initialData$padj), ]
            if(is.null(initialData)){
                output[[sprintf("famFGSEAbarplot%01d",i)]] <- renderPlot(plotOnlyText())
            }else if(nrow(initialData[!is.na(initialData$pathway),]) == 0){
                output[[sprintf("famFGSEAbarplot%01d",i)]] <- renderPlot(plotOnlyText())
            }else{
                signToShow <- unname(unlist(lapply(fgsea_table[[i]]$Pathways, 
                                                   function(x) x[1])))
                initialData <- initialData[initialData$pathway %in% 
                                               signToShow, ]
                
                output[[sprintf("famFGSEAbarplot%01d",i)]] <- 
                    renderPlot(suppressMessages(barplot(fgseaToEnrichResult(initialData,
                                                                            enrichExample = enrich_example[[1]], 
                                                                            workdir,
                                                                            as.data.frame(rowData(dds)),
                                                                            forPlot = "barplot",
                                                                            maxNgenes = maxNgenes), 
                                                        showCategory=showCateg) +
                                                    scale_fill_continuous(low="red", high="gray", 
                                                                          name="p.adjust", 
                                                                          guide=guide_colorbar(reverse=TRUE))))
            }
        })
    })
}

# Shiny App
showShiny <- function(){
    
    #!/usr/bin/env Rscript
    # Change R.home("library") to your path to R packages if you use RStudio
    #.libPaths(R.home("library"))
    
    #workdir <- "/data/cabezas/group/pavlovich/data/shared/storshow/runs/outputdir140"
    #source("/data/cabezas/group/pavlovich/data/shared/storshow/R/dev_Polina/Storshow/scripts/readYaml.R")
    #source("/data/cabezas/group/pavlovich/data/shared/storshow/R/dev_Polina/Storshow/scripts/4_bulk_RNA_seq_ShinyInternal.R")
    
    #storshow_functions <- "/data/cabezas/group/pavlovich/data/shared/storshow/R/dev_Polina/Storshow/scripts/storshow_functions.R"
    #gmt.file <- "/data/cabezas/group/pavlovich/data/shared/storshow/R/dev_Polina/Storshow/scripts/signatures/Signatures_for_GSEA_updated_20200714.gmt"
    
    stopifnot(length(ctr_cond) == length(treat_cond))
    # titlePage should be changed inside the app
    titlePage <- ""
    
    #source("/data/cabezas/group/pavlovich/data/shared/storshow/storshow_functions.R")
    # storshow_functions is specified in defaults.yaml file
    source(storshow_functions)
    
    subdir <- "for_DAG"
    
    fn_df <- read.delim(file.path(workdir, "DESeq/fn_geneClTrend.tsv"),
                        stringsAsFactors = FALSE)
    
    dds <- readRDS(file.path(workdir, "DESeq", fn_df$Name[fn_df$Type == "dds"]))
    coldata <- as.data.frame(SummarizedExperiment::colData(dds))
    
    # style_file colmes from defaults.yaml
    levels_colors <- parseStyle(style_file, coldata)
    levels_ordered <- levels_colors[[1]]
    colors_ordered <- levels_colors[[2]]
    
    # extra_data_dds is specified in JSON file
    #if(extra_data_dds != ""){
    #    dds_publ <- readRDS(extra_data_dds)
    #}
    
    padj_df <- read.delim(file.path(workdir, fn_df$Name[fn_df$Type == "padj_df"]),
                          stringsAsFactors = FALSE)

    #cond_in_right_order <- c("Ctrl_DMSO_vs_Ctrl_4oxoRA",
    #                         "Ctrl_DMSO_vs_KO_4oxoRA",
    #                         "KO_DMSO_vs_KO_4oxoRA",
    #                         "Ctrl_4oxoRA_vs_KO_4oxoRA", # 4oxoRA
    #                         "Ctrl_DMSO_vs_Ctrl_atRA",
    #                         "Ctrl_DMSO_vs_KO_atRA",
    #                         "KO_DMSO_vs_KO_atRA",
    #                         "Ctrl_atRA_vs_KO_atRA", # atRA
    #                         "Ctrl_DMSO_vs_Ctrl_retinol",
    #                         "Ctrl_DMSO_vs_KO_retinol",
    #                         "KO_DMSO_vs_KO_retinol",
    #                         "Ctrl_retinol_vs_KO_retinol", # retinol
    #                         "Ctrl_DMSO_vs_KO_DMSO")
    pairs_df <- data.frame(pairNum = c(1:length(ctr_cond)),
                           cond1 = ctr_cond,
                           cond2 = treat_cond)
    
    helpInfo_df <- data.frame(Pair = paste(pairs_df$pairNum, 
                                           pairs_df$cond1, "vs",
                                           pairs_df$cond2))
    
    inputDir <- file.path(workdir, "DESeq/DESeq2_results_sign_DE")
    fns <- reorderFns(list.files(inputDir, pattern = "tsv"), 
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    results_list <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        results_list[[i]] <- read.delim(file.path(inputDir, fns[i]), stringsAsFactors = FALSE)
        rownames(results_list[[i]]) <- 1:nrow(results_list[[i]])
    }
    
    for(list_i in 1:length(results_list)){
        results_list[[list_i]][,1:4] <- round(results_list[[list_i]][,1:4], 2)
        
        results_list[[list_i]]$pvalue <- signif(results_list[[list_i]]$pvalue, 3)
        results_list[[list_i]]$pvalue <- format(results_list[[list_i]]$pvalue, digits = 2)
        
        results_list[[list_i]]$padj <- signif(results_list[[list_i]]$padj, 3)
        results_list[[list_i]]$padj <- format(results_list[[list_i]]$padj, digits = 2)
    }
    
    #DE_df_height <- "400px"
    styling_opt <- c("striped", "condensed")
    
    geneClusters_df <- read.delim(file.path(workdir, 
                                            fn_df$Name[fn_df$Type == "Gene_cl_trend"]),
                                  stringsAsFactors = FALSE)
    
    if(any(grepl("so.so", colnames(geneClusters_df)))){
        geneClusters_df <- geneClusters_df[, c("mask", "freq", "allUp", "allDown", "so.so")]
    }else{
        geneClusters_df <- geneClusters_df[, c("mask", "freq", "allUp", "allDown")]
    }
    
    # BP
    bp_egoS_dir <- file.path(workdir, "enrichGO_simpl/BP")
    fns <- reorderFns(list.files(bp_egoS_dir, pattern = "pos.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    bp_ego_s_up <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        bp_ego_s_up[[i]] <- readRDS(file.path(bp_egoS_dir, fns[i]))
        
    }
    fns <- reorderFns(list.files(bp_egoS_dir, pattern = "neg.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    bp_ego_s_down <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        bp_ego_s_down[[i]] <- readRDS(file.path(bp_egoS_dir, fns[i]))
        
    }
    
    # MF
    mf_egoS_dir <- file.path(workdir, "enrichGO_simpl/MF")
    fns <- reorderFns(list.files(mf_egoS_dir, pattern = "pos.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    mf_ego_s_up <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        mf_ego_s_up[[i]] <- readRDS(file.path(mf_egoS_dir, fns[i]))
        
    }
    fns <- reorderFns(list.files(mf_egoS_dir, pattern = "neg.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    mf_ego_s_down <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        mf_ego_s_down[[i]] <- readRDS(file.path(mf_egoS_dir, fns[i]))
        
    }
    
    # CC
    cc_egoS_dir <- file.path(workdir, "enrichGO_simpl/CC")
    fns <- reorderFns(list.files(cc_egoS_dir, pattern = "pos.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    cc_ego_s_up <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        cc_ego_s_up[[i]] <- readRDS(file.path(cc_egoS_dir, fns[i]))
        
    }
    fns <- reorderFns(list.files(cc_egoS_dir, pattern = "neg.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    cc_ego_s_down <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        cc_ego_s_down[[i]] <- readRDS(file.path(cc_egoS_dir, fns[i]))
        
    }
    
    # KEGG
    kegg_dir <- file.path(workdir, "KEGG")
    fns <- reorderFns(list.files(kegg_dir, pattern = "pos.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    kegg_up <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        kegg_up[[i]] <- readRDS(file.path(kegg_dir, fns[i]))
        
    }
    fns <- reorderFns(list.files(kegg_dir, pattern = "neg.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    kegg_down <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        kegg_down[[i]] <- readRDS(file.path(kegg_dir, fns[i]))
        
    }
    
    #showCateg <- 25
    #showCategCnet <- 6
    #maxNgenes <- 7
    
    data(examplePathways) # pathways from fgsea package
    resRDSdir <- file.path(workdir, "DESeq/DESeq2_results_allG_rds")
    
    rowdata <- as.data.frame(rowData(dds))
    rowdata <- rowdata[!duplicated(rowdata$Entrez.Gene.ID),]
    
    fns_allG <- reorderFns(list.files(resRDSdir, pattern = ".rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    allG_filePrefix <- gsub("_all.*", "", fns_allG)
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
    
    #gmt.file <- "/data/cabezas/group/pavlovich/data/shared/storshow/signatures/Signatures_for_GSEA_updated_20200309.gmt"
    pathways <- gmtPathways(gmt.file)
    
    # FGSEA custom
    fgsea_dir <- file.path(workdir, "FGSEA/custom")
    fns <- reorderFns(list.files(fgsea_dir, pattern = "custom.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    fgsea_custom <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        fgsea_custom[[i]] <- readRDS(file.path(fgsea_dir, fns[i]))
        
    }
    
    # FGSEA allPath
    fgsea_dir <- file.path(workdir, "FGSEA/allPath")
    fns <- reorderFns(list.files(fgsea_dir, pattern = "allPath.rds"),
                      ctr_cond = ctr_cond, treat_cond = treat_cond)
    fgsea_allPath <- vector("list", length = length(fns))
    for(i in 1:length(fns)){
        fgsea_allPath[[i]] <- readRDS(file.path(fgsea_dir, fns[i]))
        
    }
    
    # Hints
    bp_upMult <- giveMult(GO_resList = bp_ego_s_up, pairs_df)
    bp_downMult <- giveMult(bp_ego_s_down, pairs_df)
    
    mf_upMult <- giveMult(mf_ego_s_up, pairs_df)
    mf_downMult <- giveMult(mf_ego_s_down, pairs_df)
    
    cc_upMult <- giveMult(cc_ego_s_up, pairs_df)
    cc_downMult <- giveMult(cc_ego_s_down, pairs_df)
    
    kegg_upMult <- giveMult(kegg_up, pairs_df)
    kegg_downMult <- giveMult(kegg_down, pairs_df)
    
    go_kegg_downList <- list(bp_downMult, mf_downMult,
                             cc_downMult, kegg_downMult)
    
    go_kegg_upList <- list(bp_upMult, mf_upMult,
                           cc_upMult, kegg_upMult)
    
    multsText <- c("BP -- Biological Process",
                   "MF -- Molecular Function",
                   "CC -- Cellular Component",
                   "KEGG -- Kyoto Encyclopedia of Genes and Genomes")
    
    fg_customMult <- giveMultFGSEA(FGSEA_resList = fgsea_custom, pairs_df)
    fg_allPathMult <- giveMultFGSEA(fgsea_allPath, pairs_df)
    
    fg_downMultList <- list(fg_customMult[[1]], fg_allPathMult[[1]])
    fg_upMultList <- list(fg_customMult[[2]], fg_allPathMult[[2]])
    
    FGSEAmultsText <- c("Custom signatures",
                        "All available pathways from fgsea r package")
    
    # Fam FGSEA custom
    dir_up <- file.path(workdir, "families/FGSEA/custom/up/tables")
    fns_up <- reorderFns(list.files(dir_up, pattern = "up.rds"),
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    famFGSEAcustom_up <- lapply(fns_up, function(x) readRDS(file.path(dir_up, x)))
    
    dir_down <- file.path(workdir, "families/FGSEA/custom/down/tables")
    fns_down <- reorderFns(list.files(dir_down, pattern = "down.rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    famFGSEAcustom_down <- lapply(fns_down, function(x) readRDS(file.path(dir_down, x)))
    
    # Fam FGSEA allPath
    dir_up <- file.path(workdir, "families/FGSEA/allPath/up/tables")
    fns_up <- reorderFns(list.files(dir_up, pattern = "up.rds"),
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    famFGSEAallPath_up <- lapply(fns_up, function(x) readRDS(file.path(dir_up, x)))
    
    dir_down <- file.path(workdir, "families/FGSEA/allPath/down/tables")
    fns_down <- reorderFns(list.files(dir_down, pattern = "down.rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    famFGSEAallPath_down <- lapply(fns_down, function(x) readRDS(file.path(dir_down, x)))
    
    # Fam BP
    dir_up <- file.path(workdir, "families/BP/up/tables")
    fns_up <- reorderFns(list.files(dir_up, pattern = "up.rds"),
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    famBP_up <- lapply(fns_up, function(x) readRDS(file.path(dir_up, x)))
    
    dir_down <- file.path(workdir, "families/BP/down/tables")
    fns_down <- reorderFns(list.files(dir_down, pattern = "down.rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    famBP_down <- lapply(fns_down, function(x) readRDS(file.path(dir_down, x)))
    
    # Fam MF
    dir_up <- file.path(workdir, "families/MF/up/tables")
    fns_up <- reorderFns(list.files(dir_up, pattern = "up.rds"),
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    famMF_up <- lapply(fns_up, function(x) readRDS(file.path(dir_up, x)))
    
    dir_down <- file.path(workdir, "families/MF/down/tables")
    fns_down <- reorderFns(list.files(dir_down, pattern = "down.rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    famMF_down <- lapply(fns_down, function(x) readRDS(file.path(dir_down, x)))
    
    # Fam CC
    dir_up <- file.path(workdir, "families/CC/up/tables")
    fns_up <- reorderFns(list.files(dir_up, pattern = "up.rds"),
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    famCC_up <- lapply(fns_up, function(x) readRDS(file.path(dir_up, x)))
    
    dir_down <- file.path(workdir, "families/CC/down/tables")
    fns_down <- reorderFns(list.files(dir_down, pattern = "down.rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    famCC_down <- lapply(fns_down, function(x) readRDS(file.path(dir_down, x)))
    
    # Fam KEGG
    dir_up <- file.path(workdir, "families/KEGG/up/tables")
    fns_up <- reorderFns(list.files(dir_up, pattern = "up.rds"),
                         ctr_cond = ctr_cond, treat_cond = treat_cond)
    famKEGG_up <- lapply(fns_up, function(x) readRDS(file.path(dir_up, x)))
    
    dir_down <- file.path(workdir, "families/KEGG/down/tables")
    fns_down <- reorderFns(list.files(dir_down, pattern = "down.rds"),
                           ctr_cond = ctr_cond, treat_cond = treat_cond)
    famKEGG_down <- lapply(fns_down, function(x) readRDS(file.path(dir_down, x)))
    
    # For plotBarPlots, levels_ordered should be provided in JSON file
    
    #colorSwitch <- function(inputCond){
    #    return(switch(inputCond, 
    #                  "Ctrl_4oxoRA" = "#A6CEE3", 
    #                  "Ctrl_atRA" = "#B2DF8A", 
    #                  "Ctrl_DMSO" = "#FB9A99",
    #                  "Ctrl_retinol" = "#FDBF6F", 
    #                  "KO_4oxoRA" = "#1F78B4",
    #                  "KO_atRA" = "#33A02C", 
    #                  "KO_DMSO" = "#E31A1C",
    #                  "KO_retinol" = "#FF7F00"))
    #}
    
    colorSwitch <- data.frame(levels_ordered, colors_ordered,
                              stringsAsFactors = F)
    
    showStates <- list()
    # one list item per pair
    for(i in 1:length(ctr_cond)){
        showStates[[i]] <- c(ctr_cond[i], treat_cond[i])
    }
    
    # height MA & Volcano
    height_MA_volcano <- paste0(335 * nrow(pairs_df), "px")
    
    if(front_target_barplot_height == "auto"){
        # x is n_rows, y is height
        x_min <- 1
        x_max <- 4
        
        y_min <- 180
        y_max <- 600
        
        slope <- round((y_max - y_min) / (x_max - x_min), 3)
        intercept <- round(0.5 * (y_min + y_max - slope * (x_min + x_max)), 3)
        
        # formula from cowplot plot_grid github
        n_cols <- ceiling(sqrt(nrow(pairs_df)))
        n_rows <- ceiling(nrow(pairs_df)/n_cols)
        
        front_target_barplot_height <- slope * n_rows + intercept
    }
    
    if(front_target_barplot_width == "auto"){
        # x is n_cols, y is width
        x_min <- 1
        x_max <- 4
        
        y_min <- 110
        y_max <- 415
        
        slope <- round((y_max - y_min) / (x_max - x_min), 3)
        intercept <- round(0.5 * (y_min + y_max - slope * (x_min + x_max)), 3)
        
        # formula from cowplot plot_grid github
        n_cols <- ceiling(sqrt(nrow(pairs_df)))
        
        front_target_barplot_width <- slope * n_cols + intercept
    }
    
    ui <- shinyUI(fluidPage(
        navbarPage(titlePage,
                   tabPanel("Gene expr",
                            sidebarLayout(
                                sidebarPanel(
                                    selectizeInput('gene', 'Select a gene', 
                                                   choices = c("", rownames(dds)),
                                                   options = list(
                                                       placeholder = 'Please select an option below',
                                                       onInitialize = 
                                                           I(paste0("function() { this.setValue('", 
                                                                    first_gene_to_show, "'); }"))
                                                   )
                                    ),
                                    actionButton("saveBarPlots", "Save plots"),
                                    width = 3),
                                
                                mainPanel("Please open the window full-size.",
                                          helpText("Thresholds used for building this app.",
                                                   "Gene filtering: at least one condition has an average of",
                                                   paste0("more than ", geneFilt, " reads."),
                                                   paste0("It left ", nrow(dds), " genes."),
                                                   paste0("DESeq2: betaPrior ", betaPrior, ","),
                                                   paste0("padj < ", DESeq2_padj_thr, ","),
                                                   paste0("abs(log2FC) > ", DESeq2_absLog2FC_thr, "."),
                                                   paste0("ORA: padj < ", 
                                                          ORA_padj_thr, "."),
                                                   paste0("GSEA: padj < ", 
                                                          GSEA_padj_thr, "."),
                                                   paste0("Families: at least ", avg_pair_overl_thr,
                                                          "% of genes overlap"),
                                                   "between two signatures."),
                                          helpText(extra_data_comment),
                                          verbatimTextOutput("savingBarPlots"),
                                          helpText("Note: 'ns' stays for non-significant,", 
                                                   "stars mean that the comparison was significant",
                                                   "with:",
                                                   "**** padj < 0.0001,",
                                                   "*** padj < 0.001,",
                                                   "** padj < 0.01,",
                                                   "* padj < 0.1."),
                                          #helpText("All saved plots from this app will be",
                                          #         "on maximus in",
                                          #         "/data/cabezas/group/pavlovich/data/ShinyApps."),
                                          uiOutput("front_target_barplot"),
                                          #uiOutput("front_extra_barplot"),
                                          width = 9)
                            )
                   ),
                   navbarMenu("DE genes",
                              tabPanel("MA & Volcano",
                                       helpText("MA plots (left) show non-significant genes as black,",
                                                "significant -- as red.",
                                                paste0("Threshold for significance: ",
                                                       "padj < ", DESeq2_padj_thr, "."),
                                                "In MA plots, no threshold for log2FC is applied.",
                                                paste0("Only top ", Volcano_nShow_geneNames, 
                                                       " up and down significant ",
                                                "DE genes are labeled.")),
                                       splitLayout(plotOutput("MAplots", height = height_MA_volcano, 
                                                              width = MA_plots_width),
                                                   plotOutput("VolcanoPlots", height = height_MA_volcano, 
                                                              width = Volcano_plots_width))
                              ),
                              tabPanel("DE tables", 
                                       helpText("Thresholds for significance: padj < 0.1 and",
                                                "no threshld for log2FC.",
                                                "Top right table contains number of significant genes",
                                                "in each comparison. Column names correspond to pairNum",
                                                "in the top left table."),
                                       #splitLayout(cellWidths = c("50%", "50%"),
                                       tableOutput("pairs_DEcounts"),
                                       #tableOutput("pairs"),
                                       #tableOutput("DE_counts")),
                                       helpText("DE tables contain only significant genes.",
                                                "Negative logFC shows up-regulated genes in cond1,",
                                                "positive logFC (scroll to the bottom of the tables)",
                                                "is for up-regulated in cond2."),
                                       uiOutput("DE_genes")
                              )
                   ),
                   navbarMenu("ORA",
                              tabPanel("BP",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("BP cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(cnetplotsUI("BPdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        cnetplotsUI("BPupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("BP tables", 
                                                            splitLayout(tablesUI("BPdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        tablesUI("BPupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("BP barplots", 
                                                            helpText(paste0("A maximum of ", 
                                                                showCateg, " categories is displayed.")),
                                                            splitLayout(barplotsUI("BPdownBarplots", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        barplotsUI("BPupBarplots", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up")))
                                       )
                              ),
                              tabPanel("MF", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("MF cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(cnetplotsUI("MFdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        cnetplotsUI("MFupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("MF tables", 
                                                            splitLayout(tablesUI("MFdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        tablesUI("MFupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("MF barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(barplotsUI("MFdownBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Down"),
                                                                        barplotsUI("MFupBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Up")))
                                       )  
                              ),
                              tabPanel("CC", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("CC cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(cnetplotsUI("CCdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        cnetplotsUI("CCupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("CC tables", 
                                                            splitLayout(tablesUI("CCdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        tablesUI("CCupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("CC barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(barplotsUI("CCdownBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Down"),
                                                                        barplotsUI("CCupBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Up")))  
                                       )    
                              ),
                              tabPanel("KEGG", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("KEGG cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(cnetplotsUI("KEGGdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        cnetplotsUI("KEGGupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("KEGG tables", 
                                                            splitLayout(tablesUI("KEGGdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        tablesUI("KEGGupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("KEGG barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(barplotsUI("KEGGdownBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Down"),
                                                                        barplotsUI("KEGGupBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Up")))
                                       )    
                              )
                   ),
                   navbarMenu("FGSEA",
                              tabPanel("Custom",
                                       sidebarPanel(
                                           selectizeInput('signat', 'Select a signature', 
                                                          choices = names(pathways),
                                                          options = list(
                                                              placeholder = 'Please select an option below',
                                                              onInitialize = I(paste0("function() { this.setValue('", 
                                                                                      first_custom_GSEA_to_show, "'); }"))
                                                          )
                                           ),
                                           actionButton("saveCustomPlots", "Save plots"),
                                           helpText("In the menu, there are signatures that are",
                                                    "significant at least in one comparison."),
                                           br(),
                                           helpText("Tables show only signatures with"),
                                           helpText(paste0("padj < ", GSEA_padj_thr, ".")),
                                           helpText("'padj' is adjusted 'pval'. The FGSEA package provides both."),
                                           br(),
                                           helpText("'nMoreExtreme' is a number of times a random gene set",
                                                    "had a more extreme enrichment score value, and",
                                                    "'size' is the number of genes from the pathway detected in the dataset."),
                                           width = 3
                                       ),
                                       mainPanel(
                                           verbatimTextOutput("savingCustomPlots"),
                                           splitLayout(uiOutput("FGSEAcustomPlots"),
                                                       uiOutput("FGSEAcustomT")),
                                           width = 9
                                       )       
                              ),
                              tabPanel("All pathways", 
                                       sidebarPanel(
                                           selectizeInput('fgseaPath', 'Select a signature', 
                                                          choices = names(examplePathways)[which((lengths(examplePathways) >= 15))],
                                                          options = list(
                                                              placeholder = 'Please select an option below',
                                                              onInitialize = 
                                                                  I(paste0("function() { this.setValue('", 
                                                                           first_allPath_GSEA_to_show, "'); }"))
                                                          )
                                           ),
                                           actionButton("saveAllPathPlots", "Save plots"),
                                           helpText("In the menu, there are signatures that are",
                                                    "significant at least in one comparison."),
                                           br(),
                                           helpText("Tables show only signatures with"),
                                           helpText(paste0("padj < ", GSEA_padj_thr, ".")),
                                           helpText("'padj' is adjusted 'pval'. The FGSEA package provides both."),
                                           br(),
                                           helpText("'nMoreExtreme' is a number of times a random gene set",
                                                    "had a more extreme enrichment score value, and",
                                                    "'size' is the number of genes from the pathway detected in the dataset."),
                                           width = 3
                                       ),
                                       mainPanel(
                                           verbatimTextOutput("savingAllPathPlots"),
                                           splitLayout(uiOutput("FGSEAallPathPlots"),
                                                       uiOutput("FGSEAallPathT")),
                                           width = 9
                                       )
                              )
                   ),
                   navbarMenu("Families",
                              tabPanel("FGSEA custom",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Fam FGSEA custom cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " genes per group is shown."),
                                                            splitLayout(famFGSEACnetplotsUI("famFGSEAcustomDownCnet", 
                                                                                       pairs_df = pairs_df,
                                                                                       upDown = "Down"),
                                                                        famFGSEACnetplotsUI("famFGSEAcustomUpCnet", 
                                                                                       pairs_df = pairs_df,
                                                                                       upDown = "Up"))),
                                                   tabPanel("Fam FGSEA custom tables", 
                                                            splitLayout(famTablesUI("famFGSEAcustomDownTables", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        famTablesUI("famFGSEAcustomUpTables", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("Fam FGSEA custom barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(famFGSEAbarplotsUI("famFGSEAcustomDownBarplots", 
                                                                                           pairs_df = pairs_df,
                                                                                           upDown = "Down"),
                                                                        famFGSEAbarplotsUI("famFGSEAcustomUpBarplots", 
                                                                                           pairs_df = pairs_df,
                                                                                           upDown = "Up")))
                                       )
                              ),
                              tabPanel("FGSEA allPath", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Fam FGSEA allPath cnetplots", 
                                                            helpText("A maximum of ", showCategCnet, 
                                                                     " categories is displayed. ",
                                                                     "A maximum of ", maxNgenes, 
                                                                     " genes per group is shown.",
                                                                     "Long names are truncated after the 50th character."),
                                                            splitLayout(famFGSEACnetplotsUI("famFGSEAallPathDownCnet", 
                                                                                            pairs_df = pairs_df,
                                                                                            upDown = "Down"),
                                                                        famFGSEACnetplotsUI("famFGSEAallPathUpCnet", 
                                                                                            pairs_df = pairs_df,
                                                                                            upDown = "Up"))),
                                                   tabPanel("Fam FGSEA allPath tables", 
                                                            splitLayout(famTablesUI("famFGSEAallPathDownTables", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        famTablesUI("famFGSEAallPathUpTables", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("Fam FGSEA allPath barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed. "),
                                                                     "Long names are truncated after the 50th character."),
                                                            splitLayout(famFGSEAbarplotsUI("famFGSEAallPathDownBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Down"),
                                                                        famFGSEAbarplotsUI("famFGSEAallPathUpBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Up")))
                                       )
                              ),
                              tabPanel("BP", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Fam BP cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(famCnetplotsUI("famBPdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        famCnetplotsUI("famBPupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("Fam BP tables", 
                                                            splitLayout(famTablesUI("famBPdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        famTablesUI("famBPupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("Fam BP barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(famBarplotsUI("famBPdownBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Down"),
                                                                        famBarplotsUI("famBPupBarplots", 
                                                                                   pairs_df = pairs_df,
                                                                                   upDown = "Up")))
                                       )
                              ),
                              tabPanel("MF", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Fam MF cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(famCnetplotsUI("famMFdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        famCnetplotsUI("famMFupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("Fam MF tables", 
                                                            splitLayout(famTablesUI("famMFdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        famTablesUI("famMFupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("Fam MF barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(famBarplotsUI("famMFdownBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Down"),
                                                                        famBarplotsUI("famMFupBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Up")))
                                       )
                              ),
                              tabPanel("CC", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Fam CC cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(famCnetplotsUI("famCCdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        famCnetplotsUI("famCCupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("Fam CC tables", 
                                                            splitLayout(famTablesUI("famCCdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        famTablesUI("famCCupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("Fam CC barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(famBarplotsUI("famCCdownBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Down"),
                                                                        famBarplotsUI("famCCupBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Up")))
                                       )
                              ),
                              tabPanel("KEGG", 
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Fam KEGG cnetplots", 
                                                            helpText(paste0("A maximum of ", showCategCnet, 
                                                                            " categories is displayed. "),
                                                                     "A maximum of ", maxNgenes, 
                                                                     " significant genes in a group is shown."),
                                                            splitLayout(famCnetplotsUI("famKEGGdownCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Down"),
                                                                        famCnetplotsUI("famKEGGupCnet", 
                                                                                    pairs_df = pairs_df,
                                                                                    upDown = "Up"))),
                                                   tabPanel("Fam KEGG tables", 
                                                            splitLayout(famTablesUI("famKEGGdownTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Down"),
                                                                        famTablesUI("famKEGGupTables", 
                                                                                 pairs_df = pairs_df,
                                                                                 upDown = "Up"))),
                                                   tabPanel("Fam KEGG barplots", 
                                                            helpText(paste0("A maximum of ", showCateg, 
                                                                            " categories is displayed.")),
                                                            splitLayout(famBarplotsUI("famKEGGdownBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Down"),
                                                                        famBarplotsUI("famKEGGupBarplots", 
                                                                                      pairs_df = pairs_df,
                                                                                      upDown = "Up")))
                                       )
                              )
                   )
        )
    ))
    
    server <- shinyServer(function(input, output, session) {
        #inputGene <- eventReactive(input$goButton, {
        #    input$gene
        #}, ignoreNULL = FALSE)
        
        observeEvent(input$gene, {
            textMessage <- NULL
        })
        
        ##### Page 1 Barplots ######
        output$plot1 <- renderCachedPlot({
            #plotSevenPlots(inputGene(), dds, padj_df)
            if(!(input$gene %in% rownames(dds))){
                plotOnlyText(textToPlot = "Select a gene")
            }else{
                plotBarPlots(input$gene, dds, padj_df,
                             levels_ordered, 
                             colorSwitch, 
                             showStates)
            }
        }, 
        #cacheKeyExpr = {inputGene()}
        cacheKeyExpr = {input$gene}
        )
        
        # wrap plotOutput in renderUI
        output$front_target_barplot <- renderUI({
            plotOutput("plot1", height = front_target_barplot_height, 
                       width = front_target_barplot_width)
        })
        
        #output$plot2 <- renderCachedPlot({
        #    #publBar(inputGene(), dds_publ)
        #    if(!(input$gene %in% rownames(dds_publ))){
        #        plotOnlyText(textToPlot = "Select a gene")
        #    }else{
        #        publBar(input$gene, dds_publ)
        #    }
        #},
        #cacheKeyExpr = {inputGene()})
        #cacheKeyExpr = {input$gene})
        
        # wrap plotOutput in renderUI
        #output$front_extra_barplot <- renderUI({
        #    plotOutput("plot2", height = front_extra_barplot_height, 
        #               width = front_extra_barplot_width)
        #})
        
        observeEvent(input$saveBarPlots, {
            textMessage <- paste0("The barplots with ", 
                                  input$gene, " expression are saved in", 
                                  " 'data/app_plots/barplots' folder.")
            output$savingBarPlots <- renderText({ textMessage })
            
            dir.create(file.path(workdir, "app_plots"), showWarnings = F)
            dir.create(file.path(workdir, "app_plots/barplots"), showWarnings = F)
            pdf(file.path(workdir, "app_plots/barplots",
                          paste0(filePrefix, input$gene, "_barplot.pdf")),
                width=1.5, height=2.8, onefile=FALSE)
            print(plotBarPlots(input$gene, dds, padj_df,
                               levels_ordered, 
                               colorSwitch, 
                               showStates))
            dev.off()
            
            #pdf(file.path(workdir, "app_plots/barplots",
            #              paste0(filePrefix, input$gene, "_publishedData_barplot.pdf")),
            #    width=7.5, height=4, onefile=FALSE)
            #print(publBar(input$gene, dds_publ))
            #dev.off()
        })
        
        ##### Page 2 DE genes ######
        output$MAplots <- renderPlot({
            giveMAplots(file.path(workdir, "DESeq/DESeq2_results_allG_rds"),
                        samp = filePrefix,
                        ylim = MA_plots_ylim,
                        DESeq2_padj_thr = DESeq2_padj_thr,
                        ctr_cond = ctr_cond, treat_cond = treat_cond)
        })
        
        output$VolcanoPlots <- renderPlot({
            plotVolcanoPlots(results_allG, pairs_df,
                             nShow = Volcano_nShow_geneNames, 
                             DESeq2_padj_thr = DESeq2_padj_thr, 
                             DESeq2_absLog2FC_thr = DESeq2_absLog2FC_thr)
        })
        
        output$pairs <- function() {
            pairs_df %>%
                knitr::kable(format = "html", row.names = TRUE) %>%
                kable_styling(styling_opt, full_width = F,
                              position = "right")
        }
        
        output$DE_counts <- function() {
            giveDEcounts(file.path(workdir, "DESeq/DESeq2_results_sign_DE"),
                         ctr_cond = ctr_cond, treat_cond = treat_cond) %>%
                knitr::kable(format = "html", row.names = TRUE) %>%
                kable_styling(styling_opt, full_width = F,
                              position = "left")
        }
        
        output$pairs_DEcounts <- function() {
            DE_counts_df <- giveDEcounts(file.path(workdir, "DESeq/DESeq2_results_sign_DE"),
                                         ctr_cond = ctr_cond, treat_cond = treat_cond)
            pairs_DEcounts_df <- merge(pairs_df, DE_counts_df,
                                       by.x="pairNum", by.y="pairNum",
                                       all.x=TRUE, all.y=TRUE, sort=FALSE)
            pairs_DEcounts_df %>%
                knitr::kable(format = "html", row.names = TRUE) %>%
                kable_styling(styling_opt, full_width = F)
        }
        
        # DE genes tables
        lapply(1:nrow(pairs_df), function(i){
            output[[sprintf("PairNum%01d",i)]] <- function() {
                results_list[[i]] %>%
                    knitr::kable(format = "html", row.names = TRUE) %>%
                    kable_styling(styling_opt, full_width = F) %>%
                    scroll_box(width = "100%", height = DE_df_height)
            }
        })
        
        output$DE_genes <- renderUI({
            lapply(1:nrow(pairs_df), function(i){
                tagList(
                    h5(paste0("Cond1 is ", pairs_df$cond1[i], ", 
                                cond2 is ", pairs_df$cond2[i], ".")),
                    tableOutput(sprintf("PairNum%01d",i))
                )
            })
        })
        
        ##### Page 3 Gene clusters ######
        output$pairs2 <- function() {
            pairs_df %>%
                knitr::kable(format = "html", row.names = TRUE) %>%
                kable_styling(styling_opt, full_width = F,
                              position = "center")
        }
        #output$geneClusters <- function() {
        #    geneClusters_df %>%
        #        knitr::kable(format = "html", row.names = TRUE) %>%
        #        kable_styling(styling_opt, full_width = F)
        #}
        
        output$geneClusters <- renderDataTable({geneClusters_df}, 
                                        options = list(scrollX = TRUE))
        
        ##### Page 4 BP enrichGO simplified #####
        # BP barplots down
        barplotsServer("BPdownBarplots", pairs_df = pairs_df, 
                       data = bp_ego_s_down, showCateg = showCateg)
        
        # BP barplots up
        barplotsServer("BPupBarplots", pairs_df = pairs_df,
                       data = bp_ego_s_up, showCateg = showCateg)
        
        # BP cnetplots Down
        cnetplotsServer("BPdownCnet", pairs_df = pairs_df, data = bp_ego_s_down, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # BP cnetplots up
        cnetplotsServer("BPupCnet", pairs_df = pairs_df, data = bp_ego_s_up, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # BP tables Down
        tablesServer("BPdownTables", pairs_df = pairs_df, data = bp_ego_s_down, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # BP tables Down
        tablesServer("BPupTables", pairs_df = pairs_df, data = bp_ego_s_up, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        ##### Page 5 MF enrichGO simplified ######
        # MF barplots down
        barplotsServer("MFdownBarplots", pairs_df = pairs_df, 
                       data = mf_ego_s_down, showCateg = showCateg)
        
        # MF barplots up
        barplotsServer("MFupBarplots", pairs_df = pairs_df,
                       data = mf_ego_s_up, showCateg = showCateg)
        
        # MF cnetplots Down
        cnetplotsServer("MFdownCnet", pairs_df = pairs_df, data = mf_ego_s_down, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # MF cnetplots up
        cnetplotsServer("MFupCnet", pairs_df = pairs_df, data = mf_ego_s_up, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # MF tables Down
        tablesServer("MFdownTables", pairs_df = pairs_df, data = mf_ego_s_down, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # MF tables Down
        tablesServer("MFupTables", pairs_df = pairs_df, data = mf_ego_s_up, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        ##### Page 6 CC enrichGO simplified #####
        # CC barplots down
        barplotsServer("CCdownBarplots", pairs_df = pairs_df, 
                       data = cc_ego_s_down, showCateg = showCateg)
        
        # CC barplots up
        barplotsServer("CCupBarplots", pairs_df = pairs_df,
                       data = cc_ego_s_up, showCateg = showCateg)
        
        # CC cnetplots Down
        cnetplotsServer("CCdownCnet", pairs_df = pairs_df, data = cc_ego_s_down, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # CC cnetplots up
        cnetplotsServer("CCupCnet", pairs_df = pairs_df, data = cc_ego_s_up, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # CC tables Down
        tablesServer("CCdownTables", pairs_df = pairs_df, data = cc_ego_s_down, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # CC tables Down
        tablesServer("CCupTables", pairs_df = pairs_df, data = cc_ego_s_up, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        ##### Page 7 KEGG pathways ######
        # KEGG barplots down
        barplotsServer("KEGGdownBarplots", pairs_df = pairs_df, 
                       data = kegg_down, showCateg = showCateg)
        
        # KEGG barplots up
        barplotsServer("KEGGupBarplots", pairs_df = pairs_df,
                       data = kegg_up, showCateg = showCateg)
        
        # KEGG cnetplots Down
        cnetplotsServer("KEGGdownCnet", pairs_df = pairs_df, data = kegg_down, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # KEGG cnetplots up
        cnetplotsServer("KEGGupCnet", pairs_df = pairs_df, data = kegg_up, 
                        maxNgenes = maxNgenes, 
                        showCategCnet = showCategCnet, 
                        geneList = geneList)
        
        # KEGG tables Down
        tablesServer("KEGGdownTables", pairs_df = pairs_df, data = kegg_down, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # KEGG tables Down
        tablesServer("KEGGupTables", pairs_df = pairs_df, data = kegg_up, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        ##### Page 8 FGSEA custom signatures #####
        # FGSEA custom signatures plots
        lapply(1:nrow(pairs_df), function(i){
            output[[sprintf("FGSEAcustomPlot%01d",i)]] <- 
                renderPlot({
                    if(!(input$signat %in% names(pathways))){
                        plotOnlyText(textToPlot = "Select a signature")
                    }else if(length(names(symbList[[i]])[names(symbList[[i]]) %in% 
                                                         pathways[[input$signat]]]) < 5){
                        plotOnlyText(textToPlot = "< 5 genes found, try another signature")
                    }else{
                        plotEnrichment(pathways[[input$signat]], symbList[[i]]) + 
                            labs(title=input$signat)
                    }
                })
        })
        
        lapply(1:nrow(pairs_df), function(i){
            output[[sprintf("FGSEAcustomPval%01d",i)]] <- 
                renderPrint({
                    if(any(fgsea_custom[[i]]$pathway %in% input$signat)){
                        paste0("padj = ", 
                               signif(fgsea_custom[[i]]$padj[fgsea_custom[[i]]$pathway %in% 
                                                                   input$signat], 2),
                               ", NES = ", 
                               signif(fgsea_custom[[i]]$NES[fgsea_custom[[i]]$pathway %in% 
                                                                input$signat], 2))
                        
                    }else{
                        paste0("padj >= ", GSEA_padj_thr, ", non-significant.")
                    }
                })
        })
        
        output$FGSEAcustomPlots <- renderUI({
            lapply(1:nrow(pairs_df), function(i){
                tagList(
                    h5(paste0("Pos ES is up in ", pairs_df$cond2[i], 
                              ", neg ES is up in ", pairs_df$cond1[i], ".")),
                    plotOutput(sprintf("FGSEAcustomPlot%01d",i), height = '360.5px'),
                    verbatimTextOutput(sprintf("FGSEAcustomPval%01d",i))
                )
            })
        })
        
        observeEvent(input$saveCustomPlots, {
            if(length(names(symbList[[i]])[names(symbList[[i]]) %in% 
                                           pathways[[input$signat]]]) < 5){
                textMessage <- "< 5 genes found, try another signature"
            }else{
                textMessage <- paste0("The enrichment plots with ", 
                                      input$signat, " signature are saved in", 
                                      " 'data/app_plots/FGSEA_custom' folder.")
                output$savingCustomPlots <- renderText({ textMessage })
                
                dir.create(file.path(workdir, "app_plots"), showWarnings = F)
                dir.create(file.path(workdir, "app_plots/FGSEA_custom"), showWarnings = F)
                dir.create(file.path(workdir, paste0("app_plots/FGSEA_custom/",
                                                           input$signat)), showWarnings = F)
                for(i in 1:nrow(pairs_df)){
                    pdf(file.path(workdir, paste0("app_plots/FGSEA_custom/",
                                                        input$signat),
                                  paste0(filePrefix, input$signat, 
                                         "_pair", i, "_FGSEA_plot.pdf")),
                        width=5, height=4, onefile=FALSE)
                    print(plotEnrichment(pathways[[input$signat]], symbList[[i]]) + 
                              labs(title=input$signat))
                    dev.off()
                }
            }
        })
        
        # FGSEA tables custom
        lapply(1:nrow(pairs_df), function(i){
            if(nrow(fgsea_custom[[i]]) == 0){
                dataToShow <- data.frame(Data = "Zero significant signatures")
            }else{
                #dataToShow <- as.data.frame(
                #    fgsea_custom[[i]][fgsea_custom[[i]]$nMoreExtreme < 
                #                            0.5 * fgsea_custom[[i]]$size, ])
                dataToShow <- as.data.frame(fgsea_custom[[i]])
            }
            output[[sprintf("FGSEAcustomTable%01d",i)]] <- function() {
                dataToShow %>%
                    knitr::kable(format = "html", row.names = TRUE) %>%
                    kable_styling(styling_opt, full_width = F) %>%
                    scroll_box(width = "100%", height = DE_df_height)
            }
        })
        
        output$FGSEAcustomT <- renderUI({
            lapply(1:nrow(pairs_df), function(i){
                tagList(
                    h5(paste0("Pos ES is up in ", pairs_df$cond2[i], 
                              ", neg ES is up in ", pairs_df$cond1[i], ".")),
                    tableOutput(sprintf("FGSEAcustomTable%01d",i))
                )
            })
        })
        
        ##### Page 9 FGSEA all signatures #####
        # FGSEA all signatures plots
        lapply(1:nrow(pairs_df), function(i){
            output[[sprintf("FGSEAallPathPlot%01d",i)]] <- 
                renderPlot({
                        if(!(input$fgseaPath %in% names(examplePathways))){
                            plotOnlyText(textToPlot = "Select a signature")
                        }else if(length(names(symbList[[i]])[names(symbList[[i]]) %in% 
                                                             pathways[[input$signat]]]) < 5){
                            plotOnlyText(textToPlot = "< 5 genes found, try another signature")
                        }else{
                            plotEnrichment(examplePathways[[input$fgseaPath]], geneList[[i]]) + 
                                labs(title=input$fgseaPath)
                        }
                    })
        })
        
        lapply(1:nrow(pairs_df), function(i){
            output[[sprintf("FGSEAallPathPval%01d",i)]] <- 
                renderPrint({
                    if(any(fgsea_allPath[[i]]$pathway %in% input$fgseaPath)){
                        paste0("padj = ", 
                               signif(fgsea_allPath[[i]]$padj[fgsea_allPath[[i]]$pathway %in% 
                                                                    input$fgseaPath], 2),
                               ", NES = ", 
                               signif(fgsea_allPath[[i]]$NES[fgsea_allPath[[i]]$pathway %in% 
                                                                 input$fgseaPath], 2))
                        
                    }else{
                        paste0("padj >= ", GSEA_padj_thr, ", non-significant.")
                    }
                })
        })
        
        output$FGSEAallPathPlots <- renderUI({
            lapply(1:nrow(pairs_df), function(i){
                tagList(
                    h5(paste0("Pos ES is up in ", pairs_df$cond2[i], 
                              ", neg ES is up in ", pairs_df$cond1[i], ".")),
                    plotOutput(sprintf("FGSEAallPathPlot%01d",i), height = '360.5px'),
                    verbatimTextOutput(sprintf("FGSEAallPathPval%01d",i))
                )
            })
        })
        
        observeEvent(input$saveAllPathPlots, {
            if(length(names(symbList[[i]])[names(symbList[[i]]) %in% 
                                           pathways[[input$signat]]]) < 5){
                textMessage <- "< 5 genes found, try another signature"
            }else{
                textMessage <- paste0("The enrichment plots with ", 
                                      input$fgseaPath, " signature are saved in", 
                                      " 'data/app_plots/FGSEA_allPath' folder.")
                output$savingAllPathPlots <- renderText({ textMessage })
                
                dir.create(file.path(workdir, "app_plots"), showWarnings = F)
                dir.create(file.path(workdir, "app_plots/FGSEA_allPath"), showWarnings = F)
                dir.create(file.path(workdir, paste0("app_plots/FGSEA_allPath/",
                                                           substr(input$fgseaPath, start = 1, stop = 40))), 
                           showWarnings = F)
                for(i in 1:nrow(pairs_df)){
                    pdf(file.path(workdir, paste0("app_plots/FGSEA_allPath/",
                                                        substr(input$fgseaPath, start = 1, stop = 40)),
                                  paste0(filePrefix, substr(input$fgseaPath, 
                                                                     start = 1, stop = 40), 
                                         "_pair", i, "_FGSEA_plot.pdf")),
                        width=5, height=4, onefile=FALSE)
                    print(plotEnrichment(examplePathways[[input$fgseaPath]], geneList[[i]]) + 
                              labs(title=input$fgseaPath))
                    dev.off()
                }
            }
        })
        
        # FGSEA tables allPath
        lapply(1:nrow(pairs_df), function(i){
            if(nrow(fgsea_allPath[[i]]) == 0){
                dataToShow <- data.frame(Data = "Zero significant signatures")
            }else{
                #dataToShow <- as.data.frame(
                #    fgsea_allPath[[i]][fgsea_allPath[[i]]$nMoreExtreme < 
                #                           0.5 * fgsea_allPath[[i]]$size, ])
                dataToShow <- as.data.frame(fgsea_allPath[[i]])
            }
            output[[sprintf("FGSEAallPathTable%01d",i)]] <- function() {
                dataToShow %>%
                    knitr::kable(format = "html", row.names = TRUE) %>%
                    kable_styling(styling_opt, full_width = F) %>%
                    scroll_box(width = "100%", height = DE_df_height)
            }
        })
        
        output$FGSEAallPathT <- renderUI({
            lapply(1:nrow(pairs_df), function(i){
                tagList(
                    h5(paste0("Pos ES is up in ", pairs_df$cond2[i], 
                              ", neg ES is up in ", pairs_df$cond1[i], ".")),
                    tableOutput(sprintf("FGSEAallPathTable%01d",i))
                )
            })
        })
        
        ##### Page Families tables #####
        # FamFGSEA tables custom down
        famTablesServer("famFGSEAcustomDownTables", pairs_df = pairs_df, data = famFGSEAcustom_down, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamFGSEA tables custom up
        famTablesServer("famFGSEAcustomUpTables", pairs_df = pairs_df, data = famFGSEAcustom_up, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamFGSEA tables allPath down
        famTablesServer("famFGSEAallPathDownTables", pairs_df = pairs_df, data = famFGSEAallPath_down, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamFGSEA tables allPath up
        famTablesServer("famFGSEAallPathUpTables", pairs_df = pairs_df, data = famFGSEAallPath_up, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamBP tables down
        famTablesServer("famBPdownTables", pairs_df = pairs_df, data = famBP_down, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamBP tables Up
        famTablesServer("famBPupTables", pairs_df = pairs_df, data = famBP_up, 
                     styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamMF tables down
        famTablesServer("famMFdownTables", pairs_df = pairs_df, data = famMF_down, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamMF tables Up
        famTablesServer("famMFupTables", pairs_df = pairs_df, data = famMF_up, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamCC tables down
        famTablesServer("famCCdownTables", pairs_df = pairs_df, data = famCC_down, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamCC tables Up
        famTablesServer("famCCupTables", pairs_df = pairs_df, data = famCC_up, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamKEGG tables down
        famTablesServer("famKEGGdownTables", pairs_df = pairs_df, data = famKEGG_down, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        # FamKEGG tables Up
        famTablesServer("famKEGGupTables", pairs_df = pairs_df, data = famKEGG_up, 
                        styling_opt = styling_opt, DE_df_height = DE_df_height)
        
        ##### Families cnetplots #####
        # Fam FGSEA custom cnetplots Down
        famFGSEACnetplotsServer("famFGSEAcustomDownCnet", pairs_df = pairs_df, 
                                fgsea_obj = fgsea_custom, 
                                fgsea_table = famFGSEAcustom_down, 
                                enrich_example = bp_ego_s_up, 
                                dds = dds, maxNgenes = maxNgenes, 
                                showCategCnet = showCategCnet, 
                                geneList = geneList, upDown = "Down")
        
        # Fam FGSEA custom cnetplots Up
        famFGSEACnetplotsServer("famFGSEAcustomUpCnet", pairs_df = pairs_df, 
                                fgsea_obj = fgsea_custom, 
                                fgsea_table = famFGSEAcustom_up, 
                                enrich_example = bp_ego_s_up, 
                                dds = dds, maxNgenes = maxNgenes, 
                                showCategCnet = showCategCnet, 
                                geneList = geneList, upDown = "Up")
            
        # Fam FGSEA allPath cnetplots Down
        famFGSEACnetplotsServer("famFGSEAallPathDownCnet", pairs_df = pairs_df, 
                                fgsea_obj = fgsea_allPath, 
                                fgsea_table = famFGSEAallPath_down, 
                                enrich_example = bp_ego_s_up, 
                                dds = dds, maxNgenes = maxNgenes, 
                                showCategCnet = showCategCnet, 
                                geneList = geneList, upDown = "Down")
        
        # Fam FGSEA allPath cnetplots Up
        famFGSEACnetplotsServer("famFGSEAallPathUpCnet", pairs_df = pairs_df, 
                                fgsea_obj = fgsea_allPath, 
                                fgsea_table = famFGSEAallPath_up, 
                                enrich_example = bp_ego_s_up, 
                                dds = dds, maxNgenes = maxNgenes, 
                                showCategCnet = showCategCnet, 
                                geneList = geneList, upDown = "Up")
        
        # FamBP cnetplots Down
        famCnetplotsServer("famBPdownCnet", pairs_df = pairs_df, 
                           data = bp_ego_s_down, famData = famBP_down,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        # FamBP cnetplots Up
        famCnetplotsServer("famBPupCnet", pairs_df = pairs_df, 
                           data = bp_ego_s_up, famData = famBP_up,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        # FamMF cnetplots Down
        famCnetplotsServer("famMFdownCnet", pairs_df = pairs_df, 
                           data = mf_ego_s_down, famData = famMF_down,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        # FamMF cnetplots Up
        famCnetplotsServer("famMFupCnet", pairs_df = pairs_df, 
                           data = mf_ego_s_up, famData = famMF_up,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        # FamCC cnetplots Down
        famCnetplotsServer("famCCdownCnet", pairs_df = pairs_df, 
                           data = cc_ego_s_down, famData = famCC_down,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        # FamCC cnetplots Up
        famCnetplotsServer("famCCupCnet", pairs_df = pairs_df, 
                           data = cc_ego_s_up, famData = famCC_up,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        # FamKEGG cnetplots Down
        famCnetplotsServer("famKEGGdownCnet", pairs_df = pairs_df, 
                           data = kegg_down, famData = famKEGG_down,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        # FamKEGG cnetplots Up
        famCnetplotsServer("famKEGGupCnet", pairs_df = pairs_df, 
                           data = kegg_up, famData = famKEGG_up,
                           maxNgenes = maxNgenes, showCategCnet = showCategCnet, 
                           geneList = geneList)
        
        ##### Family barplots #####
        # Fam FGSEA custom barplots Down
        famFGSEAbarplotsServer("famFGSEAcustomDownBarplots", pairs_df = pairs_df, 
                                fgsea_obj = fgsea_custom, 
                                fgsea_table = famFGSEAcustom_down, 
                                enrich_example = bp_ego_s_up, 
                                dds = dds, maxNgenes = maxNgenes, 
                                showCategCnet = showCategCnet, 
                                geneList = geneList, upDown = "Down")
        
        # Fam FGSEA custom barplots Up
        famFGSEAbarplotsServer("famFGSEAcustomUpBarplots", pairs_df = pairs_df, 
                                fgsea_obj = fgsea_custom, 
                                fgsea_table = famFGSEAcustom_up, 
                                enrich_example = bp_ego_s_up, 
                                dds = dds, maxNgenes = maxNgenes, 
                                showCategCnet = showCategCnet, 
                                geneList = geneList, upDown = "Up")
        
        # Fam FGSEA allPath barplots Down
        famFGSEAbarplotsServer("famFGSEAallPathDownBarplots", pairs_df = pairs_df, 
                               fgsea_obj = fgsea_allPath, 
                               fgsea_table = famFGSEAallPath_down, 
                               enrich_example = bp_ego_s_up, 
                               dds = dds, maxNgenes = maxNgenes, 
                               showCategCnet = showCategCnet, 
                               geneList = geneList, upDown = "Down")
        
        # Fam FGSEA allPath barplots Up
        famFGSEAbarplotsServer("famFGSEAallPathUpBarplots", pairs_df = pairs_df, 
                               fgsea_obj = fgsea_allPath, 
                               fgsea_table = famFGSEAallPath_up, 
                               enrich_example = bp_ego_s_up, 
                               dds = dds, maxNgenes = maxNgenes, 
                               showCategCnet = showCategCnet, 
                               geneList = geneList, upDown = "Up")
        
        # Fam BP barplots down
        famBarplotsServer("famBPdownBarplots", pairs_df = pairs_df, 
                          data = bp_ego_s_down, famData = famBP_down, 
                          showCateg = showCateg)
        
        # Fam BP barplots up
        famBarplotsServer("famBPupBarplots", pairs_df = pairs_df, 
                          data = bp_ego_s_up, famData = famBP_up, 
                          showCateg = showCateg)
        
        # Fam MF barplots down
        famBarplotsServer("famMFdownBarplots", pairs_df = pairs_df, 
                          data = mf_ego_s_down, famData = famMF_down, 
                          showCateg = showCateg)
        
        # Fam MF barplots up
        famBarplotsServer("famMFupBarplots", pairs_df = pairs_df, 
                          data = mf_ego_s_up, famData = famMF_up, 
                          showCateg = showCateg)
            
        # Fam CC barplots down
        famBarplotsServer("famCCdownBarplots", pairs_df = pairs_df, 
                          data = cc_ego_s_down, famData = famCC_down, 
                          showCateg = showCateg)
        
        # Fam CC barplots up
        famBarplotsServer("famCCupBarplots", pairs_df = pairs_df, 
                          data = cc_ego_s_up, famData = famCC_up, 
                          showCateg = showCateg)
        
        # Fam KEGG barplots down
        famBarplotsServer("famKEGGdownBarplots", pairs_df = pairs_df, 
                          data = kegg_down, famData = famKEGG_down, 
                          showCateg = showCateg)
        
        # Fam KEGG barplots up
        famBarplotsServer("famKEGGupBarplots", pairs_df = pairs_df, 
                          data = kegg_up, famData = famKEGG_up, 
                          showCateg = showCateg)
    })
    
    # Create Shiny app ----
    shinyApp(ui = ui, server = server)
}
