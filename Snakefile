localrules: preProcess,deSeq2,geneCl,saveEgoS,saveKegg,saveFgseaCustom,saveFgseaAllPath,families,geneClTrend,saveAppRscript
rule all:
    input:
        config["outdir"] + "/preProc/fn_preProc.tsv",
        config["outdir"] + "/DESeq/fn_deSeq2.tsv",
        config["outdir"] + "/DESeq/fn_geneCl.tsv",
        config["outdir"] + "/enrichGO_simpl/saveEgoS.tsv",
        config["outdir"] + "/KEGG/saveKegg.tsv",
        config["outdir"] + "/FGSEA/saveFgseaCustom.tsv",
        config["outdir"] + "/FGSEA/saveFgseaAllPath.tsv",
        config["outdir"] + "/FGSEA/families.tsv",
        config["outdir"] + "/DESeq/fn_geneClTrend.tsv",
        config["outdir"] + "/app.R"
    
rule preProcess:
    input:
        config["cmat"],
        config["coldata"]
    output: 
        config["outdir"] + "/preProc/fn_preProc.tsv",
    params:
        inc = config["cmat"],
        coldata = config["coldata"],
        outdir = config["outdir"] + "/preProc",
        organism = config["organism"],
        rm_samples = config["rm_samples"],
        geneFilt = config["geneFilt"],
        workdir = config["outdir"],
        filePrefix = config["fPrefix"]
    log:
        out = config["outdir"] + "/Log/preProcess.log"
    script: "scripts/preProc.R"

rule deSeq2:
    input:
        config["outdir"] + "/preProc/fn_preProc.tsv",
    output:
        config["outdir"] + "/DESeq/fn_deSeq2.tsv",
    params:
        indir = config["outdir"] + "/preProc/",
        outdir = config["outdir"] + "/DESeq",
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        filePrefix = config["fPrefix"],
        DESeq2_padj = config["DESeq2_padj"],
        DESeq2_absLog2FC = config["DESeq2_absLog2FC"],
        workdir = config["outdir"]
    log:
        out = config["outdir"] + "/Log/deSeq2.log"
    script: "scripts/DESeq2.R"
    
rule geneCl:
    input:
        config["outdir"] + "/DESeq/fn_deSeq2.tsv",
    output:
        config["outdir"] + "/DESeq/fn_geneCl.tsv"
    params:
        deseq2_dir = config["outdir"] + "/DESeq",
        outdir = config["outdir"],
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        filePrefix = config["fPrefix"],
        workdir = config["outdir"]
    log:
        out = config["outdir"] + "/Log/geneCl.log"
    script: "scripts/geneClusters.R"
    
rule saveEgoS:
    input:
        config["outdir"] + "/DESeq/fn_deSeq2.tsv"
    output:
        config["outdir"] + "/enrichGO_simpl/saveEgoS.tsv"
    params:
        deseq2_dir = config["outdir"] + "/DESeq",
        outdir = config["outdir"] + "/enrichGO_simpl",
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        workdir = config["outdir"],
        organism = config["organism"],
        egoS_maxCores = config["egoS_maxCores"],
        ORA_padj = config["ORA_padj"]
    log:
        out = config["outdir"] + "/Log/saveEgoS.log"
    script: "scripts/saveEgoS.R"
    
rule saveKegg:
    input:
        config["outdir"] + "/DESeq/fn_deSeq2.tsv"
    output:
        config["outdir"] + "/KEGG/saveKegg.tsv"
    params:
        deseq2_dir = config["outdir"] + "/DESeq",
        outdir = config["outdir"] + "/KEGG",
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        workdir = config["outdir"],
        organism = config["organism"],
        ORA_padj = config["ORA_padj"]
    log:
        out = config["outdir"] + "/Log/saveKegg.log"
    script: "scripts/saveKegg.R"
    
rule saveFgseaCustom:
    input:
        config["outdir"] + "/DESeq/fn_deSeq2.tsv"
    output:
        config["outdir"] + "/FGSEA/saveFgseaCustom.tsv"
    params:
        deseq2_dir = config["outdir"] + "/DESeq",
        outdir = config["outdir"] + "/FGSEA",
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        workdir = config["outdir"],
        GSEA_padj = config["GSEA_padj"]
    log:
        out = config["outdir"] + "/Log/saveFgseaCustom.log"
    script: "scripts/saveFgseaCustom.R"
    
rule saveFgseaAllPath:
    input:
        config["outdir"] + "/DESeq/fn_deSeq2.tsv"
    output:
        config["outdir"] + "/FGSEA/saveFgseaAllPath.tsv"
    params:
        deseq2_dir = config["outdir"] + "/DESeq",
        outdir = config["outdir"] + "/FGSEA",
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        workdir = config["outdir"],
        GSEA_padj = config["GSEA_padj"]
    log:
        out = config["outdir"] + "/Log/saveFgseaAllPath.log"
    script: "scripts/saveFgseaAllPath.R"
    
rule families:
    input:
        config["outdir"] + "/enrichGO_simpl/saveEgoS.tsv",
        config["outdir"] + "/KEGG/saveKegg.tsv",
        config["outdir"] + "/FGSEA/saveFgseaCustom.tsv",
        config["outdir"] + "/FGSEA/saveFgseaAllPath.tsv"
    output:
        config["outdir"] + "/FGSEA/families.tsv"
    params:
        deseq2_dir = config["outdir"] + "/DESeq",
        outdir = config["outdir"] + "/FGSEA",
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        workdir = config["outdir"],
        avg_pair_overl = config["avg_pair_overl"],
        ORA_padj = config["ORA_padj"]
    log:
        out = config["outdir"] + "/Log/families.log"
    script: "scripts/families.R"
    
rule geneClTrend:
    input:
        config["outdir"] + "/DESeq/fn_geneCl.tsv"
    output:
        config["outdir"] + "/DESeq/fn_geneClTrend.tsv"
    params:
        deseq2_dir = config["outdir"] + "/DESeq",
        outdir = config["outdir"],
        compar = config["compar"],
        storshow_dir = workflow.basedir,
        workdir = config["outdir"],
        filePrefix = config["fPrefix"]
    log:
        out = config["outdir"] + "/Log/geneClTrend.log"
    script: "scripts/geneClTrend.R"
    
rule saveAppRscript:
    input:
        config["cmat"]
    output:
        config["outdir"] + "/app.R"
    params:
        storshow_dir = workflow.basedir,
        workdir = config["outdir"],
        style = config["style"],
        compar = config["compar"]
    log:
        out = config["outdir"] + "/Log/saveAppRscript.log"
    script: "scripts/saveAppRscript.R"
