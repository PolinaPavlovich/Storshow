#!/usr/bin/env Rscript
.libPaths(R.home("library"))

# Read parameters/set up sink
workdir <- snakemake@params[['workdir']]
storshow_dir <- snakemake@params[['storshow_dir']]
#extra_data_dds <- snakemake@params[['extra_dds_dir']]
log_fn <- snakemake@log[['out']]
style_file <- snakemake@params[['style']]
compar_file <- snakemake@params[['compar']]

#workdir <- "/data/cabezas/group/pavlovich/data/shared/storshow/runs/outputdir139"
#storshow_dir <- "/data/cabezas/group/pavlovich/data/shared/storshow/R/dev_Polina/Storshow"
#style_file <- "/data/cabezas/group/pavlovich/data/shared/storshow/runs/input/CypKO_style.tsv"
#log_fn <- "/data/cabezas/group/pavlovich/data/shared/storshow/runs/outputdir139/Log/saveAppRscript.log"

logfile <- file(log_fn, open="w+")
sink(logfile, type = c("output", "message"))

# Copy necessary files to display the Shiny app
dir.create(file.path(workdir, "showShiny"), showWarnings = F)
system(paste0("cp ", file.path(storshow_dir, "scripts/readYaml.R "), 
              file.path(workdir, "showShiny/readYaml.R")))

system(paste0("cp ", file.path(storshow_dir, "scripts/storshow_functions.R "), 
              file.path(workdir, "showShiny/storshow_functions.R")))

system(paste0("cp ", file.path(storshow_dir, "scripts/signatures/Signatures_for_GSEA_updated_20200714.gmt "), 
              file.path(workdir, "showShiny/Signatures_for_GSEA_updated_20200714.gmt")))

system(paste0("cp ", file.path(storshow_dir, "scripts/bulk_RNA_seq_ShinyInternal.R "), 
              file.path(workdir, "showShiny/bulk_RNA_seq_ShinyInternal.R")))

system(paste0("cp ", style_file, " ",
              file.path(workdir, "showShiny/style.tsv")))

system(paste0("cp ", compar_file, " ",
              file.path(workdir, "showShiny/ctr_treat.tsv")))

# Write the app.R script
app_file <- file.path(workdir, "app.R")
file_con <- file(app_file, open="w+")

writeLines("#!/usr/bin/env Rscript", con = file_con)
writeLines('require(rstudioapi)', con = file_con)

writeLines('', con = file_con)

writeLines('workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)', con = file_con)
writeLines('gmt.file <- file.path(workdir, "showShiny", "Signatures_for_GSEA_updated_20200714.gmt")', 
           con = file_con)
writeLines('style_file <- file.path(workdir, "showShiny", "style.tsv")', 
           con = file_con)
writeLines('compar_file <- file.path(workdir, "showShiny", "ctr_treat.tsv")', 
           con = file_con)
writeLines('storshow_functions <- file.path(workdir, "showShiny/storshow_functions.R")', 
           con = file_con)

writeLines('', con = file_con)

writeLines(paste0('source(file.path(workdir, "showShiny/readYaml.R"))'), con = file_con)
writeLines(paste0('source(file.path(workdir, "showShiny/bulk_RNA_seq_ShinyInternal.R"))'), 
           con = file_con)

writeLines('', con = file_con)

writeLines("# start Shiny application", con = file_con)
writeLines("showShiny()", con = file_con)

close(file_con)

print("Warnings: ", warnings())
print(sessionInfo())

# Close sink
close(logfile)