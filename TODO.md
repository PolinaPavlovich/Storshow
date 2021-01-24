### Todo

   - [x] in / out dir params
   - [x] sinks for logging.
   - [x] create DB on the fly
   - [x] yaml w/ thresholds etc.
   - [x] Have conda environments built / checked (?) before executing workflow.
   - [x] parse count matrices
   - [ ] allow for multiple matrices as input ?
   - [x] rewrite for loop in preproc
   - [x] padj out
   - [x] DAVID gene symbols
   - [x] enrichGO
   - [x] KEGG results
   - [x] FGSEA

### Polina's Todo

   - [x] Create an environment with all R packages using conda installations
   - [ ] Add download buttons to the app
   - [ ] Test with a human count_mat
   - [ ] Test with a fly count_mat
   - [ ] (Frozen) Include render Rmd into the Rscripts
   - [x] Add compare2cond into the DESeq2 Rscript
   - [x] Rewrite PreProc to take coldata file as input
   - [x] Rewrite PreProc to take treat-ctr comparisons file as input
   - [ ] Add rules from StoreData
   - [ ] Add ShinyInternal
   - [x] Add functions.R that are needed for StoreData and ShinyIternal
   - [ ] Make a script that produces and saves app.R script to launch the app (conda .Rlib from Rstudio?)
   - [ ] Rewrite ShinyInternal with Shiny modules
   - [ ] Make GO/KEGG barplots height as a function of number of pathways
   - [ ] Add baseCond into DESeq2 result tables in the app
   - [ ] Move plot params from config to Shiny
   - [ ] Add lollipop plots into the app
   - [ ] Add PCA plots into the app
   - [ ] Write a func for DAVID analysis in R with gene Symbols
   - [ ] Intergrate RDAVID func into the app
   - [ ] Write documentation
   - [x] Add gene filtering options e.g. "avg_100" and "total_10"
   - [ ] Add DESeq2 analysis for betaPrior = F
   - [ ] Allow overwriting into existing output directory
   - [ ] Allow a user to provide a gmt file with GSEA signatures
