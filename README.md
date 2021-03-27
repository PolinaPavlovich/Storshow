### Storshow

Storshow takes a count matrix of bulk RNA-seq and returns a Shiny app

1) Clone Storshow directory from GitHub to storshow_dir
2) Run conda *env create -f CondaEnv.yaml* to create conda environment
3) Run *conda activate storshow*
4) Run
*python storshow_dir/Storshow -c input/counts.tsv --coldata input/coldata.tsv --compar input/ctr_treat.tsv --style input/style.tsv -o outputdir*

### Command line options

[-h] [--organism organism] [--fPrefix prefix] [--rm_samples rm_samples] [--geneFilt geneFilt]
[--DESeq2_padj DESeq2_padj] [--DESeq2_absLog2FC DESeq2_absLog2FC] [--ORA_padj ORA_padj]
[--GSEA_padj GSEA_padj] [--avg_pair_overl avg_pair_overl] [--egoS_maxScores egoS_maxScores]
[--style style]
cmat coldata compar outdir

Storshow workflow

positional arguments:
cmat    Provide path to an input count matrix.
coldata Provide path to the tsv with coldata.
compar  Provide a tsv with pairs to compare. First column is Ctr, second is Treated.
outdir  Provide an output directory.

optional arguments:
-h, --help            Show this help message and exit
--organism organism   Provide an organism. Either human, mouse or fly. Defaults to mouse.




