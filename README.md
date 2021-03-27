### Storshow

Storshow takes a count matrix of bulk RNA-seq and returns a Shiny app

1) Clone Storshow directory from GitHub to storshow_dir
2) Run conda *env create -f CondaEnv.yaml* to create conda environment
3) Run *conda activate storshow*
4) Run
*python storshow_dir/Storshow input/counts.tsv input/coldata.tsv input/ctr_treat.tsv --style style.tsv outputdir*

### Command line options

[-h] [--organism organism] [--fPrefix prefix] [--rm_samples rm_samples] [--geneFilt geneFilt]

[--DESeq2_padj DESeq2_padj] [--DESeq2_absLog2FC DESeq2_absLog2FC] [--ORA_padj ORA_padj]

[--GSEA_padj GSEA_padj] [--avg_pair_overl avg_pair_overl] [--egoS_maxScores egoS_maxScores]

[--style style]

cmat coldata compar outdir

Storshow workflow

#### positional arguments:

cmat    Provide path to an input count matrix.
coldata Provide path to the tsv with coldata.
compar  Provide a tsv with pairs to compare. First column is Ctr, second is Treated.
outdir  Provide an output directory.

#### optional arguments:

*-h, --help*                Show this help message and exit

*--organism organism*       Provide an organism. Either human, mouse or fly. Defaults to mouse.

*--fPrefix fPrefix*         String to prepend in output file names. It is usually an experiment name.

*--rm_samples rm_samples*   A string with samples to exclude, e.g. 'sampleA,sampleB,sampleC'.
 
*--geneFilt geneFilt*       A criterion to filter genes in a format 
                       'criterion_threshold,' where the 'criterion' can be
                        either 'avgCond' or 'total,' and 'threshold' is a
                        positive integer. For example, for keeping genes with
                        more than 100 reads per condition on average, specify
                        'avgCond_100'. Or for keeping genes with more than 10
                        reads per gene specify 'total_10'.

*--DESeq2_padj DESeq2_padj* The threshold for padj to select significant differentially expressed genes.
*--DESeq2_absLog2FC DESeq2_absLog2FC*
                        The threshold for absolute value of log2FC to select
                        significant differentially expressed genes.



