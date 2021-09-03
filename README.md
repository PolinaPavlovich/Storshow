## Storshow

Storshow takes a count matrix of bulk RNA-seq and returns a Shiny app

1) Clone Storshow directory from GitHub to storshow_dir
2) Run conda *env create -f CondaEnv.yaml* to create conda environment
3) Run *conda activate storshow*
4) Put input files in the *input* folder (see example in https://github.com/PolinaPavlovich/Storshow_input_example)
5) Run
*python storshow_dir/Storshow input/counts.tsv input/coldata.tsv input/ctr_treat.tsv style.tsv outputdir*

## Command line options

[-h] [--organism organism] [--fPrefix prefix] [--rm_samples rm_samples] [--geneFilt geneFilt]

[--DESeq2_padj DESeq2_padj] [--DESeq2_absLog2FC DESeq2_absLog2FC] [--ORA_padj ORA_padj]

[--GSEA_padj GSEA_padj] [--avg_pair_overl avg_pair_overl] [--egoS_maxScores egoS_maxScores]

[--style style]

cmat coldata compar outdir

## Storshow workflow

### Positional arguments:

**cmat**    Provide path to an input count matrix.

**coldata** Provide path to the tsv with coldata.

**compar**  Provide a tsv with pairs to compare. First column is Ctr, second is Treated.

**outdir**  Provide an output directory.

### Optional arguments:

**-h, --help**                Show this help message and exit

**--organism organism**       Provide an organism. Either human, mouse or fly. Defaults to mouse.

**--fPrefix fPrefix**         String to prepend in output file names. It is usually an experiment name.

**--rm_samples rm_samples**   A string with samples to exclude, e.g. 'sampleA,sampleB,sampleC'.
 
**--geneFilt geneFilt**       A criterion to filter genes in a format 
                       'criterion_threshold,' where the 'criterion' can be
                        either 'avgCond' or 'total,' and 'threshold' is a
                        positive integer. For example, for keeping genes with
                        more than 100 reads per condition on average, specify
                        'avgCond_100'. Or for keeping genes with more than 10
                        reads per gene specify 'total_10'.

**--DESeq2_padj DESeq2_padj** The threshold for padj to select significant differentially expressed genes.

**--DESeq2_absLog2FC DESeq2_absLog2FC**
                        The threshold for absolute value of log2FC to select
                        significant differentially expressed genes.
                        
**--ORA_padj ORA_padj**   P-adjusted threshold for GO term and KEGG pathways
                        enrichment.
                        
**--GSEA_padj GSEA_padj** P-adjusted threshold for FGSEA R package.

**--avg_pair_overl avg_pair_overl**
                        A percentage of genes that overlap in two significant
                        gene sets to assign them to one family.
                        
**--egoS_maxCores egoS_maxCores**
                        The maximum number of cores to calculate GO term
                        enrichment. Usually, if it is equal to the number of
                        comparisons, then the calculations of this step would
                        take 3-5 minutes.
                        
**--style style**         A path to a file containing data frame with one or two
                        columns. The first column, called 'levels', specifies
                        conditions from the most untreated to the most
                        treated. It will allow paired bar plots to have a
                        control condition always on the left and treated on
                        the right side. The second column, called 'colors', is
                        the corresponding colors. If absent, the order of
                        conditions and the colors will be assigned
                        automatically.
 
## Example data
Example data can be found here https://github.com/PolinaPavlovich/Storshow_input_example

### Example cmat counts.tsv (sep = "\t")

	Ctrl_1	Ctrl_2	Ctrl_3	Treat_1	Treat_2	Treat3	KO_Treat1 KO_Treat2 KO_Treat3
 
ENSMUSG00000000001.4	5709	5103	5516	6370	6342	4926	5319	4880	5177

ENSMUSG00000000003.15	0	3	1	1	2	0	0	0	0

ENSMUSG00000000028.14	1743	2087	1852	2842	2627	2142	1905	2631	2393

ENSMUSG00000000031.15	8	14	16	8	2	0	15	21	3

ENSMUSG00000000037.16	245	298	197	257	290	150	165	163	122

### Example coldata coldata.tsv (sep = "\t")

sampleName	state

Ctrl_1	Ctrl

Ctrl_2	Ctrl

Ctrl_3	Ctrl	

Treat_1	Treat

Treat_2	Treat

Treat_3	Treat

KO_Treat_1	KO_Treat

KO_Treat_2	KO_Treat

KO_Treat_3	KO_Treat

### Example style style.tsv (sep = "\t")

levels	colors

Ctrl	#FB9A99

Treat	#A6CEE3

KO_Treat	#B2DF8A

### Example compar ctr_treat.tsv (sep = "\t")

ctr_cond	treat_cond

Ctrl	Treat

Ctrl KO_Treat

Treat KO_Treat

                        



