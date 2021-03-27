### Storshow

Storshow takes a count matrix of bulk RNA-seq and returns a Shiny app

1) Clone Storshow directory from GitHub to storshow_dir
2) Run conda *env create -f CondaEnv.yaml* to create conda environment
3) Run *conda activate storshow*
4) Run
*python storshow_dir/Storshow -c input/counts.tsv --coldata input/coldata.tsv --compar input/ctr_treat.tsv --style input/style.tsv -o outputdir*

### Command line options

[-h] [--organism organism] [--fPrefix prefix] [

