### Storshow

What does it do.

Example data.

Example to run.

Current commands to test the pipeline (run from a directory containing 'input' folder):

cd /data/cabezas/group/pavlovich/data/shared/storshow/runs

python /data/cabezas/group/pavlovich/data/shared/storshow/R/dev_Polina/Storshow/Storshow -c input/CypKO_counts.tsv --coldata input/CypKO_coldata.tsv --compar input/CypKO_ctr_treat.tsv --rm_samples Ctrl_4oxoRA_142,Ctrl_DMSO_142,Ctrl_atRA_142,Ctrl_retinol_142 --style input/CypKO_style.tsv -o outputdir
