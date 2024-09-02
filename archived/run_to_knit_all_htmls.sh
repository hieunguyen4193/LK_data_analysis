maindir=$(pwd);
for dataset in Dataset1 Dataset2 Dataset3 Dataset4;do Rscript ${maindir}/${dataset}/knit_html.R;done
