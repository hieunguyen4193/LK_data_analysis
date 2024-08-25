maindir="/home/hieunguyen/CRC1382/src_2023/LKopplin/20231101_OFFICIAL";
logfile=${maindir}/logs_run_pipeline.txt;

echo -e "Start running pipeline for dataset 1" >> $logfile;
Rscript ${maindir}/Dataset1/run_pipeline.R;

echo -e "Start running pipeline for dataset 2" >> $logfile;
Rscript ${maindir}/Dataset2/run_pipeline.R;

echo -e "Start running pipeline for dataset 3" >> $logfile;
Rscript ${maindir}/Dataset3/run_pipeline.1st_round.R;
Rscript ${maindir}/Dataset3/run_pipeline.2nd_round.R;

echo -e "Start running pipeline for dataset 4" >> $logfile;
Rscript ${maindir}/Dataset4/run_pipeline.R;

echo -e "finish running all four pipelines" >> $logfile

