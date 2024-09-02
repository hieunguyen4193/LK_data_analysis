# Trong-Hieu Nguyen, last modified on 24.08.2022.
current_dir=$(pwd);

path_to_cellranger="/home/hieunguyen/CRC1382/src/src_pipeline/CellRanger/cellranger-6.1.2";
export PATH=${path_to_cellranger}:$PATH;

for sample_id in m331 m330 m332;do \
path_to_save_output="/media/hieunguyen/HNSD_MBPro/raw_data/230215_Kopplin_Pabst_added_NC_000001_merged_zcat";
mkdir -p ${path_to_save_output};

# M A I N - C E L L R A N G E R - C O M M A N D S
mkdir -p ${path_to_save_output}/${sample_id}_added_NC_000001_merged_zcat;

cellranger multi --id=${sample_id}_added_NC_000001_merged_zcat --csv=multi_config_Sample_${sample_id}.csv --localcores=25;

rsync -avh --progress ${current_dir}/${sample_id}_added_NC_000001_merged_zcat ${path_to_save_output}/${sample_id}_added_NC_000001_merged_zcat;

rm -rf ${current_dir}/${sample_id}_added_NC_000001_merged_zcat;
done
# EOF

