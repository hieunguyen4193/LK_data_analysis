# Trong-Hieu Nguyen, last modified on 24.08.2022.
current_dir=$(pwd);

path_to_cellranger="/home/hieunguyen/CRC1382/src/src_pipeline/CellRanger/cellranger-6.1.2";
export PATH=${path_to_cellranger}:$PATH;

for sample_id in m366 m367 m368 m369;do \
path_to_save_output="/media/hieunguyen/CRC1382_HDD2/CRC1382/raw_data/230316_Kopplin";
mkdir -p ${path_to_save_output};

# M A I N - C E L L R A N G E R - C O M M A N D S
mkdir -p ${path_to_save_output}/${sample_id}_added_NC_000001;

cellranger multi --id=${sample_id}_added_NC_000001 --csv=multi_config_Sample_${sample_id}.csv --localcores=25;

rsync -avh --progress ${current_dir}/${sample_id}_added_NC_000001 ${path_to_save_output}/${sample_id}_added_NC_000001;

rm -rf ${current_dir}/${sample_id}_added_NC_000001;
done
# EOF

