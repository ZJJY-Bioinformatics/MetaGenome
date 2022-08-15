#!/bin/sh


# created by wangjiaxuan 2022.7

# 设置工作路径
if [ ! -e 0.Input/ ]; then mkdir 0.Input  ;fi
if [ ! -e 3.Result_Sum ]; then mkdir 3.Result_Sum  ;fi
if [ ! -e 2.Humann2_Quantity ]; then mkdir 2.Humann2_Quantity  ;fi
if [ ! -e 1.Kneaddata_Clean ]; then mkdir 1.Kneaddata_Clean  ;fi
if [ ! -e shell/ ]; then mkdir shell  ;fi

func() {
    echo "Usage:"
    echo "[-i]: The input tsv without header including four column which mean group, sample ,read1 and read2 fq file path"
    echo "[-s]: which type want select, 1:st_superman,2:st"
    echo "[-h]: The help document"
    # shellcheck disable=SC2242
    # shellcheck disable=SC2242
    exit -1
}

input="0.Input/sample_input_path.tsv"

while getopts "i:s:h" opt; do
    case $opt in
      i) input="$OPTARG";;
      s) step="$OPTARG";;
      h) func;;
      ?) func;;
    esac
done


# 环境变量
fastuniq=/data/wangjiaxuan/biosoft/miniconda3/envs/meta/bin/fastuniq
kneaddata_db=/tools/db/metagenome/KneadData_db/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1
run_mem=506384m
trimmomatic_p="SLIDINGWINDOW:5:20 MINLEN:36 LEADING:3 TRAILING:3 ILLUMINACLIP:/tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10"
# shellcheck disable=SC1073


cat ${input} | while read group sample fq1 fq2
do
  cat > shell/${sample}_main_run.sh<<EOF

  # 激活conda环境
  export PATH="/home/tangwenli/miniconda3/bin:\$PATH"
  source activate humann2
  export PATH="/home/tangwenli/miniconda3/envs/humann2/bin:\$PATH"

  ## PCR Deduplication-----------
  echo -e "${fq1}\n${fq2}" > 0.Input/fq.list
  ${fastuniq} -i 0.Input/fq.list -o 0.Input/${sample}.uniq.R1.fq -p 0.Input/${sample}.uniq.R2.fq
  # Kneaddata clean fq and remove host contamination
  kneaddata \
  -i 0.Input/${sample}.uniq.R1.fq \
  -i 0.Input/${sample}.uniq.R2.fq \
  -o 1.Kneaddata_Clean \
  -db ${kneaddata_db} \
  -t 10 \
  --max-memory ${run_mem} \
  --trimmomatic-options "${trimmomatic_p}" \
  --bowtie2-options "--very-sensitive --dovetail" \
  --remove-intermediate-output \
  --trf /data/wangjiaxuan/biosoft/miniconda3/envs/meta/bin/ \
  --trimmomatic /data/wangjiaxuan/biosoft/Trimmomatic \
  --run-fastqc-start \
  --run-fastqc-end \
  --fastqc /data/wangjiaxuan/biosoft/miniconda3/envs/meta/bin/

  cat 1.Kneaddata_Clean/${sample}.uniq*_kneaddata_[pu]*.fastq > 1.Kneaddata_Clean/${sample}_humann2_in4.fq

  /home/tangwenli/miniconda3/envs/humann2/bin/humann2 \
  --threads 25 \
  --input  1.Kneaddata_Clean/${sample}_humann2_in4.fq \
  --output 2.Humann2_Quantity \
  --search-mode uniref90

  humann2_renorm_table -i 2.Humann2_Quantity/${sample}_humann2_in4_genefamilies.tsv -o 2.Humann2_Quantity/${sample}_humann2_in4_genefamilies_cpm.tsv --units cpm
  echo "${sample} have analysis finised!"
EOF
done

if [ -e 0.Input/fq.list ]; then rm 0.Input/fq.list  ;fi

# 批量运行
datetime=$(date)
echo "the workflow run at ${datetime}"
pid=()
for script in shell/*_main_run.sh
  do
    bash ${script} 2>shell/${sample}.err 1>shell/${sample}.log &
    pid+=("$!")
  done
wait ${pid[@]}
datetime=$(date)
echo "the workflow finished at ${datetime}"

# 收集bug list
/data/wangjiaxuan/biosoft/miniconda3/envs/rnaseq/bin/collect-columns \
3.Result_Sum/all.sample_buglist.tsv \
2.Humann2_Quantity/*_temp/*_bugs_list.tsv


# rename tabel

/data/wangjiaxuan/biosoft/miniconda3/envs/meta/bin/kneaddata_read_count_table \
--input 1.Kneaddata_Clean \
--output 1.Kneaddata_Clean/kneaddata_qc_result.tsv

#
/data/wangjiaxuan/biosoft/miniconda3/envs/meta/bin/multiqc \
-d 1.Kneaddata_Clean/fastqc \
-o 1.Kneaddata_Clean/multiqc_result

#
humann2_join_tables \
-i 2.Humann2_Quantity \
-o 3.Result_Sum/all.sample_genefamilies.tsv \
--file_name genefamilies

humann2_join_tables \
-i 2.Humann2_Quantity \
-o 3.Result_Sum/all.sample_pathabundance.tsv \
--file_name pathabundance

humann2_join_tables \
-i 2.Humann2_Quantity \
-o 3.Result_Sum/all.sample_pathcoverage.tsv \
--file_name pathcoverage

humann2_join_tables \
-i 2.Humann2_Quantity \
-o 3.Result_Sum/all.sample_genefamilies_cpm.tsv \
--file_name genefamilies_cpm

humann2_rename_table \
--input 3.Result_Sum/all.sample_genefamilies.tsv \
--output 3.Result_Sum/all.sample_Functionfamilies.tsv \
--names uniref90

humann2_rename_table \
--input 3.Result_Sum/all.sample_genefamilies_cpm.tsv \
--output 3.Result_Sum/all.sample_Functionfamilie_cpms.tsv \
--names uniref90

humann2_rename_table \
--input 3.Result_Sum/all.sample_genefamilies_cpm.tsv \
--output 3.Result_Sum/all.sample_KO_cpms.tsv \
--names kegg-pathway

humann2_rename_table \
--input 3.Result_Sum/all.sample_genefamilies_cpm.tsv \
--output 3.Result_Sum/all.sample_GO_cpms.tsv \
--names go

humann2_rename_table \
--input 3.Result_Sum/all.sample_genefamilies_cpm.tsv \
--output 3.Result_Sum/all.sample_metacyc_cpms.tsv \
--names metacyc-pwy

mkdir -p 4.Out2CAMP

/data/wangjiaxuan/biosoft/miniconda3/bin/Rscript \
util/out2cmap.r


mkdir -p 5.HostGene_EXP


cat ${input} | while read group sample fq1 fq2
do
  fq1=$(ls 1.Kneaddata_Clean/${sample}*_contam_1.fastq)
  fq2=$(ls 1.Kneaddata_Clean/${sample}*_contam_2.fastq)
  echo -e "${group}\t${sample}\t${fq1}\t${fq2}" >> 5.HostGene_EXP/hostrna_input.tsv
done

cd  5.HostGene_EXP
/data/wangjiaxuan/workflow/bulk-rna-seq/run_RNAseq -i hostrna_input.tsv


