#!/bin/sh


# created by wangjiaxuan 2022.7

# 设置工作路径
mkdir 3.Result_Sum
mkdir 2.Humann2_Quantity
mkdir 1.Kneaddata_Clean

func() {
    echo "Usage:"
    echo "[-i]: The input tsv without header including four column which mean group, sample ,read1 and read2 fq file path"
    echo "[-s]: which type want select, 1:st_superman,2:st"
    echo "[-h]: The help document"
    # shellcheck disable=SC2242
    # shellcheck disable=SC2242
    exit -1
}

input="/data/wangjiaxuan/workflow/meta-genome/0.Input/sample_input_path.tsv"

while getopts "s:g:c:l:r:h" opt; do
    case $opt in
      i) input="$OPTARG";;
      s) step="$OPTARG";;
      h) func;;
      ?) func;;
    esac
done

export PATH="/home/tangwenli/miniconda3/bin:$PATH"
which conda
source activate humann2
export PATH="/home/tangwenli/miniconda3/envs/humann2/bin:$PATH"

# 环境变量
fastuniq=/data/wangjiaxuan/biosoft/miniconda3/envs/meta/bin/fastuniq
kneaddata_db=/tools/db/metagenome/KneadData_db/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1
run_mem=506384m
trimmomatic_p="SLIDINGWINDOW:5:20 MINLEN:36 LEADING:3 TRAILING:3 ILLUMINACLIP:/tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10"
# shellcheck disable=SC1073
cat ${input} | while read group sample fq1 fq2
do
  ## PCR Deduplication-----------
  echo -e "${fq1}\n${fq2}" > fq.list
  ${fastuniq} -i fq.list -o 0.Input/${sample}.uniq.R1.fq -p 0.Input/${sample}.uniq.R2.fq
  ## Kneaddata clean fq and remove host contamination
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
done
rm fq.list

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