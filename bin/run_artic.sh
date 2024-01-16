#!/bin/bash
set -euo pipefail

sample_name=$1
fastq_file=$2
min_len=$3
max_len=$4
medaka_model=$5
scheme_name=$6
scheme_dir=$7
scheme_version=$8
threads=$9
max_softclip_length=${10}
normalise=${11}

if [[ "$scheme_version" == "None" ]];
then
    scheme_version="."
fi

function mock_artic {
    # write an empty VCF
    echo "Mocking artic results"
    TAB="$(echo -e '\t')"
    cat << EOF |
##fileformat=VCFv4.2
##source=Longshot v0.4.0
#CHROM${TAB}POS${TAB}ID${TAB}REF${TAB}ALT${TAB}QUAL${TAB}FILTER${TAB}INFO${TAB}FORMAT${TAB}SAMPLE
EOF
    bgzip > "${sample_name}.pass.vcf.gz"
    cp "${sample_name}.pass.vcf.gz" "${sample_name}.merged.gvcf.vcf.gz"
    # This is picked up later in process allConsensus
    echo -e ">${sample_name} Artic-Fail\nN" > "${sample_name}.consensus.fasta"
}

# name everything by the sample rather than barcode
if [[ "${fastq_file}" != "${sample_name}.fastq.gz" ]]; then
    echo "Moving input: '${fastq_file}' to '${sample_name}/${sample_name}.fastq.gz'"
    mkdir ${sample_name}
    mv ${fastq_file} ${sample_name}/${sample_name}.fastq.gz
fi

artic guppyplex --skip-quality-check \
    --min-length ${min_len} --max-length ${max_len} \
    --directory ${sample_name} --prefix ${sample_name} \
    && echo " - artic guppyplex finished"
# the output of the above will be...
READFILE="${sample_name}_${sample_name}.fastq"

artic minion --medaka --normalise ${normalise} --threads ${threads} \
    --read-file ${READFILE} \
    --medaka-model ${medaka_model} \
    --scheme-directory ${scheme_dir} \
    --scheme-version ${scheme_version} \
    --max-softclip-length ${max_softclip_length} \
    ${scheme_name} ${sample_name} &> ${sample_name}.artic.log.txt \
    || mock_artic

for vcf_set in "pass" "merged.gvcf"; do
    zcat "${sample_name}.${vcf_set}.vcf.gz" | sed "s/SAMPLE/${sample_name}/" | bgzip > "${sample_name}.${vcf_set}.named.vcf.gz"
    bcftools index -t "${sample_name}.${vcf_set}.named.vcf.gz"
done;

# rename the consensus sequence
sed -i "s/^>\S*/>${sample_name}/" "${sample_name}.consensus.fasta"

# calculate depth stats. Final output is single file annotated with primer set and sample name
for i in 1 2; do
    bam="${sample_name}.primertrimmed.${i}.sorted.bam"
    samtools view -r ${i} -b "${sample_name}.primertrimmed.rg.sorted.bam" > ${bam}
    samtools index ${bam}
    stats_from_bam ${bam} > "${bam}.stats" || echo "stats_from_bam failed, probably no alignments"
    coverage_from_bam -s 20 -p ${bam} ${bam}
    # TODO: we're assuming a single reference sequence here
    awk 'BEGIN{OFS="\t"}{if(NR==1){print $0, "sample_name", "primer_set"}else{print $0, "'${sample_name}'", "'${i}'"}}' *${bam}*".depth.txt" > "${sample_name}.depth.${i}.txt"
    rm -rf ${bam} ${bam}.bai
done
cat "${sample_name}.depth.1.txt" <(tail -n+2 "${sample_name}.depth.2.txt") > "${sample_name}.depth.txt"
