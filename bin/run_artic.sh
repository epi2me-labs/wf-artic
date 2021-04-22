#!/bin/bash
set -euo pipefail

sample_name=$1
directory=$2
min_len=$3
max_len=$4
medaka_model=$5
full_scheme_name=$6
threads=$7

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
    # This is picked up later in process allConsensus
    echo -e ">${sample_name} Artic-Fail\nN" > "${sample_name}.consensus.fasta"
}

# name everything by the sample rather than barcode
if [[ "${directory}" != "${sample_name}" ]]; then
    echo "Moving input: '${directory}' to '${sample_name}'"
    mv ${directory} ${sample_name}
fi

artic guppyplex --skip-quality-check \
    --min-length ${min_len} --max-length ${max_len} \
    --directory ${sample_name} --prefix ${sample_name} \
    && echo " - artic guppyplex finished"
# the output of the above will be...
READFILE="${sample_name}_${sample_name}.fastq"

artic minion --medaka --normalise 200 --threads ${threads} \
    --read-file ${READFILE} \
    --medaka-model ${medaka_model} \
    --scheme-directory primer_schemes \
    ${full_scheme_name} ${sample_name} \
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
