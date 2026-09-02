#!/bin/bash

FILTERED_VCF="riverine_cichlids.vcf.gz"
SITES=("AstGiIfa" "AstGiMoh" "AstRuhBKid" "AstRuhBKih" "AstRuhBMtr")


echo "Working VCF: ${FILTERED_VCF}"
echo "===================================================="

# --- Loop Through Site-Level Populations ---
for SITE in "${SITES[@]}"; do
    echo "Processing site: ${SITE}"
    
    # Extract sample names directly from the clean VCF header
    # and format as TWO columns (FID and IID) identically for PLINK/VCFtools compatibility
    bcftools query -l ${FILTERED_VCF} | grep "^${SITE}" | awk '{print $1, $1}' > ${SITE}_samples.txt
    
    NUM_SAMPLES=$(wc -l < ${SITE}_samples.txt)
    echo " -> Found ${NUM_SAMPLES} high-quality individuals."

    # Skip population if no individuals match
    if [ "$NUM_SAMPLES" -eq 0 ]; then
        echo "Warning: No individuals found for ${SITE}. Skipping."
        continue
    fi

    # 1. Site-Level Tajima's D (10kb windows)
    vcftools --gzvcf ${FILTERED_VCF} \
             --keep ${SITE}_samples.txt \
             --TajimaD 10000 \
             --out ${SITE}_tajima

    # 2. Site-Level Nucleotide Diversity (10kb windows)
    vcftools --gzvcf ${FILTERED_VCF} \
             --keep ${SITE}_samples.txt \
             --window-pi 10000 \
             --window-pi-step 5000 \
             --out ${SITE}_pi
             
    # 3. Runs of Homozygosity via PLINK
    plink --vcf ${FILTERED_VCF} \
      --keep ${SITE}_samples.txt \
      --allow-extra-chr \
      --homozyg \
      --homozyg-window-snp 50 \
      --homozyg-window-het 1 \
      --homozyg-snp 100 \
      --homozyg-kb 500 \
      --homozyg-gap 2000 \
      --homozyg-density 100 \
      --out ${SITE}_roh
done

echo "Workflow complete. Ready for R consolidation."
