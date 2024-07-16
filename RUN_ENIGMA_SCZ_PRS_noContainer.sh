#!/bin/bash

############### Please Fill in the paths below

## Path to project directory with all downloaded files
export Project_Path=/PATH/TO/DOWNLOADED/FOLDER/ENIGMA_SCZ_PRS

## Path to Your Sample's PLINK Files
export Sample_Dir=/PATH/TO/GENETIC/DATA

## Plink Sample Prefix - the plink file prefix (filename before .bed, .fam, .bim suffix)
export Prefix=YOUR_PLINK_FILES_PREFIX

## Sample name - if you plan to run this protocol for multiple samples this will be used as an identifier
export Sample_Name=YOUR_SAMPLE_NAME

## number of available cores for processing
export NCORES=8

export MEMORY=$((${NCORES} * 7700))

## If you plan to run this script for multiple samples in parallel from the same project directory
# please specify clean as no below - some intermediate files will remain to avoid conflicts
CLEAN=yes

## Path to downloaded Plink v1.9 executable file
export PLINK=/path/to/plink

## Path to downloaded Plink v2 executable file
export PLINK2=/path/to/plink2

# Ensure you have installed (or if on a server loaded like below): R v4.0.0+, Python v3.7.4+ 
# Also ensure the necessary packages installed for both R and Python
# module load R/4.0.0-foss-2020a
# module load  Python/3.7.4-GCCcore-8.3.0
# source /cluster/projects/p697/users/nadinepa/py3/bin/activate


############### END NO NEED TO CHANGE ANYTHING BEYOND THIS POINT ###################

set -euo pipefail 

export OUTDIR=${Project_Path}/output_${Sample_Name}
export RSCRIPT=Rscript
export PYTHON=python


mkdir -p ${OUTDIR}

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>${OUTDIR}/RUN_ENIGMA_SCZ_PRS_${Sample_Name}.log 2>&1

date '+%A %W %Y %X'


## Running initial data check
${RSCRIPT} ${Project_Path}/scripts/DataCheck.R \
	--target-dir ${Sample_Dir} \
	--prefix ${Prefix} \
	--project-dir ${Project_Path} \
	--out-dir ${Project_Path}/output_${Sample_Name}

cat ${Project_Path}/output_${Sample_Name}/LOG_DataCheck

### Make non-duplicated SNP plink files

OUT_DIR=${Project_Path}/${Sample_Name}_tmp
mkdir -p ${OUT_DIR}

${PLINK} --bfile ${Sample_Dir}/${Prefix} \
	--threads ${NCORES} --memory ${MEMORY} \
	--extract ${OUTDIR}/NonDup_SNPlist.txt \
	--make-bed --out ${OUT_DIR}/${Sample_Name}_rmDupSNP


### Ancestry + PCs
# set up singularity plink2

# Input parameters.
target_bfile=${OUT_DIR}/${Sample_Name}_rmDupSNP	#${Sample_Dir}/${Prefix}
target_name=${Sample_Name}
OUT_PREFIX=${Project_Path}
REF_DIR=${Project_Path}/ref
plink_prefix_1kg=${REF_DIR}/1000_Genomes_2022_from_plink2/all_hg38
ref_fasta=${REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst
deg2_relatives=${REF_DIR}/1000_Genomes_2022_from_plink2/deg2_hg38.king.cutoff.out.id


# Variables derived from input parameters.
pvar_chr_pos_ref_alt=${REF_DIR}/1000_Genomes_2022_from_plink2/all_hg38.chr_pos_ref_alt
snplist_prefix=${REF_DIR}/1000_Genomes_2022_from_plink2/all_hg38
common_snps_1kg=${REF_DIR}/1000_Genomes_2022_from_plink2/all_hg38.common_snps
target_pvar_chr_pos_ref_alt=${OUT_DIR}/${target_name}.chr_pos_ref_alt
common_snps_1kg_target=${OUT_DIR}/${target_name}.common_snps 	# rename not target specific

### Produces one file that we can comment out and add to the files delivered to the sites
## alters 1000G referece
if false; then
    # Processing of the reference files needs to be done only once (no need to repeat for each target dataset).
    # Make pvar file with chr:pos:ref:alt ids
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pfile ${plink_prefix_1kg} vzs \
        --set-all-var-ids "@:#:\$r:\$a" --new-id-max-allele-len 23 missing \
        --make-just-pvar zs cols="" --out ${pvar_chr_pos_ref_alt}

    # Identify common SNPs in each ancestry.
    POPULATIONS="EUR AFR SAS EAS AMR"
    for ANCESTRY in ${POPULATIONS}; do
    out=${snplist_prefix}.${ANCESTRY}
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pgen ${plink_prefix_1kg}.pgen --psam ${plink_prefix_1kg}.psam \
        --pvar ${pvar_chr_pos_ref_alt}.pvar.zst --autosome --remove ${deg2_relatives} --snps-only \
        --keep-if "SuperPop=${ANCESTRY}" --maf 0.01 --geno 0.01 --rm-dup 'force-first' --write-snplist --out ${out}
    done

    # Identify SNPs which are common in all ancestries
    N_POP=$(echo ${POPULATIONS} | awk '{print(NF)}')
    sort --parallel=${NCORES} ${snplist_prefix}.*.snplist | uniq -c | awk -v n_pop="${N_POP}" '$1==n_pop {print $2}' > ${common_snps_1kg}
    echo "N common SNPs = $(wc -l ${common_snps_1kg})"
fi
## above will be run before and output shared with sites -- can comment out in future


if true; then
    # See https://groups.google.com/g/plink2-users/c/ohV0F2e5onk for more details about the logic behind the following processing steps.
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --bfile ${target_bfile} --make-just-pvar cols=+info \
        --autosome --snps-only 'just-acgt' --ref-from-fa --fa ${ref_fasta} \
        --out "${OUT_DIR}/${target_name}.allign_to_ref"
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pvar "${OUT_DIR}/${target_name}.allign_to_ref.pvar" --write-snplist \
        --require-info PR --out "${OUT_DIR}/${target_name}.ref_allele_missmatch"
    ${PLINK} --threads ${NCORES} --memory ${MEMORY} --bfile ${target_bfile} --allow-extra-chr --make-just-bim \
        --flip "${OUT_DIR}/${target_name}.ref_allele_missmatch.snplist" --out "${OUT_DIR}/${target_name}.flip"
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --bed ${target_bfile}.bed --fam ${target_bfile}.fam --bim "${OUT_DIR}/${target_name}.flip.bim" \
        --make-just-pvar cols=+info --autosome --snps-only 'just-acgt' --ref-from-fa --fa ${ref_fasta} \
        --out "${OUT_DIR}/${target_name}.flip.allign_to_ref"
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pvar "${OUT_DIR}/${target_name}.flip.allign_to_ref.pvar" --write-snplist \
        --require-info PR --out "${OUT_DIR}/${target_name}.flip.ref_allele_missmatch"
    # --mac 1 filter in the following command ensures that there are no "0" alleles in the new variant IDs, otherwise it can be 1:65342:A:0
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --bed ${target_bfile}.bed --fam ${target_bfile}.fam --bim "${OUT_DIR}/${target_name}.flip.bim" \
        --make-pgen vzs 'pvar-cols=' --autosome --snps-only 'just-acgt' --mac 1 --rm-dup 'force-first' \
        --ref-from-fa --fa ${ref_fasta} --exclude "${OUT_DIR}/${target_name}.flip.ref_allele_missmatch.snplist" \
        --set-all-var-ids "@:#:\$r:\$a" --new-id-max-allele-len 23 missing --out ${target_pvar_chr_pos_ref_alt}

    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pfile ${target_pvar_chr_pos_ref_alt} vzs \
        --extract ${common_snps_1kg} --write-snplist --out ${common_snps_1kg_target}
fi


# Get independent SNPs in 1000 Genomes data taking into account only SNPs which also present in in-house data.
indep_snps="${OUT_DIR}/1kg_indep_snps"
if true; then
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pgen ${plink_prefix_1kg}.pgen --psam ${plink_prefix_1kg}.psam \
        --pvar ${pvar_chr_pos_ref_alt}.pvar.zst --extract ${common_snps_1kg_target}.snplist \
        --remove ${deg2_relatives} --nonfounders --indep-pairwise 200kb 0.5 --out ${indep_snps}
fi


## Generate Sample Related individuals
${PLINK2} --pgen ${target_pvar_chr_pos_ref_alt}.pgen --psam ${target_pvar_chr_pos_ref_alt}.psam --pvar ${target_pvar_chr_pos_ref_alt}.pvar.zst \
	--extract ${indep_snps}.prune.in \
	--threads ${NCORES} --memory ${MEMORY} \
	--make-king-table --king-table-filter 0.1 --out ${OUTDIR}/${Sample_Name}_rel


# Calculate PCA for 1000 Genomes unrelated
pca_1kg="${OUT_DIR}/pca_1kg"
n_pc=20
if true; then
   ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pgen ${plink_prefix_1kg}.pgen --psam ${plink_prefix_1kg}.psam \
        --pvar ${pvar_chr_pos_ref_alt}.pvar.zst --extract ${indep_snps}.prune.in \
        --remove ${deg2_relatives} --nonfounders --pca ${n_pc} 'allele-wts' --freq counts \
        --out ${pca_1kg}
fi

# Project 1000 Genomes PCA
pca_1kg_project="${OUT_DIR}/pca_1kg_proj"
if true; then
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pgen ${plink_prefix_1kg}.pgen --psam ${plink_prefix_1kg}.psam \
        --pvar ${pvar_chr_pos_ref_alt}.pvar.zst --read-freq ${pca_1kg}.acount \
        --score ${pca_1kg}.eigenvec.allele 2 5 'header-read' 'no-mean-imputation' 'variance-standardize' \
        --score-col-nums 6-$((${n_pc} + 5)) --out ${pca_1kg_project}
fi

# Project Sample PCA
pca_top_project="${OUT_DIR}/pca_${target_name}_proj"
if true; then
    ${PLINK2} --threads ${NCORES} --memory ${MEMORY} --pfile ${target_pvar_chr_pos_ref_alt} vzs --read-freq ${pca_1kg}.acount \
        --score ${pca_1kg}.eigenvec.allele 2 5 'header-read' 'no-mean-imputation' 'variance-standardize' \
        --score-col-nums 6-$((${n_pc} + 5)) --out ${pca_top_project}
fi


# Predict ancestry and produce PCA plots with color-coded ancestry.
${PYTHON} --version

${PYTHON} ${Project_Path}/scripts/anc_proj.py "${pca_1kg_project}.sscore" "${pca_top_project}.sscore" "${pca_top_project}"

mv ${OUT_DIR}/pca_${Sample_Name}_proj.ancestries.png ${OUTDIR}
mv ${OUT_DIR}/pca_${Sample_Name}_proj.ancestries.txt ${OUTDIR}

## Clean up tmp dir
rm ${Project_Path}/${Sample_Name}_tmp/*1kg*
rm ${Project_Path}/${Sample_Name}_tmp/*ref*
rm ${Project_Path}/${Sample_Name}_tmp/*flip*
rm ${Project_Path}/${Sample_Name}_tmp/*pca*
rm ${Project_Path}/${Sample_Name}_tmp/*common*

####### Genetic QC

## split anc_proj
${RSCRIPT} ${Project_Path}/scripts/MakeAncProjLists.R \
	--project-dir ${Project_Path} \
	--anc-dir ${OUTDIR} \
	--sample-name ${Sample_Name} \
	--out-dir ${Sample_Name}_tmp

## make plink files for each anc_proj
mkdir -p ${Project_Path}/anc_proj_files

cd ${Project_Path}/${Sample_Name}_tmp

for file in ${Sample_Name}_ancproj_AFR ${Sample_Name}_ancproj_AMR ${Sample_Name}_ancproj_EAS ${Sample_Name}_ancproj_EUR ${Sample_Name}_ancproj_SAS
do
	if [ -f "$file" ]; 
	then
		${PLINK} --bfile ${OUT_DIR}/${Sample_Name}_rmDupSNP --threads ${NCORES} --memory ${MEMORY} --keep ${file} --make-bed --out ${Project_Path}/anc_proj_files/${file}
	else
		echo "$file does not exist"
	fi
done



## anc_proj specific QC

cd ${Project_Path}/anc_proj_files

for ancFile in ${Sample_Name}_ancproj_AFR ${Sample_Name}_ancproj_AMR ${Sample_Name}_ancproj_EAS ${Sample_Name}_ancproj_EUR ${Sample_Name}_ancproj_SAS
do

	if [ -f "${ancFile}.fam" ];
        then

		count=$(wc -l ${ancFile}.fam | awk '{ print $1 }')
	
		if [ ${count} -gt 10 ];
		then
			
			## Initial SNP QC
			${PLINK} --bfile ${ancFile} --threads ${NCORES} --memory ${MEMORY} --maf 0.01 --hwe 1e-10 midp --geno 0.01 --write-snplist --out ${ancFile}.qc.snps

			## Additional Sample and SNP QC
			${PLINK} --bfile ${ancFile} --threads ${NCORES} --memory ${MEMORY} --extract ${ancFile}.qc.snps.snplist --mind 0.01 --maf 0.01 --hwe 1e-50 midp --geno 0.01 --write-snplist --make-just-fam --out ${ancFile}.qc.samples 

			## Heterozygosity
			${PLINK} --bfile ${ancFile} --threads ${NCORES} --memory ${MEMORY} --extract ${ancFile}.qc.samples.snplist --keep ${ancFile}.qc.samples.fam --het --out ${ancFile}.qc

			## R to remove outliers in heterozygosity
			$RSCRIPT ${Project_Path}/scripts/het_rm.R --het-file ${Project_Path}/anc_proj_files/${ancFile}.qc.het

			## Generate final QC files
			${PLINK} --bfile ${ancFile} --threads ${NCORES} --memory ${MEMORY} --keep ${ancFile}.valid.sample --extract ${ancFile}.qc.samples.snplist --make-bed --out ${ancFile}.FinalQC

		else
			
			## Initial SNP QC
			${PLINK} --bfile ${ancFile} --threads ${NCORES} --memory ${MEMORY} --geno 0.5 --write-snplist --out ${ancFile}.qc.snps

			## Additional Sample and SNP QC
			${PLINK} --bfile ${ancFile} --threads ${NCORES} --memory ${MEMORY} --extract ${ancFile}.qc.snps.snplist --mind 0.01 --geno 0.5 --write-snplist --make-just-fam --out ${ancFile}.qc.samples 

			## Generate final QC files
			${PLINK} --bfile ${ancFile} --threads ${NCORES} --memory ${MEMORY} --keep ${ancFile}.qc.samples.fam --extract ${ancFile}.qc.samples.snplist --make-bed --out ${ancFile}.FinalQC

		fi

	else
		echo "${ancFile}.fam does not exist"
	fi
done



cd ${Project_Path}/anc_proj_files

for ancFile in ${Sample_Name}_ancproj_AFR.FinalQC ${Sample_Name}_ancproj_AMR.FinalQC ${Sample_Name}_ancproj_EAS.FinalQC ${Sample_Name}_ancproj_EUR.FinalQC ${Sample_Name}_ancproj_SAS.FinalQC
do

	if [ -f "${ancFile}.fam" ];
        then
		echo ${ancFile} >> ${Sample_Name}_list.txt
	else
		echo "${ancFile} does not exist"
	fi
done

REF_PLINKFILE=$(head -n 1 ${Sample_Name}_list.txt)
echo "${REF_PLINKFILE} is the ref for merging"
tail -n +2 ${Sample_Name}_list.txt > ${Sample_Name}_merge_list.txt

plinkcount=$(wc -l ${Sample_Name}_merge_list.txt | awk '{ print $1 }')
echo "${plinkcount} + 1 plink files to merge"
	
if [ ${plinkcount} -gt 0 ];
then
	${PLINK} --bfile ${REF_PLINKFILE} --merge-list ${Sample_Name}_merge_list.txt --out ${Sample_Name}_QC_ENIGMA_SCZ
else
	cp ${REF_PLINKFILE}.bed ${Sample_Name}_QC_ENIGMA_SCZ.bed
	cp ${REF_PLINKFILE}.bim ${Sample_Name}_QC_ENIGMA_SCZ.bim
	cp ${REF_PLINKFILE}.fam ${Sample_Name}_QC_ENIGMA_SCZ.fam
fi

## Clean up anc_proj_files
rm ${Sample_Name}_list.txt
rm ${Sample_Name}_merge_list.txt
rm ${Sample_Name}_ancproj*

## Clean up tmp dir
rm ${Project_Path}/${Sample_Name}_tmp/*ancproj*
rm ${Project_Path}/${Sample_Name}_tmp/${Sample_Name}_rmDupSNP*

date '+%A %W %Y %X'

### Run PRS
cd ${Project_Path}

# Convert from plink format to bigSNPR .rds/.bk files
${RSCRIPT} ${Project_Path}/scripts/createBackingFile2.R \
	--file-input ${Project_Path}/anc_proj_files/${Sample_Name}_QC_ENIGMA_SCZ.bed \
	--file-output ${Project_Path}/${Sample_Name}_tmp/${Sample_Name}_QC_ENIGMA_SCZ.rds

# impute missing values
${RSCRIPT} ${Project_Path}/scripts/imputeGenotypes2.R \
	--impute-simple mean0 \
	--geno-file-rds ${Project_Path}/${Sample_Name}_tmp/${Sample_Name}_QC_ENIGMA_SCZ.rds \
	--cores $NCORES

export fileGenoRDS=${Project_Path}/${Sample_Name}_tmp/${Sample_Name}_QC_ENIGMA_SCZ.rds
export LD_FILES=${Project_Path}/ref/ldref_hm3_plus/LD_with_blocks_chr@.rds

# Generate PGS using LDPRED2-auto
cd ${Project_Path}/sumstats

for SUMSTATS in PGC_SCZ*GRCh38.sumstats.gz
do
	NEW_NAME=$(echo ${SUMSTATS} | awk -F'[.]' '{print $1}')

	${RSCRIPT} ${Project_Path}/scripts/ldpred2v2.R \
		 --ldpred-mode auto \
		 --ld-file ${LD_FILES} \
		 --ld-meta-file ${Project_Path}/ref/map_hm3_plus.rds \
		 --cores ${NCORES} \
		 --col-stat B \
		 --col-stat-se SE \
		 --col-chr CHR \
        	 --col-snp-id RSID \
		 --col-A1 EffectAllele \
		 --col-A2 OtherAllele \
		 --col-bp POS \
		 --col-pvalue P \
		 --col-n N \
		 --stat-type BETA \
		 --genomic-build hg38 \
		 --geno-file-rds ${fileGenoRDS} \
		 --sumstats ${SUMSTATS} \
		 --out ${Project_Path}/${Sample_Name}_tmp/${Sample_Name}_${NEW_NAME}_prs.txt
done

## Merge ldpred prs outputs and write to outdir
${RSCRIPT} ${Project_Path}/scripts/MergeLDpred2PRS.R \
	--sample-name ${Sample_Name} \
	--project-dir ${Project_Path} \
	--out-dir ${OUTDIR}

date '+%A %W %Y %X'

## Clean up tmp dir
rm ${Project_Path}/${Sample_Name}_tmp/*bk
rm ${Project_Path}/${Sample_Name}_tmp/*rds
rm -r ${Project_Path}/${Sample_Name}_tmp
if [ ${CLEAN} == yes ]
        then
		rm -r ${Project_Path}/anc_proj_files
	else
		echo "Clean is no - some intermediate files will remain"
	fi

date '+%A %W %Y %X'