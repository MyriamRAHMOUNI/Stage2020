#!/bin/bash
DIR_1000_genome='/media/bionasNFS2/1KG_PHASE3_OCT2014/VCF'
DIR_results_1000_genome='/media/sigridNFS2/myriam/Predixcan/1000_genomes'
BATCH="/media/sigridNFS2/myriam/Predixcan/BACH"
VCFtools="/media/sigridNFS2/myriam/Predixcan/executables/vcftools"
for i in {1..22}
	do 
		BATCHFILE="${BATCH}/commande.vcftools.chr${i}"
		echo -e " # ${VCFtools} --gzvcf ${DIR_1000_genome}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.MAC3.vcf.gz --remove-indels --recode --out ${DIR_results_1000_genome}/chr${i}_SNPs_only\n # sed '/##/d' ${DIR_results_1000_genome}/chr${i}_SNPs_only.recode.vcf > ${DIR_results_1000_genome}/chr${i}_SNPs_only_prop.recode.vcf\n #rm ${DIR_results_1000_genome}/chr${i}_SNPs_only.recode.vcf\n awk < ${DIR_results_1000_genome}/chr${i}_SNPs_only_prop.recode.vcf '(NR>1){chr=\$1;pos=\$2;rs=\$3}(NR=1) && (rs!=\".\") {print chr\":\"pos, rs}' > ${DIR_results_1000_genome}/chr${i}_SNPs_only_prop.rs\n #rm ${DIR_results_1000_genome}/chr${i}_SNPs_only_prop.recode.vcf" > ${BATCHFILE}
		qsub -q ALL ${BATCHFILE}
	done
