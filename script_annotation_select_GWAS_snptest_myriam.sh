#!/bin/bash
chromosome="1"

for chr in ${chromosome}



	do
	echo "CERIES_B37.CHR"${chr}".larnier1.MODE_1.SNPTEST.gz"
	zcat /media/sigrid/Lacie/genotype/6_annotation/"CERIES_B37.CHR"${chr}".larnier1.MODE_1.SNPTEST.gz" | grep -v '^#'> "CERIES_B37.CHR"${chr}".larnier1.MODE_1.SNPTEST"
	
	perl script_annotation_select_GWAS_plink.pl larnier_final_TWAS_10-3 "CERIES_B37.CHR"${chr}".larnier1.MODE_1.SNPTEST" 10000 "CERIES_B37.CHR"${chr}".larnier1.MODE_1.SNPTEST.select_GWAS" 0.05 ${chr}
	rm "CERIES_B37.CHR"${chr}".larnier1.MODE_1.SNPTEST"
	done

#echo -e "ensembl_id gene chr start stop p" > final_TWAS
#cat *.select_TWAS >> final_TWAS

###########################################################
models="en_Adipose_Subcutaneous en_Brain_Cortex en_Skin_Sun_Exposed_Lower_leg en_Brain_Cerebellum en_Skin_Not_Sun_Exposed_Suprapubic"
phenos="Lentigine Relachement Ride"

#Lentigine Mrahmo

../GWAS_results_CNAM/${pheno}/CERIES.chr15.mode1.${pheno}_z.snptest.gz
CNAM_final_TWAS_${model}_${pheno}_0.0005

for pheno in ${phenos}; do  for chr in {1..22} ; do echo "../GWAS_results_CNAM/${pheno}/CERIES.chr"${chr}".mode1.${pheno}_z.snptest.gz"; zcat "../GWAS_results_CNAM/${pheno}/CERIES.chr"${chr}".mode1.${pheno}_z.snptest.gz"| grep -v '^#' > CERIES.chr"${chr}".mode1.${pheno}_z.snptest; done; done


pheno="Lentigine"
for model in ${models}; do for chr in {1..22}; do perl script_annotation_select_GWAS_snptest_myriam.pl CNAM_final_TWAS_${model}_${pheno}_0.0005 CERIES.chr"${chr}".mode1.${pheno}_z.snptest 10000 ${chr}".${model}.${pheno}.select_GWAS_T0005_G0005" 0.0005 ${chr}; done; done

for model in ${models}
	do
	cat *.${model}.${pheno}.select_GWAS_T0005_G0005 > ${model}.Lentigine.select_GWAS_T0005_G0005.all
	rm *.${model}.${pheno}.select_GWAS_T0005_G0005
	done