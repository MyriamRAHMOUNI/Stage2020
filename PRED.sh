#!/bin/bash
ssh leclersi@vlad.cnam.fr #fug7a7wy
ssh sigrid@bioclustnew01.cnam.fr

#********************************METAXCAN*************************************
#Préparation des fichiers de 1000 gènome
/. qsub_1000g.sh 
#supprimer la 1ere ligne ":" 
cd /media/sigridNFS2/myriam/Predixcan/1000_genomes
for i in {2..22} 
	do 
	sed -i '1d' chr${i}_SNPs_only_prop.rs
	done 
cat chr1_SNPs_only_prop.rs chr2_SNPs_only_prop.rs chr3_SNPs_only_prop.rs chr4_SNPs_only_prop.rs chr5_SNPs_only_prop.rs chr6_SNPs_only_prop.rs chr7_SNPs_only_prop.rs chr8_SNPs_only_prop.rs chr9_SNPs_only_prop.rs chr10_SNPs_only_prop.rs chr11_SNPs_only_prop.rs chr12_SNPs_only_prop.rs chr13_SNPs_only_prop.rs chr14_SNPs_only_prop.rs chr15_SNPs_only_prop.rs chr16_SNPs_only_prop.rs chr17_SNPs_only_prop.rs chr18_SNPs_only_prop.rs chr19_SNPs_only_prop.rs chr20_SNPs_only_prop.rs chr21_SNPs_only_prop.rs chr22_SNPs_only_prop.rs > all_chr_SNPs_only_prop.rs

#preparation des fichiers pour metaxcan 

	#taizou 
	DIR_chinois="/media/cedcoulNFS2/GWAS/CERIES/data_TAIZHOU/TAIZHOU_20180909.txt"
	DIR_1000g="/media/sigridNFS2/myriam/Predixcan/1000_genomes/all_chr_SNPs_only_prop.rs" 
	DIR_results="/media/sigridNFS2/myriam/Predixcan/chinois" 
 	awk 'NR==FNR{a[$1]=$1"\t"$4"\t"$5"\t"$10"\t"$13"\t"$14"\t"$17"\t"$18"\t"$21; next} ($1 in a){print $2"\t"a[$1]}' ${DIR_chinois} ${DIR_1000g} > ${DIR_results}/TAIZOU_propre.rs
	awk '(NR>=1) {print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' < ${DIR_results}/TAIZOU_propre.rs > ${DIR_results}/TAIZOU_all_pheno.rs
echo -e "SNP\tRef\tAlt\tPigmentedSpots_Estimate\tPigmentedSpots_Pvalue\tWrinkle_Estimate\tWrinkle_Pvalue\tLaxity_Estimate\tLaxity_Pvalue"> ${DIR_results}/header 
cat  ${DIR_results}/header ${DIR_results}/TAIZOU_all_pheno.rs > ${DIR_results}/TAIZOU_all.rs #7259573
rm  ${DIR_results}/header ${DIR_results}/TAIZOU_all_pheno.rs ${DIR_results}/TAIZOU_propre.rs 

	#salia
	DIR_results="/media/sigridNFS2/myriam/Predixcan/salia" 
	DIR_salia_lentigines="/media/cedcoulNFS2/GWAS/CERIES/GWAS_IUF_Z-score/lentigines"
	#PHENOS="lentigines wrinkles sagging"
	#PHENOS="flecken falten laxity"	
	for i in {1..22}
		do
		awk '(NR>=1 && $2!="POS"){print $24, $10, $11, $3, $6}' ${DIR_salia_lentigines}/GWAS_flecken_z_imputed_genotypes_chr${i}.txt > ${DIR_results}/salia_lentigines_ch${i}_NA.rs #chr1 746208
		sed '/NA/d' <${DIR_results}/salia_lentigines_ch${i}_NA.rs > ${DIR_results}/salia_lentigines_ch${i}.rs 
		rm ${DIR_results}/salia_lentigines_ch${i}_NA.rs #chr1 739923
	done

	echo "ID REF.0. ALT.1. Estimate P_val" > header 
	cat header salia_lentigines_ch10.rs salia_lentigines_ch20.rs salia_lentigines_ch11.rs salia_lentigines_ch21.rs salia_lentigines_ch12.rs  salia_lentigines_ch22.rs salia_lentigines_ch13.rs salia_lentigines_ch2.rs salia_lentigines_ch14.rs salia_lentigines_ch3.rs salia_lentigines_ch15.rs salia_lentigines_ch4.rs salia_lentigines_ch16.rs salia_lentigines_ch5.rs salia_lentigines_ch17.rs salia_lentigines_ch6.rs salia_lentigines_ch18.rs salia_lentigines_ch7.rs salia_lentigines_ch19.rs salia_lentigines_ch8.rs salia_lentigines_ch1.rs salia_lentigines_ch9.rs > salia_lentigines.rs
	rm header

DIR_salia_wrinkles="/media/cedcoulNFS2/GWAS/CERIES/GWAS_IUF_Z-score/wrinkles"
for i in {1..22}
	do
	awk '(NR>1 && $2!="POS"){print $24, $10, $11, $3, $6}' ${DIR_salia_wrinkles}/GWAS_falten_z_imputed_genotypes_chr${i}.txt > ${DIR_results}/salia_wrinkles_ch${i}_NA.rs 
	sed '/NA/d' <${DIR_results}/salia_wrinkles_ch${i}_NA.rs > ${DIR_results}/salia_wrinkles_ch${i}.rs 
	rm ${DIR_results}/salia_wrinkles_ch${i}_NA.rs 
done	

echo "ID REF.0. ALT.1. Estimate P_val" > header 
cat header salia_wrinkles_ch10.rs salia_wrinkles_ch11.rs salia_wrinkles_ch12.rs salia_wrinkles_ch13.rs salia_wrinkles_ch14.rs salia_wrinkles_ch15.rs salia_wrinkles_ch16.rs salia_wrinkles_ch17.rs salia_wrinkles_ch18.rs salia_wrinkles_ch19.rs salia_wrinkles_ch1.rs salia_wrinkles_ch20.rs salia_wrinkles_ch21.rs salia_wrinkles_ch22.rs salia_wrinkles_ch2.rs salia_wrinkles_ch3.rs salia_wrinkles_ch4.rs salia_wrinkles_ch5.rs salia_wrinkles_ch6.rs salia_wrinkles_ch7.rs salia_wrinkles_ch8.rs salia_wrinkles_ch9.rs > salia_wrinkles.rs
rm header

DIR_salia_sagging="/media/cedcoulNFS2/GWAS/CERIES/GWAS_IUF_Z-score/sagging"	
for i in {1..22}
	do
	awk '(NR>1){print $24, $10, $11, $3, $6}' ${DIR_salia_sagging}/GWAS_laxity_z_imputed_genotypes_chr${i}.txt > ${DIR_results}/salia_sagging_ch${i}_NA.rs 
	sed '/NA/d' <${DIR_results}/salia_sagging_ch${i}_NA.rs > ${DIR_results}/salia_sagging_ch${i}.rs 
	rm ${DIR_results}/salia_sagging_ch${i}_NA.rs 
done

echo "ID REF.0. ALT.1. Estimate P_val" > header 
cat header salia_sagging_ch20.rs salia_sagging_ch21.rs salia_sagging_ch22.rs salia_sagging_ch2.rs salia_sagging_ch3.rs salia_sagging_ch4.rs salia_sagging_ch5.rs salia_sagging_ch6.rs salia_sagging_ch7.rs salia_sagging_ch8.rs salia_sagging_ch9.rs salia_sagging_ch10.rs salia_sagging_ch11.rs salia_sagging_ch12.rs salia_sagging_ch13.rs salia_sagging_ch14.rs salia_sagging_ch15.rs salia_sagging_ch16.rs salia_sagging_ch17.rs salia_sagging_ch18.rs salia_sagging_ch19.rs salia_sagging_ch1.rs > salia_sagging.rs

Spredixcan="/./media/sigridNFS2/myriam/Predixcan/executables/software/SPrediXcan.py"
Pred_model="/media/sigridNFS2/myriam/Predixcan/DGN-WB_0.5.db"
Covar_file="/media/sigridNFS2/myriam/Predixcan/covariance.DGN-WB_0.5.txt.gz"
Gwas_file="/media/sigridNFS2/myriam/Predixcan/salia/salia_wrinkles.rs"
DIR_results="/media/sigridNFS2/myriam/Predixcan/salia/results"
#lancer metaxcan 
	#use the PredictDB pipeline to generate my own prediction models :https://github.com/hakyimlab/PredictDBPipeline/wiki/Tutorial

${Spredixcan} \
--model_db_path ${Pred_model} \
--covariance ${Covar_file} \
--gwas_file ${Gwas_file} \
--snp_column ID \
--effect_allele_column REF.0. \
--non_effect_allele_column ALT.1. \
--beta_column Estimate \
--pvalue_column P_val \
--output_file ${DIR_results}/salia_wrinkles.res
#***********************************PREDIXCAN********************************************



DIR_Cnam="/media/cedcoulNFS2/GWAS/CERIES/4_IMPUTE"
${DIR_Cnam}/chr${i}.dose.maf0.01_info0.3.vcf.gz












