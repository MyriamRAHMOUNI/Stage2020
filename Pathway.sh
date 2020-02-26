#Je suis sur HG37
# seuil ld 0.8  
#info score 0.3

########################################preparer la liste des snp à garder dans l'analyse (snp presents dans le fichier snpgenemap.hg19 sans ld eur et eas)

#recuperer la liste des snp presents dans le fichier snpgenemap.hg19 => list_extract_snp
cd Myriam_2020/data_myriam/
awk < snpgenemap.hg19 '(NR>1){rs=$1} (NR>1){print rs}' > 1000genomes/list_extract_snp

#recuperer les id des individus non reliés dans les ped de 1000 genome pour eur => EUR.txt et eas => EAS.txt
cd 1000genomes
awk < integrated_call_samples_v2.20130502.ALL.ped '(NR>1){id=$2;pop=$7;related=$8} (NR>=1) && (related=="unrel") && (pop=="CEU"||pop=="TSI"||pop=="FIN"||pop=="GBR"||pop=="IBS") {print id}' > EUR.txt #317 individus
awk < integrated_call_samples_v2.20130502.ALL.ped '(NR>1){id=$2;pop=$7;related=$8} (NR>=1) && (related=="unrel") && (pop=="CHB"||pop=="JPT"||pop=="CHS"||pop=="CDX"||pop=="KHV") {print id}' > EAS.txt #437 individus 

ls *.vcf.gz > list-files-1000g
mkdir results_snp_1000g
mkdir results_eur_1000g
mkdir results_eas_1000g

for i in `cat  list-files-1000g`
	do
	#extraire list_extract_snp des fichiers de 1000g 
	vcftools --gzvcf ${i} --snps list_extract_snp --recode --out results_snp_1000g/${i::-7}_snps
	#extraire les population europenne avec individus non reliés
	vcftools --vcf results_snp_1000g/${i::-7}_snps.recode.vcf --keep EUR.txt --recode --out results_eur_1000g/${i::-7}_eur
	#extraire les population eas avec individus non reliés
	vcftools --vcf results_snp_1000g/${i::-7}_snps.recode.vcf --keep EAS.txt --recode --out results_eas_1000g/${i}_eas
	rm results_snp_1000g/${i::-7}_snps
	done

#Calcule ld garder et ne garder que les snps en ld < 0.8
cd results_eas_1000g
mkdir res_prun_eas
ls *.vcf > list-eas-1000g
for i in `cat  list-eas-1000g`
	do
	./plink2 --vcf ${i} --indep-pairwise 50 5 0.8 --out res_prun_eas/${i}_propre_eas
	done

cd ../results_eur_1000g 
ls *.vcf > list-eur-1000g
mkdir res_prun_eur
for i in `cat  list-eur-1000g`
	do
	./plink2 --vcf ${i} --indep-pairwise 50 5 0.8 --out res_prun_eur/${i}_propre_eur
	done


#Nos fichier (sauf les allemant) contiennent des ch:pos et non des rs donc il faut Retrouver les ch po de ces liste rs 
#eas
cd ../results_eas_1000g
mkdir san_##
mkdir fichiers_ch_pos_rs_propre
for i in {1..22}
	do
	sed '/##/d' ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_eas.recode.vcf > san_##/ALL.chr${i}_EAS_san_##
	awk 'NR==FNR{a[$3]=$1" "$2" "$3; next} ($1 in a){print a[$1], a[$2], a[$3]}' san_##/ALL.chr${i}_EAS_san_## res_prun_eas/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_eas.recode.vcf_propre_eas.prune.in > fichiers_ch_pos_rs_propre/chr${i}_eas_pos_rs_propre
	done
	#=> listes rs ch pos de eas
#eur
cd ../results_eur_1000g 
mkdir san_##
mkdir fichiers_ch_pos_rs_propre
for i in {1..22}
	do
	sed '/##/d' ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_eur.recode.vcf > san_##/ALL.chr${i}_EUR_san_##
	awk 'NR==FNR{a[$3]=$1" "$2" "$3; next} ($1 in a){print a[$1], a[$2], a[$3]}' san_##/ALL.chr${i}_EUR_san_## res_prun_eur/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_eur.recode.vcf_propre_eur.prune.in > fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre
	done
	#=> listes rs ch pos eur

#donc la j'ai les rs ch pos des sn que je vais prendre chez chanel et salia : /home/mrahmouni/Myriam_2020/data_myriam/1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre et chez les chinois : /home/mrahmouni/Myriam_2020/data_myriam/1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/chr${i}_eas_pos_rs_propre

############################################################CHANEL
cd Myriam_2020/data_myriam/chanel
#Déplacer les .md5 dans un dossier aprt 
mkdir md5-files 
ls *.md5 > list-md5.txt
for i in `cat list-md5.txt`
	do
	mv /home/mrahmouni/Myriam_2020/Stage2020/data_myriam_copie/chanel/${i} /home/mrahmouni/Myriam_2020/Stage2020/data_myriam_copie/chanel/md5-files
	done
#Preparations des fichiers chr par chr (gardes que la liste des snp préparée au préalable)

mkdir file_1_results 
cd file_1_results
mkdir relachement_files_1 lentigine_files_1 ride_files_1
cd lentigine 
for i in {1..22}
	do
	zcat CERIES.chr${i}.mode1.lentigine_z.snptest.gz | sed '/#/d' > CERIES.chr${i}.mode1.lentigine_z.snptest_sans#
	awk 'NR==FNR{a[$4]=$4" "$21; next} ($2 in a){print $3" "a[$2]}' CERIES.chr${i}.mode1.lentigine_z.snptest_sans# ../../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > lentigine_chr${i}_file1_pos
	rm CERIES.chr${i}.mode1.lentigine_z.snptest_sans#
	awk < lentigine_chr${i}_file1_pos '(NR>=1){snp=$1;pvalue=$3} {print snp,"\t"pvalue}' > ../file_1_results/lentigine_files_1/lentigine_chr${i}_file1	
	rm lentigine_chr${i}_file1_pos
	done
cd ..
cd relachement 

for i in {1..22}
	do
	zcat CERIES.chr${i}.mode1.relachement_z.snptest.gz | sed '/#/d' > CERIES.chr${i}.mode1.relachement_z.snptest_sans#
	awk 'NR==FNR{a[$4]=$4" "$21; next} ($2 in a){print $3" "a[$2]}' CERIES.chr${i}.mode1.relachement_z.snptest_sans# ../../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > relachement_chr${i}_file1_pos
	rm CERIES.chr${i}.mode1.relachement_z.snptest_sans#
	awk < relachement_chr${i}_file1_pos '(NR>=1){snp=$1;pvalue=$3} {print snp,"\t"pvalue}' > ../file_1_results/relachement_files_1/relachement_chr${i}_file1
	rm relachement_chr${i}_file1_pos
	done
cd .. 
cd rides

for i in {1..22}
	do
	zcat CERIES.chr${i}.mode1.ride_z.snptest.gz | sed '/#/d' > CERIES.chr${i}.mode1.ride_z.snptest_sans#
	awk 'NR==FNR{a[$4]=$4" "$21; next} ($2 in a){print $3" "a[$2]}' CERIES.chr${i}.mode1.ride_z.snptest_sans# ../../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > ride_chr${i}_file1_pos
	rm CERIES.chr${i}.mode1.ride_z.snptest_sans#
	awk < ride_chr${i}_file1_pos '(NR>=1){snp=$1;pvalue=$3} {print snp,"\t"pvalue}' > ../file_1_results/ride_files_1/ride_chr${i}_file1
	rm ride_chr${i}_file1_pos
	done
# Merger et add header pour avoir un fichier final pour chacun des 3 phènotypes 
cd 
cd Myriam_2020/data_myriam/chanel/file_1_results

cd relachement_files_1 
cd ../lentigine_files_1 
cd ../ride_files_1

echo -e "Marker\tP_value" > header
cat header lentigine_chr10_file1 lentigine_chr11_file1 lentigine_chr12_file1 lentigine_chr13_file1 lentigine_chr14_file1 lentigine_chr15_file1 lentigine_chr16_file1 lentigine_chr17_file1 lentigine_chr18_file1 lentigine_chr19_file1 lentigine_chr1_file1 lentigine_chr20_file1 lentigine_chr21_file1 lentigine_chr22_file1 lentigine_chr2_file1 lentigine_chr3_file1 lentigine_chr4_file1 lentigine_chr5_file1 lentigine_chr6_file1 lentigine_chr7_file1 lentigine_chr8_file1 lentigine_chr9_file1 > lentigine_chanel_all_file1
## add header 3694994

cd
cd Myriam_2020/data_myriam/chanel/file_1_results/ride_files_1
echo -e "Marker\tP_value" > header
cat header ride_chr10_file1 ride_chr16_file1 ride_chr21_file1 ride_chr6_file1 ride_chr11_file1  ride_chr17_file1 ride_chr22_file1 ride_chr7_file1 ride_chr12_file1 ride_chr18_file1  ride_chr2_file1 ride_chr8_file1 ride_chr13_file1 ride_chr19_file1 ride_chr3_file1 ride_chr9_file1 ride_chr14_file1 ride_chr1_file1 ride_chr4_file1 ride_chr15_file1  ride_chr20_file1 ride_chr5_file1 > rides_chanel_all_file1

cd 
cd Myriam_2020/data_myriam/chanel/file_1_results/relachement_files_1
echo -e "Marker\tP_value" > header
cat header relachement_chr10_file1 relachement_chr18_file1 relachement_chr4_file1 relachement_chr11_file1 relachement_chr19_file1 relachement_chr5_file1 relachement_chr12_file1 relachement_chr1_file1 relachement_chr6_file1 relachement_chr13_file1 relachement_chr20_file1 relachement_chr7_file1 relachement_chr14_file1 relachement_chr21_file1 relachement_chr8_file1 relachement_chr15_file1 relachement_chr22_file1 relachement_chr9_file1 relachement_chr16_file1 relachement_chr2_file1 relachement_chr17_file1 relachement_chr3_file1 > relachement_chanel_all_file1
#3694995

#PATHWAY CHANEL
cd 
cd Myriam_2020/data_myriam

	#ride kegg 
sed -i 's/ //g' chanel/file_1_results/ride_files_1/rides_chanel_all_file1
perl calculate_gsea.pl chanel/file_1_results/ride_files_1/rides_chanel_all_file1 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_chanel_Kegg --pvalue_flag --logfile log_rides_chanel_Kegg --weight 1 --large_es > Pathway_results/rides_chanel_Kegg
	#ride go
perl calculate_gsea.pl chanel/file_1_results/ride_files_1/rides_chanel_all_file1 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_chanel_go --pvalue_flag --logfile log_rides_chanel_go --weight 1 --large_es > Pathway_results/rides_chanel_go

	#lentigines kegg 
sed -i 's/ //g' chanel/file_1_results/lentigine_files_1/lentigine_chanel_all_file1
perl calculate_gsea.pl chanel/file_1_results/lentigine_files_1/lentigine_chanel_all_file1 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_chanel_Kegg --pvalue_flag --logfile log_lentigines_chanel_Kegg --weight 1 --large_es > lentigines_chanel_Kegg
	#lentigines go 
perl calculate_gsea.pl chanel/file_1_results/lentigine_files_1/lentigine_chanel_all_file1 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigine_chanel_go --pvalue_flag --logfile log_lentigines_chanel_go --weight 1 --large_es > Pathway_results/lentigines_chanel_go

	#relachement kegg 
sed -i 's/ //g' chanel/file_1_results/relachement_files_1/relachement_chanel_all_file1
perl calculate_gsea.pl chanel/file_1_results/relachement_files_1/relachement_chanel_all_file1 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_chanel_Kegg --pvalue_flag --logfile log_relachement_chanel_Kegg --weight 1 --large_es > Pathway_results/relachement_chanel_Kegg
	#relachement go 
perl calculate_gsea.pl  chanel/file_1_results/relachement_files_1/relachement_chanel_all_file1 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_chanel_go --pvalue_flag --logfile log_relachement_chanel_go --weight 1 --large_es > Pathway_results/relachement_chanel_go


################################################################SALIA 

	#Creer file 1 salia rs	pval pour les 3 phenos (rs existants mais on match quand meme par rapport à l'ld ici on peux ne prendre que le rs) 
cd 
cd Myriam_2020/data_myriam/SALIA_allemand

cd laxity
mkdir file_1_results
for i in {1..22}
	do awk 'NR==FNR{a[$26]=$26"\t"$6; next} ($3 in a){print a[$3]}' GWAS_Zlaxity_imputed_genotypes_chr${i}.txt ../../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > file_1_results/relachement_chr${i}_file1
	done
cd file_1_results
echo -e "Marker\tP_value" > header
cat header relachement_chr10_file1 relachement_chr18_file1 relachement_chr4_file1 relachement_chr11_file1 relachement_chr19_file1 relachement_chr5_file1 relachement_chr12_file1 relachement_chr1_file1 relachement_chr6_file1 relachement_chr13_file1 relachement_chr20_file1 relachement_chr7_file1 relachement_chr14_file1 relachement_chr21_file1 relachement_chr8_file1 relachement_chr15_file1 relachement_chr22_file1 relachement_chr9_file1 relachement_chr16_file1 relachement_chr2_file1 relachement_chr17_file1 relachement_chr3_file1 > relachement_SALIA_all_file1 # 3665150

cd ../../lentigines 
mkdir file_1_results
for i in {1..22}
	do awk 'NR==FNR{a[$26]=$26"\t"$6; next} ($3 in a){print a[$3]}' GWAS_Zpigment_imputed_genotypes_chr${i}.txt ../../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > file_1_results/lentigines_chr${i}_file1
	done
cd file_1_results
echo -e "Marker\tP_value" > header 
cat header lentigines_chr10_file1 lentigines_chr18_file1 lentigines_chr4_file1 lentigines_chr11_file1 lentigines_chr19_file1 lentigines_chr5_file1 lentigines_chr12_file1 lentigines_chr1_file1 lentigines_chr6_file1 lentigines_chr13_file1 lentigines_chr20_file1 lentigines_chr7_file1 lentigines_chr14_file1 lentigines_chr21_file1 lentigines_chr8_file1 lentigines_chr15_file1 lentigines_chr22_file1 lentigines_chr9_file1 lentigines_chr16_file1 lentigines_chr2_file1 lentigines_chr17_file1 lentigines_chr3_file1  > lentigines_SALIA_all_file1 #3665150

cd ../../wrinkles
mkdir file_1_results
for i in {1..22}
	do awk 'NR==FNR{a[$26]=$26"\t"$6; next} ($3 in a){print a[$3]}' GWAS_Zwrinkles_imputed_genotypes_chr${i}.txt ../../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > file_1_results/rides_chr${i}_file1
	done
cd file_1_results
echo -e "Marker\tP_value" > header 
cat header rides_chr10_file1 rides_chr16_file1 rides_chr21_file1 rides_chr6_file1 rides_chr11_file1 rides_chr17_file1 rides_chr22_file1 rides_chr7_file1 rides_chr12_file1 rides_chr18_file1 rides_chr2_file1 rides_chr8_file1 rides_chr13_file1 rides_chr19_file1 rides_chr3_file1 rides_chr9_file1 rides_chr14_file1 rides_chr1_file1 rides_chr4_file1 rides_chr15_file1 rides_chr20_file1 rides_chr5_file1 > rides_SALIA_all_file1 #3665150

cd 
cd Myriam_2020/data_myriam

##pathway salia
	# ride kegg
	perl calculate_gsea.pl SALIA_allemand/wrinkles/file_1_results/rides_SALIA_all_file1 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_salia_Kegg --pvalue_flag --logfile log_rides_salia_Kegg --weight 1 --large_es > Pathway_results/rides_salia_Kegg
	#ride go 
	perl calculate_gsea.pl SALIA_allemand/wrinkles/file_1_results/rides_SALIA_all_file1 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_salia_go --pvalue_flag --logfile log_rides_salia_go --weight 1 --large_es > Pathway_results/rides_salia_go

	#lentigines kegg 
	perl calculate_gsea.pl SALIA_allemand/lentigines/file_1_results/lentigines_SALIA_all_file1 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_salia_Kegg --pvalue_flag --logfile log_lentigines_salia_Kegg --weight 1 --large_es > Pathway_results/lentigine_salia_Kegg
	#lentigine go
	perl calculate_gsea.pl SALIA_allemand/lentigines/file_1_results/lentigines_SALIA_all_file1 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_salia_go --pvalue_flag --logfile log_lentigines_salia_go --weight 1 --large_es > Pathway_results/lentigines_salia_go
	
	#relachement kegg 
	perl calculate_gsea.pl SALIA_allemand/laxity/file_1_results/relachement_SALIA_all_file1 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_salia_Kegg --pvalue_flag --logfile log_relachement_salia_Kegg --weight 1 --large_es > Pathway_results/relachement_salia_Kegg
	#relachement go
perl calculate_gsea.pl SALIA_allemand/laxity/file_1_results/relachement_SALIA_all_file1 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_salia_go --pvalue_flag --logfile log_relachement_salia_go --weight 1 --large_es > Pathway_results/relachement_salia_go

################################################################################# TAIZHOU 

cd ../TAIZOU_chinois
	#merger les données 1000 genome eas car pour taizou les fichiers ne sont pas répartis par chr
cd /../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre
cat chr10_eas_pos_rs_propre chr13_eas_pos_rs_propre chr16_eas_pos_rs_propre chr19_eas_pos_rs_propre chr21_eas_pos_rs_propre chr3_eas_pos_rs_propre  chr6_eas_pos_rs_propre chr9_eas_pos_rs_propre chr11_eas_pos_rs_propre  chr14_eas_pos_rs_propre chr17_eas_pos_rs_propre chr1_eas_pos_rs_propre chr22_eas_pos_rs_propre chr4_eas_pos_rs_propre chr7_eas_pos_rs_propre chr12_eas_pos_rs_propre chr15_eas_pos_rs_propre chr18_eas_pos_rs_propre  chr20_eas_pos_rs_propre chr2_eas_pos_rs_propre chr5_eas_pos_rs_propre chr8_eas_pos_rs_propre > all_chr_eas_rs_propre 

cd ../../../TAIZOU_chinois
#zeydin mais bon c'est plus clair
awk ' (NR>=1) {print $1"\t"$21}' < TAIZHOU_20180909.txt > TAIZHOU_20180909_chpos_pval_relachement
awk ' (NR>=1) {print $1"\t"$17}' < TAIZHOU_20180909.txt > TAIZHOU_20180909_chpos_pval_rides
awk ' (NR>=1) {print $1"\t"$13}' < TAIZHOU_20180909.txt > TAIZHOU_20180909_chpos_pval_lentigines

#Extraire les snp qu'il faut 
awk 'NR==FNR{a[$1]=$2; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' TAIZHOU_20180909_chpos_pval_relachement ../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/all_chr_eas_rs_propre > relachement_all_TAIZOU
awk 'NR==FNR{a[$1]=$2; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' TAIZHOU_20180909_chpos_pval_rides ../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/all_chr_eas_rs_propre > rides_all_TAIZOU
awk 'NR==FNR{a[$1]=$2; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' TAIZHOU_20180909_chpos_pval_lentigines ../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/all_chr_eas_rs_propre > lentigines_all_TAIZOU

## Pathway TAIZHOU
cd 
cd Myriam_2020/data_myriam
	#ride kegg 
	perl calculate_gsea.pl TAIZOU_chinois/rides_all_TAIZOU c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_TAIZOU_Kegg --pvalue_flag --logfile log_rides_TAIZOU_Kegg --weight 1 --large_es > Pathway_results/rides_TAIZOU_Kegg
	#ride go
	perl calculate_gsea.pl TAIZOU_chinois/rides_all_TAIZOU c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_TAIZOU_go --pvalue_flag --logfile log_rides_TAIZOU_go --weight 1 --large_es > Pathway_results/rides_TAIZOU_go
	
	#lentigines kegg 
	perl calculate_gsea.pl TAIZOU_chinois/lentigines_all_TAIZOU c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_TAIZOU_Kegg --pvalue_flag --logfile log_lentigines_TAIZOU_Kegg --weight 1 --large_es > Pathway_results/lentigines_TAIZOU_Kegg
	#lentigines go 
	perl calculate_gsea.pl TAIZOU_chinois/lentigines_all_TAIZOU c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_TAIZOU_go --pvalue_flag --logfile log_lentigines_TAIZOU_go --weight 1 --large_es > Pathway_results/lentigines_TAIZOU_go

	#relachement kegg 
	perl calculate_gsea.pl TAIZOU_chinois/relachement_all_TAIZOU c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_TAIZOU_Kegg --pvalue_flag --logfile log_relachement_TAIZOU_Kegg --weight 1 --large_es > Pathway_results/relachement_TAIZOU_Kegg
	#relachement go 
	perl calculate_gsea.pl TAIZOU_chinois/relachement_all_TAIZOU c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_TAIZOU_go --pvalue_flag --logfile log_relachement_TAIZOU_go --weight 1 --large_es > Pathway_results/relachement_TAIZOU_go


##données meta analyse salia+chanel 
cd salia-chanel-meta
#Préparation des fichiers 
	#ride 
mkdir rides_results
for i in {1..22}
	do
	zcat meta_CERIES_IUF_chr${i}_ride_z_random.gz|awk '(NR>=1){snp=$1;pvalue=$9} {print snp,"\t"pvalue}' > meta_CERIES_IUF_chr${i}_ride_z_random_propre
	awk 'NR==FNR{a[$1]=$2"\t"$1; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' meta_CERIES_IUF_chr${i}_ride_z_random_propre ../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > rides_results/chr${i}_eur_propre
	rm meta_CERIES_IUF_chr${i}_ride_z_random_propre
	done

cd rides_results 
echo -e "Marker\tP_value\ttej" > header 
cat header chr10_eur_propre chr16_eur_propre chr21_eur_propre chr6_eur_propre chr11_eur_propre chr17_eur_propre chr22_eur_propre chr7_eur_propre chr12_eur_propre chr18_eur_propre chr2_eur_propre chr8_eur_propre chr13_eur_propre chr19_eur_propre chr3_eur_propre chr9_eur_propre chr14_eur_propre chr1_eur_propre chr4_eur_propre chr15_eur_propre chr20_eur_propre chr5_eur_propre > ride_all_meta2_post

awk < ride_all_meta2_post '(NR>=1){snp=$1;pvalue=$2} {print snp,"\t"pvalue}' > ride_all_meta2 #3440939
rm ride_all_meta2_post

	#lentigines
cd .. 
mkdir lentigines_results
for i in {1..22}
	do
	zcat meta_CERIES_IUF_chr${i}_lentigine_z_random.gz|awk '(NR>=1){snp=$1;pvalue=$9} {print snp,"\t"pvalue}' > meta_CERIES_IUF_chr${i}_lentigine_z_random_propre
	awk 'NR==FNR{a[$1]=$2"\t"$1; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' meta_CERIES_IUF_chr${i}_lentigine_z_random_propre ../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > lentigines_results/chr${i}_eur_propre
	rm meta_CERIES_IUF_chr${i}_lentigine_z_random_propre
	done

cd lentigines_results 
echo -e "Marker\tP_value\ttej" > header 
cat header chr10_eur_propre chr16_eur_propre chr21_eur_propre chr6_eur_propre chr11_eur_propre chr17_eur_propre chr22_eur_propre chr7_eur_propre chr12_eur_propre chr18_eur_propre chr2_eur_propre chr8_eur_propre chr13_eur_propre chr19_eur_propre chr3_eur_propre chr9_eur_propre chr14_eur_propre chr1_eur_propre chr4_eur_propre chr15_eur_propre chr20_eur_propre chr5_eur_propre > lentigines_all_meta2_post
awk < lentigines_all_meta2_post '(NR>=1){snp=$1;pvalue=$2} {print snp,"\t"pvalue}' > lentigines_all_meta2 #3440939
rm lentigines_all_meta2_post

	#relachement
cd .. 
mkdir relachement_results
for i in {1..22}
	do
	zcat meta_CERIES_IUF_chr${i}_relachement_z_random.gz|awk '(NR>=1){snp=$1;pvalue=$9} {print snp,"\t"pvalue}' > meta_CERIES_IUF_chr${i}_relachement_z_random_propre
	awk 'NR==FNR{a[$1]=$2"\t"$1; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' meta_CERIES_IUF_chr${i}_relachement_z_random_propre ../1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre > relachement_results/chr${i}_eur_propre
	rm meta_CERIES_IUF_chr${i}_relachement_z_random_propre
	done

cd relachement_results 
echo -e "Marker\tP_value\ttej" > header 
cat header chr10_eur_propre chr16_eur_propre chr21_eur_propre chr6_eur_propre chr11_eur_propre chr17_eur_propre chr22_eur_propre chr7_eur_propre chr12_eur_propre chr18_eur_propre chr2_eur_propre chr8_eur_propre chr13_eur_propre chr19_eur_propre chr3_eur_propre chr9_eur_propre chr14_eur_propre chr1_eur_propre chr4_eur_propre chr15_eur_propre chr20_eur_propre chr5_eur_propre > relachement_all_meta2_post

awk < relachement_all_meta2_post '(NR>=1){snp=$1;pvalue=$2} {print snp,"\t"pvalue}' > relachement_all_meta2 #3440939
rm relachement_all_meta2_post

##pathway meta analyse de salia + chanel 
cd 
cd Myriam_2020/data_myriam
##rides
sed -i 's/ //g' salia-chanel-meta/rides_results/ride_all_meta2
		#ride kegg 
perl calculate_gsea.pl salia-chanel-meta/rides_results/ride_all_meta2 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_salia-chanel_Kegg --pvalue_flag --logfile log_rides_salia-chanel_Kegg --weight 1 --large_es > Pathway_results/rides_salia-chanel_Kegg
		#ride go
perl calculate_gsea.pl salia-chanel-meta/rides_results/ride_all_meta2 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_salia-chanel_go --pvalue_flag --logfile log_rides_salia-chanel_go --weight 1 --large_es > Pathway_results/rides_salia-chanel_go

##lentigines
sed -i 's/ //g' salia-chanel-meta/lentigines_results/lentigines_all_meta2
		#lentigines kegg 
perl calculate_gsea.pl salia-chanel-meta/lentigines_results/lentigines_all_meta2 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_salia-chanel_Kegg --pvalue_flag --logfile log_lentigines_salia-chanel_Kegg --weight 1 --large_es > Pathway_results/lentigines_salia-chanel_Kegg
		#lentigines go 

perl calculate_gsea.pl salia-chanel-meta/lentigines_results/lentigines_all_meta2 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_salia-chanel_go --pvalue_flag --logfile log_lentigines_salia-chanel_go --weight 1 --large_es > Pathway_results/lentigines_salia-chanel_go

##Relachement 
sed -i 's/ //g' salia-chanel-meta/relachement_results/relachement_all_meta2
		#relachement kegg 
perl calculate_gsea.pl salia-chanel-meta/relachement_results/relachement_all_meta2 c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_salia-chanel_Kegg --pvalue_flag --logfile log_relachement_salia-chanel_Kegg --weight 1 --large_es > Pathway_results/relachement_salia-chanel_Kegg
		#relachement go 
perl calculate_gsea.pl salia-chanel-meta/relachement_results/relachement_all_meta2 c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_salia-chanel_go --pvalue_flag --logfile log_relachement_salia-chanel_go --weight 1 --large_es > Pathway_results/relachement_salia-chanel_go



