#serveur ssh sigrid@bioinfoclust01.cnam.fr
# cd media/bionasNFS2/GWAS/CERIES/4_IMPUTE
#bioinfoclust01.cnam.fr 
#QC (prunning/500kb de tout genes)  Je suis sur HG37
#code ? zagfan1402


# donnée meta analyse?  
# je peux utiliser les rs de chanel pour retrouver ceux des chinois ? non on va plutot utiliser 1000 genome chr par chr sans ouvrir les .gz 
# ou sont les données des chinois ? juste devant toi ;)
# fichier mapping celui du toto
# c'est ok pour le serveur ld ? population oui avant je vais garder que les snp dans le fichier de mapping 
# base de donnée pathway du tuto? non celle du lien 
# seuil ld 0.8 ? oui 
#zcat file.gz   to view the contenent
#france genomique 


#selection des logs 

#preparer la liste des snp à garder dans l'analyse 
cd Myriam_2020/data_myriam/
awk < snpgenemap.hg19 '(NR>1){rs=$1} (NR>1){print rs}' > 1000genomes/list_extract_snp
cd 1000genomes
#EUR.txt EAS.txt
awk < integrated_call_samples_v2.20130502.ALL.ped '(NR>1){id=$2;pop=$7;related=$8} (NR>=1) && (related=="unrel") && (pop=="CEU"||pop=="TSI"||pop=="FIN"||pop=="GBR"||pop=="IBS") {print id}' > EUR.txt #317 individus

awk < integrated_call_samples_v2.20130502.ALL.ped '(NR>1){id=$2;pop=$7;related=$8} (NR>=1) && (related=="unrel") && (pop=="CHB"||pop=="JPT"||pop=="CHS"||pop=="CDX"||pop=="KHV") {print id}' > EAS.txt #437 individus 

#================> 56459228 list_extract_snp

ls *.vcf.gz > list-files-1000g

mkdir results_snp_1000g
mkdir results_eur_1000g
mkdir results_eas_1000g
for i in `cat  list-files-1000g`
	do
	#extraire les snp dans le map de 1000g
	vcftools --gzvcf ${i} --snps list_extract_snp --recode --out results_snp_1000g/${i::-7}_snps
	#extraire les population europenne avec individus non relié
	vcftools --vcf results_snp_1000g/${i::-7}_snps.recode.vcf --keep EUR.txt --recode --out results_eur_1000g/${i::-7}_eur
	#extraire les population eas avec individus non relié
	vcftools --vcf results_snp_1000g/${i::-7}_snps.recode.vcf --keep EAS.txt --recode --out results_eas_1000g/${i}_eas
	rm results_snp_1000g/${i::-7}_snps
	done

#calcule ld garder que les snps en ld < 0.8
cd results_eas_1000g
mkdir res_prun_eas
ls *.vcf > list-eas-1000g
for i in `cat  list-eas-1000g`
	do
	./plink2 --vcf ${i} --indep-pairwise 50 5 0.8 --out res_prun_eas/${i}_propre_eas
	done
cd ..
cd results_eur_1000g 
ls *.vcf > list-eur-1000g
mkdir res_prun_eur
for i in `cat  list-eur-1000g`
	do
	./plink2 --vcf ${i} --indep-pairwise 50 5 0.8 --out res_prun_eur/${i}_propre_eur
	done
#======================> on prend à chaque fois les prune.in apres c'est les snps qui existent sans le fchier map avec ld < 0.8

# Retrouver les ch po de ces liste rs 
#eas
cd .. 
cd results_eas_1000g
mkdir san_##
mkdir fichiers_ch_pos_rs_propre
for i in {1..22}
	do
	sed '/##/d' ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_eas.recode.vcf > san_##/ALL.chr${i}_EAS_san_##
	awk 'NR==FNR{a[$3]=$1" "$2" "$3; next} ($1 in a){print a[$1], a[$2], a[$3]}' san_##/ALL.chr${i}_EAS_san_## res_prun_eas/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_eas.recode.vcf_propre_eas.prune.in > fichiers_ch_pos_rs_propre/chr${i}_eas_pos_rs_propre
	done
	#=> listes rs ch pos eas
#eur
cd ..
cd results_eur_1000g 
mkdir san_##
mkdir fichiers_ch_pos_rs_propre
for i in {1..22}
	do
	sed '/##/d' ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_eur.recode.vcf > san_##/ALL.chr${i}_EUR_san_##
	awk 'NR==FNR{a[$3]=$1" "$2" "$3; next} ($1 in a){print a[$1], a[$2], a[$3]}' san_##/ALL.chr${i}_EUR_san_## res_prun_eur/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_eur.recode.vcf_propre_eur.prune.in > fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre
	done
	#=> listes rs ch pos eas

#donc la j'ai les rs ch pos des sn que je vais prendre chez chanel et salia : /home/mrahmouni/Myriam_2020/data_myriam/1000genomes/results_eur_1000g/fichiers_ch_pos_rs_propre/chr${i}_eur_pos_rs_propre et chez les chinois : /home/mrahmouni/Myriam_2020/data_myriam/1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/chr${i}_eas_pos_rs_propre

###################Créer les files1 chr par chr de chanel############################################ 
cd Myriam_2020/data_myriam/chanel
#Déplacer les .md5 dans un dossier aprt 
mkdir md5-files 
ls *.md5 > list-md5.txt
for i in `cat list-md5.txt`
	do
	mv /home/mrahmouni/Myriam_2020/Stage2020/data_myriam_copie/chanel/${i} /home/mrahmouni/Myriam_2020/Stage2020/data_myriam_copie/chanel/md5-files
	done
######
	#info score = 0.3 never forget 
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
# il faut merger et add header
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

#######################PATHWAY CHANEL
cd 
cd Myriam_2020/data_myriam

#lentigines lancé ttlh 

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


###################Créer les files1 chr par chr de SALIA (rs existants mais on match quand meme par rapport à l'ld ici on peux ne prendre que le rs) 

	#Creer file 1 salia rs	pval pour les 3 phenos 
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

### TAIZHOU 

cd .. 
cd TAIZOU_chinois

	#merger les données 1000 genome eas 
cd /../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre
cat chr10_eas_pos_rs_propre chr13_eas_pos_rs_propre chr16_eas_pos_rs_propre  chr19_eas_pos_rs_propre chr21_eas_pos_rs_propre chr3_eas_pos_rs_propre  chr6_eas_pos_rs_propre chr9_eas_pos_rs_propre chr11_eas_pos_rs_propre  chr14_eas_pos_rs_propre chr17_eas_pos_rs_propre chr1_eas_pos_rs_propre   chr22_eas_pos_rs_propre chr4_eas_pos_rs_propre chr7_eas_pos_rs_propre chr12_eas_pos_rs_propre chr15_eas_pos_rs_propre chr18_eas_pos_rs_propre  chr20_eas_pos_rs_propre chr2_eas_pos_rs_propre chr5_eas_pos_rs_propre  chr8_eas_pos_rs_propre > all_chr_eas_rs_propre 

cd ../../../TAIZOU_chinois
#zeydin mais bon 
awk ' (NR>=1) {print $1"\t"$21}' <TAIZHOU_20180909.txt > TAIZHOU_20180909_chpos_pval_relachement
awk ' (NR>=1) {print $1"\t"$17}' <TAIZHOU_20180909.txt > TAIZHOU_20180909_chpos_pval_rides
awk ' (NR>=1) {print $1"\t"$13}' <TAIZHOU_20180909.txt > TAIZHOU_20180909_chpos_pval_lentigines

#haaya tawa ess7i7


awk 'NR==FNR{a[$1]=$2; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' TAIZHOU_20180909_chpos_pval_relachement ../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/all_chr_eas_rs_propre > relachement_all_TAIZOU

awk 'NR==FNR{a[$1]=$2; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' TAIZHOU_20180909_chpos_pval_rides ../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/all_chr_eas_rs_propre > rides_all_TAIZOU

awk 'NR==FNR{a[$1]=$2; next} ($1":"$2 in a){print $3"\t"a[$1":"$2]}' TAIZHOU_20180909_chpos_pval_lentigines ../1000genomes/results_eas_1000g/fichiers_ch_pos_rs_propre/all_chr_eas_rs_propre > lentigines_all_TAIZOU



### Pathway TAIZHOU

cd 
cd Myriam_2020/data_myriam

##rides
		#ride kegg 

perl calculate_gsea.pl TAIZOU_chinois/rides_all_TAIZOU c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_TAIZOU_Kegg --pvalue_flag --logfile log_rides_TAIZOU_Kegg --weight 1 --large_es > Pathway_results/rides_TAIZOU_Kegg

		#ride go
perl calculate_gsea.pl TAIZOU_chinois/rides_all_TAIZOU c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_rides_TAIZOU_go --pvalue_flag --logfile log_rides_TAIZOU_go --weight 1 --large_es > Pathway_results/rides_TAIZOU_go

##lentigines
		#lentigines kegg 
perl calculate_gsea.pl TAIZOU_chinois/lentigines_all_TAIZOU c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_TAIZOU_Kegg --pvalue_flag --logfile log_lentigines_TAIZOU_Kegg --weight 1 --large_es > Pathway_results/lentigines_TAIZOU_Kegg
		#lentigines go 

perl calculate_gsea.pl TAIZOU_chinois/lentigines_all_TAIZOU c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_lentigines_TAIZOU_go --pvalue_flag --logfile log_lentigines_TAIZOU_go --weight 1 --large_es > Pathway_results/lentigines_TAIZOU_go

##Relachement 

		#relachement kegg 
perl calculate_gsea.pl TAIZOU_chinois/relachement_all_TAIZOU c2.cp.kegg.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_TAIZOU_Kegg --pvalue_flag --logfile log_relachement_TAIZOU_Kegg --weight 1 --large_es > Pathway_results/relachement_TAIZOU_Kegg

		#relachement go 

perl calculate_gsea.pl TAIZOU_chinois/relachement_all_TAIZOU c5.bp.v7.0.symbols.gmt -map snpgenemap.hg19 -c 10000 --setstatfile stat_relachement_TAIZOU_go --pvalue_flag --logfile log_relachement_TAIZOU_go --weight 1 --large_es > Pathway_results/relachement_TAIZOU_go






