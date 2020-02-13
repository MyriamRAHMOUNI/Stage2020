#serveur ssh sigrid@bioinfoclust01.cnam.fr
# cd media/bionasNFS2/GWAS/CERIES/4_IMPUTE
#bioinfoclust01.cnam.fr 
#QC (prunning/500kb de tout genes)  Je suis sur HG37
#code ? zagfan1402


# donnée meta analyse? il viennent 
# je peux utiliser les rs de chanel pour retrouver ceux des chinois ? non on va plutot utiliser 1000 genome chr par chr sans ouvrir les .gz 
# ou sont les données des chinois ? juste devant toi ;)
# fichier mapping celui du toto
# c'est ok pour le serveur ld ? population oui avant je vais garder que les snp dans le fichier de mapping 
# base de donnée pathway du tuto? non celle du lien 
# seuil ld 0.8 ? oui 
#zcat file.gz   to view the contenent


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
mkdir results_preparation_1000g
mkdir fichiers_ch_pos_rs_propre
#eas
for i in {1..22}
	do
	sed '/##/d' | awk < results_eas_1000g/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_eas.recode.vcf ' (NR==1){print "CHR POS RS"} (NR>1){chr=$1;pos=$2;rs=$3}(NR>1){print chr, pos, rs}' > results_preparation_1000g/chr${i}_prep
	awk 'NR==FNR{a[$1]=$1; next} ($3 in a){print $1, $2, $3}' results_eas_1000g/res_prun/file1.in results_preparation_1000g/chr${i}_prep > fichiers_ch_pos_rs_propre/chr${i}_pos_rs_propre
	#=> listes rs ch pos

sed '/##/d' | awk 'NR==FNR{a[$1]=$1; next} ($3 in a){print $1, $2, $3}' ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_eas.recode.vcf ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_eas.recode.vcf_propre_eas.prune > voila

	done

###################Créer les files1 chr par chr de chanel############################################ 
cd Myriam_2020/Stage2020/data_myriam_copie/chanel
#Déplacer les .md5 dans un dossier aprt 
mkdir md5-files 
ls *.md5 > list-md5.txt
for i in `cat list-md5.txt`
	do
	mv /home/mrahmouni/Myriam_2020/Stage2020/data_myriam_copie/chanel/${i} /home/mrahmouni/Myriam_2020/Stage2020/data_myriam_copie/chanel/md5-files
	done

ls *.relachement_z.snptest.gz > list-file-relachement
ls *.lentigine_z.snptest.gz > list-file-lentigine
ls *.ride_z.snptest.gz > list-file-ride

mkdir file_1_results 
cd file_1_results
mkdir relachement_chr_pos_files_1 lentigine_chr_pos_files_1 ride_chr_pos_files_1
cd .. 

for j in list-file-relachement list-file-lentigine list-file-ride
	do 
	for i in `cat ${j}`
		do
		#gzip -d ${i}
		sed '/#/d' < ${i::-3} > ${i::-3}_sans#
		#info score = 0.3 never forget 
		awk < ${i::-3}_sans# '(NR==1) {print "CHR POS pvalue"} (NR>1){chr=$3;pos=$4;pvalue=$21}(NR>1){print chr, pos, pvalue}' > /home/mrahmouni/Myriam_2020/Stage2020/data_myriam_copie/chanel/file_1_results/${j:10}_chr_pos_files_1/${i::-11}_chr_pos_files_1
		rm ${i::-3}_sans#
		done
	done

#Retrouver les rs de mes SNP 
cd file_1_results
mkdir lentigine_rs_pval_files_1
mkdir relachement_rs_pval_files_1
mkdir ride_rs_pval_files_1
cd ..
cd ..
cd ..
for i in {1..22}
	do python3 find-rs.py data_myriam_copie/chanel/file_1_results/lentigine_chr_pos_files_1/CERIES.chr${i}.mode1.lentigine_z_chr_pos_files_1 results_1000g/chr${i}_pos_rs data_myriam_copie/chanel/file_1_results/lentigine_rs_pval_files_1
	done


#garder que ce qui'il y a dans le fichier map
#meme principe tu va voir 





	#Calcule ld http://grch37.ensembl.org/Homo_sapiens/Tools/LD à la main 
	#python recrer le fichier avec uniquement les snp sans ld *

###################Créer les files1 chr par chr de SALIA (rs existants) 

	#Creer file 1 salia rs	pval pour les 3 phenos 
cd 
cd home/mrahmouni/Myriam_2020/Stage2020/data_myriam/SALIA_allemand/laxity
ls *.txt > list-laxity
for i in `cat list-laxity` 
	do  awk < ${i} '(NR>1){snp=$1;pvalue=$26}(NR>=1){print snp,	pvalue}'> ${i::-4}-file1

cd ..
cd lentigines 
ls *.txt > list-lentigines
for i in `cat list-lentigines` 
	do  awk < ${i} '(NR>1){snp=$1;pvalue=$26}(NR>=1){print snp,	pvalue}'> ${i::-4}-file1

cd ..
cd wrinkles 
ls *.txt > list-wrinkles
for i in `cat list-wrinkles` 
	do  awk < ${i} '(NR>1){snp=$1;pvalue=$26}(NR>=1){print snp,	pvalue}'> ${i::-4}-file1

ls *-file1 > list

### TAIZHOU 

cd .. 
cd TAIZOU_chinois

awk < header ' (NR==1){print CP	chr	pos	pval} (NR >= 1) {CP=$1; chr=$2; pos=$3; Laxity_Pvalue=$20} {print CP	chr	pos	Laxity_Pvalue} '> Laxity

awk < header ' (NR==1){print CP	chr	pos	pval} (NR >= 1) {CP=$1; chr=$2; pos=$3; Wrinkle_Pvalue=$16} {print CP	chr	pos	Wrinkle_Pvalue} '> Wrinkle

awk < header ' (NR==1){print CP	chr	pos	pval} (NR >= 1) {CP=$1; chr=$2; pos=$3; PigmentedSpots_Pvalue=$20} {print CP	chr	pos	PigmentedSpots_Pvalue} '> lentigines



### Meta-analyses 






#Input file 2: SNP-gene mapping : SNP	gene	dist2gene

#Changer les rs en ch:pos 
	#fichier mapping existant en rs  en suppose qu'il est dans le bon forma 

#tr "\t" " " < file.map > file-inter
awk < file.map '(NR==1) print {SNP	gene	dist2gene} (NR>=1){snp=$1;gene=$2;dist2gene=$3}(NR>=1){print snp}' > /home/mrahmouni/Myriam_2020/data_myriam_test/chanel/liste-rs-map

#liste-rs-map ====================================> http://grch37.ensembl.org/biomart/martview/f8b7b70c58bc97e105000a02a5bb1f41 
# https://www.internationalgenome.org/faq/can-i-find-genomic-position-list-dbsnp-rs-numbers-0/



#Input file 3: gene-pathway : Gene Ontology (level 4), KEGG, Biocarta

	#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2
	# KEGG / go biological process 

#=========> Calculate_gsea.pl















#vcf-merge liste | gzip -c > multi-échantillon.vcf.gz
