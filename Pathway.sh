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


cd /home/mrahmouni/Myriam_2020/data_myriam_test/chanel

#Déplacer les .md5 dans un dossier aprt 
ls *.md5 > list-md5 
mkdir md5-files 
for i in `cat list-md5`
	do
	mv /home/mrahmouni/Myriam_2020/data_myriam_test/chanel/${i} /home/mrahmouni/Myriam_2020/data_myriam_test/chanel/md5-files/${i}
	done

#Input file 1: GWAS single-SNP results (after QC) :  rsid	Chi-2 or p-val

#Creer file1 chanel ch:po
ls *.relachement_z.snptest.gz > list-file-relachement
ls *.lentigine_z.snptest.gz > list-file-lentigine
ls *.ride_z.snptest.gz > list-file-ride

mkdir relachement_chr-pos-pval-files
mkdir lentigine_chr-pos-pval-files
mkdir ride_chr-pos-pval-files

for j in list-file-relachement list-file-lentigine list-file-ride
	do 
	for i in `cat ${j}`
		do
		gzip -d ${i}
		#info score = 0.3
		awk < ${i::-3} '(NR==1) print {rsid	pvalue} (NR>=1){snp=$1;pvalue=$21}(NR>=1){print snp,	pvalue}' > /home/mrahmouni/Myriam_2020/data_myriam_test/chanel/${j:10}_chr-pos-pval-files/${i::-11}-chr-pos
		done
	done

#Creer file1 channel filtre ld 0.8 

	#retrouver les rs de mes snp 
	#awk recuperer la liste 
	# tr changer les ':' en \t 
	#script python qui à partir du .ped de base retrouve les rs 
	#recuperer la liste des RS 
	#http://grch37.ensembl.org/Homo_sapiens/Tools/LD
	#python recrer le fichier avec uniquement les snp sans ld *




####SALIA 
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


### préparation exemple 1000genome

sed '/##/d' < exemple_1000_genomes > exemple_1000_genomes_sans#
for i in {1..22}
	do echo ${i} 
	awk < exemple_1000_genomes_sans# ' (NR==1){print "CHR POS RS"} (NR>1){chr=$1;pos=$2;rs=$3}(NR>1 && chr=='${i}'){print chr,	pos,	rs}'> resultat-${i}
	done

===> 



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
