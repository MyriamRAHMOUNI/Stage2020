#serveur ssh sigrid@bioinfoclust01.cnam.fr
# cd media/bionasNFS2/GWAS/CERIES/4_IMPUTE
#bioinfoclust01.cnam.fr 
#QC (prunning/500kb de tout genes)  Je suis sur HG37

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

#Creer file2 channel filtre ld 0.8 
#http://grch37.ensembl.org/Homo_sapiens/Tools/LD




#Input file 2: SNP-gene mapping : SNP	gene	dist2gene

#Changer les rs en ch:pos 
	#fichier mapping existant en rs  en suppose qu'il est dans le bon forma 

#tr "\t" " " < file.map > file-inter
awk < file.map '(NR==1) print {SNP	gene	dist2gene} (NR>=1){snp=$1;gene=$2;dist2gene=$3}(NR>=1){print snp}' > /home/mrahmouni/Myriam_2020/data_myriam_test/chanel/liste-rs-map

#liste-rs-map ====================================> http://grch37.ensembl.org/biomart/martview/f8b7b70c58bc97e105000a02a5bb1f41 
# https://www.internationalgenome.org/faq/can-i-find-genomic-position-list-dbsnp-rs-numbers-0/



#Input file 3: gene-pathway : Gene Ontology (level 4), KEGG, Biocarta


#=========> Calculate_gsea.pl
