#!/bin/bash
#mv lentigine_chanel_all_file1 cnam_lentigine
#mv rides_chanel_all_file1 cnam_ride
#mv relachement_chanel_all_file1 cnam_relachement

#mv lentigines_all_meta2 saliametacnam_lentigine
#mv relachement_all_meta2 saliametacnam_relachement
#mv ride_all_meta2 saliametacnam_ride

#mv relachement_all_TAIZOU taizou_lentigine
#mv lentigines_all_TAIZOU taizou_relachement
#mv rides_all_TAIZOU taizou_ride

#mv rides_SALIA_all_file1 salia_ride
#mv lentigines_SALIA_all_file1 salia_lentigines
#mv relachement_SALIA_all_file1 salia_relachement

#mkdir taizou_results_dist10_inf03
#mkdir cnam_results_dist10_inf03
#mkdir salia_results_dist10_inf03
#mkdir saliametacnam_results_dist10_inf03
#mkdir BACH

DIR="/media/sigridNFS2/myriam/GSEA/"
GSEA_PROG="/media/sigridNFS2/myriam/GSEA/calculate_gsea.pl"
BATCH="/media/sigridNFS2/myriam/GSEA/BATCH/"

PHENOS="lentigine relachement ride"
POPS="cnam salia saliametacnam taizou"
REF_PATHS="c2.cp.kegg.v7.0.symbols.gmt c5.bp.v7.0.symbols.gmt"

for POP in ${POPS} 
	do
	for PHENO in ${PHENOS}
		do
		for REF_PATH in ${REF_PATHS}
			do
			BATCHFILE="${BATCH}commande_GSEA_dist10.${POP}"
			echo "perl ${GSEA_PROG} ${DIR}${POP}_${PHENO} ${DIR}${REF_PATH} -map ${DIR}snpgenemap.hg19 -c 10000 --distance 10 --setstatfile ${DIR}stat_${POP}_${PHENO}_${REF_PATH::-17} --pvalue_flag --logfile ${DIR}log_${POP}_${PHENO}_${REF_PATH::-17} --weight 1 > ${DIR}${POP}_results_dist10_inf03/GSEA_${POP}_${PHENO}_${REF_PATH::-17}" > ${BATCHFILE}
			qsub -q ALL ${BATCHFILE}
		done	
	done
done 

#/usr/local/maui/bin/showq
#qdel
