#!/bin/bash
#

Red="$(tput setaf 1)"
Green="$(tput setaf 2)"
Bold=$(tput bold)
reset=`tput sgr0` # turns off all atribute

for file in *_gene.fasta
do
	echo "${Red}${Bold}Blastn ${reset}: "${file}"" 
    base=$(basename ${file} .fasta)
    blastn -query "${file}" -task blastn -db pecto_genomes_db/pecto_genomes_db -outfmt 6 -out ${base}_blastn_hits.tsv -evalue 1e-10 -num_threads 18 
	echo "${Green}${Bold}Blastn done ${reset}: "${file}""
    echo "${Red}${Bold}Making bed file ${reset}: "${file}"" 
    awk '{print $2,$9,$10,"Rice_"$1"_"NR}' ${base}_blastn_hits.tsv > ${base}_blastn_hits.bed
    echo "${Green}${Bold}Bed file made ${reset}: "${file}""
    echo "${Red}${Bold}Sorting bed file ${reset}: "${file}"" 
    cat  ${base}_blastn_hits.bed | awk -v OFS='\t' '{ if ($2 < $3) {print $1,$2,$3,$4} else {print $1,$3,$2,$4} }' >  ${base}_blastn_hits_sorted.bed
    echo "${Green}${Bold}Bed file sorted ${reset}: "${file}""
    echo "${Red}${Bold}Extracting blastn hit sequences ${reset}: "${file}"" 
    bedtools getfasta -fi pecto_genomes.fasta -bed ${base}_blastn_hits_sorted.bed > ${base}_blastn_hits_seqs.fasta
    echo "${Green}${Bold}Extracted blastn hit sequences ${reset}: "${file}""
    echo ""
done

