#!/bin/bash
#

Red="$(tput setaf 1)"
Green="$(tput setaf 2)"
Bold=$(tput bold)
reset=`tput sgr0` # turns off all atribute

for file in *_gene.fasta
do
	echo ""
    echo ""
    echo "${Red}${Bold}Blastn ${reset}: "${file}"" 
    base=$(basename ${file} .fasta)
    blastn -query "${file}" -task blastn -db concatenated_genomes_db/concatenated_genomes_db -outfmt 6 -out ${base}_blastn_hits.tsv -evalue 1e-10 -num_threads 18 
	echo "${Green}${Bold}Done ${reset}: "${file}""
    echo ""
    echo "${Red}${Bold}Making bed file ${reset}: "${file}"" 
    awk '{print $2,$9,$10,""$1"_"NR}' ${base}_blastn_hits.tsv > ${base}_blastn_hits.bed
    echo "${Green}${Bold}Done ${reset}: "${file}""
    echo ""
    echo "${Red}${Bold}Sorting bed file ${reset}: "${file}"" 
    cat  ${base}_blastn_hits.bed | awk -v OFS='\t' '{ if ($2 < $3) {print $1,$2,$3,$4} else {print $1,$3,$2,$4} }' >  ${base}_blastn_hits_sorted.bed
    echo "${Green}${Bold}Done ${reset}: "${file}""
    echo ""
    echo "${Red}${Bold}Extracting blastn hit sequences ${reset}: "${file}"" 
    bedtools getfasta -fi concatenated_genomes.fasta -bed ${base}_blastn_hits_sorted.bed -fullHeader > ${base}_blastn_hits.fasta
    echo "${Green}${Bold}Done ${reset}: "${file}""
    echo ""
    echo "${Red}${Bold}Naming headers with the corresponding gene and genome names ${reset}: "${file}""
    awk -v var=$(basename $file .fasta) '{ gsub(/contig/, var "_contig") } 1' ${base}_blastn_hits.fasta > ${base}_blastn_hits_seqs.fasta
    echo "${Green}${Bold}Done ${reset}: "${file}""
    echo ""
    echo "${Red}${Bold}Removing the coordinates from the headers ${reset}: "${file}""
    sed -r '/^>/s/:[0-9]+-[0-9]+//' ${base}_blastn_hits_seqs.fasta > ${base}_blastn_hits_seqs_together.fasta
    echo "${Green}${Bold}Done ${reset}: "${file}""
    echo ""
    echo "${Red}${Bold}Splitting the blastn hit sequences into individual fasta files ${reset}: "${file}""
    while read line
    do
        if [[ ${line:0:1} == '>' ]] # files starting with '>'
        then
            outfile=$(echo "${line#>}" | cut -d ':' -f1)${baseName}.fasta # '${line#>}' is the enter heading, 
                                                                       # 'cut' separates parts based on space
                                                                       # '-f1' picks up the first part
                                                           
            echo $line > "$outfile"
        else
            echo $line >> "$outfile"
        fi
    done < ${base}_blastn_hits_seqs_together.fasta
    echo "${Green}${Bold}Done ${reset}: "${file}""
    rm -r *.bed *.tsv *_blastn_hits_seqs.fasta *_gene_blastn_hits.fasta
    echo ""
    echo "${Green}${Bold}Find the blastn hits as individual fasta files. However, _together.fasta contains all of them ${reset}"
    echo ""
    echo ""
done

