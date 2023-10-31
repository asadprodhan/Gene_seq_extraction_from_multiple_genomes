# **How do I automate extracting multiple gene sequences from multiple genomes?** <br />


### **Asad Prodhan PhD** 

**https://asadprodhan.github.io/**

<br />


## **Tools**

- ### Create a conda environment as follows


```
conda create -n gene_seq_extraction
```


- ### Activate your conda environment


```
conda activate gene_seq_extraction
```


- ### Install blast


```
conda install -c bioconda blast
```


- ### Install bedtools


```
conda install -c bioconda bedtools
```


## **Steps**


- ### Concatenate all the genomes into a single fasta file


```
cat Pectbacterium* > pecto_genomes.fasta
```


- ### Convert pecto_genomes.fasta into a blast database


```
makeblastdb -in pecto_genomes.fasta -out pecto_genomes_db -dbtype 'nucl' -hash_index
```


- ### Name your genes of interest as follows


> 16s_rRNA_gene.fasta


- ### Run the following script


```
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

```


## **Output files**


<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/gene_seq_extraction_from_multiple_genomes/blob/main/Output_files.PNG"
  >
</p>
<p align = "center">
Figure 1: Output files for the 16s rRNA gene.
</p>


### The following file contains the 16s_rRNA_gene sequences from all the supplied Pectbacterium genomes  



<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/gene_seq_extraction_from_multiple_genomes/blob/main/Extracted_sequences.PNG"
  >
</p>
<p align = "center">
Figure 2. 16s rRNA gene sequences extracted from all the supplied Pectobacterium genomes
</p>


## The end


