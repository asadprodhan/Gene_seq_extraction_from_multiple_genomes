# **How do I automate extracting multiple gene sequences from multiple genomes?** <br />


### **Asad Prodhan PhD** 

**https://asadprodhan.github.io/**

<br />


## **Required tools**


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


<br />


## **Methods**


<br />


## Prepare the genome sequences


### I. Concatenate all the genomes into a single fasta file


```
cat *.fasta > concatenated_genomes.fasta
```


### II. Convert the concatenated_genomes.fasta into a blast database


```
makeblastdb -in concatenated_genomes.fasta -out concatenated_genomes_db -dbtype 'nucl' -hash_index
```


### III. Make a directory named "concatenated_genomes_db"


### IV. Move all the concatenated_genomes_db.* files into the concatenated_genomes_db directory


<br />


## Prepare the gene sequences


### I. Name your genes of interest as follows


> All gene sequence file names must have _gene.fasta extension

> rsmD_NZ_CP065044_gene.fasta


<br />


## Run the extraction


### I. Put all the sequence files along with the follwoing script in the same directory 


### II. Run the following commands one-by-one


```
chmod +x *
```


> This will make the files executable


```
dos2unix *
```


> This will make sure that all the files are in unix format



### III. Run the script as follows


```
./name-of-the-script.sh
```



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
    blastn -query "${file}" -task blastn -db concatenated_genomes_db/concatenated_genomes_db -outfmt 6 -out ${base}_blastn_hits.tsv -evalue 1e-10 -num_threads 18 
    echo "${Green}${Bold}Blastn done ${reset}: "${file}""
    echo "${Red}${Bold}Making bed file ${reset}: "${file}"" 
    awk '{print $2,$9,$10,""$1"_"NR}' ${base}_blastn_hits.tsv > ${base}_blastn_hits.bed
    echo "${Green}${Bold}Bed file made ${reset}: "${file}""
    echo "${Red}${Bold}Sorting bed file ${reset}: "${file}"" 
    cat  ${base}_blastn_hits.bed | awk -v OFS='\t' '{ if ($2 < $3) {print $1,$2,$3,$4} else {print $1,$3,$2,$4} }' >  ${base}_blastn_hits_sorted.bed
    echo "${Green}${Bold}Bed file sorted ${reset}: "${file}""
    echo "${Red}${Bold}Extracting blastn hit sequences ${reset}: "${file}"" 
    bedtools getfasta -fi concatenated_genomes.fasta -bed ${base}_blastn_hits_sorted.bed -name > ${base}_blastn_hits.fasta
    echo "${Green}${Bold}Extracted blastn hit sequences ${reset}: "${file}""
    echo "${Red}${Bold}Naming headers with the corresponding gene and genome names ${reset}: "${file}"" 
    awk -v var=$(basename $file .fasta) '{ gsub(/contig/, var "_contig") } 1' ${base}_blastn_hits.fasta > ${base}_blastn_hits_seqs.fasta
    echo "${Green}${Bold}Naming headers done ${reset}: "${file}""
    rm -r ${base}_blastn_hits.fasta *.bed *.tsv 
    echo ""
done


```


<br />



## **Output files**


<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/gene_seq_extraction_from_multiple_genomes/blob/main/Output_files.PNG"
  >
</p>
<p align = "center">
Figure 1: Output files for the rsmD_NZ_CP065044 gene.
</p>


<br />


### The following file (rsmD_NZ_CP065044_gene_blastn_hits_seqs.fasta) contains the rsmD_NZ_CP065044_gene sequences from all the supplied genomes  



<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/gene_seq_extraction_from_multiple_genomes/blob/main/Extracted_sequences.PNG"
  >
</p>
<p align = "center">
Figure 2. rsmD_NZ_CP065044 gene sequences extracted from all the supplied genomes
</p>


## The end


