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


```
mkdir concatenated_genomes_db
``` 


### IV. Move all the concatenated_genomes_db.* files into the concatenated_genomes_db directory


```
mv concatenated_genomes_db* concatenated_genomes_db 
``` 


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


[DOWNLOAD the script here](https://github.com/asadprodhan/Gene_seq_extraction_from_multiple_genomes/blob/main/gene_seq_extraction_auto_v10.sh) 


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


> Note how the headers of the blastn hits have been named using the corresponding gene and genome names


> This is for the convenience of tracking information


<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/Gene_seq_extraction_from_multiple_genomes/blob/main/Extracted_sequences_v2.png"
  >
</p>
<p align = "center">
Figure 2. rsmD_NZ_CP065044 gene sequences extracted from all the supplied genomes
</p>


## The end


