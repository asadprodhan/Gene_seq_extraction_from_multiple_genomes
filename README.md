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


Convert pecto_genomes.fasta into a blast database
makeblastdb -in pecto_genomes.fasta -out pecto_genomes_db -dbtype 'nucl' -hash_index
Name your genes of interest as follows
16s_rRNA_gene.fasta
Run the following script
The following file contains the 16s_rRNA_gene sequences from all the supplied Pectbacterium genomes  

