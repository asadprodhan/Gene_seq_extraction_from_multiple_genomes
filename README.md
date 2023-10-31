# **How do I automate extracting multiple gene sequences from multiple genomes?** <br />


Tools
create a conda environment as follows
conda create -n gene_seq_extraction
activate your conda environment
conda activate gene_seq_extraction
install blast 
conda install -c bioconda blast
install bedtools
conda install -c bioconda bedtools

Steps
concatenate all the genomes into a single fasta file
cat Pectbacterium* > pecto_genomes.fasta
Convert pecto_genomes.fasta into a blast database
makeblastdb -in pecto_genomes.fasta -out pecto_genomes_db -dbtype 'nucl' -hash_index
Name your genes of interest as follows
16s_rRNA_gene.fasta
Run the following script
The following file contains the 16s_rRNA_gene sequences from all the supplied Pectbacterium genomes  

