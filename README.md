# Code-Book-for-Final-Project
A) Obtain putative homologs for the protein of interest.
1. Create a directory for the BLAST database:

  mkdir ~/data/blast

2. Uncompress the proteomes

  gunzip proteomes/*.gz
	
3. Put all the protein sequences into a single file.

  cat  proteomes/* > ~/data/blast/allprotein.fas
	
 4. Build a BLAST database with proteomes 
 
  makeblastdb -in ~/data/blast/allprotein.fas -parse_seqids -dbtype prot
	
 5. Download protein sequence of interest for BLAST
 
  ncbi-acc-download -F fasta -m protein XP_001618798.1
	
 6. Perform a BLAST search to obtain proteins that are homologous to the query protein
 
  blastp -db ~/data/blast/allprotein.fas -query XP_001618798.1.fa -outfmt 0 -max_hsps 1 > XP_008798.blastp.typical.out
	
	
