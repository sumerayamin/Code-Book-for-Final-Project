# Code-Book-for-Final-Project
#A) Obtain putative homologs for the protein of interest.
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
  
 7. Create a more detailed process output file of the same analysis
  
  blastp -db ~/data/blast/allprotein.fas -query XP_001618798.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -     max_hsps 1 -out XP_001618798.blastp.detail.out
  
 8. Filter the output file to get the e-value that has to be less than 1e-14.:
 
  awk '{if ($6<0.00000000000001)print $1 }' XP_001618798.blastp.detail.out > XP_001618798.blastp.detail.filtered.out
   
 9. Obtain the sequences of proteins that resulted after the BLAST output using seqkit
  
   seqkit grep --pattern-file XP_001618798.blastp.detail.filtered.out ~/data/blast/allprotein.fas > XP_001618798.blastp.detail.filtered.fas
   
  #B) Perform a global multiple sequence alignment on Protein of interest putative homologs.
   
 10. Perform a global multiple sequence alignment in muscle
  
   muscle -in XP_001618798.blastp.detail.filtered.fas -out XP_001618798.blastp.detail.filtered.aligned.fas
    
 11. Provide some statistics about the alignment using t_coffee:
  
   t_coffee -other_pg seq_reformat -in XP_001618798.blastp.detail.filtered.aligned.fas -output sim
	
 12. Remove any column that contains greater than 50% gapped residues using t_coffee
 
  t_coffee -other_pg seq_reformat -in XP_001618798.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out allhomologs.aligned.r50.fa


