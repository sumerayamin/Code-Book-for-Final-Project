# Code-Book-for-Final-Project

# A) Obtain putative homologs for the protein of interest.

  1) Create a directory for the BLAST database:

    mkdir ~/data/blast

  2) Uncompress the proteomes

    gunzip proteomes/*.gz
	
  3) Put all the protein sequences into a single file.

         cat  proteomes/* > ~/data/blast/allprotein.fas
	
  4) Build a BLAST database with proteomes 
 
    makeblastdb -in ~/data/blast/allprotein.fas -parse_seqids -dbtype prot
	
  5) Download protein sequence of interest for BLAST
 
    ncbi-acc-download -F fasta -m protein XP_001618798.1
	
6)Perform a BLAST search to obtain proteins that are homologous to the query protein
  
    blastp -db ~/data/blast/allprotein.fas -query XP_001618798.1.fa -outfmt 0 -max_hsps 1 > XP_008798.blastp.typical.out
  
 7) Create a more detailed process output file of the same analysis
  
        blastp -db ~/data/blast/allprotein.fas -query XP_001618798.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -     max_hsps 1 -out XP_001618798.blastp.detail.out
  
 8) Filter the output file to get the e-value that has to be less than 1e-14.:
 
        awk '{if ($6<0.00000000000001)print $1 }' XP_001618798.blastp.detail.out > XP_001618798.blastp.detail.filtered.out
   
 9) Obtain the sequences of proteins that resulted after the BLAST output using seqkit
  
        seqkit grep --pattern-file XP_001618798.blastp.detail.filtered.out ~/data/blast/allprotein.fas > XP_001618798.blastp.detail.filtered.fas
   
  # B) Perform a global multiple sequence alignment on Protein of interest Putative Homologs.
   
 10) Perform a global multiple sequence alignment in muscle
  
    muscle -in XP_001618798.blastp.detail.filtered.fas -out XP_001618798.blastp.detail.filtered.aligned.fas
    
 11) Provide some statistics about the alignment using t_coffee:
  
    t_coffee -other_pg seq_reformat -in XP_001618798.blastp.detail.filtered.aligned.fas -output sim
	
 12) Remove any column that contains greater than 50% gapped residues using t_coffee
 
    t_coffee -other_pg seq_reformat -in XP_001618798.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out allhomologs.aligned.r50.fa
  
  # Install Newick Utilities software from GitHub
  
     Install newick on the instance: 
  
     git clone git://github.com/tjunier/newick_utils.git
   
     cd newick_utils/
   
     autoreconf -fi
   
     ./configure
   
      make
   
     sudo make install
   
   # C) Create a Phylogenic Tree: 
     
   13. Run the IQ-TREE command to create a phylogenetic tree using default parameters and arbitrary rooting
    
    iqtree -s mygene.aligned.r50.fa -nt 2
     
   14) Create a file of the unrooted tree
   
    nw_display mygene.aligned.r50.fa.treefile
    
   15) Create an image of unrooted tree
   
    gotree draw png -w 1000 -i mygene.aligned.r50.fa.treefile  -r -o  mygene.aligned.r50.fa.png
    
   16) Create a midpoint rooted tree
   
    gotree reroot midpoint -i mygene.aligned.r50.fa.treefile -o mygene.aligned.r50.fa.midpoint.treefile
    
   17) Root the Optimal phylogeny 
   
    nw_reroot mygene.aligned.r50.fa.treefile Drosophila_melanogaster_Disks_large_1_DLG1_P31007 Nematostella_vectensis_disks_large_homolog_1_XP_001638123.2 Pocillopora_damicornis_A0A3M6TL78 Homo_sapiens_Disks_large_homolog_3_DLG3_Q92796 Homo_sapiens_Disks_large_homolog_4_DLG4_P78352 Homo_sapiens_Disks_large_homolog_1_DLG1_Q12959 Homo_sapiens_Disks_large_homolog_2_DLG2_Q15700 >mygene.aligned.r50.fa.MendozaRoot.treefile

# D) Evaluating Node Support using the Bootstrap

 18) Commad to obtain the Bootstrap support
 
    iqtree -s mygene.aligned.r50.fa -bb 1000 -nt 2 --prefix mygene.r50.ufboot
  
 19) Look at the bootstrap support

    nw_display mygene.r50.ufboot.Midpointroot.treefile 

 20) Produce a graphic to see the bootstrap values

    nw_display -s mygene.alligned.r50.ufboot.Midpointroot.treefile  -w 1000 -b 'opacity:0' > mygene.alligned.r50.ufboot.Midpointroot.svg
 
 Using Notung
 
    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar --help 
 
 21) Reconcile the toy gene tree and toy species tree in notung

    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b toybatch.txt --reconcile --speciestag prefix --savepng --treestats --events --homologtabletabs -- phylogenomics 

 22) Run the following once to install a missing component

    sudo easy_install -U ete3 

 23) Now, generate a RecPhyloXML object by running the following command:

    python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g toygenetree.tre.reconciled --include.species
  
  # Reconciling gene family with the species tree
  
 24) The species tree
 
    (((Homo_sapiens,Strongylocentrotus_purpuratus)Deuterostomia,Drosophila_melanogaster)Bilateria,(Nematostella_vectensis,Pocillopora_damicornis)Cnidaria)Eumetazoa;

25) Make a batch file for Notung in nano

        nano mybatch.txt
 
 26) Perform the reconciliation:
 
    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b mybatch.txt --reconcile --speciestag prefix  --savepng --treestats --events  --phylogenomics 
  
 27) Generate a RecPhyloXML object to view the gene-within-species tree
 
    python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g mygene.aligned.r50.ufboot.Midpointroot.treefile.reconciled --include.species
  
 28) Create a visual aid that includes the internal node names
   
    nw_display  species.tre 
   
 29) In notung, replace the --reconcile command with the --root flag to re-root the tree during using the reconciliation process.
 
    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b mybatch.txt --root --speciestag prefix  --savepng --treestats --events  --phylogenomics 
  
# E) Reconciliation Based Rearrangement and Topology Tests
  
 30) Reconcile and re-arrange your tree
 
    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b mybatch.txt --rearrange --speciestag prefix  --savepng --treestats --events  --outputdir zoreconcileRearrange --edgeweights name --threshold 90 
 
 31) Visualize the re-arranged gene tree within the species tree
 
    python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g mygene.aligned.r50.ufboot.Midpointrrot.treefile.rearrange.0 --include.species
    
 32) Change the format from Notung to Newick
 
    java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -g zoreconcileRearrange/mygene.aligned.r50.ufboot.MendozaRoot.treefile.rearrange.0   -s species.tre --reconcile --speciestag prefix  --treeoutput newick --nolosses
    
 33) Unroot the Notung re-arranged tree 
  
         gotree unroot -i mygene.aligned.r50.ufboot.Midpointroot.treefile.rearrange.0.reconciled -o mygene.aligned.r50.ufboot.unrooted.treefile.rearrange
	 
 34)  Put the two trees into a single file
 
     cat mygene.aligned.r50.ufboot.treefile mygene.aligned.r50.ufboot.unrooted.treefile.rearrange > mygene.aligned.r50.alternativetrees
     
 35) Run the topology test in IQ-TREE 
 
    iqtree -s mygene.aligned.r50.fa -z mygene.aligned.r50.alternativetrees -au -zb 10000 --prefix ZO1_altTrees -m LG+F+R5 -nt 2 -te mygene.aligned.r50.ufboot.treefile
    
# F) Domain Prediction using Interproscan

 36) A one time installation of perl module
 
         sudo cpan LWP::Protocol::https
 
 37) Installation of the DataMash utility
 
          cd ~/tools
          wget http://ftp.gnu.org/gnu/datamash/datamash-1.3.tar.gz
          tar -xzf datamash-1.3.tar.gz  
          cd datamash-1.3
          ./configure
          make
          make check
          sudo make install 

 38) Iprscan5 is a perl script that will send the sequences to the inteproscan servers to be analyzed. It may take up to 15 minutes to finish analyzing the sequences.
 
    iprscan5   --email sumera.aymin@stonybrook.edu  --multifasta --useSeqId --sequence   mygene.aligned.blastp.detail.filtered.ren.fas
 
 39) After getting the domains identified for each gene use the cat command to concatenate them all into a single file
  
    cat *.tsv.txt > mygene.domains.all.tsv
  
 40) To filter domains defined by the Pfam database
 
    grep Pfam mygene.domains.all.tsv >  mygene.domains.pfam.tsv
    
 41) Fix up the output by re-arranging the interproscan output
  
    awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' mygenefamily.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > mygenefamily.domains.pfam.evol.tsv
    
  After getting the results go to https://www.evolgenius.info/evolview/

  Upload the tree by clicking on the File folder in the upper left hand corner. Give the tree a name. The upload button is at the bottom of the box that appears. After clicking submit, tree should be displayed. 

 
 
