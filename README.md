### SAPA

SAPA v1.1 (Sequence Alignment and PAML Analysis) aligns coding sequences, translates them into proteins, calculates pairwise Ka/Ks, and performs PAML analyses.

				
## Requirements:
		
		1. In addition to the BioPerl modules required for this script to work (see code header), the following programs/scripts are required 
		(make sure they are in your path variable):
			
			- ucsc utilities (faOneRecord & faSomeRecords, found here -> http://hgdownload.soe.ucsc.edu/admin/exe/)
			- PHAST (tree_doctor, found here -> http://compgen.cshl.edu/phast/)
			- fastaSortByName.pl, found here -> http://raven.iab.alaska.edu/~ntakebay/teaching/programming/perl-scripts/fastaSortByName.pl
			- paml_prep.pl  (ask yazahmed@gmail.com)

		2. CDS file should contain all sequences to be analyzed. Header should have follwoing format:
			>gene_id .......... species=<species_id>

		3. If the orthologous genes to be analyzed have different gene_id's, then must use "--include_outgroups" option and provide "orthology_file"

		
 ## Required arguments:

--CDS_file|c <string>				sequence file containing CDS from multiple species.

--gene_list|g <string>				file containing list of genes to be analyzed (one gene_id/line).
	OR
--gene_id|i <string>				individual gene_id to be analyzed
  
--output|o					name of directory for outputs 


  Options:

--show_samples					(optional) display species IDs in CDS_file and exit (only requires "-c <CDS_file>" argument).

--include_outgroups				(must be specifies if orthologous sequecnes have different gene_id's)

	--orthology_file <string>		file with orthology information. format: <transcript_id>\t<orthologue_id>


--view_DNA_alignment				(optional. Requires Geneious, obviously) open DNA alignment in Geneious.

--view_Protein_alignment			(optional. Requires Geneious, obviously) open Protein alignment in Geneious.

--save_alignments				(optional) retain alignment files

--calculate_KaKs				Calculate pairwise Ka/Ks.

--run_PAML					Perform PAML codeml analyses

	--tree_file				(required) newick format phylogeny
	
	*** Running PAML will execute the default NULL, i.e. model = 0, NSsites = 0. To execute the branch model and/or the branch-site model(s), specify as follows:

	--branch	(optional)
	--branchSite	(optional)	
	--restrict_samples			File with species IDs to use in branch-type PAML tests. Required if --branch or --branchSite is specificied. 
								One ID per line. must match species ID in *.tree file and FASTA file.

--help						

					
