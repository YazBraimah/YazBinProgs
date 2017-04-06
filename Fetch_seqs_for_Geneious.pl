#!/bin/bash

if [ $# == 0 ]; then
    
    echo "Usage: Fetch_seqs_for_Geneious.pl <TCONS_CDS_ID>"
fi 

faOneRecord merged.CDS.ALL.uniq.fasta $1 | sed '/amr/ s/>/>D.amr /g' | sed '/lum/ s/>/>D.lum /g' | sed '/vir/ s/>/>D.vir /g' | sed '/nov/ s/>/>D.nov /g'| sed 's/ species=.*//g' > $1.seqs.fa

