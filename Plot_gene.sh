#!/bin/bash

if [ $# == 0 ]; then
    echo "Your command line contains no arguments"
    echo "Usage: Plot_gene.sh <matrix_file> <samples_file> <gene_name> <log_mode>"
    exit
fi

if [ -n "$4" ]; then
	PtR --matrix $1 --samples $2 --output $3 --no_reuse --per_gene_plots --gene_grep $3 --log2
    rm $3.R $3.minCol10.minRow10.log2.dat
	mv $3.minCol10.minRow10.log2.per_gene_plots.pdf $3.TPM_log2_plot.pdf
else 
    PtR --matrix $1 --samples $2 --output $3 --no_reuse --per_gene_plots --gene_grep $3
    rm $3.R $3.minCol10.minRow10.dat
	mv $3.minCol10.minRow10.per_gene_plots.pdf $3.TPM_plot.pdf   
fi

