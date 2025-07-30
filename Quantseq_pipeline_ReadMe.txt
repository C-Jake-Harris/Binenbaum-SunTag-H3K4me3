##################################################################################
#### Overview of the pipeline used to quantify data from Quant-seq libraries #####
##################################################################################

# As Quantseq (QuantSeq 3‘ mRNA-Seq Library Prep Kit for Illumina (FWD) only recommends using read 1, move to dir with only the read 1 fastq files to be processed
# see https://www.lexogen.com/wp-content/uploads/2015/11/015UG009V0211_QuantSeq-Illumina.pdf

#########################################
### Step 1: align the reads with STAR ###
#########################################
# ensure that the STAR aligner is set up and paths are correct
# run the loop script below (it generates an output dir called Aligned)
# note for Quantseq, the gene counts that are from colum 3 of the *ReadsPerGene.out.tab file

./Quant_pipe_loop.sh ./ /pathtoyour/STAR/genomeIndex_2 Aligned

#########################################
### Step 2: extract count table data ####
#########################################
# the above creates a gene level counts table called Merged_CountsTable.txt
# Use the counts data for downstream analysis in R e.g. DEseq2

./Quantseq_extractcounts.sh


