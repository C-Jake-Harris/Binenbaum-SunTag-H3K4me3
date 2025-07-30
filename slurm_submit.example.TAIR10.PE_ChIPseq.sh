#!/bin/bash

#SBATCH -J X # Name of the job:
#SBATCH -A X # Which project should be charged:
#SBATCH -p X # choose partition
#SBATCH --nodes=1 # How many whole nodes should be allocated?
#SBATCH --ntasks=56
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-ccl              # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
module load bowtie/2.5.0
module load samtools/1.10


#### USER ARGUMENTS START ###

infileR1="sample1_1.fq.gz"
infileR2="sample1_2.fq.gz"
name="sample1"

#### USER ARGUMENTS END ###


#! bowtie2 alignment (no unaligned, no discordants)
bowtie2 -p 56 --no-unal --no-discordant -x ~/pathtoyour/Bowtie2Index/genome -1 ./"${infileR1}" -2 ./"${infileR2}" 2> "${name}".bowtie_summary.txt 1> "${name}".sam

#! sam to bam
samtools view -@ 56 -bS "${name}".sam > "${name}".bam

#! remove sam
rm "${name}".sam

####### De-duplication #######
#! sort by n
samtools sort -@ 56 -n -o "${name}".nsort.bam -O BAM "${name}".bam
rm "${name}".bam

#! this will fill in mate coordinates and insert size fields
samtools fixmate -@ 56 -m "${name}".nsort.bam "${name}".fixmate.nsort.bam
rm "${name}".nsort.bam

#! This will sort based on chromosome number and coordinates
samtools sort -@ 56 -o "${name}".chrsorted.fixmate.nsort.bam "${name}".fixmate.nsort.bam
rm "${name}".fixmate.nsort.bam

#! This will remove all the duplicates and also print some basic stats about the result file (-s)
samtools markdup -@ 56 -r "${name}".chrsorted.fixmate.nsort.bam "${name}".final.bam -s 2> "${name}".markdup.txt
rm "${name}".chrsorted.fixmate.nsort.bam
####### De-duplication end #######

#! index and tracks
samtools index "${name}".final.bam
bamCoverage -b "${name}".final.bam --numberOfProcessors 56 --normalizeUsing RPGC --blackListFileName /pathtoyour/Reference/blacklist_greenscreen_Klasfeld2022.bed --effectiveGenomeSize 135000000 --binSize 10 -o "${name}".RPGC_10bp.bw





