#!/bin/bash

#SBATCH -J X # Name of the job:
#SBATCH -A X # Which project should be charged:
#SBATCH -p X # choose partition
#SBATCH --nodes=1 # How many whole nodes should be allocated?
#SBATCH --ntasks=56
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

#! Assign user-supplied arguments
infileR1="Sample1_1.fastq.gz"
infileR2="Sample1_2.fastq.gz"
name="Sample1"

#! sensitive alignment with up to 200 alignments per read (no unaligned, no discordants) | sam to bam | sort bam | exclude secondary alignments (-F 256), quality score above 5 (-q 5), retain header (-h)
bowtie2 -p 56 --very-sensitive -k 200 --no-unal --no-discordant -x /pathtoyour/Reference/Col_CEN/ColCEN_bowtie2_index/ColCEN_index -1 ./"${infileR1}" -2 ./"${infileR2}" 2> "${name}".bowtie_summary.txt | samtools view -@ 56 -bS - | samtools sort -@ 56 - | samtools view -@ 56 -F 256 -q 5 -h -b - > "${name}".primary.sort.bam

#! sort by n
samtools sort -@ 56 -n -o "${name}".nsorted.primary.sort.bam -O BAM "${name}".primary.sort.bam
rm "${name}".primary.sort.bam

#! this will fill in mate coordinates and insert size fields
samtools fixmate -@ 56 -m "${name}".nsorted.primary.sort.bam "${name}".fixmate.nsorted.primary.sort.bam
rm "${name}".nsorted.primary.sort.bam

#! This will sort based on chromosome number and coordinates
samtools sort -@ 56 -o "${name}".chrsorted.fixmate.nsorted.primary.sort.bam "${name}".fixmate.nsorted.primary.sort.bam
rm "${name}".fixmate.nsorted.primary.sort.bam

#! This will remove all the duplicates and also print some basic stats about the result file (-s)
samtools markdup -@ 56 -r "${name}".chrsorted.fixmate.nsorted.primary.sort.bam "${name}".final.bam -s 2> "${name}".markdup.txt
rm "${name}".chrsorted.fixmate.nsorted.primary.sort.bam

#! index and tracks
samtools index "${name}".final.bam
bamCoverage -b "${name}".final.bam --numberOfProcessors 56 --normalizeUsing RPGC --effectiveGenomeSize 135000000 --binSize 10 -o "${name}".RPGC_10bp.bw
bamCoverage -b "${name}".final.bam --numberOfProcessors 56 --normalizeUsing RPGC --effectiveGenomeSize 135000000 --binSize 10 --outFileFormat bedgraph -o "${name}".RPGC_10bp.bg
