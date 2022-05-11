# Scripts That Involve Loops

Often times you'll need to run several files through some loop to do the same command. A common NGS example is to loop through fastq files to align them to a genome. Let's create a script to align our fastq files in the `data` folder to the Mouse Genome! 

```
[tutln01@c1cmp047 ~]$ cd introHPC
[tutln01@c1cmp047 introHPC]$ nano starLoop.sh
```
```
#!/bin/bash
#SBATCH --job-name=starLoop                   #
#SBATCH --time=03-00:00:00                    #
#SBATCH --partition=preempt                   #
#SBATCH --nodes=1                             #
#SBATCH --mem=32G                             #
#SBATCH --output=%j.out                       #
#SBATCH --error=%j.err                        
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YourEmail@tufts.edu
 
## Load STAR aligner
module load STAR/2.7.0a

#make directories
mkdir /cluster/home/tutln01/introHPC/STAR
mkdir /cluster/home/tutln01/introHPC/STAR/output

## Defing reference genome, gtf, raw data and output directories
REF_DIR=/cluster/tufts/bio/data/genomes/Mus_musculus/UCSC/mm10/Sequence/STAR
GTF_DIR=/cluster/tufts/bio/data/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf
DATA_DIR=/cluster/home/tutln01/introHPC/data
OUT_DIR=/cluster/home/tutln01/introHPC/STAR/output

for i in ${DATA_DIR}/*fastq.gz
do
    #create output file name
    f=$(basename -- "$i")
   
    #run STAR
    STAR --genomeDir ${REF_DIR} \
    --readFilesIn $i\
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUT_DIR}/$f \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 4 \
    --alignIntronMin 1 \
    --alignIntronMax 2500 \
    --sjdbGTFfile ${GTF_DIR} \
    --sjdbOverhang 99
done
```

