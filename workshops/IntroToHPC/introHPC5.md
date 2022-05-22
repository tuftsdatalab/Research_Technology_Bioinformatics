# Optional: STAR Loop Script

Often times you'll need to run several files through some loop to do the same command. A common NGS example is to loop through fastq files to align them to a genome. Let's create a script to align our fastq files in the `data` folder to the Mouse Genome! 

```
[tutln01@c1cmp047 ~]$ cd introHPC
[tutln01@c1cmp047 introHPC]$ nano starLoop.sh
```

Please copy and paste this script into starLoop.sh and then hit `Esc`+`Cntl`+`x`+`y`+`Enter` to save and exit:

```
#!/bin/bash
#SBATCH --job-name=starLoop              # name your job
#SBATCH --time=03-00:00:00               # how long your job might take
#SBATCH --partition=preempt              # which partition you want to run it on
#SBATCH --nodes=1                        # how many nodes do you want
#SBATCH --mem=32Gb                       # how much memory do you want
#SBATCH --output=%j.out                  # name of output file
#SBATCH --error=%j.err                   # name of error file
#SBATCH --mail-type=ALL                  # request to be emailed when job begins and ends
#SBATCH --mail-user=YourEmail@tufts.edu  # provide your email
 
## Load STAR aligner
module load STAR/2.7.0a

#make directories
mkdir /cluster/home/tutln01/introHPC/STAR

## Defing reference genome, gtf, raw data and output directories
REF_DIR=/cluster/tufts/bio/data/genomes/Mus_musculus/UCSC/mm10/Sequence/STAR
GTF_DIR=/cluster/tufts/bio/data/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf
DATA_DIR=/cluster/home/tutln01/introHPC/data
OUT_DIR=/cluster/home/tutln01/introHPC/STAR

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
Quite the script! Let's break it down:
- First we have our `#SBATCH` headers which tell SLURM how to run our job
- Next we load the STAR module with `module load STAR/2.7.0a`
- Now we will make an output directory for our results with `mkdir /cluster/home/tutln01/introHPC/STAR`
- Next we define variables that point to where:
 - `$REF_DIR` our mouse reference data is 
 - `$GTF_DIR` where the gtf file is
 - `$DATA_DIR` where our input fastq files are
 - `$OUT_DIR` where our output directory is
- Then we loop through all files in `$DATA_DIR` that end in `fastq.gz`
- These files are stripped to get the "base name" (SRR accession without "fastq.gz")
- Now the file is run through the STAR command which will align it to the reference genome in `$REF_DIR`, use the GTF file in `$GTF_DIR` for coordinates, and output a bam file in the `$OUT_DIR` folder. 

We can now submit this script with sbatch like so:

```
[tutln01@c1cmp047 introHPC]$ sbatch starLoop.sh
```

Once the script has run you'll note that output bam files are located in the STAR directory that we created in the script:

```
[tutln01@c1cmp047 introHPC]$ cd STAR
[tutln01@c1cmp047 introHPC]$ ls
```
NOTE: For a more in-depth look at RNA-seq workflows, check out our [RNA-seq tutorial](https://huoww07.github.io/Bioinformatics-for-RNA-Seq/)

_____________________________________________________________________________________________________________________________________________________

[Next](./introHPC6.md)

[Previous](./introHPC4.md)
