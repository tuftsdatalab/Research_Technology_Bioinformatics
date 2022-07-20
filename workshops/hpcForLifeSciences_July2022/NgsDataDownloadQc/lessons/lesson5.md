# Process Reads

## Learning Objectives

 - Understand when to trim NGS data
 - Trim NGS data with Trim-Galore
 - Understand Trim-Galore Output
 
## Trimming NGS Data

NGS data can contain reads with poor base pair quality and adapters still present. These poor bases and adapters can interfere with the accuracy of downstream analyses. To get rid of these poor quality bases and adapters we will need a tool that can perform both. Here we will use the tool Trim Galore to do just that. Please navigate back to the terminal tab and use the following command to load the Trim Galore module:

```
module load trim-galore/0.6.4_dev
```

To trim our data we will enter:

```
trim_galore --illumina --paired --fastqc -o trim_galore/ fastq/SRR15607266_pass_1.fastq.gz fastq/SRR15607266_pass_2.fastq.gz
```

- `--illumina` Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter
                        'AGATCGGAAGAGC' instead of the default auto-detection of adapter sequence.
- `--paired` This option performs length trimming of quality/adapter/RRBS trimmed reads for
                        paired-end files
- `-o` specifies that our output directory will be a folder called trim_galore

## Trim Galore Output

When you run the command above you will note several lines of output to the screen. We will focus on the summary section.

_________________________________________________________________________________________________________________________________________________________

[Main Page](../README.md)

[Previous](lesson4.md)
