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
trim_galore --illumina --paired -o trim_galore/ fastq/SRR15607266_pass_1.fastq.gz fastq/SRR15607266_pass_2.fastq.gz
```

- `--illumina` Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter
                        'AGATCGGAAGAGC' instead of the default auto-detection of adapter sequence.
- `--paired` This option performs length trimming of quality/adapter trimmed reads for
                        paired-end files
- `-o` specifies that our output directory will be a folder called trim_galore

## Trim Galore Output

When you run the command above you will note several lines of output to the screen. We will focus on the summary section:

```
=== Summary ===

Total reads processed:               1,264,290
Reads with adapters:                   451,277 (35.7%)
Reads written (passing filters):     1,264,290 (100.0%)

Total basepairs processed:    96,086,040 bp
Quality-trimmed:                 257,279 bp (0.3%)
Total written (filtered):     95,110,835 bp (99.0%)
```
Here we note that only a small fraction of bases were removed, 0.3%, due to quality control issues. You will also note that the adapter sequence was detected in 451,277 reads while FastQC did not detect any adapters. This is due to the fact that FastQC requires the entire adapter to be present while Trim Galore does not. And as we can see while the adapter sequence was detected it did not trigger Trim Galore to remove any reads with those partial adapter sequences.

_________________________________________________________________________________________________________________________________________________________

Next Lesson: [Ngs Data Alignment and Viewing](../../NgsDataAlignmentViewing/lessons/lesson1.md)

[Previous](lesson4.md)
