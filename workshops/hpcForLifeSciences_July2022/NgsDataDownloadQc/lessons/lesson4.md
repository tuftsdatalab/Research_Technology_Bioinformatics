# Quality Control

## FastQC

Next Generation Sequencing can produce a large number of reads in each experiment, giving low-cost and in-depth information about the underlying RNA or DNA sample. However, every platform will produce errors (incorrect nucleotides in the sequence). Hence, quality control is an important step in data analysis. FastQC provides several modules to asses the quality of sequencing data:

- Sequence Quality
- GC content
- Per base sequence content
- Adapters in Sequence

To run FastQC we will need to load the FastQC module:

```
module load fastqc/0.11.9
```

Now that we have it loaded we can run FastQC on our sequencing data:

```
fastqc SRR15607266.fastq.gz
```

### Sequence Quality

### GC Content

### Per base sequence content

### Adapters in Sequence

