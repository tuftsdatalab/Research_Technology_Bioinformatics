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

After running FastQC, you will notice several files in your directory:

```
SOdaskfsjdfdfasadsklfjdlsafj
```

To investigate the quality of our sequence, we will need to view the `.html` file that was produced. Navigate to the OnDemand Tab and click on `Files > Home Directory`:

![](../images/files_home.png)

Now navigate to this workshop's directory and click on the `.html` file and then click View at the top of the screen:

![]()

## FastQC Output

Here we notice several quality control plots, let's take a minute to discuss what some of these plots mean.

### Sequence Quality

![](../seq_qual_hist.png)

### GC Content

### Per base sequence content

### Adapters in Sequence

_________________________________________________________________________________________________________________________________________________________

[Next](lesson5.md)

[Previous](lesson3.md)
