# Obtain Public NGS Data

## Learning Objectives

- Retrieve reference data from [US National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/)
- Download SRA Next Generation Sequencing Data
- Understand Fasta, GFF, and Fastq file formats

## SARS-Cov-2 Reference Data

The [US National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/) (NCBI) hosts repositories for many types of biomedical and genomics data. Today we'll retrieve reference data from the [Genomes Database FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/) as well as the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). 

### Downloading Reference Data

We will enter our `data` directory and use the `wget` command to download our reference data from the NCBI repository:

```
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz
```
You will notice that both these files end in `.gz` - this indicates they are compressed. Compressing files is useful when storing data but to use it you will often need to decompress it with `gunzip -d`:

```
gunzip -d GCF_009858895.2_ASM985889v3_genomic.fna.gz
gunzip -d GCF_009858895.2_ASM985889v3_genomic.gff.gz
```

You will also notice that one file ends in `.fna` and the other ends in `.gff`. These are Fasta and GFF files, respectively.

### Fasta Format
The virus genome is in fasta format. Fasta format has two parts, a sequence identifier preceeded by a ">" symbol, followed by the sequence on subsequent lines.

<p align="center">
<img src="../images/genome_view_2,PNG.png" width="900">
</p>


### GFF Format
The gene annotation file is in Generic Feature Format (GFF). This format tells us where genes are located in the reference genome. Note that we must always be sure that our gene information and genome come from the same source. 

<p align="center">
<img src="../images/gff.PNG">
</p>

## SARS-Cov-2 NGS Sequencing Data

We are interested in obtaining reads from the sample [Viral genomic RNA sequencing of a B.1.617.2/Delta isolate; Severe acute respiratory syndrome coronavirus 2; RNA-Seq](https://www.ncbi.nlm.nih.gov/sra/?term=SRR15607266)
<p align="center">
<img src="../images/srr.png" width="900">
</p>

### Download NGS Sequencing Data

To download our NGS Data we will need to load the sratoolkit module. We can access this module by running:

```
module load sra/2.10.8
```
Now that we have this module loaded we can use it to pull our SARS-Cov-2 NGS sequencing data from the [SRA database](https://www.ncbi.nlm.nih.gov/sra):

```
fastq-dump --outdir fastq --gzip --skip-technical --read-filter pass --dumpbase --split-files --clip SRR15607266
```
- `fastq-dump` is the command to pull SRA data
- `--outdir` specifies where you want your SRA data deposited
- `--gzip` compress the output using gzip
- `--skip-technical` only download biological reads and skip technical reads
- `--read-filter pass` return reads that pass filtering, so reads without N's
- `--dumpbase` formats sequences by base space
- `--split-files` splits files into left and right ends, if these reads are not paired-end they are combined into a single file
- `--clip` this option will remove reads that SRA tags to be removed

We can enter our `fastq` directory and see that we have two files for this sample, seeing as this is paired-end data:

```
cd fastq
ls
```

>```
>SRR15607266_pass_1.fastq.gz  SRR15607266_pass_2.fastq.gz
>```

### Fastq format
Fastq format is a way to store both sequence data and information about the quality of each sequenced position.Each block of 4 lines contains one sequencing reads, for example:
```
@SRR15607266.1 1 length=76
NTTATCTACTTTTATTTCAGCAGCTCGGCAAGGGTTTGTTGATTCAGATGTAGAAACTAAAGATGTTGTTGAATGT
+SRR15607266.1 1 length=76
#8ACCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
```

1. Sequence identifier
2. Sequence
3. \+ (optionally lists the sequence identifier again)
4. Quality string

Paired end sequencing data will typically be stored as two fastq files, one for the forward and one for the reverse.  Each file should contain the same number of reads, with the same labels, in the same order. If this convention is not followed, it could cause errors with downstream tools. Fortunately there are tools such as [BBTools Repair](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/repair-guide/) that can help restore pairing information.

### Base Quality Scores

The symbols we see in the read quality string are an encoding of the quality score:

<img src="../images/base_qual_0.png" width="500">

A quality score is a prediction of the probability of an error in base calling: 

<img src="../images/base_qual_1.png" width="750">

Going back to our read, we can see that for most of our read the quality score is "G" –> "Q" =  38 -> Probability < 1/1000 of an error.

_________________________________________________________________________________________________________________________________________________________

[Next](lesson4.md)

[Previous](lesson2.md)
