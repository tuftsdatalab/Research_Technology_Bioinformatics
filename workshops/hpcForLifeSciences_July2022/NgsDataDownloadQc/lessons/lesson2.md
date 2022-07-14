# Obtain Public NGS Data

The [US National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/) (NCBI) hosts repositories for many types of biomedical and genomics data. Today we'll retrieve reference data from the [Genomes Database FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/) as well as the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). 

We will use the `wget` command to download our reference data from the NCBI repository:

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz
```
You will notice that both these files end in `.gz` - this indicates they are compressed. Compressing files is useful when storing data but to use it you will often need to decompress it with `gunzip -d`:

```
gunzip -d https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
gunzip -d https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz
```

You will also notice that one file ends in `.fna` and the other ends in `.gff`. These are Fasta and GFF files, respectively.

### Fasta Format
The virus genome is in fasta format. Fasta format has two parts, a sequence identifier preceeded by a ">" symbol, followed by the sequence on subsequent lines.

<p align="center">
<img src="../images/genome_view_2,PNG.png" width="900">
</p>


### GFF Format
The gene annotation file is in Generic Feature Format (GFF). This formet tells us where genes are located in the reference genome.
Note that we must always be sure that our gene information and genome come from the same source. 

<p align="center">
<img src="../images/gff.PNG">
</p>
