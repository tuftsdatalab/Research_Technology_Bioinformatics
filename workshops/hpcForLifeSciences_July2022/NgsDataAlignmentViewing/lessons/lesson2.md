Approximate time: 10 minutes

## Learning Objectives
- View alignment in IGV

## BAM Visualization with IGV

1. With a Chrome web browser, visit [https://ondemand.pax.tufts.edu/](https://ondemand.pax.tufts.edu/)
2. Login with your Tufts credentials
3. On the top grey menu bar, choose `Interactive Apps->IGV`.
<img src="../images/igv_dropdown.pdf" width="300">

4. Choose the following compute resource parameters: Number of house = 1, Number of cores = 1s, Amount of Memory = 8GB, Resevation = Bioinformatics Workshop, No Directory (leave blank)

<img src="../images/igv_params.png" width="500">

5. Click the blue button `Launch IGV` when it appears

After this the IGV window will appear, probably as a small window on a grey background.
Click the square icon in the top right corner to maximize the window.

<img src="../images/igv_start.pdf" width="500">

### Load reference genome 

1. Choose reference genome by clicking the `Genomes` menu and selecting `Load Genome from File...`

<img src="../images/igv_ref_1.png" width="300">

2. Navigate from your home directory to `hpcDay2/ref/` and select the fasta file `GCF_009858895.2_ASM985889v3_genomic.fna`

<img src="../images/igv_ref_2.png" width="600">

4. Click `Open`

You will see the name of the file and name of the sequence populate.


### Load the GFF file

1. Choose the GFF file by clicking on `File` menu and selecting `Load from File...`
<img src="../images/igv_from_file.png" width="300">

2. As before, navigate to `hpcDay2/ref/` and select the GFF file `GCF_009858895.2_ASM985889v3_genomic.gff` and click `Open`

### Load the BAM file

1. Choose the BAM file by clicking on `File` menu and selecting `Load from File...`

2. Navigate to `hpcDay2/results/` and select the sorted BAM file `sarscov2.srt.bam`

3. Click `Open`

4. It will take a minute or so to load all the reads. You can view progress in the lower right hand corner.


5. When it's done, you will have the following view. Each row of data is called a track. There are four tracks visible: the top track shows the reference genome coordinates, followed by two tracks of our alignment (coverage and reads) followed by the GFF track showing the gene locations on our refererence genome.
<img src="../images/igv_result.pdf" width="800">

### Examining a Gene

1. Zoom in on the Gene that encodes the Spike protein by hovering with your mouse on the genome track and clicking and dragging over the portion that is directly above the "S" protein.
2. We should see the coordinate box show roughly basepairs 22,000-25,000 `NC_045512.2:22,000-25,000`.

<img src="../images/igv_coordinates.png" width="500">

3. Looking at the Coverage track we see there are five positions that are colored, which indicates a basepair mismatch with respect to the references sequence in over 20% of the reads a that position. we will examine one more closely.


### Examining a Variant

1. This region contains one of the 4 mutations that differentiate the delta variant from the originally characterized sequence. 

2. Select region around the variant at 22,995.

3. View  amino acid change is a T>K change at protein position 478, which corresponds to a C>A SNP at nucleotide position 22,995. 

2. Explain hovering on a variant.

<img src="../images/igv_variant_info.png" width="200">

## Summary of Workshop

[Previous: Alignment ](lesson1.md)
[Home ](lesson1.md)

