# Background

## SARS-CoV-2 Spike Protein

- SARS-CoV-2 is a virus with an RNA genome, contained within an envelope that is studded with spike protiens
- Spike proteins bind with host cell receptor angiotensin-converting enzyme 2 (ACE2)129 to enter cells
- As such spike protein is under selective pressure to change the receptor binding domain in order to more readily infect cells

![](../images/spike.png)

<details>
<summary><b>More Info</b></summary>
<br>
Severe acute respiratory syndrome-coronavirus (SARS-CoV) has an RNA genome, which encodes proteins necessary for it’s replication inside cells. The genome is contained in an envelope, which is studded with envelope and Spike proteins. It enters host cells when the receptor binding domain of the Spike protein binds to the host cell receptor angiotensin-converting enzyme 2 (ACE2)129, and the ACE2–virus complex enters the cell where downstream steps in the virus lifecycle take place. Viruses have a high mutation rate, and not surprisingly the Spike protein is under selective pressure to change the receptor binding domain in order to more readily infect cells.
</details>

## Variants of Concern (VOC)

- Mutations in the spike protein have been noted to result in higher ACE2 affinity/transmissibility
- Here we note five residues and in the spike protein and how they've changed over different SARS-CoV-2 variants
- Today we are studying the delta variant

![](../images/voc.png)

<details>
<summary><b>More Info</b></summary>
<br>
There are several varaints of the originally characterized SARS-CoV-2 sequence. They contain 5 major mutations with respect to the original sequence, shown here. For example take T478, which is spike protein residue number 478 which is threonine (T) in the original sequence, and changes to lysine (K) in the delta variant. These mutations in the spike proteins have been shown to result in higher ACE2 affinity, transmissibility among other phenotypes that make them an increased thread to public health. Today, we’ll study a delta variant sample and compare these positions to the originally characterized sequence using bioinformatics methods.
</details>

## SARS-CoV-2 Resources

- The US National Center for Biotechnology Information (NCBI) hosts repositories for many types of biomedical and genomics data. 
- We will use it to download our SARS-CoV-2 reference genome and the raw next generation sequencing data

![](../images/sars_resources.png)

<details>
<summary><b>More Info</b></summary>
<br>
The US National Center for Biotechnology Information hosts repositories for many types of biomedical and genomics data. Today we'll retrieve reference data from the nucleotide repository, which contains among other things sequence records for SARS-CoV-2 genomes as well as the Sequence Read Archive (SRA) where raw next generation sequencing data can be easily obtained for reanalysis.
</details>

## Viral Genome Next Generation Sequencing (NGS)

- RNA is extracted from the sample, contains virus and host RNA
- Virus specific primers are used to capture SARS-CoV-2 and not host RNA
- Transcribed into complementary DNA, and amplified
- Adapters are added to attach to flowcell
- Sequence of the fragments is determined using fluorescently tagged nucleotides 
  - They emit a colored light signal when they attach to a base on the fragment

![](../images/viral_ngs.png)

<details>
<summary><b>More Info</b></summary>
<br>
Let's take a moment ot go through the Viral NGS Workflow. Total RNA is extracted from the sample, this contains both virus and host RNA. At this point a decision is made about downstream protocol based on whether the virus is known. In the case we will study, the authors are interested in sequencing SARS-CoV-2 which has a known sequence, so they can design virus specific primers which will capture SARS-CoV2 and not the host RNA, and transcribe it into complementary DNA which can be prepped for sequencing. Importantly, the primers they design will limit what they capture, so they must design the primers to bind in regions that they expect are conserved between viral variants. An amplification procedure has to be performed, because the current NGS technologies require a high input DNA amount and the viral genome amount is several orders of magnitude lower. Then, the preparation for NGS sequencing begins. We start with viral cRNA fragments that can be sequenced the same as DNA on an NGS sequencer. Sequencing adapters are added (these are the blue rectangles) that will allow the cDNA fragments to attach to the sequencing flowcell. The sequence of the fragments is determined using fluorescently tagged nucleotides that emit a colored light signal when they attach to a base on the fragment being sequences. It is the optical signal that the instrument reads to determine which base in the sequence was read. 
</details>

## Types of Read Data

- We often have two kinds of data, single-end and paired-end
- **single-end** sequence each DNA fragement from one end only
- **paired-end** sequence each DNA fragement from both sides
  - paired-end data is useful when sequencing highly repetitive sequences. Today we will be working with paired end data

![](../images/single_paired.png)

_______________________________________________________________________________________________________________________________________________

[Next](lesson2.md)

[Previous](../README.md)
