# Shell Scripting and Variables

## Variables and Modules

Now we are going to get to the fun part - scripting! Here we will download some NGS data using the SRA toolkit and then align these reads to a genome using STAR. So let's start by discussing what a variable is. A variable is a word that you assign value to:

```
[tutln01@lc1cmp047 data]$ testVar='/cluster/home/tutln01/introHPC/data'
[tutln01@lc1cmp047 data]$ echo $testVar
/cluster/home/tutln01/introHPC/data
```
Here we create a variable `testVar` and assign it to the directory `/cluster/home/tutln01/introHPC/data` we can then call that variable with `echo` by adding a `$` to the front of the word! Now to write a script often times you will need something called a **module**, basically a software package. For our NGS download script we will need to do load the sratoolkit module. We do that by:

```
[tutln01@lc1cmp047 data]$ module av sra
-------------------------------------------------- /opt/shared/Modules/modulefiles-rhel6 ---------------------------------------------------
sra/2.10.8 sra/2.5.0  sra/2.9.2
[tutln01@lc1cmp047 data]$ module load sra/2.10.8
[tutln01@lc1cmp047 data]$ module list
sra/2.10.8 
[tutln01@lc1cmp047 data]$ module purge
```
Here we did a few things: 
- we checked which modules were available that matched the pattern `sra` with `module av sra`
- we loaded the newest sra module with `module load sra/2.10.8`
- we checked which modules were loaded with `module list`
- we purged all loaded modules with `module purge`

Now that we understand variables/modules let's create a script to download NGS data:

```
[tutln01@lc1cmp047 data]$ nano download.sh
```
Now copy the following script and paste in download.sh, then save the file with `Esc`+`Cntl`+`x`+`y`+`Enter`:

```
#!/bin/bash
#SBATCH --job-name=sraPull               # name your job
#SBATCH --time=03-00:00:00               # how long your job might take
#SBATCH --partition=preempt              # which partition you want to run it on
#SBATCH --nodes=1                        # how many nodes do you want
#SBATCH --mem=5Gb                        # how much memory do you want
#SBATCH --output=%j.out                  # name of output file
#SBATCH --error=%j.err                   # name of error file
#SBATCH --mail-type=ALL                  # request to be emailed when job begins and ends
#SBATCH --mail-user=YourEmail@tufts.edu  # provide your email

# load the sra module and configure it
module load sra/2.10.8 
vdb-config --interactive

# assign the acc variable and use it to download fastq data
acc='/cluster/home/tutln01/introHPC/sraAccList.txt'

# run the fastq-dump command
fastq-dump --split-e --gzip $(<$acc)
```

Here you will notice a few things:
- `#!/bin/bash`: this signifies that this is a batch script
- `#SBATCH` headers: these headers are going to be used by SLURM to dictate how your job is scheduled.
- `module load sra/2.10.8` This loads the sratoolkit module
- `vdb-config --interactive` configures this sratoolkit
- we assign the `$acc` variable to our accession list file with SRA accession numbers
- we run the `fastq-dump` command (with the `--split-e` option in case to pull paired data if applicable) to grab our fastq files 



## loops

what is a bash loop

## SLURM headers and how to submit batch jobs

slurm/ slurm headers 
sbatch
check job
cancel job
