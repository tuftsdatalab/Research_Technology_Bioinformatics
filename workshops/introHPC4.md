# Shell Scripting and Variables

## Variables and Modules

Now we are going to get to the fun part - scripting! Here we will download some NGS data using the SRA toolkit and then align these reads to a genome using STAR. So let's start by discussing what a variable is. A variable is a word that you assign value to:

```
[tutln01@lc1cmp047 data]$ testVar='/cluster/home/tutln01/introHPC/data'
[tutln01@lc1cmp047 data]$ echo $testVar
/cluster/home/tutln01/introHPC/data
```
Here we create a variable `testVar` and assign it to the directory `/cluster/home/tutln01/introHPC/data` we can then call that variable with `echo` by adding a `$` to the front of the word! Now to write a script often times you will need something called a **module** or software package 


Now let's create a script to download NGS data using variables:

```
[tutln01@lc1cmp047 data]$ nano download.sh
```
Now copy the following script:

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
- `module load sra/2.10.8` This loads the sratoolkit modu



## loops

what is a bash loop

## SLURM headers and how to submit batch jobs

slurm/ slurm headers 
sbatch
check job
cancel job
