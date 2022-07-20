### Batch Job ###

Write a batch submission script e.g. **sbatch.sh**

```bash
#!/bin/bash
#SBATCH --job-name=job            # job name is "job"
#SBATCH --nodes=1                 # 1 nodes #for many shared-memory programs,please leave -N as 1.
#SBATCH -n 2                      # 2 tasks total and 1 cpu per task, that gives you 2 cpu cores for this job
#SBATCH --partition=batch         # running on "batch" partition/queue
#SBATCH --reservation=bioworkshop # running on a reservation, named "bioworkshop", if no access to reservation, omit this line
#SBATCH --mem=8Gb                 # requesting 8GB of RAM total for the number of cpus you requested
#SBATCH --time=0-24:00:00         # requested time (DD-HH:MM:SS) 24 hours
#SBATCH --output=%j.out           # saving standard output to file, %j=JOBID
#SBATCH --error=%j.err            # saving standard error to file, %j=JOBID
#SBATCH --mail-type=ALL           # email optitions
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu  # use your own Tufts email address

# [this is a comment]
# The order of the "#SBATCH" options doesn't matter

#[commands_you_would_like_to_exe_on_the_compute_nodes]
# for example, running blast 
# load the module so the correct version of blast is available to you

module load blast-plus/2.11.0

# running blast
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6

```

**Submit** the job using the following command from command line interface:

`$ sbatch sbatch.sh`

**Sample Scripts** including R, conda, matlab, gaussian, .etc

`/cluster/tufts/hpc/tools/slurm_scripts`

**Useful Link**

https://tufts.box.com/v/HPC-New-User

---

NEXT - [Job Status](Job_Status.md)
