### Batch Job ###

Write a batch submission script e.g. **myjob.sh**

```bash
#!/bin/sh
#SBATCH -J My_Job_Name   #job name
#SBATCH --time=00-00:20:00  #requested time (DD-HH:MM:SS)
#SBATCH -p preempt    #running on "preempt" partition/queue
#SBATCH -N 1    #1 nodes #for many shared-memory programs,please leave -N as 1.
#SBATCH -n 2   #2 tasks total and 1 cpu per task, that gives you 2 cpu cores for this job
#SBATCH --mem=2g  #requesting 2GB of RAM total for the number of cpus you requested
#SBATCH --output=MyJob.%j.%N.out  #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=MyJob.%j.%N.err   #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL    #email optitions
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]
# for example, running a python script 
# load the module so the correct version python is available to you
module load anaconda/2021.05
# If you have a conda env that you would like to use, activate it here using "source activate xxx"
# run python script
python myscript.py #make sure myscript.py exists in the current directory
# make sure you save all plots, data, outputs generated to files in your script
# Don't forget to deactivate your conda env if you are using one
```

**Submit** the job using the following command from command line interface:

`$ sbatch myjob.sh`

**Sample Scripts** including R, conda, matlab, gaussian

`/cluster/tufts/hpc/tools/slurm_scripts`

**Useful Link**

https://tufts.box.com/v/HPC-New-User

