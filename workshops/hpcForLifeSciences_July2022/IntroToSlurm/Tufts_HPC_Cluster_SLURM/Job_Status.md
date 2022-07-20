### Job Status

- Checking your **active** jobs

```bash
[ymalon01@c1cmp044 LS]$ squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
          24063163     batch      job ymalon01  R       0:17      1 c1cmp044 
            
[ymalon01@c1cmp044 LS]$ squeue -u ymalon01
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
          24063163     batch      job ymalon01  R       0:54      1 c1cmp044 
```

To check your active jobs in the queue:

`$ squeue -u $USER` or `$ squeue -u your_utln`

To cancel a specific job:

`$ scancel JOBID`

To cancel all of your jobs:

`$ scancel -u $USER` or `$ scancel -u your_utln`

To check details of your active jobs (running "R" or pending "PD"):

`$ scontrol show jobid -dd JOBID`

```bash
[ymalon01@c1cmp044 LS]$ scontrol show jobid -dd 24063163
JobId=24063163 JobName=job
   UserId=ymalon01(31003) GroupId=ymalon01(5343) MCS_label=N/A
   Priority=12833 Nice=0 Account=normal QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=0 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   DerivedExitCode=0:0
   RunTime=00:01:31 TimeLimit=1-00:00:00 TimeMin=N/A
   SubmitTime=2022-07-20T12:33:14 EligibleTime=2022-07-20T12:33:14
   AccrueTime=2022-07-20T12:33:14
   StartTime=2022-07-20T12:33:15 EndTime=2022-07-21T12:33:15 Deadline=N/A
   PreemptEligibleTime=2022-07-20T12:33:15 PreemptTime=None
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2022-07-20T12:33:15
   Partition=batch AllocNode:Sid=c1cmp044:27677
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=c1cmp044
   BatchHost=c1cmp044
   NumNodes=1 NumCPUs=2 NumTasks=2 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=2,mem=8G,node=1,billing=2
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   JOB_GRES=(null)
     Nodes=c1cmp044 CPU_IDs=2-3 Mem=8192 GRES=
   MinCPUsNode=1 MinMemoryNode=8G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Reservation=bioworkshop
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/cluster/home/ymalon01/workshop/LS/sbatch.sh
   WorkDir=/cluster/home/ymalon01/workshop/LS
   StdErr=/cluster/home/ymalon01/workshop/LS/24063163.err
   StdIn=/dev/null
   StdOut=/cluster/home/ymalon01/workshop/LS/24063163.out
   Power=
   MailUser=ymalon01 MailType=NONE
```

- Checking your **finished** jobs

*You can no longer see these jobs in `squeue` command output.*

**Querying finished jobs helps users make better decisions on requesting resources for future jobs. **

Display job CPU and memory usage:

`$ seff JOBID`

```bash
[ymalon01@c1cmp044 LS]$ seff 24063163
Job ID: 24063163
Cluster: pax
User/Group: ymalon01/ymalon01
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 2
CPU Utilized: 00:03:00
CPU Efficiency: 50.00% of 00:06:00 core-walltime
Job Wall-clock time: 00:03:00
Memory Utilized: 54.16 MB
Memory Efficiency: 0.66% of 8.00 GB

```

Display job detailed accounting data:

`$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j JOBID`

```bash
[ymalon01@c1cmp044 LS]$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j 24063163
 Partition      State  Timelimit               Start                 End    Elapsed     MaxRSS     ReqMem  MaxVMSize   NNodes      NCPUS        NodeList 
---------- ---------- ---------- ------------------- ------------------- ---------- ---------- ---------- ---------- -------- ---------- --------------- 
     batch  COMPLETED 1-00:00:00 2022-07-20T12:33:15 2022-07-20T12:36:15   00:03:00                   8Gn                   1          2        c1cmp044 
            COMPLETED            2022-07-20T12:33:15 2022-07-20T12:36:15   00:03:00     55464K        8Gn    198364K        1          2        c1cmp044 
           OUT_OF_ME+            2022-07-20T12:33:15 2022-07-20T12:36:15   00:03:00          0        8Gn    108052K        1          2        c1cmp044 
```

NOTE: there are more format options, see [sacct](https://slurm.schedmd.com/sacct.html)

---
If you have any HPC related questions, please feel free to contact us at **tts-research@tufts.edu**.

*End of Day 1*


