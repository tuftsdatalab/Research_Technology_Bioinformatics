### Job Status

- Checking your **active** jobs

```bash
[ymalon01@login-prod-01 ~]$ squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
            296794   preempt     bash ymalon01  R       5:12      1 cc1gpu001 
            
[ymalon01@login-prod-01 ~]$ squeue -u ymalon01
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
            296794   preempt     bash ymalon01  R       5:21      1 cc1gpu001 
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
[ymalon01@cc1gpu001 ~]$ scontrol show jobid -dd 296794
JobId=296794 JobName=bash
   UserId=ymalon01(31003) GroupId=ymalon01(5343) MCS_label=N/A
   Priority=10833 Nice=0 Account=(null) QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=0 Restarts=0 BatchFlag=0 Reboot=0 ExitCode=0:0
   DerivedExitCode=0:0
   RunTime=00:10:33 TimeLimit=1-02:30:00 TimeMin=N/A
   SubmitTime=2021-03-22T22:18:50 EligibleTime=2021-03-22T22:18:50
   AccrueTime=2021-03-22T22:18:50
   StartTime=2021-03-22T22:18:55 EndTime=2021-03-24T00:48:55 Deadline=N/A
   PreemptEligibleTime=2021-03-22T22:18:55 PreemptTime=None
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2021-03-22T22:18:55
   Partition=preempt AllocNode:Sid=login-prod-01:34458
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=cc1gpu001
   BatchHost=cc1gpu001
   NumNodes=1 NumCPUs=2 NumTasks=2 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=2,mem=2G,node=1,billing=2
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   JOB_GRES=(null)
     Nodes=cc1gpu001 CPU_IDs=30-31 Mem=2048 GRES=
   MinCPUsNode=1 MinMemoryNode=2G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=bash
   WorkDir=/cluster/home/ymalon01
   Power=
   MailUser=ymalon01 MailType=NONE
```

- Checking your **finished** jobs

*You can no longer see these jobs in `squeue` command output.*

**Querying finished jobs helps users make better decisions on requesting resources for future jobs. **

Display job CPU and memory usage:

`$ seff JOBID`

```bash
[ymalon01@login-prod-01 ~]$ seff 296794
Job ID: 296794
Cluster: pax
Use of uninitialized value $user in concatenation (.) or string at /usr/bin/seff line 154, <DATA> line 602.
User/Group: /ymalon01
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 2
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:22:12 core-walltime
Job Wall-clock time: 00:11:06
Memory Utilized: 1.16 MB (estimated maximum)
Memory Efficiency: 0.06% of 2.00 GB (2.00 GB/node)
```

Display job detailed accounting data:

`$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j JOBID`

```bash
[ymalon01@login-prod-01 ~]$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j  296794
 Partition      State  Timelimit               Start                 End    Elapsed     MaxRSS     ReqMem  MaxVMSize   NNodes      NCPUS        NodeList 
---------- ---------- ---------- ------------------- ------------------- ---------- ---------- ---------- ---------- -------- ---------- --------------- 
   preempt  COMPLETED 1-02:30:00 2021-03-22T22:18:55 2021-03-22T22:30:01   00:11:06                   2Gn                   1          2       cc1gpu001 
           OUT_OF_ME+            2021-03-22T22:18:55 2021-03-22T22:30:01   00:11:06         8K        2Gn    135100K        1          2       cc1gpu001 
            COMPLETED            2021-03-22T22:18:56 2021-03-22T22:30:01   00:11:05       592K        2Gn    351672K        1          2       cc1gpu001 
```

NOTE: there are more format options, see [sacct](https://slurm.schedmd.com/sacct.html)

