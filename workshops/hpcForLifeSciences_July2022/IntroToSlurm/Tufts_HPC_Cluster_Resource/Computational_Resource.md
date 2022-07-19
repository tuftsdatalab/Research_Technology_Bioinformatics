
### CPUs
Resources are orgnized into **partitions** on the cluster based on functionality and priority.
After logging in on the HPC cluster, you can use command `sinfo` to check the `partition` you have access to (all partitions listed in the `sinfo` output).

```
[ymalon01@login-prod-01 ~]$ sinfo
PARTITION    AVAIL  TIMELIMIT  NODES  STATE NODELIST 
interactive     up    4:00:00      1    mix c1cmp064 
interactive     up    4:00:00      1   idle c1cmp063 
batch*          up 7-00:00:00      1  down* p1cmp005 
batch*          up 7-00:00:00      1  drain p1cmp056 
batch*          up 7-00:00:00     16   resv c1cmp[009,033,035-039,044-049],p1cmp[004,009,054] 
batch*          up 7-00:00:00     34    mix c1cmp[003-008,010-020,023-024,034,040-043,051-052,054],i2cmp001,p1cmp[003,012,015,018,020-021] 
batch*          up 7-00:00:00     17  alloc c1cmp[021-022,053],i2cmp003,p1cmp[001,006-008,010-011,013-014,019,022-024,055] 
batch*          up 7-00:00:00      2   idle p1cmp[016-017] 
mpi             up 7-00:00:00      1  down* p1cmp005 
mpi             up 7-00:00:00      1  drain p1cmp056 
mpi             up 7-00:00:00     16   resv c1cmp[009,033,035-039,044-049],p1cmp[004,009,054] 
mpi             up 7-00:00:00     34    mix c1cmp[003-008,010-020,023-024,034,040-043,051-052,054],i2cmp001,p1cmp[003,012,015,018,020-021] 
mpi             up 7-00:00:00     16  alloc c1cmp[021-022,053],p1cmp[001,006-008,010-011,013-014,019,022-024,055] 
mpi             up 7-00:00:00      2   idle p1cmp[016-017] 
gpu             up 7-00:00:00      1    mix p1cmp073 
gpu             up 7-00:00:00      2  alloc c1cmp[025-026] 
largemem        up 7-00:00:00      7    mix c1cmp[027-028,030,057,061-062],i2cmp055 
largemem        up 7-00:00:00      2  alloc p1cmp[049-050] 
largemem        up 7-00:00:00      3   idle c1cmp[032,058-059] 
preempt         up 7-00:00:00      2   mix$ p1cmp[094-095] 
preempt         up 7-00:00:00      4  maint p1cmp[090,092,103,109] 
preempt         up 7-00:00:00      1  down* p1cmp005 
preempt         up 7-00:00:00      2  drain p1cmp[038,056] 
preempt         up 7-00:00:00      3   resv p1cmp[004,009,054] 
preempt         up 7-00:00:00     71    mix cc1gpu[001-005],i2cmp[010-032,038-043,045-051],p1cmp[003,012,015,018,020-021,070-077,079-080,091,093,096,098-102,104-108,110] 
preempt         up 7-00:00:00     25  alloc c1cmp[025-026],i2cmp[004-006,008-009,033-035,037,052-053],p1cmp[006-008,010-011,013-014,019,022-024,055] 
preempt         up 7-00:00:00     20   idle p1cmp[016-017,031-037,039-042,081-086,097] 

```

*  [**OnDemand**](https://ondemand.pax.tufts.edu) `Misc`-->`Inventory ` shows more node details (core count & memory)

<img src="https://github.com/tuftsdatalab/Research_Technology_Bioinformatics/blob/3a85fd4bbc3b5b2ebe2f03321835aa5b834febef/workshops/hpcForLifeSciences_July2022/IntroToSlurm/images/Misc2.png" alt="Misc" width=60%>

<img src="https://github.com/tuftsdatalab/Research_Technology_Bioinformatics/blob/3a85fd4bbc3b5b2ebe2f03321835aa5b834febef/workshops/hpcForLifeSciences_July2022/IntroToSlurm/images/Inventory.png" alt="Inventory" width=60%>


### GPUs

__NVIDIA GPUs__ are available in `gpu` and `preempt` partitions

- Request GPU resources with `--gres`. See details below.
- Please **DO NOT** manually set `CUDA_VISIBLE_DEVICES`. 
- Users can ONLY see GPU devices that are assigned to them with `$ nvidia-smi`.
- `gpu` partition`-p gpu`:
  - NVIDIA P100s
    - In "gpu" partition
    - Request with: `--gres=gpu:p100:1`(one P100 GPU, can request up to 6 on one node)
  - NVIDIA Tesla K20xm
    - In "gpu" partition
    - Request with: `--gres=gpu:k20xm:1`(one Tesla K20xm GPU, can request up to 1 on one node)
- `preempt` partition `-p preempt`:
  - `a100`, `v100`, `p100`, ` rtx_6000`, `t4`
  - NVIDIA T4
    - In "preempt" partition
    - Request with: `--gres=gpu:t4:1`(one T4 GPU, can request up to 4 on one node)
  - NVIDIA P100
    - In "preempt" partition
    - Request with: `--gres=gpu:p100:1`(one P100 GPU, can request up to 4 on one node)
  - NVIDIA rtx_6000
    - In "preempt" partition
    - Request with: `--gres=gpu:rtx_6000:1`(one RTX_6000 GPU, can request up to 8 on one node)
  - NVIDIA V100
    - In "preempt" partition
    - Request with: `--gres=gpu:v100:1`(one V100 GPU, can request up to 4 on one node)
  - NVIDIA A100
    - In "preempt" partition
    - Request with: `--gres=gpu:a100:1`(one A100 GPU, can request up to 8 on one node)
