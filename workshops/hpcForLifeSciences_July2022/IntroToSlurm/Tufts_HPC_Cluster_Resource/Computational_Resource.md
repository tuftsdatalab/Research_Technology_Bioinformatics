
### CPUs

```
PARTITION       TIMELIMIT      
batch*          7-00:00:00          
gpu             7-00:00:00        
interactive     4:00:00        
largemem        7-00:00:00        
mpi             7-00:00:00         
preempt         7-00:00:00     
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
