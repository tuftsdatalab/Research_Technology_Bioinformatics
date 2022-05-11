

## What is the Cluster?

Before getting to the cluster it is worth discussing what a cluster is and some of the terminology. First, let's discuss the difference between a CPU and a GPU.

### CPU -- Central Processing Unit
  - A CPU can never be fully replaced by a GPU
  - Can be thought of as the taskmaster of the entire system, coordinating a wide range of general-purpose computing tasks
### GPU -- Graphics Processing Unit
  - GPUs were originally designed to create images for computer graphics and video game consoles
  - Performing a narrower range of more specialized tasks

![](../images/cpuGpu.png)

You'll notice that in the picture above the CPU is composed of a smaller unit, a **core**. A core is the computing unit in a CPU. You'll also note that the whole system (including CPUs, GPUs and Storage) is a single computer in the system called a **node**.

![](../images/coreNode.png)

When a CPU performs some computation they use a storage hierarchy. This hierarchy places small/fast storage options close to the CPU and slower/larger options away from the CPU. These small/fast options are called **memory/RAM** while the slower/larger options are simply called **storage**.

![](../images/memStore.png)

Now that we now the components we can put together an image of what a computer cluster is. A **computer cluster** is a group of loosely or tightly connected computers that work together as a single system. A **HPC (High Performance Compute) cluster** is a computer cluster capable of performing computations at high speeds.

![](../images/hpcImage.png)



what is a cluster
getting to onDemand
what is linux

## Navigate to the Cluster

To get the Tufts HPC cluster you'll first need an account. If you haven't already done so please [request an account](https://tufts.qualtrics.com/jfe/form/SV_5bUmpFT0IXeyEfj) with Research Technology. You will also need to either be on a Tufts network or be connected to the [VPN](https://access.tufts.edu/vpn). Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu) and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

![](../images/ondemandLayout.png)



## Open an Interactive Session

different node architechtures
slurm

## Orienting Yourself in Terminal

ls,pwd,cd,clear
