# Logging In

## Navigate to the Cluster

To get the Tufts HPC cluster you'll first need an account. If you haven't already done so please [request an account](https://tufts.qualtrics.com/jfe/form/SV_5bUmpFT0IXeyEfj) with Research Technology. You will also need to either be on a Tufts network or be connected to the [VPN](https://access.tufts.edu/vpn). Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu) and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

<img src="../images/ondemandLayout.png" height=250px/>


## Open an Interactive Session

Now that we are logged in let's use the cluster! To start click on `Clusters > Tufts HPC Shell Access`. You'll notice the following: 

<img src="../images/cli.png"/>

Where:
- tutln01 is your username
- @login-prod-01 is the node you are on

Now it is **IMPORTANT** to note that when you log in you are on the login node. This is a shared node, sort of like a waiting room. You can't run anything from this login node. For that you'll need to request compute resources so type and enter this into your terminal:

```
srun --time=0-24:00:00 --partition=batch -n 1 --mem=8Gb --reservation=bioworkshop --pty bash
```
So what did you do? Well you just used what is called a SLURM command. SLURM is what is known as a job scheduler and it is used to organize how jobs are run on the HPC. Let's break down what you did above:

- `srun` runs a parallel job on the cluster
- `--time` How long do we want to use this resource? The format is in day-hour:minute:second, so here we requested 0 day, 24 hours, 0 minutes and 0 seconds
- `--partition` identifies the partition you want to use - here we use the batch parition
- `-n` How many CPU cores do we want to use? Here we asked for 1
- `--mem` How much memory do we want to use? Here we asked for 8 gigabytes
- `--pty` What kind of terminal do we want? Here we asked for a bash terminal
- `--reservation` what reservation do we want to use? sometimes nodes are "reserved" for a special use case like this workshop. Here we ask for the `bioworkshop` reservation

Now you'll notice that the node has changed:

```
[tutln01@i2cmp008 ~]$ 
```

Now let's create a directory for our workshop with `mkdir` and change into it with `cd `: 

```
mkdir hpcDay2
cd hpcDay2
```

In our workshop directory we will set up a directory structure to organize our inputs/outputs:

```
mkdir data
mkdir results
mkdir scripts
```

_________________________________________________________________________________________________________________________________________________________

[Next](lesson3.md)

[Previous](lesson1.md)
