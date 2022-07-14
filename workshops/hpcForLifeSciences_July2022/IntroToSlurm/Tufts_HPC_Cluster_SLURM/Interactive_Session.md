
### Interactive Session

- Particularly good for debugging and working with software GUI. 

  `$ srun [options] --pty [command]`

- Command 

  - command to run an application, given the module is already loaded.
  - `bash` for a bash shell

- Options

  - Pseudo terminal `--pty`
  - Partition `-p` 
    - Default batch if not specified
  - Time `-t` or `--time=`
    - Default 15 minutes if not specified on non-interactive partition
  - Number of CPU cores `-n` 
    - Default 1 if not specified
  - Memory `--mem=`
    - Default 2GB if not specified
  - GPU `--gres=`
    - Default none
  - X Window `--x11=first`
    - Default none	

  Starting an interactive session of bash shell on preempt partition with 2 CPU cores and 2GB of RAM, with X11 forwarding for 1 day, 2 hours, and 30 minutes (use `exit` to end session and release resources).

```bash
[ymalon01@login-prod-01 ~]$ srun -p preempt -t 1-2:30:00 -n 2 --mem=2g --x11=first --pty bash
srun: job 296794 queued and waiting for resources
srun: job 296794 has been allocated resources
[ymalon01@cc1gpu001 ~]$ 
```

