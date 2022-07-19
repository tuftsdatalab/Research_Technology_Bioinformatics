> Hostname: `xfer.cluster.tufts.edu`
>
> Protocol: SCP or SFTP
>
> Use port 22 for SFTP
> 
**Hostname for file transfer: xfer.cluster.tufts.edu**

    > NOTE:
    >
    > Local_Path is the path to your files or directory on your local computer
    > Cluster_Path is the path to your files or directory on the cluster
    > Home Directory: /cluster/home/your_utln/your_folder
    > Research Project Storage Space Directory: /cluster/tufts/yourlabname/your_utln/your_folder

---
### Command Line

***Execute from your local machine terminal.***

- General Format:

    `$ scp From_Path To_Path`

    `$ rsync From_Path To_Pat`

***NOTE: If you are transfering very large files that could take hours to finish, we would suggest using `rsync` as it has ability to restart from where it left if interrupted.***

**File** Transfer with `scp`or `rsync`:

- Download from cluster

    `$ scp your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path  `

    `$ rsync your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path`

- Upload to cluster

    `$ scp Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path `

    `$ rsync Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

**Directory** Transfer with `scp` or `rsync`:

- Download from cluster

    `$ scp -r your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path  `

    `$ rsync -azP your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path`

- Upload to cluster

    `$ scp -r Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

    `$ rsync -azP Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`
---    
### OnDemand

***Only for transfering files size less than 976MB per file.***

Go to **[OnDemand]( https://ondemand.pax.tufts.edu/)** 

Under **`Files`**, using the **`Upload`** or **`Download`** buttons to transfer. 

# Add Screenshots

---
### File Transfer Client

> Hostname: `xfer.cluster.tufts.edu`
>
> Protocol: SCP or SFTP
>
> Use port 22 for SFTP

- **[WinSCP](https://winscp.net/eng/index.php)**<img src="https://miro.medium.com/max/500/1*Of7JOwV0wZgDIjgaS4qKlQ.png" alt="WinSCP" width=20%>

- **[FileZilla](https://filezilla-project.org/)**  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/0/01/FileZilla_logo.svg/1200px-FileZilla_logo.svg.png" alt="FileZilla" width=10%>

- **[Cyberduck<img src="https://cdn.cyberduck.io/img/cyberduck-icon-384.png" alt="CyberDuck" width=10%>](https://cyberduck.io/)**

  ---

