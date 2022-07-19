### Modules

- What are modules

  - A tool that **simplify** shell initialization and lets users easily modify their environment during the session with modulefiles
  - Each modulefile contains the **information** needed to configure the shell for an application. (PATH, LD_LIBRARY_PATH, CPATH, etc.)
  - Modules are useful in managing **different versions** of applications. 
  - Modules can also be bundled into metamodules that will load an entire **set of different applications (dependencies)**. 

  

  To check available modules installed on the cluster:

  ```
  [ymalon01@login-prod-01 ~]$ module av
  ```

  Upon login, environment `PATH` is set for the system to search executables:

  ```
  [ymalon01@login-prod-01 ~]$ echo $PATH
  /usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/ymalon01/bin:/cluster/home/ymalon01/.local/bin
  ```

  For example, I would like to use `gcc` compiler, to check what versions of gcc compiler is available, load the version I would like to use, and use it:

  ```
  [ymalon01@login-prod-01 ~]$ module av gcc
  
  ----------------------------------------------------------- /opt/shared/Modules/modulefiles-rhel6 ------------------------------------------------------------
  gcc/4.7.0 gcc/4.9.2 gcc/5.3.0 gcc/7.3.0
  
  -------------------------------------------------------------- /cluster/tufts/hpc/tools/module ---------------------------------------------------------------
  gcc/8.4.0 gcc/9.3.0
  ```

  ```
  [ymalon01@login-prod-01 ~]$ module load gcc/7.3.0
  [ymalon01@login-prod-01 ~]$ module list
  Currently Loaded Modulefiles:
    1) use.own     2) gcc/7.3.0
  ```

  ```
  [ymalon01@login-prod-01 ~]$ which gcc
  /opt/shared/gcc/7.3.0/bin/gcc
  [ymalon01@login-prod-01 ~]$ echo $PATH
  /opt/shared/gcc/7.3.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/ymalon01/bin:/cluster/home/ymalon01/.local/bin
  [ymalon01@login-prod-01 ~]$ gcc --version
  gcc (GCC) 7.3.0
  Copyright (C) 2017 Free Software Foundation, Inc.
  This is free software; see the source for copying conditions.  There is NO
  warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  ```

  I can swap a module for another (doesn't have to be the same software):

  ```
  [ymalon01@login-prod-01 ~]$ module swap gcc/7.3.0 gcc/9.3.0 
  [ymalon01@login-prod-01 ~]$ module list
  Currently Loaded Modulefiles:
    1) use.own     2) gcc/9.3.0
  ```

  I can also unload loaded modules:

  ```
  [ymalon01@login-prod-01 ~]$ module unload gcc
  [ymalon01@login-prod-01 ~]$ echo $PATH
  /usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/ymalon01/bin:/cluster/home/ymalon01/.local/bin
  ```

  I can unload ALL of the loaded modules:

  ```
  [ymalon01@login-prod-01 ~]$ module purge
  ```

  

- Install Software/Packages

  - [R](https://tufts.box.com/s/qximkv5ke2y4k0vbg6m04m6fc6exh88h) (R command line recommanded)
    - R/4.0.0
    - gcc 
    - gdal
    - curl
  - [Python](https://tufts.box.com/v/CondaEnvonHPC) (Conda env recommanded)
    - anaconda/3 (older version, source activate)
    - anaconda/2021.05 (newer version, source activate)
    - Use the same version of conda on one conda env every time
  - Other software compiled from source
    - gcc
    - cmake
    - ... any dependencies, load if available, install if not.
    - Follow instructions (read it through)
    - Use "--prefix=" to install in non-standard locations
    - Modify the environment variables !!! (such as PATH, LD_LIBRARY_PATH, CPATH, .etc)


---

