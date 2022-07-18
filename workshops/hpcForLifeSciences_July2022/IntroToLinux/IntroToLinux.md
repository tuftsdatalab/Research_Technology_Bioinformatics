# The Shell
=========


## Objectives
----------

- What is the shell?
- How do you access it?
- How do you use it and what is it good for?

  * Running commands
  * File Directory Structure
  * Manipulating files
  * Simple Bash Scripts

What is the shell?
------------------

The *shell* is a program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard combination.

There are many reasons to learn about the shell.  A few specific ones:

* For most bioinformatics tools, you have to use the shell. There is no
  graphical interface. If you want to work in metagenomics or genomics you're
  going to need to use the shell.

* The shell gives you *power*. The command line gives you the power to
  do your work more efficiently and more quickly. Shell allows users to automate repetitive tasks.

* To use remote computers or cloud computing, you need to use the shell.


Automation


The most important reason to learn the shell is to learn about
**automation**.  Any time you find yourself doing roughly the same
computational task more than few times, it may be worth automating it;
the shell is often the best way to automate anything to do with files.

Today we're going to go through how to access Unix/Linux and some of the basic
shell commands.

### Information on the shell
------------------------

The challenge with UNIX is that it's not particularly simple - it's a
power tool, with its own deep internal logic with lots of details.
The joke is that Unix is user-friendly - it's just very selective
about who its friends are!

shell cheat sheets:

* https://files.fosswire.com/2007/08/fwunixref.pdf
* https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md

Explain shell - a web site where you can see what the different
components of a shell command are doing.

* http://explainshell.com


(We'll look at these at the bottom.)

## How to access the shell
--------------------------

https://ondemand.pax.tufts.edu/pun/sys/dashboard

put the image here


Mac
~~~

On Mac the shell is available through Terminal  
Applications -> Utilities -> Terminal  
Go ahead and drag the Terminal application to your Dock for easy access.

Windows
~~~~~~~

For Windows, we're going to be using gitbash.  
Download and install `gitbash <https://gitforwindows.org/>`__;
Open up the program.

Other options: 

* https://docs.microsoft.com/en-us/windows/terminal/install
* https://conemu.github.io/

Linux


You probably already know how to find the shell prompt.

### Starting with the shell
-----------------------

We will spend most of our time learning about the basics of the shell
by manipulating some experimental data.

Now we're going to download the data for the tutorial. For this you'll need
internet access, because you're going to get it off the web.

Open up the shell and type the command::

   whoami

and then hit ENTER 

This is a good question for Mondays ....


### Running Commands
--------------------------

Let's try another.

Much like text shortcuts, shell commands often use abbreviations to get their point across.

For example, the command *pwd* is short for "pass working directory."

Now type the command

```

   pwd

```

You should see something similar to this:

```

/cluster/home/username01/

```

Try this command


```

ls


```

It may be empty for the moment, let's circle back to the command shortly.

===============
Key Takeaway
===============

'pwd' and 'ls' are examples of commands - programs you run at the shell
prompt that do stuff. pwd stands for 'print working directory', while
'ls' stands for 'list files'.

================


### Navigating in the Shell
--------------------------

We are going to make a place to work for this workshop.

The following command makes a new directory.

```

mkdir JulyWorkshop


```

===================

Smart Tips:

Avoid spaces and special characters in names.

Spelling and Capitalization are literal in unix, be careful when making and using files to remember your convention.

====================

You can check that the new directory was created by repeating the list command.

```

ls

```

A directory is like a desk drawer. We create them to store files that relate to each other mostly.

When creating directories and filenames it is helpful to put some information about the project and the date of activity.

Let's go into our directory and look around.

Another command you'll find yourself using a lot is 'cd', which stands
for 'change directory'.  Try typing::

```

cd JulyWorkshop


```


and then


```

   pwd

```

You should now see something like this:

```

/cluster/home/username01/JulyWorkshop

```


Let's make an empty file here to start creating our file structure.


```

touch emptyfile.txt


```


This is an empty file, to demonstrate file structure.


```

ls


```


You should just see the file name.

```

emptyfile.txt


```


If you want to go back to the directory that is in the level above our current file, another common shortcut used in bahs is `..`.


```

cd ..

```

`..` is a reference to a **RELATIVE PATH**


```

ls JulyWorkshop


```


You should just see the file name.



```

emptyfile.txt


```

A Relative Path means that the command only works from the relative location that you are in.

This can get confusing if you are moving around a lot in your directories or sending commands to SLURM, so the alternative method to navigating around the cluster is using an *ABSOLUTE PATH*.


Let's go back into the JulyWorkshop directory, but this time use your ABSOLUTE path by changing *username01** to your username. If you forget your username, try *whoami*


```

cd /cluster/home/username01/JulyWorkshop


```


Many commands in bash can be used with the ABSOLUTE PATH.

```

ls /cluster/home/username01/JulyWorkshop


```


This is helpful for checking for outputs from SLURM jobs when they are running.


### Going Home


Sometimes we get lost, so it is useful to know a few ways to get back to where you started.

```

cd


```

This command returns you to your home directory.

```

pwd

```

Other options for this command are

```

cd ~


```


```

cd $HOME

```



### Parameters for Bash commands

Many bash commands have special **parameters**, sometimes referred to as **flags** that open up a lot more possibilities.

Let's start by going to your home directory (you choose the command)


As you start using bash more and more, you will find a mix of files and directories/folders. If we want to know which is which, we can type::

```

    ls -F

```

Anything with a "/" after it is a directory.  Things with a "*" after
them are programs.  It there's nothing there it's an otherwise
unremarkable file (e.g. a data file).

Depending on which terminal you are using, some of the file types may have different colors. 

In our ondemand shell:

Files are white
Directories are blue
Programs are green
Compressed files are red (e.g. files that end in .zip or .gzip or .tar)


You can also use the command::

    ls -l

to see whether items in a directory are files or directories. `ls -l`
gives a lot more information too, such as the size of the file.


It also shows the permissions of who can read, write or execute a file.


```

drwxrwx--- 2 username05 username05     4096 Jul 18 09:57 JulyWorkshop


```

The first 10 letters in this line indicates the permission settings.









Example for Class


