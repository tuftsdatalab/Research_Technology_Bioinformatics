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


<img width="692" alt="Ondemand_Shell" src="https://user-images.githubusercontent.com/8632603/179539946-5d4fa52d-95ae-4215-ab16-24c912879aeb.png">


Mac
---

On Mac the shell is available through Terminal  
Applications -> Utilities -> Terminal  
Go ahead and drag the Terminal application to your Dock for easy access.

Windows
-------

For Windows, an easy one to install and use right away is  gitbash.  
Download and install `gitbash <https://gitforwindows.org/>`__;
Open up the program.

Other options: 

* https://docs.microsoft.com/en-us/windows/terminal/install
* https://conemu.github.io/

Linux
-----


You probably already know how to find the shell prompt.

### Starting with the shell
---------------------------

We will spend most of our time learning about the basics of the shell
by manipulating some experimental data.

Now we're going to download the data for the tutorial. For this you'll need
internet access, because you're going to get it off the web.

Open up the shell and type the command::

   whoami

and then hit ENTER 

(This is a good question for Mondays ....)


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


Key Takeaway
===============

`pwd` and `ls` are examples of commands - programs you run at the shell
prompt that do stuff. `pwd` stands for 'print working directory', while
`ls` stands for 'list files'. It is similar to the abbreviations used in texting, it takes less time to get the point across (lol, tbh, imho, afaik, ftw ....)

================


### Navigating in the Shell
--------------------------

We are going to make a place to work for this workshop.

The following command makes a new directory.

```

mkdir JulyWorkshop


```

===================

### Pro Tips:

Avoid spaces and special characters in names.

Spelling and Capitalization are literal in unix, be careful when making and using files to remember your convention.

====================

You can check that the new directory was created by repeating the list command.

```

ls

```

A directory is like a desk drawer. We create them to store files that relate to each other mostly.

When creating directories and filenames it is helpful to put some information about the project and the date of activity.


<img width="786" alt="File_Folder_Structure" src="https://user-images.githubusercontent.com/8632603/179539866-ecd6e880-f468-4151-bbaa-149f52c328b4.png">


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


Let's make a file here using a common command "echo" to start creating our file structure.


```

echo "Hello World " > helloworld.txt


```


The ">" in this command tells the command to place the output into the place it is pointing.


To open the file and check the contents, there are a few options. "cat" is a useful command for many reasons, let's see a demo here. Typing this command will print the contents of the file to the screen.

```
cat helloworld.txt
```

"cat" will open the entire file, so this is not the best command for long files.

In that case "head" is a good option. Head pulls the top ten lines of the file and prints them to the screen.

```
head helloworld.txt"

```

It does not look any different from cat in this case because there is only one line in the file.

A third way to check file contents is by using a program called "less" (or "more").

"less" will open the file interactively, then you can scroll through it and when you are done, push "q" on your keyboard to close the file.

```
less helloworld.txt
```

Press "q" to close the file.

There are many versions of these tools on command line, but "cat", "head" and "less" are very common.


### Absolute and Relative Paths
===============================



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


Let's go back into the JulyWorkshop directory, but this time use your ABSOLUTE path by changing *username01* to your username. If you forget your username, try *whoami*


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


<img width="523" alt="File_Permissions" src="https://user-images.githubusercontent.com/8632603/179539739-75f4edf9-5f5d-4de9-b20c-97abc7869be6.png">







Example for Class (Credit to https://angus.readthedocs.io/en/2019/)




## What is BLAST?
BLAST is the **B**asic **L**ocal **A**lignment **S**earch **T**ool.
It uses an index to rapdily search large sequence databases;
it starts by finding small matches between the two sequences and extending those matches.
For more information on how BLAST works and the different BLAST functionality,
check out the summary on [Wikipedia](https://en.wikipedia.org/wiki/BLAST) or
the NCBI's list of [BLAST resources](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs).

BLAST can be helpful for identifying the source of a sequence,
or finding a similar sequence in another organism.
In this lesson, we will use BLAST to find zebrafish proteins that
are similar to a small set of mouse proteins.

## Why use the command line?
BLAST has a very nice graphical interface for searching sequences in NCBI's database.
However, running BLAST through the commmand line has many benefits:
  * It's much easier to run many BLAST queries using the command line than the GUI
  * Running BLAST with the command line is reproducible and can be documented in a script
  * The results can be saved in a machine-readable format that can be analyzed later on
  * You can create your own databases to search rather than using NCBI's pre-built databases
  * It allows the queries to be automated
  * It allows you to use a remote computer to run the BLAST queries
  
Later on in the workshop we will talk more about these advantages and have a more in-depth explanation of the shell.

## Running BLAST

Many common programs are pre-loaded into the Tufts HPC using a system called "modules".

To see whether blast is available as a module, try running this command.

```

module av blast

```


As of July 2022, these are the modules you might see displayed.

<img width="711" alt="Blast_modules" src="https://user-images.githubusercontent.com/8632603/179539551-1d0c8933-30f2-43d5-957c-f4216d849ca6.png">


Choose the blast-plus version of the module and load it.

```

module load blast-plus/2.11.0

```

Confirm that the module is loaded.

```
module list

```





We need some data!  Let's grab the mouse and zebrafish RefSeq
protein data sets from NCBI, and put them in our home directory.

Now, we'll use `curl` to download the files from a Web site onto our
computer; note, these files originally came from the
[NCBI FTP site](ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot)

```
curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
curl -o mouse.2.protein.faa.gz -L https://osf.io/j2qxk/download
curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
```

If you look at the files in the current directory:

```
ls -l
```

You should now see these 3:

```
total 29908
-rw-rw-r-- 1 username01 username01 12553742 Jun 29 08:41 mouse.1.protein.faa.gz
-rw-rw-r-- 1 username01 username01  4074490 Jun 29 08:41 mouse.2.protein.faa.gz
-rw-rw-r-- 1 username01 username01 13963093 Jun 29 08:42 zebrafish.1.protein.faa.gz
```

The three files you just downloaded are the last three on the list - the
`.faa.gz` files.

All three of the files are FASTA protein files (that's what the .faa
suggests) that are compressed with `gzip` (that's what the .gz means).

Uncompress them:

```
gunzip *.faa.gz
```

and let's look at the first few sequences in the file:

```
head mouse.1.protein.faa 
```

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's pretty
ubiquitous.  It's a text file, containing records; each record
starts with a line beginning with a '>', and then contains one or more
lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with '>', which says "take
all the output and put it into this file here."

```
head -n 11 mouse.1.protein.faa > mm-first.faa
```

So now, for example, you can do `cat mm-first.faa` to see the contents of
that file (or `less mm-first.faa`). TIP: if you try `less mm-first.faa` you will need to exit by pressing the `q` key in your keyboard.

Now let's BLAST these two sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done
by calling 'makeblastdb':

```
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```

Next, we call BLAST to do the search:

```
blastp -query mm-first.faa -db zebrafish.1.protein.faa
```

This should run pretty quickly, but you're going to get a lot of output!!
To save it to a file instead of watching it go past on the screen,
ask BLAST to save the output to a file that we'll name `mm-first.x.zebrafish.txt`:

```
blastp -query mm-first.faa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt
```

and then you can 'page' through this file at your leisure by typing:

```
less mm-first.x.zebrafish.txt
```

(Type spacebar to move down, and 'q' to get out of paging mode.)

-----

Let's do some more sequences (this one will take a little longer to run):

```
head -n 498 mouse.1.protein.faa > mm-second.faa
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.txt
```

will compare the first 96 sequences.  You can look at the output file with:

```
less mm-second.x.zebrafish.txt
```

(and again, type 'q' to get out of paging mode.)

Notes:

* you can copy/paste multiple commands at a time, and they will execute in order;

* why did it take longer to BLAST ``mm-second.faa`` than ``mm-first.faa``?

Things to mention and discuss:

* `blastp` options and -help.
* command line options, more generally - why so many?
* automation rocks!

----

Last, but not least, let's generate a more machine-readable version of that
last file --

```
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6
```

You can open the file with `less mm-second.x.zebrafish.tsv` to see how the file looks like.

See [this link](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) for a description of the possible BLAST output formats.


