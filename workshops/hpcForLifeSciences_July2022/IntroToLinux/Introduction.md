# Intro to Linux
=======================================

* TOC
{:toc}

This short workshop provides some basic training on bash and shell scripting on the command line on the Linux-based Tufts HPC cluster.

This course is not meant to be comprehensive, but provides some insights into how the command line works as well as some strategic resources for studying and understanding command line on the HPC cluster.

#### Helpful Tip: Some Vocabulary
===================================

"Command line" is a more general term to indicate that you are using text commands on a terminal (linux bash shell or similar). Command line differs from "Graphical User Interface (GUI)" because all commands are texts instead of drag-and-drop or interactive formats such as the Windows or Mac Operating Sytems provide.

"HPC" stands for High Performance Computing, "cluster" refers to a shared computer resource to enable more powerful computation than regularly available on an individual machine.

"Linux" can refer to any of the free open source version of "Unix" from AT&T Bell labs who pioneered the language in 1965. There are a number of Linux operating systems installed on HPC clusters (Ubuntu, Debian, RedHat Enterprise License (RHEL), CentOs, Fedora, etc.) Each of these systems have slight differences that may impact the commands demoed here. Tufts University Research Cluster is currently using RHEL7.

"Bash" is one type of languages used in a "shell", the text interface on the Linux system. This lesson introduces a few objectives to help users understand how to use bash commands on the Linux RHEL shell of our HPC. Other shell languages have slight differences that affect how commands are run (e.g. new MacOSX ship with "zsh" as the default shell language on their installed terminal programs).


## Learning Objectives
-----------------------

- What is the shell?
- How do you access it?
- How do you use it and what is it good for?

  * Running commands
  * File Directory Structure
  * Manipulating files
  * Simple Bash Scripts

## What is the shell?
------------------

The **shell** is a program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard combination.

There are many reasons to learn about the shell.  A few specific ones:

* For most bioinformatics tools, you have to use the shell. There is no
  graphical interface. If you want to work in metagenomics or genomics you're
  going to need to use the shell.

* The shell gives you **power**. The command line gives you the power to
  do your work more efficiently and more quickly. Shell allows users to automate repetitive tasks.

* To use remote computers or cloud computing, you need to use the shell.


### Knowing Shell Increases Speed and Efficiency Through Automation

The most important reason to learn the shell is to learn about
**automation**.  Any time you find yourself doing roughly the same
computational task more than few times, it may be worth automating it;
the shell is often the best way to automate anything to do with files.

In this lesson, we're going to go through how to access Unix/Linux and some of the basic
shell commands. We will finish with a demonstration of how to run programs interactively as well by submitting a job to SLURM (https://it.tufts.edu/sites/default/files/uploaded-files/2020-03/QuickStart%20for%20Slurm.pdf). Slurm is a scalable cluster management and job scheduling system for Linux clusters. Other job scheduling systems you may be familiar with from other universities are "PBS" and "SGE_Batch".

### Where to learn shell commands
-----------------------------------------

The challenge with bash for the command line is that it's not particularly simple - it's a
power tool, with its own deep internal logic with lots of details.

Practice is the best way to learn, but here are some helpful shell command resources:

* [Fun With Unix Cheat Sheet](https://files.fosswire.com/2007/08/fwunixref.pdf)
* [Shell Cheatsheet - Software Carpentry](https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md)
* [Explain shell](http://explainshell.com) - a web site where you can see what the different
components of a shell command are doing.
