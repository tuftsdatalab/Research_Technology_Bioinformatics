# Searching and Redirection

## Searching Files

Often times it we will need to answer quick questions about our file like does it contain a certain pattern? We can leverage a powerful searching command called `grep`. Say for instance we want to know how many sequences are in our fastq file? Let's first review the anatomy of a fastq file:

<img src="../images/fastqFile.jpeg" />

So if we wanted to see how many sequences there were we need to count patterns that only occur when new sequence information is started. We will use the unique label for each sequence. In our data each sequence has a label that starts with `@SRRxxxxxxx.x`. So we will just see how many lines  start with `@SRR`:

```
[tutln01@c1cmp047 data]$ zcat SRR1552453.fastq.gz | grep ^@SRR -c
1000
```
Here we pipe our files content into `grep` and ask grep to find all instances of `@SRR` at the beginning of a line with `^` and then count those instances with `-c`. So we find that we have 1000 sequences in our file! When finding patterns you can specify more complex patterns with regular expressions:

|regular expression|description|
|-|-|
|```.``` | a single character.|
|```?``` | the preceding character matches 0 or 1 times only.|
|```*``` | the preceding character matches 0 or more times.|
|```+``` | the preceding character matches 1 or more times.|
|```{n}``` | the preceding character matches exactly n times.|
|```{n,m}``` | the preceding character matches at least n times and not more than m times.|
|```[agd]``` | the character is one of those included within the square brackets.|
|```[^agd]``` | the character is not one of those included within the square brackets.|
|```[c-f]``` | the dash within the square brackets operates as a range. In this case it means either the letters c, d, e or f.|
|```()``` | allows us to group several characters to behave as one.|
|```|``` | (pipe symbol) the logical OR operation.|
|```^``` | matches the beginning of the line.|
|```$``` | matches the end of the line.|

We won't go through every test case, but for a thorough review take a look at [Ryan's Tutorials](https://ryanstutorials.net/linuxtutorial/grep.php)

## Redirection

Now that we are a bit familiar with piping we should also mention another powerfull tool, redirection. We can use redirection to add content to files or overwrite them completely. Let's use our test file as an example:

```
[tutln01@lc1cmp047 data]$ cat testfile.txt
Hello World
[tutln01@lc1cmp047 data]$ echo Hello Again >> testfile.txt
[tutln01@lc1cmp047 data]$ cat testfile.txt
Hello World
Hello Again
```
Here we used the `echo` command to call some text (`Hello Again`) and append it to the file. Now if we used `>` we would completely overwrite the data:

```
[tutln01@lc1cmp047 data]$ cat testfile.txt
Hello World
Hello Again
[tutln01@lc1cmp047 data]$ echo New Content > testfile.txt
[tutln01@lc1cmp047 data]$ cat testfile.txt
New Content
```

_______________________________________________________________________________________________________________________________________________________

[Next](./introHPC4.md)

[Previous](./introHPC2.md)
