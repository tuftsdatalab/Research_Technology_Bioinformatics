# Manipulating Data

So now that we have downloaded and inspected our data we can get to manipulating it! So to start, let's talk about accessing parts of your data. To grab the first column in a data frame/matrix you can do so like:

    iris[,1]

    5.1 4.9 4.7 4.6 5.0 5.4 ...

To grab the first row:

    iris[1,]

    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          5.1         3.5          1.4         0.2  setosa
    

Now if your data is a data frame you have a special way of accessing coluns with the ```$``` operator:

    iris$Sepal.Length

    5.1 4.9 4.7 4.6 5.0 ...

This comes in handy for readability. While you can grab your data by column number, it is much easier to read that you are grabbing Sepal Length. To grab mulitple columns/rows, you can do the following for both data frames and matrices:

    iris[,c(1,3,4)]``` (grabbing the first, third, and fourth columns):

    Sepal.Length Petal.Length Petal.Width
    1          5.1          1.4         0.2
    2          4.9          1.4         0.2
    3          4.7          1.3         0.2
    4          4.6          1.5         0.2
    5          5.0          1.4         0.2
    6          5.4          1.7         0.4

In a data frame, to access columns you can be more specific and specify by column name:

    iris[c("Petal.Length","Species")]

    Petal.Length Species
    1          1.4  setosa
    2          1.4  setosa
    3          1.3  setosa
    4          1.5  setosa
    5          1.4  setosa
    6          1.7  setosa

## Subsetting Data

To subset our data we need to know a little bit about the different logical operators:

| Operator | Description |
:-------|:-----|
| > | greater than | 
| >= | greater than or equal |
| < | less than |
| <= | less than or equal |
| == | equals | 
| != | not equal |
| & | and |
| \| | or|

Let's go through a few of these!

Subsetting so that we only have rows where the Petal Length is greater than 1:

    iris[iris$Petal.Length > 1,]

Subsetting so that we only have rows where Petal Length is less than 1.5:

    iris[iris$Petal.Length < 1.5,]

Subsetting so that we only have rows where the Species is setosa:

    iris[iris$Species == "setosa",]

Subsetting so that we only have rows where the Species is not setosa:

    iris[iris$Species != "setosa",]

Subsetting so that we only have rows where the Species is setosa or versicolor:

    iris[iris$Species != "setosa" | iris$Species != "versicolor" ,]

## Using Dplyr

When subsetting data we should also mention the R package dplyr. This package has functionality to neatly modify data frames. Let's go through a quick example:

    filtered <- iris %>%
        filter(Species == "setosa") %>%
        select(-Petal.Length) %>%
        mutate(SpeciesWithPetalWidth = paste(Species,Petal.Width,sep="")) %>%
        arrange(Petal.Width)

Here we go through a few commonly used functions:

  * ```filter()``` will filter by some logical experssion
  * ```select()``` will select a column - here we remove a column by saying ```select(-Petal.Length)```
  * ```mutate()``` will create a new column
  * ```arrange()``` will arrange your rows by some variable - here we choose ```Petal.Width```
    
    > NOTE:  For more dplyr data wrangling tips check out the [Data Wrangling with dplyr and tidyr Cheat Sheet](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
 
 
## References:

1. [HBC Training](https://hbctraining.github.io/Intro-to-R-flipped/lessons/05_introR-data-wrangling.html)
2. [Dplyr](https://dplyr.tidyverse.org/)

_________________________________________________________________________________________________________________________________________________________________________________


[Back To Introduction to R](../IntroToR.md)

[Back To The Main Page](../../index.md)
