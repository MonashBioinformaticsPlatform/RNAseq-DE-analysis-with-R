---
layout: page
title: Programming with R
subtitle: Instructor's Guide
minutes: 0
---



## Legend

We are using a dataset with records on inflammation from patients following an
arthritis treatment. With it we explain `R` data structure, basic data
manipulation and plotting, writing functions and loops.

## Overall

This lesson is written as an introduction to R, but its real purpose is to
introduce the single most important idea in programming: how to solve problems
by building functions, each of which can fit in a programmer's working memory.
In order to teach that, we must teach people a little about the mechanics of
manipulating data with lists and file I/O so that their functions can do things
they actually care about.
Our teaching order tries to show practical uses of every idea as soon as it is
introduced; instructors should resist the temptation to explain the "other 90%"
of the language as well.

The secondary goal of this lesson is to give them a usable mental model of how
programs run (what computer science educators call a
[notional machine](reference.html#notional-machine) so that they can debug
things when they go wrong.
In particular, they must understand how function call stacks work.

The final example asks them to build a command-line tool that works with the
Unix pipe-and-filter model.
We do this because it is a useful skill and because it helps learners see that
the software they use isn't magical.
Tools like `grep` might be more sophisticated than the programs our learners can
write at this point in their careers, but it's crucial they realize this is a
difference of scale rather than kind.

The `R` novice inflammation contains a lot of material to cover.
Remember this lesson does not spend a lot of time on data types, data
structure, etc.
It is also on par with the similar lesson on Python.
The objective is to explain modular programming with the concepts of functions,
loops, flow control, and defensive programming (i.e. SWC best practices).
Supplementary material is available for R specifics
([Addressing Data](01-supp-addressing-data.html),
[Data Types and Structure](01-supp-data-structures.html),
[Understanding Factors](01-supp-factors.html),
[Introduction to RStudio](01-supp-intro-rstudio.html),
[Reading and Writing .csv](01-supp-read-write-csv.html),
[Loops in R](03-supp-loops-in-depth.html),
[Best Practices for Using R and Designing Programs](06-best-practices-R.html),
[Dynamic Reports with knitr](07-knitr-R.html),
[Making Packages in R](08-making-packages-R.html)).

A typical, half-day, lesson would use the first three lessons:

1. [Analyzing Patient Data](01-starting-with-data.html)
2. [Creating Functions](02-func-R.html)
3. [Analyzing Multiple Data Sets](03-loops-R.html)

An additional half-day could add the next two lessons:

4.  [Making choices](04-cond.html)
5.  [Command-Line Programs](05-cmdline.html)

Time-permitting, you can fit in one of these shorter lessons that cover bigger picture ideas like best practices for organizing code, reproducible research, and creating packages:

6.  [Best practices for using R and designing programs](06-best-practices-R.html)
7.  [Dynamic reports with knitr](07-knitr-R.html)
8.  [Making packages in R](08-making-packages-R.html)

## [Analyzing Patient Data](01-starting-with-data.html)

* Check learners are reading files from the correct location (set working
  directory); remind them of the shell lesson

* Provide shortcut for the assignment operator (`<-`) (RStudio: Alt+- on
  Windows/Linux; Option+- on Mac)


~~~{.r}
dat <- read.csv("data/inflammation-01.csv", header = FALSE)
~~~



~~~{.output}
Warning in file(file, "rt"): cannot open file 'data/inflammation-01.csv':
No such file or directory

~~~



~~~{.output}
Error in file(file, "rt"): cannot open the connection

~~~



~~~{.r}
animal <- c("m", "o", "n", "k", "e", "y")
# Challenge - Slicing (subsetting data)
animal[4:1]  # first 4 characters in reverse order
~~~



~~~{.output}
[1] "k" "n" "o" "m"

~~~



~~~{.r}
animal[-1]  # remove first character
~~~



~~~{.output}
[1] "o" "n" "k" "e" "y"

~~~



~~~{.r}
animal[-4]  # remove fourth character
~~~



~~~{.output}
[1] "m" "o" "n" "e" "y"

~~~



~~~{.r}
animal[-1:-4]  # remove first to fourth characters
~~~



~~~{.output}
[1] "e" "y"

~~~



~~~{.r}
animal[c(5, 2, 3)]  # new character vector
~~~



~~~{.output}
[1] "e" "o" "n"

~~~



~~~{.r}
# Challenge - Subsetting data
max(dat[5, 3:7])
~~~



~~~{.output}
Error in eval(expr, envir, enclos): object 'dat' not found

~~~


~~~{.r}
sd_day_inflammation <- apply(dat, 2, sd)
plot(sd_day_inflammation)
~~~

## [Addressing Data](01-supp-addressing-data.html)

* Note that the data frame `dat` is not the same set of data as in other lessons

## [Data Types and Structure](01-supp-data-structures.html)

* Lesson on data types and structures

## [Understanding Factors](01-supp-factors.html)

## [Introduction to RStudio](01-supp-intro-rstudio.html)

## [Reading and Writing .csv](01-supp-read-write-csv.html)























