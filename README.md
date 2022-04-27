# MAna - R functions for Metagenomic Analyses and visualization

---------

MAna is a package-like R script including frequently used in-house functions for data analyses and visualization of "typical" metagenomic studies.

## Installation

MAna is compilation-free (up to now). You may simply download the file ***func.R*** and load the functions into your working enviroment by running `source("func.R")` command in R.

However, you may still need to install labraries used by the function(s) you would like to use manually.

## Usage

Standard input of MAna is a table of profile and another of factors/phenotypes. Both of them should be presented in type of ***data.frame***. which is formatted following rules:

1. Samples in row, while variables in column;
1. With no header nor column of sample ID;
1. Rows and columns should be named properly;
1. Categorical variable should be presented in type of ***factor*** or ***ordered***.

You may check the declaration of spesific function(s) (or even the codes) for more details on usage. Particularly, the parameter ***df*** will always wait for a argument of the profile table metioned above.

When you have multiple profiles for the same set of samples, you may also combine them into a ***list*** to make your working enviroment neat. You may further embed the list of profiles and the table of factors together into a ***nested list*** (that is, a ***dataset***), and perform some simple operation by using `dataset.subset()` and `dataset.trim()`.
