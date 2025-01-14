This repository contains the R code for reproducing the computational results in the paper "Geometric ergodicity of trans-dimensional Markov chain Monte Carlo algorithms."

The core code is contained in the R Markdown file `Reproduce.Rmd` in the `code` folder.
It takes several minutes to knit.
Once knitted, it produces an `.html` file.

The data sets used in the paper and code are stored in the `data` folder.

Some of the code chunks in the R Markdown file are computationally intensive to run, so their results are preloaded into `.RData` files in the `output` folder.

The code in this repository is written in R, version 4.3.1. The package `renv` is used to cache the computing environment.
