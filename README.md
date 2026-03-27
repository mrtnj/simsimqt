# simsimqt


## Introduction

`simsimqt` is a simulation package for quantitative traits controlled by a finite
number of unlinked loci. It is meant for simulating scenarios with directional
or stabilising selection while allowing flexibility about the breeding structure.
One of the `sim`s stands for "simple", and the other for "simulation". For
simulating the details of breeding programs or simulations that need to be
genomically explicit, you might turn to a package like `AlphaSimR` or `SLiM`.


## Note about status of the package

This package being developed in conjunction with research projects where it is
being used. While we will strive to maintain the interface and behaviour of
extant functions, changes may occur as new needs arise.

If you have questions about its uses for any particular purpose, feel free to
contact the author via email at _martin full stop johnsson full stop slu full stop se_.


## Installing

Installation from this repository with the `remotes` package.

The package contains some C++ code, which means that you need a compiler set up
for package installation. On Windows, this means that you need `RTools` installed.
On a *nix system, this should be nothing out of the ordinary.

```
remotes::install_github("mrtnj/simsimqt", build_vignettes = TRUE)
```


## Reading the vignette

There is an introductory vignette that walks through an example simulation.
Find it in R after installation by running:

```
vignette("simsimqt-introduction")
```
