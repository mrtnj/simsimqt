# simsimqt

## Note about status of the package

This package is in development, and while published, not recommended for use
outside the research group yet.


## Introduction

`simsimqt` is a simulation package for quantitative traits controlled by a finite number of unlinked loci. It is meant for simulating scenarios with directional or stabilising selection while allowing flexibility about the breeding structure. One of the `sim`s stands for "simple", and the other for "simulation". For simulating the details of breeding programs or simulations that need to be genomically explicit, you might turn to a package like `AlphaSimR` or `SLIM`.


## Installing

Installation from this repository with the `remotes` package.

The package contains some C++ code, which means that you need a compiler set up for package installation. On Windows, this means that you need `RTools` installed. On a *nix system, this should be nothing out of the ordinary.

```
remotes::install_github("mrtnj/simsimqt")
```

**Note: Again, this package is not yet recommended to anyone outside of the research group.**


## Reading the vignette

There is an introductory vignette that walks through an example simulation. Find it in R after installation by running:

```
vignette("simsimqt-introduction")
```
