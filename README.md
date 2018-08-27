## VirVarDP

### Introduction

Virus variant detection pipeline
A pipeline based on `python` to detect synonnymous or non-synonnymous substitutions within homologous gene of Virus complete genome, Developed by Dechun Lin.

This lists the basic information for using [`VirVarDP`](https://github.com/lindechun/VirVarDP).

### VirVarDP Installation Notes

#### Required dependencies

* A UNIX based operating system.

* python 2.7 installed.

* megacc

#### Installation

Download VirVarDP from GitHub:

```
git clone https://github.com/lindechun/VirVarDP
```

And then, install required modules.

```
$ cd /your/path/to/VirVarDP/bin/
$ conda install  dependencies.r
```
* you will need to add `/your/path/to/VirVarDP/bin/` to your PATH environment variable.

### Basic test data set

See `/your/path/to/VirVarDP/data/` for test data set.

### Usage

```VirVarDP [options]```

See `/your/path/to/VirVarDP/example/` for example of variant calling of 100 Zika sequences.

