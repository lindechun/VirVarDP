## VirVarDP

Virus variant detection pipeline

## Author

Dechun Lin

### Introduction

A pipeline based on `python` to detect synonnymous or non-synonnymous substitutions within homologous gene of Virus complete genome

This lists the basic information for using [`VirVarDP`](https://github.com/lindechun/VirVarDP).

#### Requirements

* A UNIX based operating system.

* python 2.7

* python packages: Biopython, scipy, pandas, statsmodels

* [megacc](https://www.megasoftware.net/)

#### Installation

Download VirVarDP from GitHub. You'll need to add Prokka's bin directory to your $PATH.:

```
git clone https://github.com/lindechun/VirVarDP.git
/your/path/to/VirVarDP/VirVarDP.py -h
```

You will nedd to install all the dependencies packages by `pip` or `conda`.


### Basic test data set

See `/your/path/to/VirVarDP/data/` for test data set.

### Usage

```VirVarDP [options]```

See `/your/path/to/VirVarDP/example/` for example of variant calling of 100 Zika sequences.

