## VirVarDP

Virus variant detection pipeline

## Author

Dechun Lin

## Introduction

A pipeline based on `python` to detect synonnymous or non-synonnymous substitutions within homologous gene of Virus complete genome

This lists the basic information for using [`VirVarDP`](https://github.com/lindechun/VirVarDP).

## Requirements

* A UNIX based operating system.

* python 2.7

* python packages: biopython, scipy, pandas, statsmodels

* [megacc](https://www.megasoftware.net/)

## Installation

Download VirVarDP from GitHub. You'll need to add Prokka's bin directory to your $PATH.:

```
git clone https://github.com/lindechun/VirVarDP.git
/your/path/to/VirVarDP/VirVarDP.py -h
```

You will nedd to install all the dependencies packages by `pip` or `conda`.


## Basic test data set

See `/your/path/to/VirVarDP/data/` for test data set.

## Usage

```
$ python ../bin/VirVarDP.py -h
usage: VirVarDP.py [-h] [-v] [-i FALIST] [-r REFERENCE] [-I INFO]
                   [-l GENELENGTH] [-t SET_TABLE] [-g GROUPCOLUMNS] [-G GAP]
                   [-j ORDERGENES] [-p PREFIX] [-o OPATH]

The flow of Virus variant detection pipeline

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show program's version number and exit
  -i FALIST        a list file consists of individual Genes, based codon align
                   (megacc)
  -r REFERENCE     which sample id regard as reference
  -I INFO          sample info table
  -l GENELENGTH    eachGene.length.txt
  -t SET_TABLE     Genetic Code Table, default is 1
  -g GROUPCOLUMNS  the columns of info table, which is the groups of sample
  -G GAP           The gaps before specific gene. eg: '**:40,**:50'
  -j ORDERGENES    The order of genes
  -p PREFIX        Output file prefix
  -o OPATH         Path of Output File
```

See `/your/path/to/VirVarDP/example/` for example of variant calling of 100 Zika sequences.