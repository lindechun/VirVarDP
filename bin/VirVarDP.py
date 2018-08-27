#!/usr/bin/env python
# coding=utf-8
from __future__ import division

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2016, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "1.0.0"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

from argparse import ArgumentParser
import re,sys,os
# import pandas as pd
# import numpy as np
# from Bio import SeqIO,AlignIO
# from Bio.Seq import Seq,translate,reverse_complement
# from Bio.Data import CodonTable
# from subprocess import Popen,PIPE
# from scipy.stats import fisher_exact,chi2_contingency
# import statsmodels.stats.multitest as ssm
# from collections import Counter

def parseCommand():
    parser = ArgumentParser(description='The main of Virus variant detection pipeline',version='1.0.0')
    parser.add_argument('-i', action='store', dest='inFa', help='Input is a Align\'s fasta file, based codon way')
    parser.add_argument('-r', action='store', dest='reference',default='the first sequence in inFa', help='which sample id regard as reference')
    parser.add_argument('-I', action='store', dest='info', help='sample info table')
    parser.add_argument('-t', action='store', dest='set_table', default=1, type=int, help='Genetic Code Table, default is 1')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample')
    parser.add_argument('-o', action='store', dest='Opath', help='Path of Output File')
    return parser.parse_args()

def main():
    para=parseCommand()
    inFa=para.inFa
    info=para.info
    reference=para.reference
    groupColumns=para.groupColumns
    set_table=para.set_table
    Opath=para.Opath


if __name__  == "__main__":
    main()