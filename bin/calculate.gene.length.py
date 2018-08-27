#!/usr/bin/env python
# coding=utf-8

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2016, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "0.0.1"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

import os,sys
if len(sys.argv) != 2:
    sys.exit("Usage python {0} <path of fasta>".format(sys.argv[0]))

import pandas as pd
from glob import glob
from Bio import SeqIO

data=pd.DataFrame()
for i in glob(sys.argv[1]+'/*.fa'):
    for line in SeqIO.parse(i,'fasta'):
        data=data.append([[line.id,os.path.basename(i)[:-3],len(line.seq)]])

data.columns=['ID','Gene','Length']
data=data.pivot_table(values='Length',index='ID',columns='Gene')
data.to_csv("eachGene.length.txt",sep="\t")
