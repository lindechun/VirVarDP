#!/usr/bin/env python
# coding=utf-8

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2016, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "0.0.1"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

from glob import glob
from subprocess import Popen,PIPE
from Bio import SeqIO
import sys,os

if len(sys.argv) != 3:
    sys.exit("Usage: python {0} < the path of raw fasta, dir: *.fa> <output path>".format(sys.argv[0]))

for i in glob(sys.argv[1]+'/*.fa'):
    filename=os.path.basename(i)
    Ofile=sys.argv[2]+'/'+filename

    os.system('megacc -a '+os.path.dirname(os.path.realpath(__file__))+'/clustal_align_coding.mao -d '+i+' -o '+Ofile+' -f Fasta -n')
    os.system('mv '+Ofile.replace('.fa','.fasta')+' '+Ofile)

    record=[]
    for line in SeqIO.parse(Ofile,'fasta'):
        record.append(line)
    SeqIO.write(record,open(Ofile,'w'),"fasta")
