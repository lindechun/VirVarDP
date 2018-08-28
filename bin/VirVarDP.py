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

def parseCommand():
    parser = ArgumentParser(description='The flow of Virus variant detection pipeline',version='1.0.0')
    parser.add_argument('-i', action='store', dest='falist', help='a list file consists of individual Genes, based codon align (megacc)')
    parser.add_argument('-r', action='store', dest='reference',default='the first sequence in inFa', help='which sample id regard as reference')
    parser.add_argument('-I', action='store', dest='info', help='sample info table')
    parser.add_argument('-l', action='store', dest='GeneLength', help='eachGene.length.txt')
    parser.add_argument('-t', action='store', dest='set_table', default=1, type=int, help='Genetic Code Table, default is 1')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample')
    parser.add_argument('-G', action='store', dest='Gap', default="0", help='The gaps before specific gene. eg: \'**:40,**:50\'')
    parser.add_argument('-j', action='store', dest='OrderGenes',default='', help='The order of genes')
    parser.add_argument('-p', action='store', dest='Prefix', help='Output file prefix')
    parser.add_argument('-o', action='store', dest='Opath', help='Path of Output File')
    return parser.parse_args()

def main():
    para=parseCommand()
    falist=para.falist
    info=para.info
    reference=para.reference
    groupColumns=para.groupColumns
    GeneLength=para.GeneLength
    Gap=para.Gap
    set_table=para.set_table
    orderGenes=para.OrderGenes
    prefix=para.Prefix
    Opath=para.Opath

    if not falist:
        sys.exit("Please set -i.")

    ## Step00
    print("# Step00. create results path")
    print('mkdir '+Opath)
    os.system('mkdir '+Opath)

    ## Step01
    print("\n# Step01. variant calling for individual genes")

    for i in open(falist,'r'):
        i=i.strip()
        geneID=os.path.basename(i).split('.')[0]
        dirpath=os.path.dirname(os.path.realpath(__file__))

        variantCalling="python {0}/varintCalling.py -i {1} -r {2} -I {3} -g {4} -o {5}".format(dirpath,i,reference,info,groupColumns,Opath+'/'+geneID)
        print(variantCalling)
        os.system(variantCalling)

    ## Step02
    print("\n# Step02. Gather every gene's variant results")
    print("## synonymous substitutions")
    syno_sub="python {0}/synonymous.substitutions.gather.py -i {1} -g {2} -I {3} -l {4} -o '{5}' -r {6} -G '{7}' -p {8}".format(dirpath,info,groupColumns,Opath,GeneLength,orderGenes,reference,Gap,Opath+'/'+prefix)
    print(syno_sub)
    os.system(syno_sub)

    print("\n## non-synonymous codon substitutions")
    nonsyno_codon="python {0}/nonsynonymous.codon.substitutions.gather.py -i {1} -g {2} -I {3} -l {4} -o '{5}' -r {6} -G '{7}' -p {8}".format(dirpath,info,groupColumns,Opath,GeneLength,orderGenes,reference,Gap,Opath+'/'+prefix)
    print(nonsyno_codon)
    os.system(nonsyno_codon)

    print("\n## non-synonymous AminoAcid substitutions")
    nonsyno_aminoAcid="python {0}/nonsynonymous.AminoAcid.substitutions.gather.py -i {1} -g {2} -I {3} -l {4} -o '{5}' -r {6} -G '{7}' -p {8}".format(dirpath,info,groupColumns,Opath,GeneLength,orderGenes,reference,Gap,Opath+'/'+prefix)
    print(nonsyno_aminoAcid)
    os.system(nonsyno_aminoAcid)

if __name__  == "__main__":
    main()
