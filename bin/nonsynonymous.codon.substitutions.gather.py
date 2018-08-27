#!/usr/bin/env python
# coding=utf-8

__author__ = "Lin Dechun"
__copyright__ = "Copyright 2016, BGI Research."
__credits__ = ["Lin Dechun"]
__version__ = "0.0.1"
__maintainer__ = "Lin Dechun"
__email__ = "lindechun@genomics.cn"

import os,sys
from argparse import ArgumentParser
import pandas as pd
from pandas import ExcelWriter
import xlsxwriter
from glob import glob
import re
from collections import Counter

def parseCommand():
    parser = ArgumentParser(description='This script use to gather synonymous substitutions of each genes',version='1.0.0')
    parser.add_argument('-i', action='store', dest='Info', help='sample info table')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample')
    parser.add_argument('-I', action='store', dest='Ipath', help='The path of CodonM\'s results')
    parser.add_argument('-l', action='store', dest='GeneLength', help='eachGene.length.txt')
    parser.add_argument('-o', action='store', dest='Order',default='', help='The order of genes')
    parser.add_argument('-r', action='store', dest='reference',default='the first sequence in info table', help='which sample id regard as reference')
    parser.add_argument('-G', action='store', dest='Gap', default="0", help='The gaps before specific gene. eg: \'**:40,**:50\'')
    parser.add_argument('-p', action='store', dest='Prefix', help='Output file prefix')
    return(parser.parse_args())

def calculateLocus(genelength,orderGene,Gap_list):
    ama={}
    for i,k in genelength.iteritems():
        ama[i]=[1,k.max()/3]

    for k,i in enumerate(orderGene):
        if k==0:
            continue

        if i in Gap_list:
            ama[orderGene[k]][0]=sum(ama[orderGene[k-1]])+Gap_list[i]/3
        else:
            ama[orderGene[k]][0]=sum(ama[orderGene[k-1]])

    for k,i in enumerate(orderGene):
        if k==0:
            pass
        else:
            ama[orderGene[k]][1]=sum(ama[orderGene[k]])-1
    return(ama)


def gatherSubstitution(Ipath,ama,orderGene):
    if not len(glob(Ipath+"/*/*.synonymous.substitution.eachSample.txt")):
        sys.exit('There isn\'t exist anyone syonymous substitution results of every gene')

    for i in orderGene:
        Gene_syno_file=Ipath+'/'+i+'/'+i+'.nonsynonymous.codon.substitution.eachSample.txt'

        if os.path.exists(Gene_syno_file):
            table=pd.read_table(Gene_syno_file,sep="\t",index_col=0,header=[0,1])
            columns_temp=zip(*[x for x in table.columns])
            Locus=list(columns_temp[1])
            locus_temp=[int(x)+ama[i][0]-1 for x in Locus]

            table.columns=pd.MultiIndex.from_tuples([tuple([i,columns_temp[0][k],x]) for k,x in enumerate(locus_temp)])
            if 'df' not in dir():
                df=table
            else:
                df=df.join(table)
    return(df)

def WriteExcel(Prefix,df,Ref):
    writer=ExcelWriter(Prefix+'.nonsynonymous.codon.substitution.xlsx',engine='xlsxwriter')
    # Add a format. Light red fill with dark red text.
    format1 = writer.book.add_format({'bg_color': '#FFC7CE',
                                   'font_color': '#9C0006'})
    # Add a format. Green fill with dark green text.
    format2 = writer.book.add_format({'bg_color': '#C6EFCE',
                                   'font_color': '#006100'})

    df.to_excel(writer,sheet_name="Sheet1")
    worksheet=writer.sheets["Sheet1"]

    count=0
    for i,k in df.iteritems():
        RefCodon=k.loc[Ref]

        count+=1
        worksheet.conditional_format(4,count,len(df)+3,count,
            {'type':'text',
            'criteria': 'not containing',
            'value':RefCodon,
            'format':format1}
            )
        worksheet.conditional_format(4,count,len(df)+3,count,
            {'type':'text',
            'criteria':'containing',
            'value':RefCodon,
            'format':format2}
            )

    writer.close()

def main():
    para=parseCommand()
    Info=para.Info
    groupColumns=para.groupColumns
    Ipath=para.Ipath
    GeneLength=para.GeneLength
    Order=para.Order
    Ref=para.reference
    Gap=para.Gap
    Prefix=para.Prefix

    if not Order:
      sys.exit("Please set -o")

    info=pd.read_table(Info,index_col=0)
    genelength=pd.read_table(GeneLength,sep="\t",index_col=0)
    orderGene=Order.split(',')

    Gap_list={}
    for i in Gap.split(','):
        gene,gap=i.split(':')
        Gap_list[gene]=int(gap)

    ama=calculateLocus(genelength, orderGene, Gap_list)

    dat=gatherSubstitution(Ipath, ama, orderGene)
    dat=dat.loc[info.index,]
    dat.insert(0,"Group", info[groupColumns])

    WriteExcel(Prefix,dat,Ref)  ## syno substitutions writed to excel

if __name__ == '__main__':
    main()
