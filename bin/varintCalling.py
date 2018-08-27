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
import pandas as pd
import numpy as np
from Bio import SeqIO,AlignIO
from Bio.Seq import Seq,translate,reverse_complement
from Bio.Data import CodonTable
from subprocess import Popen,PIPE
from scipy.stats import fisher_exact,chi2_contingency
import statsmodels.stats.multitest as ssm
from collections import Counter

def parseCommand():
    parser = ArgumentParser(description='Virus variant detection pipeline',version='1.0.0')
    parser.add_argument('-i', action='store', dest='inFa', help='Input is a Align\'s fasta file, based codon way')
    parser.add_argument('-r', action='store', dest='reference',default='the first sequence in inFa', help='which sample id regard as reference')
    parser.add_argument('-I', action='store', dest='info', help='sample info table')
    parser.add_argument('-t', action='store', dest='set_table', default=1, type=int, help='Genetic Code Table, default is 1')
    parser.add_argument('-g', action='store', dest='groupColumns',help='the columns of info table, which is the groups of sample')
    parser.add_argument('-o', action='store', dest='Opath', help='Path of Output File')
    return parser.parse_args()

def tranTable():
    print '''
    Genetic Code Table

        ID        List of Names
        1    ->    ['Standard', 'SGC0']
        2    ->    ['Vertebrate Mitochondrial', 'SGC1']
        3    ->    ['Yeast Mitochondrial', 'SGC2']
        4    ->    ['Mold Mitochondrial', 'Protozoan Mitochondrial', 'Coelenterate Mitochondrial', 'Mycoplasma', 'Spiroplasma', 'SGC3']
        5    ->    ['Invertebrate Mitochondrial', 'SGC4']
        6    ->    ['Ciliate Nuclear', 'Dasycladacean Nuclear', 'Hexamita Nuclear', 'SGC5']
        9    ->    ['Echinoderm Mitochondrial', 'Flatworm Mitochondrial', 'SGC8']
        10    ->    ['Euplotid Nuclear', 'SGC9']
        11    ->    ['Bacterial', 'Plant Plastid']
        12    ->    ['Alternative Yeast Nuclear']
        13    ->    ['Ascidian Mitochondrial']
        14    ->    ['Alternative Flatworm Mitochondrial']
        15    ->    ['Blepharisma Macronuclear']
        16    ->    ['Chlorophycean Mitochondrial']
        21    ->    ['Trematode Mitochondrial']
        22    ->    ['Scenedesmus obliquus Mitochondrial']
        23    ->    ['Thraustochytrium Mitochondrial']
        24    ->    ['Pterobranchia Mitochondrial']
    '''

    """ ###氨基酸三字母和一字母对应表
    amacid_dict={'F':'Phe','S':'Ser','Y':'Tyr','C':'Cys','L':'Leu','END':'Ter','W':'Trp',
                 'P':'Pro','H':'His','R':'Arg','Q':'Gln','I':'Ile','T':'Thr','N':'Asn',
                 'K':'Lys','M':'Met','V':'Val','A':'Ala','D':'Asp','G':'Gly','E':'Glu'}
    """

def translation(records,Opath,set_table):
    "translation nucleotide"
    amacid_record=[]
    for fasta in records:
        fasta.seq=translate(fasta.seq,table=set_table)
        amacid_record.append(fasta)
    SeqIO.write(records,Opath+".faa","fasta")   ###保存翻译结果
    return amacid_record

def statSameSite(same_virus,Opath,inFa,mu_type,ref):
    "一致位点的密码子在特定氨基酸下的使用频率统计"
    if ref == 'the first sequence in inFa':
        same_virus=same_virus.iloc[0]
    else:
        same_virus=same_virus.loc[ref,:]

    same_virus=pd.DataFrame(same_virus).reset_index()
    same_virus.columns=['refLocus','refAmAcid']

    seq_id=[]
    for fasta in SeqIO.parse(inFa,'fasta'):
        seq_id.append(fasta.id)
        same_virus[fasta.id]=['{0}'.format(fasta.seq[3*i-3:3*i]) for i in same_virus.refLocus]

    same_virus.set_index(['refAmAcid','refLocus'],inplace=True)
    same_virus.columns=seq_id
    same_virus.sort_index(level='refLocus',inplace=True)
    # same_virus.to_csv(Opath+mu_type+".substitution.eachSample.txt",sep='\t')
    return same_virus

def muscleAndDiff(records,Opath):
    '''specific locus calling'''

    order_id=[]
    for i in records:
        order_id.append(i.id)

    child=Popen('muscle',stdin=PIPE,stdout=PIPE,stderr=PIPE)
    SeqIO.write(records,child.stdin,"fasta")
    child.stdin.close()
    child.stderr.close()

    AmAcid_align=dict()

    for line in SeqIO.parse(child.stdout,"fasta"):
        AmAcid_align[line.id]=str(line.seq)

    virus=pd.DataFrame()
    for i in order_id:
        temp=list(AmAcid_align[i])
        temp.insert(0,i)
        virus=virus.append([temp])

    handle=open(Opath+".align.faa","w")

    for seq_id in order_id:
        handle.write(">{0}\n{1}\n".format(seq_id,AmAcid_align[seq_id]))
    handle.close()

    ##将序列每个字符都分开
    virus=virus.fillna('None')
    virus.set_index(0,inplace=True)

    diff_record=[]          #一致位点和不一致位点
    same_record=[]
    for k,v in virus.iteritems():
        if (len(v.drop_duplicates())>1):
            temp10=v.drop_duplicates()
            if ('-' in temp10.values and len(temp10)==2):  ##序列gap的处理
                same_record.append(k)
            else:
                diff_record.append(k)
        else:
            same_record.append(k)

    diff_virus=virus.loc[:,diff_record]
    diff_virus_bak=diff_virus.apply(lambda x:''.join(x),axis=1)

    # diff_locus=','.join(map(str,map(lambda x:x,diff_record)))  #差异位点的坐标
    diff_locus=','.join(map(str,diff_record))  #差异位点的坐标

    print >> open(Opath+'.aminoAcid.substitution.coord','w'),diff_locus   ###保存差异位点坐标
    same_virus=virus.loc[:,same_record] #一致位点

    return diff_virus_bak,diff_locus,same_virus,diff_virus

def marker_more50(series,compare_list):
    """
    各个分组的优势密码子必须达到各自的50%，其优势密码子不同，才保留这个位点
    """
    more_stats=[]
    more_codons=[]
    for eachGroup in compare_list:
        groupCodon=series.loc[eachGroup]
        group_conpose=Counter(groupCodon).most_common(1)

        more_codons.append(group_conpose[0][0])
        more_codon=group_conpose[0][1]/float(series.shape[0])>0.5
        more_stats.append(more_codon)

    if sum(more_stats) and len(set(more_codons))>1:
        return True
    else:
        False

def main():
    para=parseCommand()
    inFa=para.inFa
    info=para.info
    reference=para.reference
    groupColumns=para.groupColumns
    set_table=para.set_table
    Opath=para.Opath

    if not para.inFa:
        print "Please set -i.\nsee -h for help"
        tranTable()
        sys.exit(0)

    if not os.path.exists(Opath):
        os.system('mkdir -p '+Opath)

    Opath=Opath+'/'+Opath.split('/')[-1]

    ################### Step01.pre-process #####################
    standard_table = CodonTable.unambiguous_dna_by_id[para.set_table]  #密码表
    codon=pd.DataFrame(data=standard_table.forward_table,index=['AmAcid']).T  #建密码子数据框

    codon.index.name='Codon'
    codon.reset_index(inplace=True)
    codon['Number']=0

    seq_record=[]
    for fasta in SeqIO.parse(inFa, 'fasta'):
        fasta.seq=Seq(str(fasta.seq).replace('-',''))
        seq_record.append(fasta)

    amacid_records=translation(seq_record,Opath,set_table)

    ## use muscle to Align amacid and find diff and same locus
    diff_virus,diff_locus,same_virus,diff_virus_2=muscleAndDiff(amacid_records,Opath)

    if diff_locus != '':
        fasta='\n'.join(['>'+k+"\n"+v for k,v in diff_virus.iteritems()])
        print >> open(Opath+'.aminoAcid.substitution.faa','w'),fasta

    fit_record_codon=statSameSite(same_virus,Opath,inFa,".synonymous",reference)
    diff_record_codon=statSameSite(diff_virus_2,Opath,inFa,".nonsynonymous.codon",reference)

    ################### Step01.pre-process #####################

    ############### Step02.comparative analysis of raw variant calling(Step01) ###############

    mapinfo=pd.read_table(info,index_col=0,sep="\t") ## read sample info table

    ### 样本分组
    compare_list=[]
    for i in mapinfo[groupColumns].unique():
        compare_list.append(mapinfo.loc[mapinfo[groupColumns]==i].index)
    ### 样本分组

    ########  synonymous substitution  ########
    # fit_record_codon=pd.read_table(Ipath+'.synonymous.substitution.eachSample.txt',sep='\t',header=0,index_col=[0,1])
    fit_record_codon=fit_record_codon.T
    fit_record_codon=fit_record_codon.loc[mapinfo.index]

    ## 去除完全一致位点
    fit_columns=[i for i,k in fit_record_codon.iteritems() if len(k.unique())>1]
    fit_record_codon=fit_record_codon.loc[:,fit_columns]
    ## 去除完全一致位点

    ## 根据分组统计每个位点的密码子分布
    temp=fit_record_codon.T
    temp=temp.reset_index()
    temp2=pd.melt(temp,id_vars=['refAmAcid','refLocus'])
    temp2['Lineage']=list(mapinfo.loc[temp2.ID][groupColumns])

    temp2.rename(columns={'value':"Codon"},inplace=True)
    temp2['value']=1
    temp2=temp2.sort_values(by=['refLocus','refAmAcid'])
    del temp2['ID']

    temp3=temp2.groupby(['refAmAcid','refLocus','Codon','Lineage'])['value'].count()

    temp3=temp3.reset_index()
    temp4=temp3.sort_values(by="refLocus")
    ## 根据分组统计每个位点的密码子分布

    def fi_exact(series):
        temp=pd.pivot_table(series,values='value',index=['Lineage'],columns=['Codon'],aggfunc=np.sum)
        temp=temp.fillna(0)
        if len(temp.columns)==1:
            return 1
        else:
            return chi2_contingency(temp)[1]

    temp5=temp4.groupby(['refAmAcid','refLocus']).apply(fi_exact)
    temp5=temp5.reset_index()
    temp5=temp5.sort_values(by="refLocus")
    temp5=temp5.rename(columns={0:'pvalue'})

    temp5.loc[:,"pvalue_corr"]=ssm.multipletests(temp5.pvalue,alpha=0.05,method="fdr_bh")[1]
    temp5.loc[:,'sign']=temp5.pvalue_corr<0.05

    temp5=temp5.loc[temp5.pvalue_corr<0.05]
    temp6=temp5.set_index(['refAmAcid','refLocus']).T

    fit_record_codon=fit_record_codon.loc[:,temp6.columns]

    column_record=[]
    for i,k in fit_record_codon.iteritems():
       if marker_more50(k,compare_list):
           column_record.append(i)

    fit_record_codon=fit_record_codon.loc[:,column_record]

    if fit_record_codon.shape[1] != 0:
        fit_record_codon.to_csv(Opath+".synonymous.substitution.eachSample.txt",sep="\t")
    ########  synonymous substitution  ########

    ########  non-synonymous substitution  ########

    ### non-sysnonymous AminoAcid
    if os.path.exists(Opath+".aminoAcid.substitution.faa"):
        nonsyno=pd.DataFrame()
        for line in SeqIO.parse(Opath+".aminoAcid.substitution.faa","fasta"):
            temp=list(line.seq)
            temp.insert(0,line.id)
            nonsyno=nonsyno.append([temp])

        nonsyno=nonsyno.set_index(0)
        nonsyno.index.name='ID'
        nonsyno.columns=open(Opath+".aminoAcid.substitution.coord",'r').readline().rstrip().split(',')
        nonsyno.columns.name="Locus"
        nonsyno.columns=nonsyno.columns.astype(int)
        nonsyno=nonsyno.loc[mapinfo.index]

        #### filter those locus that only have - and X
        nonsyno_columns=[i for i,k in nonsyno.iteritems() if len(k.unique())>1]
        nonsyno=nonsyno.loc[:,nonsyno_columns]

        nonsyno_columns=[]
        for k,v in nonsyno.iteritems():
            if (len(v.drop_duplicates())>1):
                temp10=v.drop_duplicates()
                ##序列gap的处理
                if ('-' in temp10.values and 'X' in temp10.values and len(temp10)==3):
                    pass
                elif ('-' in temp10.values and 'X' not in temp10.values and len(temp10)==2):
                    pass
                elif ('-' not in temp10.values and 'X' in temp10.values and len(temp10)==2):
                    pass
                else:
                    nonsyno_columns.append(k)
        nonsyno=nonsyno.loc[:,nonsyno_columns]
        #### filter those locus that only have - and X

        def fi_exact_amacid(series):
            temp=pd.pivot_table(series,values='value',index=['Lineage'],columns=['AmAcid'],aggfunc=np.sum)
            temp=temp.fillna(0)
            if len(temp.columns)==1:
                return 1
            else:
                return chi2_contingency(temp)[1]

        nonsyno_temp=nonsyno.T
        nonsyno_temp=nonsyno_temp.reset_index()

        temp11=pd.melt(nonsyno_temp,id_vars=['Locus'])
        temp11['Lineage']=list(mapinfo.loc[temp11.ID][groupColumns])
        temp11.rename(columns={'value':'AmAcid'},inplace=True)
        temp11=temp11.loc[temp11.AmAcid !='X']
        temp11['value']=1
        temp11=temp11.sort_values(by=['Locus'])
        del temp11['ID']

        temp22=temp11.groupby(['Locus','AmAcid','Lineage'])['value'].count()
        temp22=temp22.reset_index()
        temp33=temp22.sort_values(by='Locus')
        
        temp44=temp33.groupby(['Locus']).apply(fi_exact_amacid)
        temp44=temp44.reset_index()
        temp44=temp44.sort_values(by="Locus")
        temp44.rename(columns={0:'pvalue'},inplace=True)

        temp44.loc[:,"pvalue_corr"]=ssm.multipletests(temp44.pvalue,alpha=0.05,method="fdr_bh")[1]
        temp55=temp44.loc[temp44.pvalue_corr<0.05]
        temp66=temp55.set_index(['Locus']).T
        nonsyno=nonsyno.loc[:,temp66.columns]

        if nonsyno.shape[1] != 0:
            nonsyno.to_csv(Opath+".nonsynonymous.AmAcid.substitution.txt",sep="\t")

        ### non-sysnonymous codon
        # diff_record_codon=pd.read_table(Ipath+'.nonsynonymous.codon.substitution.eachSample.txt',sep='\t',header=0,index_col=[0,1])
        diff_record_codon.reset_index(inplace=True)
        diff_record_codon=diff_record_codon.loc[diff_record_codon.refLocus.isin(nonsyno.columns)]

        if diff_record_codon.shape[0] != 0:
            diff_record_codon.T.to_csv(Opath+'.nonsynonymous.codon.substitution.eachSample.txt',sep="\t",header=None)
        ########  non-synonymous substitution  ########

    ############### Step02.comparative analysis of raw variant calling(Step01) ###############

    ############ Step03. remove intermidiate file ############
        os.system('rm '+Opath+".aminoAcid.substitution.faa")
        os.system('rm '+Opath+".aminoAcid.substitution.coord")
    ############ Step03. remove intermidiate file ############

if __name__  == "__main__":
    main()