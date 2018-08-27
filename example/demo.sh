#!/bin/sh
dirpath=$(cd `dirname $0`; pwd)
prefix='Zika'

mkdir -p results

# Step01. variant calling for individual genes
while read lines
do
    DataId=`basename $lines .fa`
    DataDir=`dirname $lines`
    python $dirpath/../bin/VirVarDP.py -i $lines -r LC002520 -I $dirpath/../data/SampleInfo.txt -g Lineage -o results/$DataId
done < data.list

# Step02. Gather every gene's variant results
## synonymous substitutions
python $dirpath/../bin/synonymous.substitutions.gather.py  -i $dirpath/../data/SampleInfo.txt -g Lineage -I results -l $dirpath/../data/eachGene.length.txt -o "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" -r LC002520 -G 'NS4B:69' -p results/$prefix

## non-synonymous codon substitutions
python $dirpath/../bin/nonsynonymous.codon.substitutions.gather.py  -i $dirpath/../data/SampleInfo.txt -g Lineage -I results -l $dirpath/../data/eachGene.length.txt -o "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" -r LC002520 -G 'NS4B:69' -p results/$prefix

## non-synonymous AminoAcid substitutions
python $dirpath/../bin/nonsynonymous.AminoAcid.substitutions.gather.py  -i $dirpath/../data/SampleInfo.txt -g Lineage -I results -l $dirpath/../data/eachGene.length.txt -o "C,prM,E,NS1,NS2A,NS2B,NS3,NS4A,NS4B,NS5" -r LC002520 -G 'NS4B:69' -p results/$prefix
