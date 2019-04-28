#!/bin/bash

# give path to nodes.dmp as first argument

if [ $# -eq 0 ]
then
	echo "Give path to nodes.dmp as argument"
	exit
fi
if [ ! -r "$1" ]
then
	echo "Cannot read file $1"
	exit
fi

echo Testing lca
echo Testing -m lca
cmp ./lca_out.tsv <(../bin/lca -t $1 -m lca -i lca_in.tsv)
echo Testing -m lowest
cmp ./lowest_out.tsv <(../bin/lca -t $1 -m lowest -i lowest_in.tsv)
echo Testing -m path
cmp ./path_out.tsv <(../bin/lca -t $1 -m path -i path_in.tsv)

echo Testing -m lca + exclusion
cmp ./exclusion_lca_out.tsv <(../bin/lca -t $1 -m lca -e exclusion_list.txt -i exclusion_lca_in.tsv)
echo Testing -m lowest + exclusion
cmp ./exclusion_lowest_out.tsv <(../bin/lca -t $1 -m lowest -e exclusion_list.txt -i exclusion_lowest_in.tsv)
echo Testing -m path + exclusion
cmp ./exclusion_path_out.tsv <(../bin/lca -t $1 -m path -e exclusion_list.txt -i exclusion_path_in.tsv)

