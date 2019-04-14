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

echo Testing lca -m lca
cmp ./lca_out.tsv <(../bin/lca -t $1 -m lca -i lca_in.tsv)
echo Testing lca -m lowest
cmp ./lowest_out.tsv <(../bin/lca -t $1 -m lowest -i lowest_in.tsv)
echo Testing lca -m path
cmp ./path_out.tsv <(../bin/lca -t $1 -m path -i path_in.tsv)

