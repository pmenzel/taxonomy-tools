# Taxonomy Tools

A collection of tools for specific use cases that deal with the NCBI taxonomy.
For all programs, first download and extract `ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`.

## blasthits2krona

This program requires [Krona](https://github.com/marbl/Krona) and
`ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz`.

Usage:
```
blastn ... -outfmt 6 -out blast.tab.out

blasthits2krona -t nodes.dmp -n names.dmp -a nucl_gb.accession2taxid -i blast.tab.out -o blast.krona

Krona/KronaTools/scripts/ImportText.pl -o blast.html blast.krona
```
Open `blast.html`.

## acc2name

This program reads a text file with GenBank nucleotide accession numbers,
appends the corresponding taxon id and the taxon name, and prints the output in
a 3-column tab-separated file.

It requires
`ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz`.

Usage:
```
acc2name -t nodes.dmp -n names.dmp -a nucl_gb.accession2taxid -i acc_list.txt -o acc_name.tsv
```

## taxonid2name

This program reads a text file with taxon IDs and appends the corresponding taxon
name, separated by a tab.

Usage:
```
taxonid2name -t nodes.dmp -n names.dmp -i taxon_id_list.txt [-o out.tsv]
```

Optional arguments:
```
-p    Print full taxon path.
-r    Print taxon path containing only ranks specified by a comma-separated list,
		  for example: superkingdom,phylum,class,order,family,genus,species
```


## lca

This program reads a file with a list of taxon IDs (separated by tabs) per line and
outputs the least common ancestor (LCA) of the taxa.

Usage:
```
lca -t nodes.dmp -i input.tsv -o out
```

By default, the program returns the standard LCA for the given list of taxon IDs.
When setting the option `-m lowest`, it instead returns that taxon that is lowest
in a path and above any bifurcation created by the input IDs.

Example:

![LCA modes](img/lca_modes.png?raw=true "LCA modes")

When given the list of nodes `a`, `b`, and `c` as input, the standard LCA of these three nodes is `c`,
whereas the option `-m lowest` would output the node just above `a` and `b`.

When setting `-m path`, it will return the leaf node that has the highest weighted path to the root.

### Exclusion list
It is possible to specify a file containing a list with taxon IDs (one per
line) using option `-e`.  These IDs are excluded from the LCA calculation in all
three modes, unless there are no other taxon IDs in an input row.

This is for example useful, to exclude taxa assigned to environmental or
unclassified nodes in the NCBI taxonomy.

## subtree

This program reads a list of taxon IDs (one taxon ID per line) and
outputs a list with all taxon IDs that are descendants of the taxon IDs
in the input list.

Usage:
```
subtree -t nodes.dmp -i in.txt -o out.txt
```

One use case for this program is to create a list of taxon IDs that
is used as an exclusion list for `lca` (see above).  
Similarly, the taxon ID list can be used for including or excluding taxa in
the NCBI BLAST+ programs, using the options `-taxidlist` or
`-negative_taxidlist` respectively, which are available with the v5 database
format from BLAST+ version >=2.8.0.

Example:
```
tar xf taxdump.tar.gz names.dmp nodes.dmp

grep "scientific name" names.dmp| grep -w -P "uncultured|environmental samples|metagenome|unclassified|unidentified" | cut -f1 >in.txt

subtree -t nodes.dmp -i in.txt > exclusion_list.txt

blastn -negative_taxidlist exclusion_list.txt -db ...
```

## taxonR

The file `taxonR.cpp` provides functions for quickly accessing the `nodes.dmp`
file from within R.

Example usage:
```{r}
> Rcpp::sourceCpp("/home/ptr/software/taxonomy-tools/src/taxonR.cpp")
> n <- loadNodesDmp("/home/ptr/temp/nodes.dmp")
> is_ancestor(n,2,122)
[1] TRUE
> parent(n,222)
[1] 506
```

## License

Copyright (c) 2018,2019 Peter Menzel

Taxonomy-tools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Taxonomy-tools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the file LICENSE for more details.

You should have received a copy of the GNU General Public License
along with the source code.  If not, see <http://www.gnu.org/licenses/>.


