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

## lca

This program reads a file with a list of taxon IDs (separated by tabs) per line and
outputs the least common ancestor (LCA) of the taxa.

Usage:
```
lca -t nodes.dmp -i input.tsv -o out
```

By default, the program returns the standard LCA for the given list of taxon IDs.
When setting the option `-m lowest`, it instead returns that taxon which is lowest
in a path and above any bifurcation creates by the input IDs.

Example:

![LCA modes](img/lca_modes.png?raw=true "LCA modes")

When given the list of nodes `a`, `b`, and `c` as input, the standard LCA of these three nodes is `c`,
whereas the option `-m lowest` would output the node just above `a` and `b`.

When setting `-m path`, it will return the leaf node that has the heighest weighted path to the root.

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


