# Taxonomy Tools

A collection of tools for specific use cases that deal with the NCBI taxonomy.
For all programs, first download and extract [taxdump.tgz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).

## blasthits2krona

This program requires [Krona](https://github.com/marbl/Krona) and
[nucl_gb.accession2taxid.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz).

Usage:
```
blastn ... -outfmt 6 -out blast.tab.out

blasthits2krona -t nodes.dmp -n names.dmp -a nucl_gb.accession2taxid -i blast.tab.out -o blast.krona

Krona/KronaTools/scripts/ImportText.pl -o blast.html blast.krona
```
Open `blast.html`.


## License

Copyright (c) 2018 Peter Menzel

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


