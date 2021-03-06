#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <stdexcept>

#include "util.hpp"

void usage(char *progname);

int main(int argc, char** argv) {

	std::unordered_map<Accession,TaxonId> acc2taxid;
	TaxTree nodes;
	std::unordered_map<TaxonId, TaxonName> node2name;

	std::string nodes_filename = "";
	std::string names_filename = "";
	std::string acc2taxid_filename = "";
	std::string in1_filename = "";
	std::string out_filename;

	bool verbose = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hva:n:t:i:o:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'a':
				acc2taxid_filename = optarg; break;
			case 'o':
				out_filename = optarg; break;
			case 'n':
				names_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'i':
				in1_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}

	if(out_filename.length() == 0) { error("Error: Please specify the name of the output file, using the -o option."); usage(argv[0]); }
	if(names_filename.length() == 0) { error("Please specify the location of the names.dmp file with the -n option."); usage(argv[0]); }
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(acc2taxid_filename.length() == 0) { error("Please specify the location of the nucl_gb.accession2taxid file, using the -a option."); usage(argv[0]); }
	if(in1_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }

	std::ifstream acc2taxid_file;
	acc2taxid_file.open(acc2taxid_filename);
	if(!acc2taxid_file.is_open()) { error("Could not open file " + acc2taxid_filename); exit(EXIT_FAILURE); }
	if(verbose) std::cerr << "Reading accession to taxon id map from file " << acc2taxid_filename << std::endl;
	parse_accession2taxid(acc2taxid,acc2taxid_file);
	acc2taxid_file.close();

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
	if(verbose) std::cerr << "Reading taxonomic tree from file " << nodes_filename << std::endl;
	parseNodesDmp(nodes,nodes_file);
	nodes_file.close();

	std::ifstream names_file;
	names_file.open(names_filename);
	if(!names_file.is_open()) { error("Could not open file " + names_filename); usage(argv[0]); }
	if(verbose) std::cerr << "Reading taxon names from file " << names_filename << std::endl;
	parseNamesDmp(node2name,names_file);
	names_file.close();


	std::ifstream in1_file;
	in1_file.open(in1_filename);
	if(!in1_file.is_open()) {  error("Could not open file " + in1_filename); exit(EXIT_FAILURE); }

	std::ofstream out_file;
	out_file.open(out_filename);
	if(!out_file.is_open()) {  error("Could not open file " + out_filename + " for writing"); exit(EXIT_FAILURE); }
	if(verbose) std::cerr << "Writing to file " << out_filename << std::endl;

	if(verbose) std::cerr << "Processing " << in1_filename <<"..." << "\n";
	std::string line;
	while(getline(in1_file,line)) {
		if(line.length() == 0) { continue; }
		Accession acc = line;
		auto it_id = acc2taxid.find(acc);
		if(it_id == acc2taxid.end()) {
			std::cerr << "Warning: Accession " << acc << " is not found in "<< acc2taxid_filename << ". Skipping.\n";
			continue;
		}
		TaxonId id = it_id->second;

		auto it_name = node2name.find(id);
		if(it_name == node2name.end()) {
			std::cerr << "Warning: Taxon ID " << id << " for accession " << acc << " is not contained in names.dmp file "<< names_filename << ". Skipping.\n";
			continue;
		}
		std::string name = it_name->second;

		out_file << acc << "\t" << id << "\t" << name << "\n";
	}

	in1_file.close();
	out_file.close();

}

void usage(char *progname) {
	fprintf(stderr, "Copyright 2018 Peter Menzel\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -n names.dmp -a nucl_gb.accession2taxid -i blast.out -o blast2krona.out\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file\n");
	fprintf(stderr, "   -n FILENAME   Name of names.dmp file\n");
	fprintf(stderr, "   -a FILENAME   Name of nucl_gb.accession2taxid file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	exit(EXIT_FAILURE);
}

