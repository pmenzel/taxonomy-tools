#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <climits>

#include "util.hpp"

void usage(char *progname);

int main(int argc, char **argv) {

	TaxTree nodes;
	std::unordered_map<TaxonId,unsigned int> node2depth;

	std::string nodes_filename;
	std::string in_filename;
	std::string out_filename;

	bool verbose = false;
	bool debug = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "ahdvrl:g:t:i:o:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'i':
				in_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'o':
				out_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(in_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }
	if(out_filename.length() == 0) { error("Please specify the name of the output file, using the -o option."); usage(argv[0]); }

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
	parseNodesDmp(nodes,nodes_file);
	nodes_file.close();

	std::ifstream inputfile;
	inputfile.open(in_filename);
	if(!inputfile.is_open()) { error("Could not open file " + in_filename); exit(EXIT_FAILURE); }

	std::ofstream out_file;
	out_file.open(out_filename);
	if(!out_file.is_open()) {  error("Could not open file " + out_filename + " for writing."); exit(EXIT_FAILURE); }

	std::set<TaxonId> ids;
	std::string line;
	while(getline(inputfile, line)){
		if(line.length() == 0) { continue; }
		line += "\t";
		ids.clear();
		if(debug) std::cerr << "Processing line: " << line << std::endl;
		size_t start = 0, end = 0;
		while((end = line.find_first_of("\t",start)) != std::string::npos) {
			std::string curr_id = line.substr(start, end - start);
			if(debug) std::cerr << "Start: " << start << " End: " << end << std::endl;
			if(debug) std::cerr << "Current field: " << ">" << curr_id << "<" << std::endl;
			TaxonId id;
			try {
				id = stoul(curr_id);
			}
			catch(const std::invalid_argument& ia) {
				std::cerr << "Error: Bad number in taxon id " << curr_id << std::endl;
				break;
			}
			catch (const std::out_of_range& oor) {
				std::cerr << "Error: Bad number (out of range error) in taxon id " << curr_id << std::endl;
				break;
			}

			if(nodes.count(id)>0) {
					ids.insert(id);
			}
			else {
				std::cerr << "Taxon ID " << id << " not found in nodes.dmp" << std::endl;
			}
			start = end+1;
		}

		if(!ids.empty()) {
			TaxonId lca = (ids.size()==1) ?  *(ids.begin()) : lca_from_ids(nodes, ids);
			if(debug) std::cerr << "LCA=" << lca << std::endl;
			out_file << lca << "\n";
		}
	}

	inputfile.close();
	out_file.close();

}

void usage(char *progname) {
	fprintf(stderr, "Copyright 2018 Peter Menzel\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -i input.tsv -o output.txt\n", progname);
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file.\n");
	fprintf(stderr, "   -i FILENAME   Name of tab-delimited input file.\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	exit(EXIT_FAILURE);
}

