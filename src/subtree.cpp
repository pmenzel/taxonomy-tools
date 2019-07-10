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
#include <map>
#include <climits>
#include <assert.h>

#include "util.hpp"

void usage(char *progname);

bool dfs(const TaxonId parent, const std::multimap<TaxonId, TaxonId> & parent2children, TaxonSet & output_ids) {

	auto range = parent2children.equal_range(parent);
	for(auto i = range.first; i != range.second; ++i) {
		output_ids.insert(i->second);
		if(parent2children.count(i->second) > 0) {
			dfs(i->second, parent2children, output_ids);
		}
	}

}

int main(int argc, char **argv) {

	TaxTree nodes;
	std::multimap<TaxonId, TaxonId> parent2children;
	TaxonSet input_ids;
	TaxonSet output_ids;

	std::string nodes_filename;
	std::string in_filename;
	std::string out_filename;

	bool verbose = false;
	bool debug = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hdvm:t:i:o:e:")) != -1) {
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

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
	parseNodesDmp(nodes,nodes_file);
	nodes_file.close();

	nodes_file.open(nodes_filename);
	if(!nodes_file) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
	parseNodesDmpTopDown(parent2children, nodes_file);
	nodes_file.close();

	std::ifstream list_file;
	list_file.open(in_filename);
	if(!list_file.is_open()) { error("Could not open file " + in_filename); exit(EXIT_FAILURE); }
	if(verbose) std::cerr << "Reading taxon IDs from file " << in_filename << std::endl;
	std::string line;
	while(getline(list_file, line)) {
		if(line.length() == 0) { continue; }
		size_t start = line.find_first_of("0123456789");
		if(start == std::string::npos) {
			continue;
		}
		size_t end = line.find_first_not_of("0123456789",start+1);
		if(end == std::string::npos) {
			end = start + line.length()-start;
		}
		try {
			uint64_t taxid = stoul(line.substr(start,end-start));
			if(debug)	std::cerr << "Found taxon id " << taxid << ", start=" << start << " end = " <<end<< std::endl;
			if(nodes.count(taxid) > 0) {
				input_ids.insert(taxid);
			}
			else {
				std::cerr << "Warning: Taxon ID " << taxid << " was not found in taxonomic tree. Skipping." << std::endl;
			}
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad taxon id in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
	list_file.close();

	std::ostream * out_stream;
	if(out_filename.length()>0) {
    std::ofstream * ofs = new std::ofstream();
    ofs->open(out_filename);
    if(!(*ofs)) {  error("Could not open file " + out_filename + " for writing"); exit(EXIT_FAILURE); }
    out_stream = ofs;
  }
  else {
    out_stream = &std::cout;
  }


	for(const TaxonId n : input_ids) {
		output_ids.insert(n);
		if(parent2children.count(n) > 0) {
			dfs(n, parent2children, output_ids);
		}
	}

	for(const TaxonId id : output_ids) {
		*out_stream << id << "\n";
	}
	out_stream->flush();
  if(out_filename.length()>0) {
    ((std::ofstream*)out_stream)->close();
    delete out_stream;
  }


}

void usage(char *progname) {
	fprintf(stderr, "Copyright 2018,2019 Peter Menzel\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -i input.txt [-o output.txt]\n", progname);
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file.\n");
	fprintf(stderr, "   -i FILENAME   Name of input file.\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -v            Enable verbose mode.\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	exit(EXIT_FAILURE);
}

