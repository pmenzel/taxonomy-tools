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
#include <assert.h>

#include "util.hpp"

void usage(char *progname);

int main(int argc, char **argv) {

	TaxTree nodes;
	std::unordered_map<TaxonId,unsigned int> node2depth;
	TaxonSet excluded_ids;

	std::string nodes_filename;
	std::string in_filename;
	std::string out_filename;
	std::string exclusion_filename;
	std::string mode = "lca";

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
			case 'm':
				mode = optarg; break;
			case 'i':
				in_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'o':
				out_filename = optarg; break;
			case 'e':
				exclusion_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(in_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }
	if(mode != "lca" && mode != "lowest" && mode != "path") { error("Error: Mode must either be 'lca', 'lowest', or 'path'."); usage(argv[0]); }

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
	parseNodesDmp(nodes,nodes_file);
	nodes_file.close();

	if(exclusion_filename.length() > 0) {
		std::ifstream exclusion_file;
		exclusion_file.open(exclusion_filename);
		if(!exclusion_file) { error("Could not open file " + exclusion_filename); exit(EXIT_FAILURE); }
		parseExclusionFile(nodes, excluded_ids, exclusion_file);
		exclusion_file.close();
	}

	std::ifstream inputfile;
	inputfile.open(in_filename);
	if(!inputfile.is_open()) { error("Could not open file " + in_filename); exit(EXIT_FAILURE); }

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

	std::set<TaxonId> ids;
	TaxonId2Count id2count;
	std::string line;
	while(getline(inputfile, line)){
		if(line.length() == 0) { continue; }
		line += "\t";
		ids.clear();
		id2count.clear();
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

			if(nodes.find(id) != nodes.end()) {
				if(mode=="path") {
					id2count[id]++;
				}
				else {
					ids.insert(id);
				}
			}
			else {
				std::cerr << "Warning: Taxon ID " << id << " not found in nodes.dmp, line: " << line << std::endl;
			}
			start = end+1;
		}

		if(ids.size()>0 || id2count.size() > 0) {
			TaxonId result = 0;
			if(mode=="lca") {
				if(ids.size()==1) {
					result = *(ids.begin());
				}
				else if(excluded_ids.size()>0) {
					TaxonSet filtered_ids;
					// go through ids and add all ids that are not in exclusion list to new set
					for(auto const & it_id : ids) {
						bool keep = true;
						for(auto const & it_excl_id : excluded_ids) {
							if(is_ancestor(nodes, it_excl_id, it_id)) {
								keep = false;
								break;
							}
						}
						if(keep) {
							filtered_ids.emplace(it_id);
						}
					}
					// if filtered set is empty, then do lca on original ids
					if(filtered_ids.size()==0) {
						result = lca_from_ids(nodes, ids);
					}
					else if(filtered_ids.size()==1) {
						result = *(filtered_ids.begin());
					}
					else if(filtered_ids.size()==2) {
						result = lca_two(nodes, *(filtered_ids.begin()), *(++filtered_ids.begin()));
					}
					else { // else do lca on new set
						result = lca_from_ids(nodes, filtered_ids);
					}
				}
				else if(ids.size()==2) {
						result = lca_two(nodes, *(ids.begin()), *(++ids.begin()));
				}
				else {
					result = lca_from_ids(nodes, ids);
				}
			}
			else if(mode=="path") {
				if(id2count.size()==1) {
					result = id2count.begin()->second ;
				}
				else if(excluded_ids.size()>0) {
					TaxonId2Count filtered_id2count;
					// go through ids and add all ids that are not in exclusion list to new set
					for(auto const & it_id : id2count) {
						bool keep = true;
						for(auto const & it_excl_id : excluded_ids) {
							if(is_ancestor(nodes, it_excl_id, it_id.first)) {
								keep = false;
								break;
							}
						}
						if(keep) {
							filtered_id2count.emplace(it_id.first,it_id.second);
						}
					}
					// if filtered set is empty, then do lca on original ids
					if(filtered_id2count.size()==0) {
						result = heaviest_path(nodes, id2count);
					}
					else { // else do lca on new set
						result = heaviest_path(nodes, filtered_id2count);
					}
				}
				else {
						result = heaviest_path(nodes, id2count);
				}
			}
			else {
				assert(conflict=="lowest");
				if(ids.size()==1) {
					result = *(ids.begin());
				}
				else if(excluded_ids.size()>0) {
					TaxonSet filtered_ids;
					// go through ids and add all ids that are not in exclusion list to new set
					for(auto const & it_id : ids) {
						bool keep = true;
						for(auto const & it_excl_id : excluded_ids) {
							if(is_ancestor(nodes, it_excl_id, it_id)) {
								keep = false;
								break;
							}
						}
						if(keep) {
							filtered_ids.emplace(it_id);
						}
					}
					// if filtered set is empty, then do lca on original ids
					if(filtered_ids.size()==0) {
						result = lowest_from_ids(nodes, ids);
					}
					else { // else do lca on new set
						result = lowest_from_ids(nodes, filtered_ids);
					}
				}
				else {
					result = lowest_from_ids(nodes, ids);
				}
			}

			if(result == 0) { std::cerr << "Warning: Could not determine LCA in line " << line << std::endl; }
			(*out_stream) << result << "\n";
		}
		else {
			std::cerr << "No taxon IDs found in line " << line << std::endl;
		}
	}

	inputfile.close();
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
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -i input.tsv -o output.txt\n", progname);
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file.\n");
	fprintf(stderr, "   -i FILENAME   Name of tab-delimited input file.\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -e FILENAME   Name of file with exlusion taxon IDs.\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -m STRING     Mode, must be either 'lca' (default), 'lowest', or 'path'.\n");
	fprintf(stderr, "   -v            Enable verbose mode.\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	exit(EXIT_FAILURE);
}

