#include "util.hpp"

void error(const std::string e) {
	std::cerr << "Error: " << e << std::endl << std::endl;
}

inline bool isalpha(const char & c) {
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

void strip(std::string & s) {
	for(auto it = s.begin(); it!=s.end(); ++it) {
		if(!isalpha(*it)) {
			s.erase(it);
			it--;
		}
	}
}

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const TaxTree & nodes, const std::string & id1, const std::string & id2) {
	TaxonId node1;
	TaxonId node2;
	try {
		node1 = stoul(id1);
		node2 = stoul(id2);
	}
	catch(const std::invalid_argument& ia) {
		std::cerr << "Error: Bad number in taxon id" << std::endl;
		return false;
	}
	catch (const std::out_of_range& oor) {
		std::cerr << "Error: Bad number (out of range error) in taxon id" << std::endl;
		return false;
	}
	return is_ancestor(nodes,node1,node2);
}

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const TaxTree & nodes, TaxonId node1, TaxonId node2) {
	if(nodes.count(node1)==0) { std::cerr << "Taxon ID " << node1 << " not found in taxonomy!" << std::endl; return false; }
	if(nodes.count(node2)==0) { std::cerr << "Taxon ID " << node2 << " not found in taxonomy!" << std::endl; return false; }
	/* climb up from node 2 and return true if encountering node 1 */
	while(nodes.count(node2)>0 && node2 != nodes.at(node2)) {
		if(node2==node1) {
			return true;
		}
		node2 = nodes.at(node2);
	}
	return false;
}

void parse_accession2taxid(std::unordered_map<Accession,TaxonId> & acc2taxid, std::ifstream & acc2taxid_file) {
	acc2taxid.reserve(150e6);
	std::string line;
	getline(acc2taxid_file, line); // skip header line
  while(getline(acc2taxid_file, line)) {
		if(line.length() == 0) { continue; }
    try {
      size_t start = line.find('\t',0);
      size_t end = line.find("\t",start+1);
      //std::string acc = line.substr(start+1,end-start-1);
      TaxonId taxid = strtoul(line.c_str() + end + 1,NULL,10);
      if(taxid == ULONG_MAX) {
        std::cerr << "Found bad taxid number (out of range error) in line: " << line << std::endl;
        continue;
      }
      //std::cerr << acc << "\t" << taxid << "\n";
      acc2taxid.emplace(line.substr(start+1,end-start-1),taxid);
    }
    catch(const std::invalid_argument& ia) {
      std::cerr << "Found bad identifier in line: " << line << std::endl;
    }
    catch (const std::out_of_range& oor) {
      std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
    }
  }

}

void parseNodesDmp(TaxTree & nodes, std::ifstream & nodes_file) {
	nodes.reserve(2e6);
	std::string line;
	while(std::getline(nodes_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			size_t end = line.find_first_not_of("0123456789");
			TaxonId node = stoul(line.substr(0,end));
			size_t start = line.find_first_of("0123456789",end);
			end = line.find_first_not_of("0123456789",start+1);
			TaxonId parent = (TaxonId)stoul(line.substr(start,end-start));
			nodes.emplace(node,parent);
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad number in line: " << line << std::endl;
		}
		catch(const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
}

void parseNodesDmpWithRank(TaxTree & nodes, std::unordered_map<TaxonId,Rank> & node2rank, std::ifstream & nodes_file) {
	nodes.reserve(2e6);
	std::string line;
	while(std::getline(nodes_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			size_t end = line.find_first_not_of("0123456789");
			//cerr << "end=" << end << "\t";
			TaxonId node = stoul(line.substr(0,end));
			size_t start = line.find_first_of("0123456789",end);
			//cerr << "start=" << start <<"\t";
			end = line.find_first_not_of("0123456789",start+1);
			//cerr << "end=" << end <<"\t";
			TaxonId parent = stoul(line.substr(start,end-start));
			start = line.find_first_of("abcdefghijklmnopqrstuvwxyz",end);
			//cerr << "start=" << start <<","<< line[start] << "\t";
			end = line.find_first_not_of("abcdefghijklmnopqrstuvwxyz ",start);
			//cerr << "end=" << end << "\t";
			std::string rank = line.substr(start,end-start);
			nodes.emplace(node,parent);
			node2rank.emplace(node,rank);
			//cerr << node << "\t" << parent << "\t" <<rank << "\n";

		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad number in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
}

void parseNamesDmp(std::unordered_map<TaxonId,TaxonName> & names, std::ifstream & names_file) {
	std::string line;
	while(std::getline(names_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			if(line.find("scientific name")==std::string::npos) continue;
			size_t start = line.find_first_of("0123456789");
			size_t end = line.find_first_not_of("0123456789",start);
			TaxonId node_id = stoul(line.substr(start,end-start));
			start = line.find_first_not_of("\t|",end);
			end = line.find_first_of("\t|",start+1);
			TaxonName name = line.substr(start,end-start);
			names.emplace(node_id,name);
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad number in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
}


TaxonName getTaxonNameFromId(const std::unordered_map<TaxonId,TaxonName> & node2name, const TaxonId & id, const std::string & names_filename) {
	TaxonName taxon_name;
	if(node2name.count(id)==0) {
		std::cerr << "Warning: Taxon ID " << id << " is not found in file "<< names_filename << "." << std::endl;
		taxon_name = "taxonid:"; taxon_name += std::to_string(id);
	}
	else {
		taxon_name = node2name.at(id);
	}
	return taxon_name;
}

// not thread-safe!
TaxonId lca_from_ids(const TaxTree & nodes, const std::set<TaxonId> & ids) {

	static std::unordered_map<TaxonId,unsigned int> node2depth;

	size_t num_ids = ids.size();
	if(num_ids == 1) {
		return *(ids.begin());
	}
	TaxonId * leafs = (TaxonId *) calloc(num_ids,sizeof(TaxonId));
	unsigned int shallowest_depth = 100000;
	unsigned int index = 0;
	for(auto it : ids) {
		// check if this id was already seen, then skip it
		if(nodes.count(it)==0) {
			std::cerr << "Warning: Taxon ID " << it << " is not contained in taxonomic tree.\n";
			num_ids--;
			continue;
		}

		leafs[index++] = it;

		//if id is alrady in the depth map then do not add it.
		auto pos = node2depth.find(it);
		if(pos == node2depth.end()) {
			unsigned int depth = 1;
			TaxonId id = it;
			while(nodes.count(id)>0 && id != nodes.at(id)) {
				depth++;
				id = nodes.at(id);
			}
			node2depth.emplace(it,depth);
			//cerr << "Inserting to depth map: " << *it <<" -> " << depth << endl;
			if(depth < shallowest_depth) { shallowest_depth = depth; }
		}
		else if(pos->second < shallowest_depth) {
			shallowest_depth = pos->second;
		}
	}

	if(num_ids<=0) {
		free(leafs);
		return 0;
	}

	//cerr << "shallowest depth = " << shallowest_depth << endl;

	for(int index = 0; index < num_ids; ++index) {
		for(int i = node2depth.at(leafs[index]) - shallowest_depth; i > 0; i--) {
			leafs[index] = nodes.at(leafs[index]);
		}
	}

	while(true) {
		//foreach element in the list, check if id is the same, otherwise go one level up in tree, i.e. one more iteration
		TaxonId first = leafs[0];
		bool found = true;
		for(size_t index = 0; index < num_ids; ++index) {
			if(first != leafs[index]) {
				found = false;
			}
			leafs[index] = nodes.at(leafs[index]);
		}
		if(found) {
			free(leafs);
			return first;
		}
	}
	free(leafs);

}



TaxonId lowest_from_ids(const TaxTree & nodes, const std::set<TaxonId> & ids) {
	size_t num_ids = ids.size();
	if(num_ids == 1) {
		if(nodes.find(*(ids.begin()))==nodes.end()) {
			std::cerr << "Warning: Taxon ID " << *(ids.begin()) << " is not contained in taxonomic tree.\n";
			return 0;
		}
		else {
			return *(ids.begin());
		}
	}
	std::set<TaxonId> s;
	for(auto it : ids) {
		if(nodes.count(it)==0) {
			std::cerr << "Warning: Taxon ID " << it << " is not contained in taxonomic tree.\n";
		}
		else {
			s.insert(it);
		}
	}

	TaxonId lca = lca_from_ids(nodes,s);
	while(s.find(lca) != s.end()) {
		s.erase(lca);
		if(s.size() > 1)
			lca = lca_from_ids(nodes,s);
		else
			return *(s.begin());
	}

	return lca;
}

TaxonId lca_two(const TaxTree & nodes, TaxonId node1, TaxonId node2) {
	std::set<TaxonId> lineage1;
	lineage1.emplace(node1);
	while(nodes.count(node1)>0 && node1 != nodes.at(node1)) {
		lineage1.emplace(nodes.at(node1));
		node1 = nodes.at(node1);
	}

	TaxonId lca = node2;
	do {
		lca = nodes.at(lca);
	} while(lineage1.count(lca)==0 && lca != nodes.at(lca));

	return lca;
}

