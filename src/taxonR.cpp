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
#include <map>
#include <vector>
#include <unordered_set>
#include <climits>
#include <assert.h>
#include <Rcpp.h>

#include "util.hpp"


// [[Rcpp::export]]
SEXP loadNodesDmp(std::string nodes_filename) {

	if(nodes_filename.length() == 0) { std::cerr << "Give file name\n"; exit(0); }
	pTaxTree nodes = new TaxTree();

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
	parseNodesDmp(*nodes,nodes_file);
	nodes_file.close();
	Rcpp::XPtr<TaxTree> p(nodes);

	return p;

}

// [[Rcpp::export]]
TaxonId parent(SEXP nodes_, TaxonId node){
	Rcpp::XPtr<TaxTree> nodes(nodes_);
	if((*nodes).count(node)==0) { throw std::range_error("Taxon id not found in nodes.dmp: "+std::to_string(node)); }
  TaxonId b = (*nodes).at(node);
  return b;
}

// [[Rcpp::export]]
bool is_ancestor(SEXP nodes_, TaxonId node1, TaxonId node2){
	Rcpp::XPtr<TaxTree> nodes(nodes_);
	if((*nodes).count(node1)==0) { throw std::range_error("Taxon id not found in nodes.dmp: "+std::to_string(node1)); }
	if((*nodes).count(node2)==0) { throw std::range_error("Taxon id not found in nodes.dmp: "+std::to_string(node2)); }

  bool r = is_ancestor(*nodes, node1, node2);
  return r;
}



