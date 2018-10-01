#ifndef TAXONOMY_TOOLS_UTIL_H
#define TAXONOMY_TOOLS_UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <set>
#include <time.h>
#include <unordered_map>
#include <fstream>
#include <limits.h>

void error(const std::string e);

void strip(std::string &s);

bool isalpha(const char & c);

void parse_accession2taxid(std::unordered_map<std::string,uint64_t> &, std::ifstream &);

void parseNodesDmp(std::unordered_map<uint64_t,uint64_t> &, std::ifstream &);

void parseNodesDmpWithRank(std::unordered_map<uint64_t,uint64_t> &, std::unordered_map<uint64_t,std::string> &, std::ifstream &);

void parseNamesDmp(std::unordered_map<uint64_t,std::string> &, std::ifstream &);

std::string getTaxonNameFromId(const std::unordered_map<uint64_t,std::string> &, uint64_t, const std::string &);

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const std::unordered_map<uint64_t,uint64_t> &, const std::string &, const std::string &);

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const std::unordered_map<uint64_t,uint64_t> &, uint64_t, uint64_t);

#endif
