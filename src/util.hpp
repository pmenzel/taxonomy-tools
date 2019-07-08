#ifndef TAXONOMY_TOOLS_UTIL_H
#define TAXONOMY_TOOLS_UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <set>
#include <time.h>
#include <unordered_map>
#include <map>
#include <fstream>
#include <limits.h>

using TaxonId = uint64_t;
using Count = uint32_t;
using Accession = std::string;
using Rank = std::string;
using TaxonName = std::string;

using TaxonId2Count = std::unordered_map<TaxonId, Count>;
using TaxonSet = std::set<TaxonId>;
using TaxTree = std::unordered_map<TaxonId,TaxonId>;
using pTaxTree = std::unordered_map<TaxonId,TaxonId> *;

void error(const std::string e);

void strip(std::string &s);

bool isalpha(const char & c);

void parse_accession2taxid(std::unordered_map<Accession,TaxonId> &, std::ifstream &);

void parseNodesDmp(TaxTree &, std::ifstream &);

void parseNodesDmpWithRank(TaxTree &, std::unordered_map<TaxonId,Rank> &, std::ifstream &);

void parseNodesDmpTopDown(std::multimap<TaxonId, TaxonId> &, std::ifstream &);

void parseNamesDmp(std::unordered_map<TaxonId, TaxonName> &, std::ifstream &);

void parseExclusionFile(const TaxTree &, TaxonSet &, std::ifstream &);

TaxonName getTaxonNameFromId(const std::unordered_map<TaxonId,TaxonName> &, const TaxonId &, const std::string &);

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const TaxTree &, const TaxonName &, const TaxonName &);

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const TaxTree &, TaxonId, TaxonId);

TaxonId lca_from_ids(const TaxTree &, const std::set<TaxonId> &);

TaxonId lca_two(const TaxTree & nodes, TaxonId node1, TaxonId node2);

TaxonId lowest_from_ids(const TaxTree &, const std::set<TaxonId> &);

TaxonId heaviest_path(const TaxTree &, const TaxonId2Count &);

#endif
