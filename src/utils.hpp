#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <vector>
#include <iostream>
#include <filesystem>

#include "kseq.h"

#define GAP_CHARACTER '-'
#define SEPARATOR_CHARACTER '#'

namespace fbg::utils {

using std::cerr, std::endl, std::string, std::filesystem::path, std::vector;

// kseq setup
KSEQ_INIT(gzFile, gzread);
path currentp("");
unsigned long long rows = 0;
gzFile fp = NULL;
kseq_t *seq = NULL;

bool contains_chars(const string &s, const string &chars) {
	for (char c : chars) {
		if (s.find(c) != std::string::npos)
			return true;
	}
	return false;
}

unsigned long long count_ones(const vector<bool> v)
{
	unsigned long long res = 0;
	for (const bool b : v) res += b;
	return res;
}

void open_msa_file(const path &p) {
	fp = gzopen(p.c_str(), "r");
	if (fp == NULL) { cerr << "failed to read input MSA file " << p << endl; exit(1); }
	seq = kseq_init(fp);
	currentp = p;
	rows = 0;
}

/* requires: open_msa_file to have been called 
 * modifies: input strings s and id
 * returns: 1 if success, 0 if file has (probably) ended */
int get_msa_line(string &s, string &id) {
	assert(fp != NULL and seq != NULL);
	long long len;
	if ((len = kseq_read(seq)) >= 0) {
		s = string(seq->seq.s);
		id = string(seq->name.s);
		rows += 1;
		return 1;
	}
	return 0;
}

unsigned long long get_rows() {
	return rows;
}

path get_path() {
	return currentp;
}

void close_msa_file() {
	gzclose(fp);
	currentp = path("");
}

} // namespace fbg::utils
#endif // ifndef UTILS_HPP
