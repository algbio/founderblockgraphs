#include <iostream>
#include <filesystem>
#include <fstream> // std::ofstream
#include <sstream> // std::istringstream
#include <string>
#include <tuple>
#include <zlib.h>  
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <sdsl/util.hpp> // init_support for rank queries
#include <sdsl/config.hpp> // util things

#include "command-line-parsing/locate_patterns.h" // gengetopt-generated parser
#include "gfakluge.hpp"
#include "kseq.h"
#include "gafanchor.hpp"

#define UNIQUE_CHAR '$'
#define CONCAT_SEPARATOR '#'
#define SUPERSOURCE_ID "supersource"

//#define LOCATE_PATTERNS_DEBUG

// GFAKluge setup
#define GFAK_LINK_TYPE 1
using gfak::GFAKluge;

// kseq setup
KSEQ_INIT(gzFile, gzread)

using std::cerr, std::cout, std::endl, std::flush;
using gafanchor::GAFAnchor;
using std::vector, std::istringstream, std::string, std::pair, std::tuple, std::tie, std::get;

typedef sdsl::csa_wt<> csa_type;
typedef csa_type::size_type csa_size_type;
struct repeat_free_index {
	//TODO handle reverse complement nodes
	vector<string>              node_ids;
	vector<int>                 node_lengths;
	vector<pair<int,int>>       edges;
	vector<int>                 overlaps;
        csa_type                    csa;
	sdsl::bit_vector            edge_leader_positions;
        sdsl::bit_vector            b_positions;
        sdsl::bit_vector            e_positions;
	sdsl::rank_support_v5<>     edge_leader_rank1_support;
        sdsl::rank_support_v5<>     b_rank1_support;
        sdsl::rank_support_v5<>     e_rank1_support;
        sdsl::select_support_mcl<>  edge_leader_select1_support;
        sdsl::select_support_mcl<>  b_select1_support;
        sdsl::select_support_mcl<>  e_select1_support;
};

// add supersource (id: SUPERSOURCE_ID, sequence: "UNIQUE_CHAR") to the graph and connect it to all sources
void normalize_graph(GFAKluge &graph);

void build_index(GFAKluge &graph, struct repeat_free_index &index);

// locate ALL occurrences of s in graph and store them as GAFAnchors in res; initial data in res is destroyed; output is 0 is no occurrence, >0 if pattern spans 3+ nodes and is equal to count, <0 if pattern spans <=2 nodes and is upper bound of count
int locate(const repeat_free_index &index, const string &s_id, const string &s, vector<GAFAnchor> &res);

int main(int argc, char* argv[])
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);

	if (argsinfo.inputs_num == 0)
		{cmdline_parser_print_help(); exit(1);};
	if (argsinfo.inputs_num == 1)
		{cerr << argv[0] << ": missing patterns file" << endl; exit(1);};
	if (argsinfo.inputs_num == 2)
		{cerr << argv[0] << ": missing output file" << endl; exit(1);};
	if (argsinfo.inputs_num > 3)
		{cerr << argv[0] << ": too many arguments" << endl; exit(1);};

	// check and open output file
	std::filesystem::path outputpath {argsinfo.inputs[2]};
	std::ofstream outputfs;
	if (std::filesystem::exists(outputpath)) {
		if (argsinfo.overwrite_flag) {
			outputfs = std::ofstream(outputpath, std::ios::out | std::ios::trunc);
		} else {
			cerr << argv[0] << ": output file already exists" << endl;
			exit(1);
		}
	} else {
		outputfs = std::ofstream(outputpath);
	}
	if (!outputfs) {cerr << argv[0] << ": error opening output file " << outputpath << endl; exit(1);};

	repeat_free_index index;
	{ // scope for gg
		cerr << "Reading the graph..." << flush;
		GFAKluge gg;
		gg.parse_gfa_file(argsinfo.inputs[0]);
		cerr << " done." << endl;
#ifdef LOCATE_PATTERNS_DEBUG 
		int debug_edges = 0;
		for (auto const& [key,neighbors] : gg.get_seq_to_edges())
			debug_edges += neighbors.size();
		cerr << "DEBUG: read graph with " << gg.get_name_to_seq().size() << " nodes and " << debug_edges << " edges\n";
#endif

		cerr << "Adding a supersource to the graph..." << flush;
		normalize_graph(gg);
		cerr << " done." << endl;
#ifdef LOCATE_PATTERNS_DEBUG 
		cerr << "DEBUG: normalized graph is\n" << gg;
#endif

		cerr << "Indexing the graph...";
		build_index(gg, index);
		cerr << " done.\n";
	}

	// check and open patterns file
	gzFile gzfp = gzopen(argsinfo.inputs[1], "r");
	if (gzfp == NULL) {cerr << argv[0] << ": error opening patterns file " << argsinfo.inputs[1] << endl; exit(1);};
	kseq_t *seq = kseq_init(gzfp);

	// locate
	cerr << "Locate\n";
	while (kseq_read(seq) >= 0) {
		const string id(seq->name.s);
		const string sequence(seq->seq.s);

		int res;
		vector<GAFAnchor> occs;
		if ((res = locate(index, id, sequence, occs)) != 0) {
			cout << id << ": occurs " << (std::abs(res)) << " times" << ((res > 0) ? "" : " (at most)") << "\n";
			for (auto &o : occs)
				outputfs << o.to_string(index.node_ids) << "\n";
		} else {
			cout << id << ": does not occur!\n";
		}
	}

	kseq_destroy(seq);  
	gzclose(gzfp); 
	return 0;
}

void normalize_graph(GFAKluge &graph)
{
	//TODO handle reverse complement nodes
	// find sources
	std::map<string,bool> is_source;
	for (auto const& [id,_] : graph.get_name_to_seq())
		is_source[id] = true;
	for (auto const& [id,edges] : graph.get_seq_to_edges())
		for (auto const &e : edges) {
			assert(e.type == GFAK_LINK_TYPE);
			is_source[e.sink_name] = false;
		}

	//  add dummy supersource
	gfak::sequence_elem s;
	s.sequence = UNIQUE_CHAR;
	s.name = SUPERSOURCE_ID;
	graph.add_sequence(s);

	for (auto const &[id,b] : is_source) {
		if (b) {
			gfak::edge_elem e;
			e.type = GFAK_LINK_TYPE;
			e.source_name = SUPERSOURCE_ID;
			e.sink_name = id;
			e.source_orientation_forward = true;
			e.sink_orientation_forward = true;
			e.alignment = "0M";
			graph.add_edge(e.source_name, e);
		}
	}
}

void build_index(GFAKluge &g, struct repeat_free_index &res)
{
	// TODO check performance of this
	auto get_name_to_seq = g.get_name_to_seq();

	std::map<string,int> node_index;
	res.node_ids.reserve(get_name_to_seq.size());
	res.node_lengths.reserve(get_name_to_seq.size());
	for (int index = 0; auto const& [id,s_elem] : get_name_to_seq) {
		res.node_ids.push_back(id);
		res.node_lengths.push_back(s_elem.sequence.size());
		node_index[id] = index++;
	}

	for (auto const & [id,edges] : g.get_seq_to_edges()) {
		for (auto const &e : edges) {
			assert(e.source_orientation_forward and e.sink_orientation_forward);
			assert(e.type == GFAK_LINK_TYPE);
			assert(node_index.contains(e.source_name) and node_index.contains(e.sink_name));
			res.edges.push_back({ node_index[e.source_name], node_index[e.sink_name] });

			assert(e.alignment.ends_with("M"));
			res.overlaps.emplace_back(std::stoi(e.alignment.substr(0, e.alignment.size() - 1)));
		}
	}

#ifdef LOCATE_PATTERNS_DEBUG 
	cerr << "DEBUG: edges are ";
	for (auto const [u,v] : res.edges)
		cerr << u << "," << v << " ";
	cerr << "\n";
	cerr << "DEBUG: overlaps are ";
	for (auto const ov : res.overlaps)
		cerr << ov << " ";
	cerr << "\n";
#endif

	// build csa of edge concatenation
	string edge_concat;
	edge_concat += CONCAT_SEPARATOR;
	for (unsigned long int i = 0; i < res.edges.size(); i++) {
		const int u = res.edges[i].first;
		const int v = res.edges[i].second;
		const int ov = res.overlaps[i];

		const string first = get_name_to_seq[res.node_ids[u]].sequence;
		edge_concat += first;

		const string second = get_name_to_seq[res.node_ids[v]].sequence.substr(ov);
		edge_concat += second;

		edge_concat += CONCAT_SEPARATOR;
	}
#ifdef LOCATE_PATTERNS_DEBUG 
	cerr << "DEBUG: edge_concat and edge_leader_positions are\n" << edge_concat << "\n";
#endif
	sdsl::construct_im(res.csa, edge_concat, 1);
	edge_concat.clear();

	// mark every separator character that starts an edge in edge_concat
	res.edge_leader_positions = sdsl::bit_vector(res.csa.size(), 0);
	for (unsigned long int i = 0, consumed = 0; i < res.edges.size(); i++) {
		const int u = res.edges[i].first;
		const int v = res.edges[i].second;
		const int ov = res.overlaps[i];

		res.edge_leader_positions[consumed] = 1;

		// TODO check performance of this
		consumed += 1 + get_name_to_seq[res.node_ids[u]].sequence.size() + get_name_to_seq[res.node_ids[v]].sequence.size() - ov;
	}
#ifdef LOCATE_PATTERNS_DEBUG 
	cerr << res.edge_leader_positions << "\n";
#endif
	res.edge_leader_rank1_support = sdsl::rank_support_v5<>(&res.edge_leader_positions);
	res.edge_leader_select1_support = sdsl::select_support_mcl<>(&res.edge_leader_positions);

	// find the lex range of every node label
	res.b_positions = sdsl::bit_vector(res.csa.size(), 0);
	res.e_positions = sdsl::bit_vector(res.csa.size(), 0);
	for (unsigned long int i = 0; i < res.node_ids.size(); i++) {
		if (SUPERSOURCE_ID == res.node_ids[i]) continue;

		csa_size_type l,r;
		const string &s = get_name_to_seq[res.node_ids[i]].sequence;
		sdsl::backward_search(res.csa, 0, res.csa.size()-1, s.begin(), s.end(), l, r);
		assert(l <= r and res.b_positions[l] == 0 and res.e_positions[r] == 0);
		res.b_positions[l] = 1;
		res.e_positions[r] = 1;

#ifdef LOCATE_PATTERNS_DEBUG 
		cerr << "DEBUG: lex interval of " << s << " is " << l << ".." << r << "\n";
#endif
	}
	res.b_rank1_support = sdsl::rank_support_v5<>(&res.b_positions);
	res.e_rank1_support = sdsl::rank_support_v5<>(&res.e_positions);
	res.b_select1_support = sdsl::select_support_mcl<>(&res.b_positions);
	res.e_select1_support = sdsl::select_support_mcl<>(&res.e_positions);

#ifdef LOCATE_PATTERNS_DEBUG 
	cerr << "DEBUG: compressed suffix array is " << std::endl;
	cerr << sdsl::extract(res.csa, 0, res.csa.size()-1) << std::endl;
	cerr << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << std::endl;
	sdsl::csXprintf(cerr, "%2I %2S %3s %3P %2p %3B   %:1T", res.csa);

	cerr << "DEBUG: b_positions and e_positions are\n";
	for (auto b : res.b_positions) cerr << (b ? "b" : " ");
	cerr << "\n";
	for (auto e : res.e_positions) cerr << (e ? "e" : " ");
	cerr << "\n";
#endif
}

int locate_edge(const repeat_free_index &index, const csa_size_type l)
{
	const int pos = index.csa[l]; // position in the text via the suffix array
	const int edge = index.edge_leader_rank1_support(pos+1)-1;
	return edge;
}

tuple<int,int> locate_edge_and_position(const repeat_free_index &index, const csa_size_type l)
{
	const int pos = index.csa[l]; // position in the text via the suffix array
	const int edge = index.edge_leader_rank1_support(pos)-1;
	const int edgestartpos = index.edge_leader_select1_support(edge+1)+1;
	return tuple(edge, pos - edgestartpos);
}

int locate(const repeat_free_index &index, const string &s_id, const string &s, vector<GAFAnchor> &res)
{
	res.clear();

	csa_size_type l = 0, r = index.csa.size() - 1, ll, rr;
	bool expanded = false;
	int firstcount = 0, count = 0;
	vector<tuple<csa_size_type,csa_size_type,unsigned long int>> expanded_ranges; // l, r, i
	for (unsigned long int i = s.size() - 1; i != (unsigned long int)0-1; i--) {
		const int newcount = sdsl::backward_search(index.csa, l, r, s[i], ll, rr);

		if (newcount > 0) {
			// continue consuming the string
			l = ll;
			r = rr;
			count = newcount;
		} else {
			// 1a. no separator character can be read
			if (backward_search(index.csa, l, r, CONCAT_SEPARATOR, ll, rr) == 0)
				return 0;

			csa_size_type const rank(index.b_rank1_support(l + 1));

			// 1b. string is not prefix of any edge label
			if (rank == 0)
				return 0;

			// expand range
			ll = index.b_select1_support(rank);
			rr = index.e_select1_support(rank);
			assert(ll <= rr);

			// 1c. string is not prefix of any edge label
			if (!(ll <= l and r <= rr))
				return 0;

			if (!expanded) {
				firstcount = count;
				expanded = true;
				expanded_ranges.push_back({ l, r, i+1 });
			} else {
				expanded_ranges.push_back({ ll, rr, i+1 });
			}
			l = ll;
			r = rr;

			// 1d. expansion was successful but we cannot continue with next character
			count = backward_search(index.csa, l, r, s[i], ll, rr);
			if (count == 0)
				return 0;

			l = ll;
			r = rr;
		}
	}

	// 2. whole string was consumed
	if (expanded) {
		assert(expanded_ranges.size() > 0);
		std::reverse(expanded_ranges.begin(), expanded_ranges.end());
		vector<int> path;
		int plength = 0;
		path.push_back(-1);
		for (auto [_,r,__] : expanded_ranges) {
			// tricky choice!
			const auto edge = locate_edge(index, r);
			const auto u = index.edges[edge].first;
			path.push_back(u);
			plength += index.node_lengths[u] - index.overlaps[edge];
		}
		path.push_back(-1);

		for (auto kstart = l; kstart <= r; kstart++) {
			const auto [edgestart,posstart] = locate_edge_and_position(index, kstart);
			const auto ovstart = index.overlaps[edgestart];
			const auto ustart = index.edges[edgestart].first;
			const auto ustartlength = index.node_lengths[ustart];
			assert(posstart < ustartlength);
			plength += ustartlength - ovstart;
			path[0] = ustart;
			for (auto kend = get<0>(expanded_ranges.back()); kend <= get<1>(expanded_ranges.back()); kend++) {
				const auto edgeend = locate_edge(index, kend);
				//const auto uend = index.edges[edgeend].first;
				const auto vend = index.edges[edgeend].second;
				//const auto ovend = index.overlaps[edgeend];
				path[path.size()-1] = vend;
				res.push_back(GAFAnchor(
						s_id, //qname
						s.size(), //qlength
						0, //qstart
						s.size(), //qend
						path, //path
						plength + index.node_lengths[vend], //plength
						posstart, //pstart
						posstart + s.size())); //pend
			}
			plength -= ustartlength - ovstart;
		}
		return firstcount * count;
	} else {
		// all occurrences are inside edge labels
		for (auto k = l; k <= r; k++) {
			const auto [edge, pos] = locate_edge_and_position(index, k);
			const auto u = index.edges[edge].first;
			const auto v = index.edges[edge].second;
			const auto ov = index.overlaps[edge];
			const auto ulength = index.node_lengths[u];

			if (pos + s.size() - 1 < (unsigned long int)ulength) {
				// A. fully inside first node u
				res.push_back(GAFAnchor(
						s_id, //qname
						s.size(), //qlength
						0, //qstart
						s.size(), //qend
						vector({u}), //path
						ulength, //plength
						pos, //pstart
						pos + s.size())); //pend
			} else if (pos >= ulength - ov) {
				// B. fully inside second node
				res.push_back(GAFAnchor(
						s_id, //qname
						s.size(), //qlength
						0, //qstart
						s.size(), //qend
						vector({v}), //path
						index.node_lengths[v], //plength
						pos - ulength + ov, //pstart
						pos - ulength + ov + s.size())); //pend
			} else {
				// C. fully inside edge
				res.push_back(GAFAnchor(
						s_id, //qname
						s.size(), //qlength
						0, //qstart
						s.size(), //qend
						vector({u,v}), //path
						ulength + index.node_lengths[v], //plength
						pos, //pstart
						pos + s.size())); //pend
			}
		}
		return -count;
	}
}

