#include <string>
#include <vector>

namespace gafanchor {
using std::string, std::vector;
class GAFAnchor {
	private:
		// query info
		string qname;
		int qlength, qstart, qend;
		// path info
		bool pstrand; // true for +, false for -
		vector<int> path; // ids are the 0-based node index in the graph
		vector<bool> orientations; // true for +, false for -
		int plength, pstart, pend;
		// we ignore residue matches, alignment block length, mapping quality

	public:
		GAFAnchor(string qname, int qlength, int qstart, int qend, vector<int> path, int plength, int pstart, int pend, bool reverse = false)
		{
			this->qname = qname;
			this->qlength = qlength;
			this->qstart = qstart;
			this->qend = qend;
			this->pstrand = true;
			this->path = path;
			if (reverse)
				this->orientations = vector<bool>(path.size(), false);
			else
				this->orientations = vector<bool>(path.size(), true);
			this->plength = plength;
			this->pstart = pstart;
			this->pend = pend;

			assert(pend <= plength);
		}

		// empty constructor
		GAFAnchor()
		{
		}

		bool operator<(const GAFAnchor &a) const
		{
			// TODO: does this order matter?
			return tie(qstart, pstart, qend, pend, path) <
				tie(a.qstart, a.pstart, a.qend, a.pend, a.path);
		}

		bool operator==(const GAFAnchor &a) const
		{
			// TODO: does this order matter?
			return tie(qstart, pstart, qend, pend, path) ==
				tie(a.qstart, a.pstart, a.qend, a.pend, a.path);
		}

		int get_query_start() const
		{
			return qstart;
		}

		string get_query_id() const
		{
			return qname;
		}

		int get_path_length() const
		{
			return path.size();
		}

		int get_length() const
		{
			return qend - qstart;
		}

		int start_distance_query(GAFAnchor &a) const
		{
			if (qstart <= a.qstart)
				return a.qstart - qstart - 1;
			else
				return qstart - a.qstart - 1;
		}

		static int gap_query(GAFAnchor &a1, GAFAnchor &a2)
		{
			assert(a1.qend <= a2.qstart);
			return (a2.qstart - a1.qend);
		}

		void reverse()
		{
			int newqstart = qlength - qend;
			int newqend = qlength - qstart;

			qstart = newqstart;
			qend = newqend;

			int newpstart = plength - pend;
			int newpend = plength - pstart;

			assert(newpstart >= 0);
			pstart = newpstart;
			pend = newpend;

			std::reverse(path.begin(), path.end());
			for (long unsigned int i = 0; i < orientations.size(); i++)
				orientations[i] = !orientations[i];
		}

		string to_string(const vector<string> &node_ids)
		{
			string out;
			// GAF format
			out += qname + "\t";                                   // query name
			out += std::to_string(qlength) + "\t";                 // query length
			out += std::to_string(qstart) + "\t";                  // query start (0-based, closed)
			out += std::to_string(qend) + "\t";                    // query end (open)
			out += std::string() + ((pstrand) ? "+" : "-") + "\t"; // strand

			for (long unsigned int i = 0; i < path.size(); i++) { // path
				out += std::string() + ((orientations[i]) ? ">" : "<"); 
				out += node_ids[path[i]];
			}
			out += "\t";

			out += std::to_string(plength) + "\t"; // path length
			out += std::to_string(pstart) + "\t";  // start position on the path
			out += std::to_string(pend) + "\t";    // end position on the path
			out += "0\t0\t255";    // FIXME

			return out;
		}
};
}
