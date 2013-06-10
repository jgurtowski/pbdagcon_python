#include <cstdint>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"

// Takes an alignment file (currently blasr -m 5) containing a single target, 
// generates a consensus and prints it to stdout.
void alnFileConsensus(const std::string alnFile) {
    std::ifstream file(alnFile);
    Alignment aln;
    file >> aln;
    AlnGraphBoost ag(aln.tlen); 
    std::string targetId = aln.tid;
    uint32_t len = aln.tlen;
    while (file >> aln) {
        Alignment a = normalizeGaps(aln);
        assert(aln.tlen = len);
        ag.addAln(a);
    }
    ag.mergeNodes();
    std::cout << ">" << targetId << std::endl;
    std::cout << ag.consensus() << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        alnFileConsensus(argv[1]);
        
    return 0;
}
