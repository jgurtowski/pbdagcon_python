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

// Takes an alignment file (currently blasr -m 5) sorted by target, 
// generates a consensus for each target that has > 8x coverage and 
// prints it to stdout in fasta format.
void alnFileConsensus(const std::string alnFile) {
    std::ifstream file(alnFile);
    int coverage = 0;

    // Initialize the graph with the first target
    Alignment aln;
    file >> aln;
    std::string targetId = aln.tid;
    AlnGraphBoost ag(aln.tlen);
    uint32_t len = aln.tlen;
    if (aln.qid != targetId) {
        aln = normalizeGaps(aln);
        ag.addAln(aln);
        coverage++;
    }

    while (file >> aln) {
        if (aln.tid != targetId){
            // we're on to a new target, generate consensus ...
            if (coverage > 8) {
                ag.mergeNodes();
                std::cout << ">" << targetId << std::endl;
                std::cout << ag.consensus() << std::endl;
            }

            // ... and then move on to the next target
            ag = AlnGraphBoost(aln.tlen);
            len = aln.tlen; 
            targetId = aln.tid;
            coverage = 0;
        }
        // sanity check
        assert(aln.tlen = len);

        // skip self hits
        if (aln.qid == targetId) continue;
        aln = normalizeGaps(aln);
        ag.addAln(aln);
        coverage++;
    }
    // print out final target consensus
    if (coverage > 8) {
        ag.mergeNodes();
        std::cout << ">" << targetId << std::endl;
        std::cout << ag.consensus() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc == 2)
        alnFileConsensus(argv[1]);
        
    return 0;
}
