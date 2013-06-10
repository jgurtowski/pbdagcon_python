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

// Prints out the merged graph
void mergeGraph(Alignment *algs, const std::string backbone) {
    AlnGraphBoost ag(backbone); 
    for (int i=0; i < 5; i++) {
        normalizeGaps(algs[i]);
        ag.addAln(algs[i]);
    }
    ag.mergeNodes();
    ag.printGraph();
}

// Checks for a consistent best path, pipe output to sort | uniq -c
void 
consistencyCheck(Alignment *algs, std::string backbone, bool normalize) {
    for (int i = 0; i < 10000; i++) {
        std::string backbone = "ATATTAGGC";
        //AlnGraphBoost ag(backbone); 
        AlnGraphBoost ag(backbone.length()); 
        for (int i=0;i < 5; i++) {
            Alignment a = normalize ? normalizeGaps(algs[i]) : algs[i];
            ag.addAln(a);
        }

        ag.mergeNodes();
        std::cout << ag.consensus() << std::endl;
    }
}

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
    //ag.printGraph();
    std::cout << ">" << targetId << std::endl;
    std::cout << ag.consensus() << std::endl;
}

int main(int argc, char* argv[]) {
    const std::string backbone = "ATATTAGGC";
    Alignment *algs = new Alignment[5];
    
    algs[0].tstr = "ATATTA---GGC";
    algs[0].qstr = "ATAT-AGCCGGC";

    algs[1].tstr = "ATATTA-GGC";
    algs[1].qstr = "ATAT-ACGGC";

    algs[2].tstr = "AT-ATTA--GGC";
    algs[2].qstr = "ATCAT--CCGGC";

    algs[3].tstr = "ATATTA--G-GC";
    algs[3].qstr = "ATAT-ACCGAG-";

    algs[4].tstr = "ATATTA---GGC";
    algs[4].qstr = "ATAT-AGCCGGC";

    for(int i=0; i < 5; i++) {
        Alignment& ra = algs[i];
        ra.tstart = 1;
    }

    //mergeGraph(algs, backbone);
    //consistencyCheck(algs, backbone, false);
    //consistencyCheck(algs, backbone, true);
    if (argc == 2)
        alnFileConsensus(argv[1]);
        
    return 0;
}
