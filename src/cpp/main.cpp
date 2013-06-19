#include <cstdint>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <log4cpp/Category.hh>
#include <log4cpp/Appender.hh>
#include <log4cpp/FileAppender.hh>
#include <log4cpp/Layout.hh>
#include <log4cpp/PatternLayout.hh>
#include <log4cpp/Priority.hh>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"
#include "BlasrM5AlnProvider.hpp"

void toFasta(std::string& targetId, AlnGraphBoost& g) {
    g.mergeNodes();
    std::string cns = g.consensus(8);
    if (cns.length() == 0) return;

    std::cout << ">" << targetId << std::endl;
    std::cout << cns << std::endl;
}

// Takes an alignment file (currently blasr -m 5) sorted by target, 
// generates a consensus for each target that has > 8x coverage and 
// prints it to stdout in fasta format.
void alnFileConsensus(const std::string alnFile, int minCov=8) {
    log4cpp::Category& logger = log4cpp::Category::getRoot();
    std::ifstream file(alnFile);
    int coverage = 0;

    // Initialize the graph with the first target
    Alignment aln;
    file >> aln;
    std::string targetId = aln.tid;
    AlnGraphBoost ag(aln.tlen);
    uint32_t len = aln.tlen;
    if (aln.qid != targetId) {
        logger.info("Processing target: %s, length: %d", targetId.c_str(), len);
        aln = normalizeGaps(aln);
        ag.addAln(aln);
        coverage++;
    }

    while (file >> aln) {
        if (aln.tid != targetId){
            // we're on to a new target, generate consensus ...
            if (coverage > minCov)
                toFasta(targetId, ag);

            // ... and then move on to the next target
            ag = AlnGraphBoost(aln.tlen);
            len = aln.tlen; 
            targetId = aln.tid;
            coverage = 0;
            logger.info("Processing target: %s, length: %d", targetId.c_str(), len);
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
    if (coverage > minCov) 
        toFasta(targetId, ag);
}

void alnFileConsensus2(const std::string file, size_t minCov=8) {
    std::vector<Alignment> alns;
    BlasrM5AlnProvider ap(file);
    bool hasNext = true;
    while (hasNext) {
        hasNext = ap.nextTarget(alns);
        if (alns.size() < minCov) continue;
        AlnGraphBoost ag(alns[0].tlen);
        for (auto it = alns.begin(); it != alns.end(); ++it) {
            Alignment aln = normalizeGaps(*it);
            ag.addAln(aln);
        }
    
        ag.mergeNodes();
        std::string cns = ag.consensus(minCov);
        std::cout << ">" << alns[0].tid << std::endl;
        std::cout << cns << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // Setup the root logger to a file
    log4cpp::Appender *fapp = new log4cpp::FileAppender("default", "gcon.log", false);
    log4cpp::PatternLayout *layout = new log4cpp::PatternLayout();
    layout->setConversionPattern("%d [%p] %m%n");
    fapp->setLayout(layout); 
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.setPriority(log4cpp::Priority::INFO);
    root.addAppender(fapp);

    if (argc == 2)
        alnFileConsensus2(argv[1]);
        
    return 0;
}
