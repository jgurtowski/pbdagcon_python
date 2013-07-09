#include <string>
#include <map>
#include <vector>
#include <gtest/gtest.h>
#include <boost/graph/adjacency_list.hpp>
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"

TEST(AlnGraphBoostTest, RawConsensus) {
    std::string backbone = "ATATTAGGC";
    AlnGraphBoost ag(backbone); 
    dagcon::Alignment *algs = new dagcon::Alignment[5];
    
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
        dagcon::Alignment& ra = algs[i];
        ra.id = "target";
        ra.len = 9;
        ra.start = 1;
    }
    ag.addAln(algs[0]);
    ag.addAln(algs[1]);
    ag.addAln(algs[2]);
    ag.addAln(algs[3]);
    ag.addAln(algs[4]);

    ag.mergeNodes();
    std::string expected = "ATATAGCCGGC";
    const std::string actual = ag.consensus();
    EXPECT_EQ(expected, actual);
}
