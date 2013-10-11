#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <gtest/gtest.h>
#include <boost/graph/adjacency_list.hpp>
#include <log4cpp/Appender.hh>
#include <log4cpp/OstreamAppender.hh>
#include <log4cpp/Layout.hh>
#include <log4cpp/PatternLayout.hh>
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
        ra.tlen = 9;
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

TEST(AlnGraphBoostTest, DanglingNodes) {
    log4cpp::Appender *fapp = new log4cpp::OstreamAppender("console", &std::cout);
    log4cpp::PatternLayout *layout = new log4cpp::PatternLayout();
    layout->setConversionPattern("%d %p [%c] %m%n");
    fapp->setLayout(layout); 
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.addAppender(fapp);
    AlnGraphBoost ag(12); 
    dagcon::Alignment a;
    a.tstr = "C-GCGGA-T-G-";
    a.qstr = "CCGCGG-G-A-T";

    ag.addAln(a);
    EXPECT_FALSE(ag.danglingNodes());
}
