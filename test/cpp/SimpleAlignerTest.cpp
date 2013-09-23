#include <cstring>
#include <string>
#include <algorithm>
#include <gtest/gtest.h>
#include "Alignment.hpp"
#include "SimpleAligner.hpp"

TEST(SimpleAligner, align) {
    SimpleAligner sa;
    dagcon::Alignment a;
    a.id = "test";
    a.start = 765;
    a.end = 1897;
    a.tlen = 2092;
    a.strand = '-';
    a.tstr = "ACAGAGATGCAAGGTAAAGTACAATTGAAAAACTAACCTCTTCCAGCGAGACTTATAGCGA";
    a.qstr = "ACAGAAGATGAAGGTAAATACAATGAAAAAACTACCTCGGTTCCAGCGAGAACTATAGCGA";
    sa.align(a);
    EXPECT_EQ("TCGCTATAAGT-CTCGCTGGAA--GAGGTTAGTTTTT-CAATTGTACTTTACCTTGCATCT-CTGT", a.tstr);
}
