#include <cstring>
#include <string>
#include <algorithm>
#include <gtest/gtest.h>
#include "Types.h"
#include "PlatformId.h"
#include "Enumerations.h"
#include "DNASequence.hpp"
#include "tuples/TupleList.hpp"
#include "tuples/DNATuple.hpp"
#include "tuples/TupleMetrics.hpp"
#include "datastructures/alignment/Path.h"
#include "datastructures/alignment/AlignmentStats.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "algorithms/alignment/AlignmentUtils.hpp"
#include "FASTASequence.hpp"
#include "FASTQSequence.hpp"
#include "algorithms/alignment/DistanceMatrixScoreFunction.hpp"
#include "Alignment.hpp"
#include "algorithms/alignment/sdp/SDPFragment.hpp"
#include "algorithms/alignment/SDPAlign.hpp"
#include "SimpleAligner.hpp"

TEST(SimpleAligner, align) {
    SimpleAligner sa;
    dagcon::Alignment a;
    a.id = "test";
    a.start = 765;
    a.end = 1897;
    a.len = 2092;
    a.strand = '-';
    a.tstr = "ACAGAGATGCAAGGTAAAGTACAATTGAAAAACTAACCTCTTCCAGCGAGACTTATAGCGA";
    a.qstr = "ACAGAAGATGAAGGTAAATACAATGAAAAAACTACCTCGGTTCCAGCGAGAACTATAGCGA";
    sa.align(a);
    EXPECT_EQ("ACAGA-GATGCAAGGTAAAGTACAATTG-AAAAACTAACCTC--TTCCAGCGAG-ACTTATAGCGA", a.tstr);
    EXPECT_EQ(195, a.start);
}
