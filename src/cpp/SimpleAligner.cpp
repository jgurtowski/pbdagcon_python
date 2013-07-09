#include <stdint.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "Types.h"
#include "PlatformId.h"
#include "Enumerations.h"
#include "DNASequence.hpp"
#include "tuples/TupleMetrics.hpp"
#include "tuples/TupleList.hpp"
#include "tuples/DNATuple.hpp"
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

SimpleAligner::SimpleAligner() {
    config_.indelRate = 0.3;
    config_.indel = 5;
    config_.match = 0;
    config_.sdpIndel = 5;
    config_.sdpIns = 5;
    config_.sdpDel = 10;
    config_.kmer = 11;
    tupleMetrics_.Initialize(config_.kmer);
    distScoreFn_.del = config_.indel;
    distScoreFn_.ins = config_.indel;
    distScoreFn_.InitializeScoreMatrix(SMRTDistanceMatrix);
}

void SimpleAligner::align(dagcon::Alignment& aln) {
    // This alignment type defined in blasr code base
    blasr::Alignment blasrAln;
    DNASequence query, target;
    query.seq = (Nucleotide*) aln.qstr.c_str();
    query.length = aln.qstr.length();

    target.seq = (Nucleotide*)aln.tstr.c_str();
    target.length = aln.tstr.length();
    SDPAlign(query, target, distScoreFn_, tupleMetrics_.tupleSize,
             config_.sdpIndel, config_.sdpIndel, config_.indelRate,
             blasrAln, Global);

    std::string queryStr, alignStr, targetStr;
    CreateAlignmentStrings(blasrAln, query.seq, target.seq, 
            targetStr, alignStr, queryStr, query.length, target.length);

    if (aln.strand == '-') {
        aln.start = aln.len - (aln.end + blasrAln.tPos);
    } else {
        aln.start += blasrAln.tPos;
    }

    aln.qstr = queryStr;
    aln.tstr = targetStr;
}
