#include <stdint.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "Types.h"
#include "PlatformId.h"
#include "Enumerations.h"
#include "algorithms/alignment/ScoreMatrices.h"
#include "DNASequence.hpp"
#include "tuples/TupleMetrics.hpp"
#include "tuples/TupleList.hpp"
#include "tuples/BaseTuple.hpp"
#include "tuples/DNATuple.hpp"
#include "tuples/TupleMetrics.hpp"
#include "datastructures/alignment/Path.h"
#include "datastructures/alignment/AlignmentStats.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "algorithms/alignment/AlignmentUtils.hpp"
#include "qvs/QualityValue.hpp"
#include "qvs/QualityValueVector.hpp"
#include "datastructures/reads/ZMWGroupEntry.hpp"
#include "FASTASequence.hpp"
#include "FASTQSequence.hpp"
#include "algorithms/alignment/BaseScoreFunction.hpp"
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
    int alignScore;
    Nucleotide* qs = new Nucleotide[aln.qstr.length()];
    Nucleotide* ts = new Nucleotide[aln.qstr.length()];
    std::strcpy((char*)qs, aln.qstr.c_str());
    std::strcpy((char*)ts, aln.tstr.c_str());
    DNASequence query, target;
    query.seq = qs;
    query.length = aln.qstr.length();

    target.seq = ts;
    target.length = aln.tstr.length();
    alignScore = SDPAlign(query, target,
                          distScoreFn_, tupleMetrics_.tupleSize,
                          config_.sdpIndel, config_.sdpIndel, config_.indelRate,
                          blasrAln,
                          Global);

    std::string queryStr, alignStr, targetStr;
    CreateAlignmentStrings(blasrAln, query.seq, target.seq, 
            targetStr, alignStr, queryStr, query.length, target.length);

    aln.qstr = queryStr;
    aln.tstr = targetStr;
}
