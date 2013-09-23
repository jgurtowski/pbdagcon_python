#include <vector>
#include <stdint.h>
#include <cstring>
#include <string>
#include <cassert>
#include <algorithm>
#include "Alignment.hpp"
#include "SimpleAligner.hpp"

SimpleAligner::SimpleAligner() {
    config_.indelRate = 0.3;
    config_.indel = 5;
    config_.match = 0;
    config_.sdpIndel = 5;
    config_.sdpIns = 5;
    config_.sdpDel = 10;
    config_.kmer = 11;
    config_.bandSize = 10;
    tupleMetrics_.Initialize(config_.kmer);
    distScoreFn_.del = config_.indel;
    distScoreFn_.ins = 4;
    distScoreFn_.InitializeScoreMatrix(SMRTDistanceMatrix);
}

void SimpleAligner::align(dagcon::Alignment& aln) {
    // This alignment type defined in blasr code base
    blasr::Alignment initialAln, refinedAln;
    FASTQSequence query;
    query.seq = (Nucleotide*)aln.qstr.c_str();
    query.length = aln.qstr.length();

    DNASequence target;
    target.seq = (Nucleotide*)aln.tstr.c_str();
    target.length = aln.tstr.length();
    SDPAlign(query, target, distScoreFn_, tupleMetrics_.tupleSize,
             config_.sdpIndel, config_.sdpIndel, config_.indelRate*2,
             initialAln, Local);

    GuidedAlign(query, target, initialAln, distScoreFn_, 
        config_.bandSize, refinedAln);

    std::string queryStr, alignStr, targetStr;

    //StickPrintAlignment(initialAln, query, target, std::cout);
    //StickPrintAlignment(refinedAln, query, target, std::cout);

    CreateAlignmentStrings(refinedAln, query.seq, target.seq, 
            targetStr, alignStr, queryStr, query.length, target.length);

    // alignment coordinates may change, update alignment object
    aln.start += refinedAln.GenomicTBegin();
    aln.end = aln.start + refinedAln.GenomicTEnd();

    if (aln.strand == '-') {
        aln.start = aln.tlen - aln.end;
        aln.qstr = revComp(queryStr);
        aln.tstr = revComp(targetStr);
    } else {
        aln.qstr = queryStr;
        aln.tstr = targetStr;
    }
}

void SimpleAligner::operator() (dagcon::Alignment& aln) {
    align(aln);
}
