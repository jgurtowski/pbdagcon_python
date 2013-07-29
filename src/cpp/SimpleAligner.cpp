#include <vector>
#include <stdint.h>
#include <cstring>
#include <string>
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
    tupleMetrics_.Initialize(config_.kmer);
    distScoreFn_.del = config_.indel;
    distScoreFn_.ins = config_.indel;
    distScoreFn_.InitializeScoreMatrix(SMRTDistanceMatrix);
}

void SimpleAligner::align(dagcon::Alignment& aln) {
    // This alignment type defined in blasr code base
    blasr::Alignment blasrAln;
    DNASequence query, target;
    query.seq = (Nucleotide*)aln.qstr.c_str();
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
        aln.qstr = revComp(queryStr);
        aln.tstr = revComp(targetStr);
    } else {
        aln.start += blasrAln.tPos;
        aln.qstr = queryStr;
        aln.tstr = targetStr;
    }
    aln.start++;
}

void SimpleAligner::operator() (dagcon::Alignment& aln) {
    align(aln);
}
