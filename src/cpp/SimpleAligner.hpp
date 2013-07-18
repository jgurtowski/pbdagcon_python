#ifndef __GCON_SIMPLE_ALIGNER__
#define __GCON_SIMPLE_ALIGNER__
#include "Types.h"
#include "PlatformId.h"
#include "Enumerations.h"
#include "DNASequence.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "algorithms/alignment/AlignmentUtils.hpp"
#include "algorithms/alignment/SDPAlign.hpp"

namespace Aligner {
struct Config {
    float indelRate;
    int indel;
    int match;
    int sdpIndel;
    int sdpIns;
    int sdpDel;
    int kmer;
};
}

class SimpleAligner {
public:
    SimpleAligner();
    void align(dagcon::Alignment& aln); 
    void operator() (dagcon::Alignment& aln);
private:
    Aligner::Config config_;
    TupleMetrics tupleMetrics_;
    DistanceMatrixScoreFunction<DNASequence, DNASequence> distScoreFn_;
};

#endif // __GCON_SIMPLE_ALIGNER__
