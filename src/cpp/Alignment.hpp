#ifndef __GCON_ALIGNMENT_HPP__
#define __GCON_ALIGNMENT_HPP__

// Super-simple alignment representation.
// XXX: Accomodate CIGAR strings? Still need target sequence for gap pushing.
struct Alignment {
    uint32_t tlen;
    
    // comforming offsets are 1-based
    uint32_t tstart;

    // query id
    std::string qid;

    // target id
    std::string tid;

    // query and target strings must be equal length
    std::string qstr;
    std::string tstr;

    // XXX: currently blasr m5 output parser.  make more flexible?
    void parse(std::istream& instrm);
};

std::istream& operator>>(std::istream& instrm, Alignment& data);

/// Simplifies the alignment by normalizing gaps.  Converts mismatches into
/// indels ... 
///      query: CAC        query:  C-AC
///             | |  --->          |  |
///     target: CGC       target:  CG-C
///
/// Shifts equivalent gaps to the right in the reference ...
///      query: CAACAT        query: CAACAT
///             | | ||  --->         |||  |
///     target: C-A-AT       target: CAA--T
///
/// Shifts equivalent gaps to the right in the read ...
///      query: -C--CGT       query: CCG--T
///              |  | |  --->        |||  |
///     target: CCGAC-T      target: CCGACT  
Alignment normalizeGaps(Alignment& aln);

#endif // __GCON_ALIGNMENT_HPP__
