#include <vector>
#include <fstream>
#include "Alignment.hpp"
#include "BlasrM5AlnProvider.hpp"

BlasrM5AlnProvider::BlasrM5AlnProvider(std::string fpath) {
    // XXX: initialize logger 
    fpath_ = fpath;

    // XXX: check format
    
    // XXX: check failbit
    fstream_.open(fpath);
    currTargetId_ = "";
    firstAln_ = true;
}

bool BlasrM5AlnProvider::nextTarget(std::vector<Alignment>& dest) {
    // first clear any previous alignments
    dest.clear();

    // process up to EOF or next target
    // need to maintain state in between calls
    if (! firstAln_)
        dest.push_back(prevAln_); 

    Alignment aln;
    while (fstream_ >> aln) {
        if (aln.tid != currTargetId_) {
            firstAln_ = false;
            prevAln_ = aln;
            currTargetId_ = aln.tid;
            break;
        }

        // skip self hits
        if (aln.qid == currTargetId_) continue;
        dest.push_back(aln);
    }

    return fstream_.good() ? true : false;
}

void BlasrM5AlnProvider::checkFormat() {
    // XXX: implement
}
