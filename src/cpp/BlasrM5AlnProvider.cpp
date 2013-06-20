#include <vector>
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include "Alignment.hpp"
#include "BlasrM5AlnProvider.hpp"


BlasrM5AlnProvider::BlasrM5AlnProvider(std::string fpath) {
    // XXX: initialize logger 
    fpath_ = fpath;

    fstream_.open(fpath);
    if (! fstream_ || fstream_.fail()) {
        throw M5Exception::FileOpenError();
    }

    checkFormat();
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
    // parse the first line and run some field checks
    std::string line;
    std::getline(fstream_, line);
    std::stringstream row(line);
    std::string col;
    std::vector<std::string> fields;

    while(std::getline(row, col, ' ')) {
        if (col == "") continue;
        fields.push_back(col);
    }

    if (fields.size() < 19) {
        boost::format msg("Expected 19 fields, found %d");
        msg % fields.size();
        throw M5Exception::FormatError(msg.str());
    }

    fstream_.seekg(0, fstream_.beg);
}
