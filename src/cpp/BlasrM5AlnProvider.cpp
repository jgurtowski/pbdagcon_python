#include <vector>
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include "Alignment.hpp"
#include "BlasrM5AlnProvider.hpp"

BlasrM5AlnProvider::BlasrM5AlnProvider(std::string fpath) {
    fpath_ = fpath;

    fstream_.open(fpath);
    if (! fstream_.is_open() || fstream_.fail()) {
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

    // check that the file is sorted by target 
    Alignment aln;
    std::vector<std::string> raw, sorted;
    int max = 50, count = 0;
    while(fstream_ >> aln && count++ < max) 
        raw.push_back(aln.tid);

    sorted = raw;
    std::sort(sorted.begin(), sorted.end());
    if (raw != sorted)
        throw M5Exception::SortError();

    // all is well, rewind stream in prep for real work
    if (fstream_.eof()) {
        fstream_.close();
        fstream_.open(fpath_);
    } else
        fstream_.seekg(0, fstream_.beg);
}
