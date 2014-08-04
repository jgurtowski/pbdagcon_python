#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <log4cpp/Category.hh>
#include <boost/format.hpp>
#include "Alignment.hpp"
#include "BlasrM5AlnProvider.hpp"


BlasrM5AlnProvider::BlasrM5AlnProvider(const std::string& fpath) :
    fpath_(fpath),
    currId_(""),
    firstAln_(true),
    fs_() {

    //checkFormat();
    fs_.open(fpath_);
    is_ = &fs_;
}

BlasrM5AlnProvider::BlasrM5AlnProvider(std::istream* stream) :
    fpath_(""),
    currId_(""),
    firstAln_(true),
    fs_(),
    is_(stream) {
}

BlasrM5AlnProvider::~BlasrM5AlnProvider() {
    delete is_;
}

bool BlasrM5AlnProvider::nextTarget(std::vector<dagcon::Alignment>& dest) {
    // first clear any previous alignments
    dest.clear();

    // process up to EOF or next target
    // need to maintain state in between calls
    if (! firstAln_)
        dest.push_back(prevAln_); 

    dagcon::Alignment aln;
    while (*is_ >> aln) {
        if (aln.id != currId_) {
            firstAln_ = false;
            prevAln_ = aln;
            currId_ = aln.id;
            break;
        }
        dest.push_back(aln);
    }

    return (*is_);
}

void BlasrM5AlnProvider::checkFormat() {
    log4cpp::Category& logger = 
        log4cpp::Category::getInstance("BlasrM5AlnProvider");
    std::ifstream ifs(fpath_);
    if (! ifs.is_open() || ifs.fail()) {
        throw M5Exception::FileOpenError();
    }
    // parse the first line and run some field checks
    std::string line;
    std::getline(ifs, line);
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

    // check how the alignments are grouped
    dagcon::Alignment aln;
    std::vector<std::string> raw, sorted;
    int max = 50, count = 0;
    while(ifs >> aln && count++ < max) 
        raw.push_back(aln.id);

    sorted = raw;
    std::sort(sorted.begin(), sorted.end());

    std::string logl = "dagcon::Alignments appear to be grouped by %s";
    if (raw != sorted) {
        logger.info(logl.c_str(), "query");
        dagcon::Alignment::groupByTarget = false;
    } else {
        logger.info(logl.c_str(), "target");
    }
    ifs.close();
}
