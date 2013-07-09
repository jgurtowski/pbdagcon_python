#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <cassert>
#include "Alignment.hpp"

using namespace dagcon;

///
/// Simple method to reverse complement a sequence.
///
std::string revComp(std::string& seq) {
    const std::string bases = "ACTG";
    std::string::iterator curr = seq.begin();
    for (; curr != seq.end(); ++curr) {
        char& c = *curr;
        c = c == 'T' ? bases[0] : 
            c == 'G' ? bases[1] :
            c == 'A' ? bases[2] : 
            c == 'C' ? bases[3] : c;
    }
    return std::string(seq.rbegin(), seq.rend());
}

// Set this to false if the alignments are grouped by query.  The parse
// routine will be adjusted to build the alignment graph based on the 
// queries. 
bool Alignment::groupByTarget = true;

Alignment::Alignment() : 
    len(0), 
    start(0), 
    id(""), 
    qstr(""), 
    tstr("") { }

// Parses blasr m5 output grouped either by target or query.
void Alignment::parse(std::istream& instrm) {
    std::string line;
    std::getline(instrm, line);
    std::stringstream row(line);
    std::string col;
    std::vector<std::string> fields;
    while(std::getline(row, col, ' ')) {
        if (col == "") continue;
        fields.push_back(col);
    }

    // avoids *some* empty lines
    if (fields.size() == 0) return;

    // base query id (without the last '/<coordinates>'), allows us to 
    // group properly by query when asked.
    std::string baseQid = fields[0].substr(0,fields[0].find_last_of("/"));
    id = groupByTarget ? fields[5] : baseQid; 

    std::istringstream ssLen(groupByTarget ? fields[6] : fields[1]);
    ssLen >> len;
    std::istringstream ssStart(groupByTarget ? fields[7] : fields[2]);
    ssStart >> start;
    start++;

    // the target is always reversed.
    char tStrand = fields[9][0];
    if (tStrand == '-' && groupByTarget) {
        // only need to reverse complement when correcting targets
        qstr = revComp(fields[16]);
        tstr = revComp(fields[18]);
    } else {
        qstr = groupByTarget ? fields[16] : fields[18];
        tstr = groupByTarget ? fields[18] : fields[16];
    }
}

std::istream& operator>>(std::istream& instrm, Alignment& data) {
    data.parse(instrm);
    return instrm;
}

Alignment normalizeGaps(Alignment& aln) {
    size_t qlen = aln.qstr.length(), tlen = aln.tstr.length();
    assert(qlen == tlen);
    std::string qNorm, tNorm;

    // convert mismatches to indels
    for (size_t i=0; i < qlen; i++) {
        char qb = aln.qstr[i], tb = aln.tstr[i];
        if (qb != tb && qb != '-' && tb != '-') {
            qNorm += '-';
            qNorm += qb;
            tNorm += tb;
            tNorm += '-';
        } else {
            qNorm += qb;
            tNorm += tb;
        }
    }

    // update lengths
    qlen = qNorm.length();
    tlen = tNorm.length();

    // push gaps to the right, but not past the end
    for (size_t i=0; i < qlen-1; i++) {
        // pushing target gaps
        if (tNorm[i] == '-') {
            size_t j = i;
            while (true) {
                char c = tNorm[++j];
                if (c != '-' || j > qlen - 1) {
                    if (c == qNorm[i]) {
                        tNorm[i] = c;
                        tNorm[j] = '-';
                    }
                    break;
                }
            }
        }

        // pushing query gaps
        if (qNorm[i] == '-') {
            size_t j = i;
            while (true) {
                char c = qNorm[++j];
                if (c != '-' || j > tlen - 1) {
                    if (c == tNorm[i]) {
                        qNorm[i] = c;
                        qNorm[j] = '-';
                    }
                    break;
                }
            }
        }
    }

    // generate the final, normalized alignment strings
    Alignment finalNorm;
    finalNorm.id = aln.id;
    finalNorm.start = aln.start;
    finalNorm.len = aln.len;
    for (size_t i=0; i < qlen; i++) {
        if (qNorm[i] != '-' || tNorm[i] != '-') {
            finalNorm.qstr += qNorm[i];
            finalNorm.tstr += tNorm[i];
        }
    }

    return finalNorm;
}
