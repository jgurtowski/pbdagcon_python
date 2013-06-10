#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <cassert>
#include "Alignment.hpp"

std::string revComp(std::string& seq) {
    const std::string bases = "ACTG";
    for (char& c : seq) {
        c = c == 'T' ? bases[0] : 
            c == 'G' ? bases[1] :
            c == 'A' ? bases[2] : 
            c == 'C' ? bases[3] : c;
    }
    return std::string(seq.rbegin(), seq.rend());
}

// parses blasr m5 output
void Alignment::parse(std::istream& instrm) {
    std::string line;
    std::getline(instrm, line);
    std::stringstream row(line);
    std::string col;
    int nCol = 0;
    char tStrand = '+';
    while(std::getline(row, col, ' ')) {
        if (col == "") continue;
        nCol++;
        if (nCol == 6) {
            tid = col;
        } else if (nCol == 7) {
            std::istringstream iss(col);
            iss >> tlen;
        } else if (nCol == 8) {
            std::istringstream iss(col);
            iss >> tstart;
            tstart++;
        } else if (nCol == 10) {
            tStrand = col[0];
        } else if (nCol == 17) {
            qstr = tStrand == '-' ? revComp(col) : col;
        } else if (nCol == 19) {
            tstr = tStrand == '-' ? revComp(col) : col;
        }
    }
}

std::istream& operator>>(std::istream& instrm, Alignment& data) {
    data.parse(instrm);
    return instrm;
}

Alignment normalizeGaps(Alignment& aln) {
    size_t qlen = aln.qstr.length(), tlen = aln.qstr.length();
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
    finalNorm.tstart = aln.tstart;
    for (size_t i=0; i < qlen; i++) {
        if (qNorm[i] != '-' || tNorm[i] != '-') {
            finalNorm.qstr += qNorm[i];
            finalNorm.tstr += tNorm[i];
        }
    }

    return finalNorm;
}
