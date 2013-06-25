#include <cstdint>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <log4cpp/Category.hh>
#include <log4cpp/Appender.hh>
#include <log4cpp/FileAppender.hh>
#include <log4cpp/Layout.hh>
#include <log4cpp/PatternLayout.hh>
#include <log4cpp/Priority.hh>
#include <boost/format.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/thread.hpp>
#include <boost/call_traits.hpp>
#include <boost/progress.hpp>
#include <boost/bind.hpp>
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"
#include "BlasrM5AlnProvider.hpp"
#include "BoundedBuffer.hpp"

///
/// Single-threaded consensus execution.
///
void alnFileConsensus(const std::string file, size_t minCov=8) {
    log4cpp::Category& logger = log4cpp::Category::getInstance("consensus");
    std::vector<Alignment> alns;
    try {
        BlasrM5AlnProvider* ap;
        if (file == "-") { 
            ap = new BlasrM5AlnProvider(&std::cin);
        } else {
            ap = new BlasrM5AlnProvider(file);
        }

        bool hasNext = true;
        while (hasNext) {
            hasNext = ap->nextTarget(alns);
            if (alns.size() < minCov) continue;
            AlnGraphBoost ag(alns[0].len);
            for (auto it = alns.begin(); it != alns.end(); ++it) {
                Alignment aln = normalizeGaps(*it);
                ag.addAln(aln);
            }
        
            ag.mergeNodes();
            std::string cns = ag.consensus(minCov);
            if (cns == "") continue;
            std::cout << ">" << alns[0].id << std::endl;
            std::cout << cns << std::endl;
        }
    } 
    catch (M5Exception::FileOpenError) {
        logger.error("Error opening file: %s", file.c_str());
        return;
    }
    catch (M5Exception::FormatError err) {
        logger.error("Format error. Input: %s, Error: %s", 
            file.c_str(), err.msg.c_str());
    }
    catch (M5Exception::SortError err) {
        logger.error("Input file is not sorted by either target or query.");
    }
}

typedef std::vector<Alignment> AlnVec;
typedef BoundedBuffer<AlnVec> AlnBuf;
typedef BoundedBuffer<std::string> CnsBuf;

class Reader {
    AlnBuf* alnBuf_;
    const std::string fpath_;
    size_t minCov_;
    int nCnsThreads_;
public:
    Reader(AlnBuf* b, const std::string fpath) : alnBuf_(b), fpath_(fpath) {
        minCov_ = 8;
    }

    void setNumCnsThreads(int n) {
        nCnsThreads_ = n;
    }

    void operator()() {
        log4cpp::Category& logger = 
            log4cpp::Category::getInstance("Reader");
        try {
            BlasrM5AlnProvider* ap;
            if (fpath_ == "-") { 
                ap = new BlasrM5AlnProvider(&std::cin);
            } else {
                ap = new BlasrM5AlnProvider(fpath_);
            }
            AlnVec alns;
            bool hasNext = true;
            while (hasNext) {
                hasNext = ap->nextTarget(alns);
                if (alns.size() < minCov_) continue;
                alnBuf_->push(alns);
            }
        } 
        catch (M5Exception::FileOpenError) {
            logger.error("Error opening file: %s", fpath_.c_str());
        }
        catch (M5Exception::FormatError err) {
            logger.error("Format error. Input: %s, Error: %s", 
                fpath_.c_str(), err.msg.c_str());
        }
        catch (M5Exception::SortError err) {
            logger.error("Input file is not sorted by either target or query.");
        }

        // write out sentinals, one per consensus thread
        AlnVec sentinel;
        for (int i=0; i < nCnsThreads_; i++)
            alnBuf_->push(sentinel);
    }
};

class Consensus {
    AlnBuf* alnBuf_;
    CnsBuf* cnsBuf_;
    int minWeight_;
public:
    Consensus(AlnBuf* ab, CnsBuf* cb) : alnBuf_(ab), cnsBuf_(cb) {
        minWeight_ = 8;
    }

    void operator()() {
        AlnVec alns;
        alnBuf_->pop(&alns);

        while (alns.size() > 0) {
            AlnGraphBoost ag(alns[0].len);
            for (auto it = alns.begin(); it != alns.end(); ++it) {
                Alignment aln = normalizeGaps(*it);
                ag.addAln(aln);
            }
            ag.mergeNodes();
            std::ostringstream fasta;
            std::string cns = ag.consensus(minWeight_);
            if (cns != "") {
                fasta << ">" << alns[0].id << std::endl;
                fasta << cns << std::endl;
                cnsBuf_->push(fasta.str()); 
            }

            alnBuf_->pop(&alns);
        }
        // write out a sentinal
        cnsBuf_->push("");
    }
};

class Writer {
    CnsBuf* cnsBuf_;
    int nCnsThreads_;
public:
    Writer(CnsBuf* cb) : cnsBuf_(cb) {}
    
    void setNumCnsThreads(int n) {
        nCnsThreads_ = n;
    }

    void operator()() {
        std::string cns;
        cnsBuf_->pop(&cns);
        int sentinelCount = 0;
        while (true) {
            if (cns != "") 
                std::cout << cns;
            if (cns == "" && ++sentinelCount == nCnsThreads_) 
                break;

            cnsBuf_->pop(&cns);
        }
    }
};

void setupLogger() {
    // Setup the root logger to a file
    log4cpp::Appender *fapp = new log4cpp::FileAppender("default", "gcon.log", false);
    log4cpp::PatternLayout *layout = new log4cpp::PatternLayout();
    layout->setConversionPattern("%d %p [%c] %m%n");
    fapp->setLayout(layout); 
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.setPriority(log4cpp::Priority::INFO);
    root.addAppender(fapp);
}

int main(int argc, char* argv[]) {
    setupLogger();
    log4cpp::Category& logger = log4cpp::Category::getInstance("main");

    // XXX: Add command line options
    if (argc == 2) {
        logger.info("Single-threaded. Input: %s", argv[1]);
        alnFileConsensus(argv[1]);
    } else if (argc == 3) {
        logger.info("Multi-threaded. Input: %s, Threads: %s", 
            argv[1], argv[2]);
    
        AlnBuf alnBuf(30);
        CnsBuf cnsBuf(30);

        int nthreads;
        std::istringstream threadArg(argv[2]);
        threadArg >> nthreads;

        if (nthreads == 0) {
            boost::format msg("Invalid argument for num threads: %s");
            msg % argv[2];
            std::cerr << msg << std::endl;
            exit(-1);
        }
        
        Writer writer(&cnsBuf);
        writer.setNumCnsThreads(nthreads);
        boost::thread threadedWriter(writer);

        std::vector<boost::thread> threadedConsensus;
        for (int i=0; i < nthreads; i++) {
            Consensus c(&alnBuf, &cnsBuf);
            threadedConsensus.push_back(boost::thread(c));
        }

        Reader reader(&alnBuf, argv[1]);
        reader.setNumCnsThreads(nthreads);
        boost::thread threadedReader(reader);

        threadedWriter.join();
        std::vector<boost::thread>::iterator it;
        for (it = threadedConsensus.begin(); it != threadedConsensus.end(); ++it)
            it->join();
    
        threadedReader.join();
    }
        
    return 0;
}
