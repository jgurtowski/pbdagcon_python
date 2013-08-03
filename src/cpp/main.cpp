#include <cstdint>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
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
#include <boost/program_options.hpp>
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"
#include "BlasrM5AlnProvider.hpp"
#include "BoundedBuffer.hpp"
#include "tuples/TupleMetrics.hpp"
#include "SimpleAligner.hpp"

namespace opts = boost::program_options;

struct FilterOpts {
    /// Minimum alignment coverage for consensus
    size_t minCov;
    /// Minimum consensus length to output
    size_t minLen;
} fopts;

bool AlignFirst = false;

///
/// Single-threaded consensus execution.
///
void alnFileConsensus(AlnProvider* ap, const FilterOpts& fopts) {
    log4cpp::Category& logger = 
        log4cpp::Category::getInstance("afc");
    std::vector<dagcon::Alignment> alns;
    std::vector<std::string> seqs;
    SimpleAligner aligner;
    bool hasNext = true;
    while (hasNext) {
        hasNext = ap->nextTarget(alns);
        if (alns.size() < fopts.minCov) continue;

        if (AlignFirst) 
            for_each(alns.begin(), alns.end(), aligner); 

        AlnGraphBoost ag(alns[0].len);
        for (auto it = alns.begin(); it != alns.end(); ++it) {
            dagcon::Alignment aln = normalizeGaps(*it);
            trimAln(aln);
            boost::format msg("%s %s %s %d %d");
            msg % aln.sid;
            msg % aln.strand;
            msg % aln.qstr.substr(aln.qstr.length()-50, 50);
            msg % aln.start;
            msg % aln.qstr.length();
            logger.debugStream() << msg.str();
            ag.addAln(aln);
        }
    
        ag.danglingNodes();
        ag.mergeNodes();
        ag.consensus(seqs, fopts.minCov);
        int i = 0;
        for (auto it = seqs.begin(); it != seqs.end(); ++it) {
            std::string cns = *it;
            if (cns.length() < fopts.minLen) continue; 
            std::cout << ">" << alns[0].id << "/" << i++ << std::endl;
            std::cout << cns << std::endl;
        }
    }
}

typedef std::vector<dagcon::Alignment> AlnVec;
typedef BoundedBuffer<AlnVec> AlnBuf;
typedef BoundedBuffer<std::string> CnsBuf;

class Reader {
    AlnBuf* alnBuf_;
    const std::string fpath_;
    size_t minCov_;
    int nCnsThreads_;
public:
    Reader(AlnBuf* b, const std::string fpath, size_t minCov) : 
        alnBuf_(b), 
        fpath_(fpath),
        minCov_(minCov)
    { }

    void setNumCnsThreads(int n) {
        nCnsThreads_ = n;
    }

    void operator()() {
        log4cpp::Category& logger = 
            log4cpp::Category::getInstance("Reader");
        try {
            AlnProvider* ap;
            if (fpath_ == "-") { 
                ap = new BlasrM5AlnProvider(&std::cin);
            } else {
                ap = new BlasrM5AlnProvider(fpath_);
            }
            AlnVec alns;
            bool hasNext = true;
            while (hasNext) {
                hasNext = ap->nextTarget(alns);
                size_t cov = alns.size();
                if (cov > 0 && cov < minCov_) {
                    logger.debug("Coverage requirement not met for %s, coverage: %d", 
                        alns[0].id.c_str(), alns.size());
                    continue;
                }
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
    size_t minLen_;
    int minWeight_;
    SimpleAligner aligner;
public:
    Consensus(AlnBuf* ab, CnsBuf* cb, size_t minLen) : 
        alnBuf_(ab), 
        cnsBuf_(cb),
        minLen_(minLen),
        minWeight_(8)
    { }

    void operator()() {
        AlnVec alns;
        alnBuf_->pop(&alns);
        std::vector<std::string> seqs;

        while (alns.size() > 0) {
            if (alns.size() < fopts.minCov) {
                alnBuf_->pop(&alns);
                continue;
            }

            if (AlignFirst) 
                for_each(alns.begin(), alns.end(), aligner); 

            AlnGraphBoost ag(alns[0].len);
            for (auto it = alns.begin(); it != alns.end(); ++it) {
                dagcon::Alignment aln = normalizeGaps(*it);
                ag.addAln(aln);
            }
            ag.mergeNodes();
            std::ostringstream fasta;
            ag.consensus(seqs, minWeight_, minLen_);
            int i = 0;
            for (auto it = seqs.begin(); it != seqs.end(); ++it) {
                std::string cns = *it;
                fasta << ">" << alns[0].id << "/" << i++ << std::endl;
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
            std::cout << cns;
            if (cns == "" && ++sentinelCount == nCnsThreads_) 
                break;

            cnsBuf_->pop(&cns);
        }
    }
};

void setupLogger(log4cpp::Priority::Value priority) {
    // Setup the root logger to a file
    log4cpp::Appender *fapp = new log4cpp::FileAppender("default", "pbdagcon.log", false);
    log4cpp::PatternLayout *layout = new log4cpp::PatternLayout();
    layout->setConversionPattern("%d %p [%c] %m%n");
    fapp->setLayout(layout); 
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.setPriority(priority);
    root.addAppender(fapp);
}

bool parseOpts(int ac, char* av[], opts::variables_map& vm) {
    opts::options_description 
        odesc("PacBio read-on-read error correction via consensus");
    odesc.add_options()
        ("help,h", "Display this help")
        ("verbose,v", "Increase logging verbosity")
        ("align,a", "Align sequences before adding to consensus")
        ("min-length,m", opts::value<int>()->default_value(500), 
            "Filter corrected reads less than length")
        ("min-coverage,c", opts::value<int>()->default_value(8),
            "Minimum coverage required to correct")
        ("correct-query,q", "Correct the queries instead of targets. "
            "If not set, and the input is a file, the code will correct the side "
            "that is grouped.")
        ("threads,j", opts::value<int>(), "Number of consensus threads to use")
        ("rbuf,r", opts::value<int>()->default_value(30), "Size of the read buffer")
        ("wbuf,w", opts::value<int>()->default_value(30), "Size of the write buffer")
        ("input", opts::value<std::string>()->default_value("-"), "Input (flag is optional)")
    ;

    opts::positional_options_description pdesc; 
    pdesc.add("input", 1);
    opts::store(opts::command_line_parser(ac, av)
                .options(odesc).positional(pdesc).run(), vm);

    opts::notify(vm);

    if (vm.count("help") || ! vm.count("input")) {
        std::cout << "Usage: " << av[0] << " [options] <input>"<< "\n\n";
        std::cout << odesc << "\n";
        return false;
    }


    return true;
}

int main(int argc, char* argv[]) {
    opts::variables_map vm;
    if (! parseOpts(argc, argv, vm)) exit(1);

    // http://log4cpp.sourceforge.net/api/classlog4cpp_1_1Priority.html
    // defaults to INFO.
    setupLogger(vm.count("verbose") ? 700 : 600);
    log4cpp::Category& logger = log4cpp::Category::getInstance("main");

    if (vm.count("align")) {
        dagcon::Alignment::parse = parsePre;
        AlignFirst = true;
    }

    if (vm.count("correct-query")) {
        dagcon::Alignment::groupByTarget = false;
        logger.info("Configured to correct queries");
    } else {
        logger.info("Configured to correct targets");
    }

    fopts.minCov = vm["min-coverage"].as<int>();
    fopts.minLen = vm["min-length"].as<int>();

    std::string input = vm["input"].as<std::string>(); 
    if (vm.count("threads")) {
        int nthreads = vm["threads"].as<int>();
        logger.info("Multi-threaded. Input: %s, Threads: %d", 
            input.c_str(), nthreads);
    
        AlnBuf alnBuf(vm["rbuf"].as<int>());
        CnsBuf cnsBuf(vm["wbuf"].as<int>());

        Writer writer(&cnsBuf);
        writer.setNumCnsThreads(nthreads);
        boost::thread writerThread(writer);

        std::vector<boost::thread> cnsThreads;
        for (int i=0; i < nthreads; i++) {
            Consensus c(&alnBuf, &cnsBuf, fopts.minLen);
            cnsThreads.push_back(boost::thread(c));
        }

        Reader reader(&alnBuf, input, fopts.minCov);
        reader.setNumCnsThreads(nthreads);
        boost::thread readerThread(reader);

        writerThread.join();
        std::vector<boost::thread>::iterator it;
        for (it = cnsThreads.begin(); it != cnsThreads.end(); ++it)
            it->join();
    
        readerThread.join();
    } else {
        AlnProvider* ap;
        try {
            if (input == "-") { 
                ap = new BlasrM5AlnProvider(&std::cin);
            } else {
                ap = new BlasrM5AlnProvider(input);
            }
            logger.info("Single-threaded. Input: %s", input.c_str());
            alnFileConsensus(ap, fopts);
        } 
        catch (M5Exception::FileOpenError) {
            logger.error("Error opening file: %s", input.c_str());
            exit(1);
        }
        catch (M5Exception::FormatError err) {
            logger.error("Format error. Input: %s, Error: %s", 
                input.c_str(), err.msg.c_str());
        }
        catch (M5Exception::SortError err) {
            logger.error("Input file is not sorted by either target or query.");
        }
    }
        
    return 0;
}
