#ifndef __GCON_BLASRM5_ALN_PROVIDER__
#define __GCON_BLASRM5_ALN_PROVIDER__

#include "AlnProvider.hpp"

///
/// Exceptions thrown by this class
///
namespace M5Exception {
    struct FormatError {
        std::string msg;
        FormatError(std::string m) { msg = m; }
    };

    struct SortError {};

    struct FileOpenError {};
}

///
/// Provides sets of alignments for a given target sequence from a blasr M5 
/// file.  File may be grouped by target or query.  The grouping determines
/// which set gets corrected.  Earlier, pre-assembly reads were corrected as
/// targets.  However, we can avoid the sort step if we can correct the reads 
/// as queries, since blasr groups alignments by query.
///
class BlasrM5AlnProvider : public AlnProvider {
public:
    /// Constructs a new alignment provider.  Checks the format of the file and
    /// throws an exception if it's malformed.
    /// \param fpath Path to the file containing alignments.
    BlasrM5AlnProvider(const std::string& fpath);

    /// Constructs a provider based on the given stream.  Note that no checks
    /// are actually made on the validity of the format, caveat emptor. This
    /// can be used to take a piped stream of alignments straight from blasr.
    BlasrM5AlnProvider(std::istream* stream);
    
    /// Cleans up some stuff.
    ~BlasrM5AlnProvider();

    /// Gets the set of alignments for the next target and puts them into the
    /// given vector.  Note this function will clear the contents of the vector
    /// prior to adding the next set of alignments.
    /// \param dest reference to a vector to hold the alignments.
    /// \return True if there are more targets, otherwise false.
    bool nextTarget(std::vector<dagcon::Alignment>& dest);
    
    /// Called during constructor, checks that the file is formatted correctly.
    /// Also determines if the input is grouped by query or target.
    void checkFormat();

private:
    /// Path to the input file
    const std::string fpath_;

    /// State variables 
    std::string currId_;
    dagcon::Alignment prevAln_;
    bool firstAln_;
    
    /// Represents an input stream to the alignments.
    std::ifstream fs_;
    std::istream* is_;
};

#endif //__GCON_BLASRM5_ALN_PROVIDER__
