#ifndef __GCON_ALN_PROVIDER__
#define __GCON_ALN_PROVIDER__

///
/// Provides sets of alignments for a given target sequence from a blasr M5 
/// file.
///
class BlasrM5AlnProvider {
public:
    /// Constructs a new alignment provider.  Checks the format of the file and
    /// throws an exception if it's malformed.
    /// \param fpath Path to the file containing alignments.
    BlasrM5AlnProvider(const std::string fpath);

    /// Gets the set of alignments for the next target and puts them into the
    /// given vector.  Note this function will clear the contents of the vector
    /// prior to adding the next set of alignments.
    /// \param dest reference to a vector to hold the alignments.
    /// \return True if there are more targets, otherwise false.
    bool nextTarget(std::vector<Alignment>& dest);
    
    /// Called during constructor, checks that the file is formatted correctly.
    void checkFormat();

private:
    /// Represents an input stream to the alignment file.
    std::ifstream fstream_;

    /// Path to the input file
    std::string fpath_;

    /// State variables 
    std::string currTargetId_;
    Alignment prevAln_;
    bool firstAln_;
};

#endif //__GCON_ALN_PROVIDER__
