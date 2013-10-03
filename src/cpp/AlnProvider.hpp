#ifndef __GCON_ALN_PROVIDER__
#define __GCON_ALN_PROVIDER__

///
/// Generic alignment provider interface.
///
class AlnProvider {
public:
    /// Gets the set of alignments for the next target and puts them into the
    /// given vector.  Note this function will clear the contents of the vector
    /// prior to adding the next set of alignments.
    /// \param dest reference to a vector to hold the alignments.
    /// \return True if there are more targets, otherwise false.
    virtual bool nextTarget(std::vector<dagcon::Alignment>& dest) = 0;

    virtual ~AlnProvider() {};
};

#endif //__GCON_ALN_PROVIDER__
