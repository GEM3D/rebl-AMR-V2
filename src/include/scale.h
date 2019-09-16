#ifndef _SCALE_H_
#define _SCALE_H_
#include "definitions.h"
#include "typedefs.h"

template <size_t N>
class FullOctreeTop
{

    private:
    uint fullOctreeLevel;
    morton<N> rootKey;
    vector<int> Nbr;  

    public:
    FullOctreeTop();
    void convertRank2Bits( int myrank );
    void constructNbrProcs();
    bool isBoundary( uint &direction );
    void checkGraphConsistency(int myrank);  
    void readRoot(morton<N> &key); 
    void readNbrs(vector<int>&Nbrs );
    ~FullOctreeTop() {};

};

#endif

