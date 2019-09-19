#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <vector>
#include <list>
#include <unistd.h>
#include <mpi.h>
#include <memory>
#include <time.h>
#include <stdexcept>
#include "zoltan.h"
#include <cstddef>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <iomanip> 	
#include <cmath>
#include <assert.h>
#include "params.h"
#include <locale>

#define hash 0
#define nonnative 1 /*!\ default STD hash and comparison functions  */


using real = double;
using uint = unsigned int;
using integer = int;
using namespace std;

// Major data structure used in this research

template <size_t N>
using morton = std::bitset<N>; /*!\var typedef template of bitset class used for Morton encoding */

#if(nonnaitive)
#if ( hash )
template <size_t N>
class morton_hash
{
    public:
    size_t operator()( const morton<N> key ) const // <-- don't forget const, means this function is not allowed to modify the object
    {
        morton<N> kt;
        /*
               for(uint i=0;i<N;i++)
               {
                kt[i]=key[N-i-1];
                }
        */
        //  kt=key;
        size_t hashval = 0;
        hashval = key.to_ulong();
        return hashval;
    }
};

template <size_t N, typename value>
using bitmap = std::unordered_map<morton<N>, value *, morton_hash<N>>; /*! \var typedef template unordered_map used in amrMesh Class*/


#else

template <size_t N>
class compare
{
    public:
    bool operator()( const morton<N> &a, const morton<N> &b ) const
    {
        return a.to_ullong() < b.to_ullong();
    }
};

template <size_t N, typename value>
using bitmap = std::map<morton<N>, value *, compare<N>>; /*! \var typedef template map (red-black tree) used in amrMesh Class*/

#endif
// using bitmapgeom=std::unordered_map<morton , uint >; /*! \var  typedef unordered_map used in amrMesh Class*/

#else
template <size_t N, typename value>
using bitmap = std::unordered_map<morton<N>, value *>; /*! \var typedef template unordered_map used in amrMesh Class*/


#endif

template <size_t M>
using bitlist = std::list<morton<M>>; /*! \var typedef template list is used in amrMesh Class*/

template <size_t N>
using bitvector = std::vector<morton<N>>;

template <size_t N>
using bitunorderedset = std::unordered_set<morton<N>>;

/*
 * \brief function to put commas in the output, e.g.  1000 >> 1,000
 *  
 * */
template<class T>
std::string FormatWithCommas(T value)
{
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
}


typedef struct sv
{
real f; /*!< source terms for the poisson equation*/
real p; /*!< cell-centered pressure values */  

sv& operator=(const sv& a)
{
    f = a.f;
    p = a.p;

    return *this;
}

sv& operator*(const double &number)
{
		f = f*number;
		p = p*number;
	return *this;		
}

sv& operator+(const sv& a)
{
		f = f+a.f;
		p = p+a.p;

	return *this;		
}




} Q;
 
template<size_t N>
bool compare(morton<N> a, morton<N> b) 
{ 
    return (a.to_ulong() < b.to_ulong()); 
} 



template<size_t N>
void getSortedIndex(vector<morton<N>> &nbrs,int *index )
{
 vector<morton<N>> nbrs_old;

    for ( int i = 0; i < nbrs.size(); i++ )
    {
        nbrs_old.push_back( nbrs.at( i ) );
    }

    std::sort( nbrs.begin(), nbrs.end(), compare<N> );

    for ( int i = 0; i < nbrs.size(); i++ )
    {
        index[i] = 0;
        for ( int j = 0; j < 4; j++ )
        {
            if ( nbrs_old.at( j ) == nbrs.at( i ) )
            {
                index[i] = j;
                break;
            }
        }
       // cout << " INDICES " << index[i] << endl;
       //     }
 }
}

#define RED "\033[01;31m"
#define GREEN "\033[22;32m"
#define YELLOW "\033[22;33m"
#define BLUE "\033[22;34m"
#define MAGENTA "\033[22;35m"
#define CYAN "\033[22;36m"
#define RESET "\033[22;0m"

#endif
