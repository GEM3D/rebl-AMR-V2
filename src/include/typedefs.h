#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_
#include "definitions.h"
#include "zoltan.h"


/*!
 * \struct CenterCoords 
 * \brief Stores the coordinate of the centroid of the elments
 */ 


typedef struct{
	real x;
	real y;
	real z;	
} CenterCoords;


typedef std::vector<CenterCoords>Center_coords;


/*!
 * \struct Zoltan_Out
 * \brief This structure is an interface to store the output from Zoltan
 *
 *
 */
typedef struct { 
    int changes;
    int numGidEntries;
    int numLidEntries;
    int numImport;
    int numExport;
    unsigned int* importGlobalGids;
    unsigned int* importLocalGids;
    unsigned int* exportGlobalGids;
    unsigned int* exportLocalGids; 
    int *importProcs;
    int *importToPart;
    int *exportProcs;
    int *exportToPart;
    int *parts;  
} Zoltan_Out;

bool zoltanGeometricPartitioner(const uint size,const uint ncube_total,const uint offset,const int method,struct Zoltan_Struct *zz,const Center_coords &XYZ,real *weight,Zoltan_Out *zoltan_out);

bool zoltanGeometricPartitionerSerial(const uint size,const uint ncube_total,const uint offset,const int method,struct Zoltan_Struct *zz,const Center_coords &XYZ,real *weight,Zoltan_Out *zoltan_out,int comsize);

void treeProcessorTopology( int argcs, char *pArgs[] );
//void fTreeProcessorTopology( int argcs, char *pArgs[] );



void readSTLGeom( int argc, char *argv[], real **triangle_center, int *nn, const real *xyz );
void delete_ship(real *geom_xyz, int *geom_nn);

inline void TwoPowN( uint b, real *result )
{
    uint two = 1;
    two      = two << b;
    *result  = (real)two;
};








#endif
