#include "typedefs.h"
#include "communicate.h"
#include "datatype.h"
#include "forest.h"
#include "phdf5.h"
#include "tree.h"
#include "scale.h"
#define PROCSIZE 64
#define TREESIZE 64


void fTreeProessorTopology( int argcs, char *pArgs[] )
{
   unsigned int i, j, k, l;
    int my_rank, comsize;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &comsize );
    int Dim[3] = {0, 0, 0};
    MPI_Comm Comm_cart;


    real ancestorlength[3] = {2.0f, 2.0f, 2.0f};
    real ancestorcoords[3] = {0.0f, 0.0f, 0.0f};
    real xyz2[3];
    real *geom_xyz = NULL;
    int geom_nn;

    real xyz1[6] = {-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f};

    // this needs to be fixed, we need to open the damn file in parallel
    // and everybody should read only the portion that resides in its portion
    // this will be done whenever I can ..., HAS TO BE DONE THOUGH

   uint meshlevel=0;

   const uint proclevel=WSIZE/3;



 if ( my_rank == 0 )
    {
//        cout << GREEN << "enter the refinement level for processor topology" << RESET << endl;
        cout << GREEN << "enter the refinement level for each tree" << RESET << endl;

  //      cin >> proclevel;

        cin >> meshlevel;
    }


    MPI_Bcast( &meshlevel, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    const uint treesize = proclevel*3;
    const uint weaksize = treesize;
    const uint forsize = TREESIZE;

    uint npx = 2, npy = 2, npz = 2;

  for(uint i=0;i<3;i++)
    {
      ancestorlength[i]=pow(2.,proclevel+1);
      xyz1[2*i]= -ancestorlength[i];
      xyz1[2*i+1]= ancestorlength[i];
    }

    Tree<treesize, uint> proc( ancestorlength, ancestorcoords, npx, npy, npz );

    if(comsize!=pow(8,proclevel))
    {
     throw std::runtime_error(RED "inconsistent full tree toplogy" RESET);
    }

    readSTLGeom( argcs, pArgs, &geom_xyz, &geom_nn, xyz1 );

    // cout << "geometrysize" << geom_nn << endl;

    FullOctreeTop<weaksize> myscale;

    myscale.convertRank2Bits( my_rank );

    myscale.constructNbrProcs();

//    myscale.checkGraphConsistency( my_rank );

//*****************************************************
//
//              assign rootKey to proc
//
// ****************************************************

   morton<weaksize> ktmp;
   myscale.readRoot(ktmp);
   proc.insertSeed(ktmp);

/*
/*
   FullTree<treesize,real> fullTree(ancestorlength, ancestorcoords);
   fullTree.insertSeed(ktmp);
   fullTree.setLevel(proclevel); 
*/
   real xyz0[6];
   proc.enclosingBoxFixedLevel(ktmp,proclevel,xyz0);

 /* 
   cout<<my_rank<<" "<< xyz0[0]<<" "<<  xyz0[1]<<" "<< xyz0[2]<<" "<<  xyz0[3]<<" "<< xyz0[4]<<" "<<  xyz0[5]<<endl;

  for(uint i=0;i<6;i++)
  {
   xyz0[i]=0.0;
  }

   fullTree.enclosingBox(ktmp,xyz0);

   cout<<my_rank<<" "<< xyz0[0]<<" "<<  xyz0[1]<<" "<< xyz0[2]<<" "<<  xyz0[3]<<" "<< xyz0[4]<<" "<<  xyz0[5]<<endl;
*/
   auto it=proc.begin();

   it->second=new uint[1];
   it->second[0]=my_rank;

//***********************************
//
//           construct forest       
//
//********************************** 

   vector<int>Nbrs;
   myscale.readNbrs(Nbrs);

   morton<treesize>key;

   proc.insertNbrs(Nbrs);

   for(auto it=proc.begin();it!=proc.end();it++)
  {
      it->second=new uint[1];
      it->second[0]=it->first.to_ulong();
 //   cout<<my_rank<<" "<<RED<<it->first<<"  "<<it->second[0]<<RESET<<endl;
}

 // to remove the singularity from 

   uint disc = 2;
   Forest<forsize, real, treesize, uint> forest( ancestorlength, ancestorcoords, proc,proclevel, disc, disc, disc );

   forest.getMaxSeedsLevel( proc );

/* 
for(uint i=0;i<Nbrs.size();i++)
{
cout<<my_rank<<"  "<<Nbrs.at(i)<<endl;
} 
*/

 forest.constructCommWeak(Nbrs);
// bottle-neck in move geometry, grows rapidly with meshlevel


    real xx[3] = {0.0, 0, 0};

   const double denom=1./pow(2.0,proclevel-1);
 // cout<<" denom "<<denom<<endl;

  real xc[3]={0.0,0.0,0.0};

/*
for(int i=0;i<geom_nn;i++)

{
 geom_xyz[3*i]=geom_xyz[3*i]*denom;
 geom_xyz[3*i+1]=geom_xyz[3*i+1]*denom;
 geom_xyz[3*i+2]=geom_xyz[3*i+2]*denom;
//cout<<geom_xyz[3*i]<<endl;
}
*/
myscale.readRoot(ktmp);


proc.centroidFixedLevel(ktmp,proclevel,xc);
//cout<<" "<<xc[0]<<" "<<xc[1]<<" "<<xc[2]<<endl;

#if(1)
for(int i=0;i<geom_nn;i++)
{
 geom_xyz[3*i]=geom_xyz[3*i]+xc[0];
 geom_xyz[3*i+1]=geom_xyz[3*i+1]+xc[1];
 geom_xyz[3*i+2]=geom_xyz[3*i+2]+xc[2];

}

    forest.moveGeom( proc,proclevel, geom_xyz, geom_nn, xx );
    Phdf5<forsize, real, treesize, uint> IO;
    uint index = 0;
    forest.refineForestBalanced( meshlevel ,proc);
    forest.getTotalMeshSize();

    IO.writeMultiBlock( forest, 1 );

  // IO.writePolyvertex( forest, 1 );

#endif

}
