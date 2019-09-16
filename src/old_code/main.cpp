#include "typedefs.h"
#include "communicate.h"
#include "datatype.h"
#include "forest.h"
#include "templateForest.h"
#include "phdf5.h"
#include "templatePhdf5.h"
#include "tree.h"
#include "definitions.h"
#include "scale.h"

/** \mainpage
 *
 *
 *  <STRONG>
 *   "rebl-AMR: An Open-Source Software for Binarized-Octree Mesh Generation around Immersed Geometries"
 * </STRONG>
 *
 * \details
 * This software generates forest of octrees for an immersed geometry. <br>
 * reble-AMR uses red-black tree data strcture to store the nodes in a Z-curve compliant fashion. <br>
 * This project is part of the <STRONG>NSF-GEM3D </STRONG>Award No. 1440638. <br>
 *
 *	Required Libraries:
 *
 *      (1)  CMAKE
 *      (2)  MPI (MPI 3.0 standard compliant version)
 *      (3)  Zoltan
 *      (4)  ParMetis
 *      (5)  HDF5
 *
 *      Usage:
 *      mpirun -np <number of processes> progName input/myGeometry.stl <params.txt
 *
 * @authors     Jaber J. Hasbestan and Inanc Senocak (PI) <br>
 *
 * @date        Feb 2018
 *
 * @copyright  (c)
 * Department of Mechanical Engineering and Materials Science <br>
 * University of Pittsburgh, Pittsburgh, PA <br>
 *
 *
 *  \image html f16_bunny.png
 *
 *  <STRONG>
 *	Subject to GPL 3.0 license
 *  </STRONG>
 *
*/

int main( int argcs, char *pArgs[] )
{
 

 if ( argcs != 4 )
    {
        printf( RED "Geometry Not Supplied\n" RESET );
        cout << RED "argc=" << argcs << endl;
    }



 int proclevel = atoi (pArgs[2]);
 int meshlevel = atoi (pArgs[3]);


    MPI_Init( &argcs, &pArgs );

if ( WEAK == 0 )
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
    real xyz1[6] = {-2.0f, 2.0f, -2.0f, 2.0f, -2.0f, 2.0f};

    // this needs to be fixed, we need to open the damn file in parallel
    // and everybody should read only the portion that resides in its portion
    // this will be done whenever I can ..., HAS TO BE DONE THOUGH
  //  const uint treesize = PROCSIZE;
  //  const uint forsize = TREESIZE;
    readSTLGeom( argcs, pArgs, &geom_xyz, &geom_nn, xyz1 );


  // delete_ship(geom_xyz, &geom_nn);

    // cout << "geometrysize" << geom_nn << endl;

//    uint npx = 8, npy = 8, npz = 8;
 //   Tree<64, real> myTree( ancestorlength, ancestorcoords, npx, npy, npz );
 //   Tree<64, uint> myTree1( ancestorlength, ancestorcoords, npx, npy, npz );
 //   Tree<32, real> myTree2( ancestorlength, ancestorcoords, npx, npy, npz );
    Tree<PROCSIZE, uint> proc( ancestorlength, ancestorcoords, npx, npy, npz );

    // cout << "size of my tree " << proc.size() << endl;
    uint mylevel;

    // remove unnecessary files
    if ( my_rank == 0 )
    {
        system( "exec rm -r /sol/*" );

    }

/*    uint proclevel;
    uint meshlevel;
    uint levels[2];



    if ( my_rank == 0 )
    {
        cout << GREEN << "enter the refinement level for processor topology" << RESET << endl;
        cout << GREEN << "enter the refinement level for each tree" << RESET << endl;

//        cin >> proclevel;

  //      cin >> meshlevel;
        levels[0] = proclevel;
        levels[1] = meshlevel;

        assert( proclevel > 0 );
        assert( meshlevel >= 0 );

        cout << " proclevel " << proclevel << " refinelevel " << meshlevel << endl;
    }
    // will be in config file
    //
    //

    MPI_Bcast( levels, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    proclevel = levels[0];
    meshlevel = levels[1];
   */
    //    cout<<"proclevel "<<proclevel<<endl;
    //    cout<<"meshlevel "<<proclevel<<endl;

    if ( proclevel > ( PROCSIZE / 3 ) || meshlevel > ( TREESIZE / 3 ) )
    {
        cout << proclevel << " " << PROCSIZE / 3 << endl;
        throw std::runtime_error( RED "template size is not big to fit all the levels" RESET );
    }

    Center_coords XYZ;
    double t1;
    MPI_Barrier( MPI_COMM_WORLD );
    t1 = MPI_Wtime();

    //**************************************************************
    //
    //                 generate processor topology
    //
    //**************************************************************

    proc.refine( 0 );

    proc.convertStl2Morton( geom_nn, geom_xyz );

    for ( uint j = 0; j < proclevel - 1; j++ )
    {
        proc.pushToRefinelist( j + 1 );
        proc.fourToOne();
        proc.refineRefineList();
        //    cout<<j<<" : "<<proc.size()<<endl;
    }
  //  cout << "proc done \n" << endl;
    double t2;
    MPI_Barrier( MPI_COMM_WORLD );
    t2 = MPI_Wtime();

    //****************************************************************
    // clock_t t1 = clock();
    // cout<<"rank "<<my_rank<<" topology creation time "<<float(t1-t0)/CLOCKS_PER_SEC<<endl;

    // find the weights for each cube


    real xyz[6];
    uint count;
    morton<PROCSIZE> key;
// this maybe iomproved by searching the encoded geometry insetead of brute force search
#if ( 1 )

    //    if ( my_rank == 0 )

    {
        for ( auto it = proc.begin(); it != proc.end(); it++ )
        {
            // count is set to 1 such that an empty box has a weighting of one
            count = 1;

            key = it->first;
            proc.enclosingBox( key, xyz );
            it->second = (uint *)malloc( sizeof( uint ) );

            for ( uint j = 0; j < geom_nn; j++ )
            {
                if ( geom_xyz[3 * j + 0] >= xyz[0] && geom_xyz[3 * j + 0] <= xyz[1] )
                {
                    if ( geom_xyz[3 * j + 1] >= xyz[2] && geom_xyz[3 * j + 1] <= xyz[3] )
                    {
                        if ( geom_xyz[3 * j + 2] >= xyz[4] && geom_xyz[3 * j + 2] <= xyz[5] )
                        {
                            count++;
                        }
                    }
                }
            }
            // cout<<(it->second)<<endl;

            it->second[0] = count;
        }

        /*
           for(auto it=proc.begin();it!=proc.end();it++)
           {
           count=0;
           key=it->first;
           cout<<"size "<<it->second[0]<<endl;

           }
          */
    }

#else

// note that the 000,000, ..., 000 most probabaly will not happen unless we put the geometry way in the corner
//

#endif

    //    clock_t t2 = clock();
    // cout << "rank " << my_rank << " " << "weight distribution " << float( t2 - t1 ) / CLOCKS_PER_SEC << endl;

    if ( proc.size() < comsize )
    {
        // cout<<proclevel<<" "<<treesize/3<<endl;
        throw std::runtime_error( RED "please increase topology level : Topolgy level is not distributable among the processes\n" RESET );
    }

    int myvalue = proc.size();
    int offset, totalvalue;
    // cout<<proc.size()<<endl;
    // mpi.getOffset( myvalue, 0, &offset );
    // mpi.getTotalNumber( myvalue, 0, offset, &totalvalue );
    // cout << "totalvalue " << totalvalue << endl;

    /*======================================================================

                                 prepare for Zoltan

     =====================================================================*/

    double t3;
    double t4;

    real xyzc[3];
    uint co = 0;
    real *weight = nullptr;

    MPI_Barrier( MPI_COMM_WORLD );
    t3 = MPI_Wtime();

    if ( ZOLTAN_ON )
    {
 
    weight = new real[proc.size()];
    for ( auto it = proc.begin(); it != proc.end(); it++ )
    {
        proc.centroid( it->first, xyzc );
        XYZ.push_back( CenterCoords() );
        XYZ.at( co ).x = xyzc[0];
        XYZ.at( co ).y = xyzc[1];
        XYZ.at( co ).z = xyzc[2];
        if ( WEIGHT )
        {
            weight[co] = it->second[0];
        }
        else
        {
            weight[co] = 1.0;
        }
        // cout<<co<<"\t" <<weight[co]<<endl;
        // specify rank
        it->second[0] = 0;
        co++;
    }

    //******************************************************
    //
    //       This part every proc generates at tree
    //        and has all the info for proc topology
    //
    //******************************************************

       float ver;
        int rc;

        totalvalue = proc.size();

        rc = Zoltan_Initialize( argcs, pArgs, &ver );

        if ( rc != ZOLTAN_OK )
        {
            printf( RED "Problem in Initialzing Zoltan\n" RESET );
            MPI_Finalize();
            exit( 0 );
        }

        struct Zoltan_Struct *zz = NULL;
        //        zz = Zoltan_Create( MPI_COMM_WORLD );
        //        MPI_COMM_SELF is used to be able to utilize parallel function in serial
        zz = Zoltan_Create( MPI_COMM_SELF );
        // I define this sturcuture to mask sending all these arrays separately
        // and send only one structure instead

        Zoltan_Out zoltan_out;

        if ( ZOLTAN_GEOMETRIC_PARTITION )
        {
            zoltanGeometricPartitionerSerial( myvalue, totalvalue, offset, 1, zz, XYZ, weight, &zoltan_out, comsize );
        }
        //    if ( my_rank == 0 )
        // non scalable part for zoltan partitioning

        uint *ID = new unsigned int[proc.size()];

        for ( uint i = 0; i < proc.size(); i++ )
        {
            ID[i] = 0;
        }
        // elements to be exported are given,

        for ( int i = 0; i < zoltan_out.numExport; i++ )
        {
            ID[zoltan_out.exportLocalGids[i]] = zoltan_out.exportToPart[i];
        }
        /*
                    for ( int i = 0; i < proc.size(); i++ )
                    {
                     cout<<" mryank "<< my_rank<<" elem_id "<<i<<" procId "<<tmp[i]<<endl;
                    }
        */

        co = 0;

        for ( auto it = proc.begin(); it != proc.end(); it++ )
        {

            it->second[0] = ID[co];
            //                      cout <<"tag  "<< it->second[0] << endl;
            co++;
        }

        delete[] ID;
        /*
                        co=0;
                        for ( auto it = proc.begin(); it != proc.end(); it++ )
                        {

                                cout<<"my_rank "<<my_rank<<"   "<<" elem_id "<<co<<"  "<<it->second[0] <<endl;
                         //                      cout <<"tag  "<< it->second[0] << endl;
                               co++;
                        }
         */

        Zoltan_LB_Free_Part( &zoltan_out.importGlobalGids, &zoltan_out.importLocalGids, &zoltan_out.importProcs, &zoltan_out.importToPart );

        Zoltan_LB_Free_Part( &zoltan_out.exportGlobalGids, &zoltan_out.exportLocalGids, &zoltan_out.exportProcs, &zoltan_out.exportToPart );

        Zoltan_Destroy( &zz );
 
       delete[] weight;

    }

    //
    //
    //  this is only to check to see if the distribution was OK
    // note that Zoltan might give zero elements to one processes
    // leading to idling or possibly dead-locks
    //
    //
    //

    MPI_Barrier( MPI_COMM_WORLD );
    t4 = MPI_Wtime();
    double t5;
    double t6;

//    uint disc = 8;

  //  cout << "tree size " << treesize << endl;
  //  cout << "for size " << forsize << endl;
  //  cout << "ancestorlength " << ancestorlength[0] << endl;

    TemplateForest<TREESIZE, real, PROCSIZE, uint, Tree<PROCSIZE, uint>> forest( proc, ancestorlength, ancestorcoords, npx, npy, npz );

if(ZOLTAN_ON)
{
    forest.checkZoltanPartConsistency( proc );
}
    //  Forest<forsize, real, treesize, uint> forest( ancestorlength, ancestorcoords, proc,0, disc, disc, disc );

    //    Forest<forsize, real, treesize, uint > forest(ancestorlength, ancestorcoords,proc,fixedlevel, disc, disc, disc );

    /*
      if(my_rank==0)
       {
          ofstream myfile0;
            myfile0.open( "proc.txt" );

        for(uint j=0;j<comsize;j++)
        {
         for ( auto it = proc.begin(); it != proc.end(); it++ )
        {
            // cout<<co<<"\t" <<weight[co]<<endl;
            // specify rank
            //
        if(j==it->second[0])
        {
            myfile0<<it->first<<" "<<it->second[0] <<endl;;
        }
        }
         myfile0<<"-------------------------------"<<endl;
       }
    //   myfile0.close();
       }
    */

    forest.getMaxSeedsLevel( proc );


    forest.setMaxProcLevel( proclevel );

    forest.etomPatternConstruct( proc );

    forest.createCommGraph( 0 );

    forest.checkGraphConsistency();

    uint refinelevel = meshlevel;// this coordinate changes bunny
    real xx[3] = {0.8, 0.50, 0};

    MPI_Barrier( MPI_COMM_WORLD );
    t5 = MPI_Wtime();


    forest.moveGeom( proc, geom_xyz, geom_nn, xx );

    MPI_Barrier( MPI_COMM_WORLD );
    t6 = MPI_Wtime();

    forest.refineForestBalanced( refinelevel, proc );

    MPI_Barrier( MPI_COMM_WORLD );

    double t7 = MPI_Wtime();

    if ( my_rank == 0 )
    {
        ofstream myfile0;
        myfile0.open( "proc.txt" );
        //

        myfile0 << "-------------------------------" << endl;
        myfile0 << "-------------------------------" << endl;
        for ( uint j = 0; j < comsize; j++ )
        {
            for ( auto it = proc.begin(); it != proc.end(); it++ )
            {
                // cout<<co<<"\t" <<weight[co]<<endl;
                // specify rank
                //
                if ( j == it->second[0] )
                {
                    myfile0 << it->first << " " << it->second[0] << endl;
                    ;
                }
            }
            myfile0 << "-------------------------------" << endl;
        }
        myfile0.close();
    }

    /*

    //================================================
    // derefine section

      forest.pushToDerefineEachTree(refinelevel, proc );
      forest.pushToDerefineEachTree(refinelevel-1, proc );
      forest.pushToDerefineEachTree(refinelevel-2, proc );
      forest.pushToDerefineEachTree(refinelevel-3, proc );
//  forest.pushToDerefineEachTree(refinelevel-4, proc );

    //================================================




        xx[0]=0.5;

        forest.moveGeom( proc, geom_xyz, geom_nn, xx );

        forest.refineForestBalanced( refinelevel, proc );


   
        xx[0]=0.25;

        forest.moveGeom( proc, geom_xyz, geom_nn, xx );

        forest.refineForestBalanced( refinelevel, proc );

    */

    forest.getTotalMeshSize();

    MPI_Barrier( MPI_COMM_WORLD );

    double t8 = MPI_Wtime();

    /*
        std::bitset<64> test (std::string("1110001001000000000000000000000000000000000000000000000000000000"));

        uint nbrseedlevel;
        forest.getNbrSeedLevel( test, 4,&nbrseedlevel,proc );

         cout<<" ??????????????????????????????"<<nbrseedlevel<<endl;
    */
    // float ver;
    // int rc = Zoltan_Initialize( argcs, pArgs, &ver );
    //   struct Zoltan_Struct *z1;
    //   z1=Zoltan_Create( MPI_COMM_WORLD );

    //     forest.nonCollectiveNbrComm();

    //    templatePhdf5<forsize,TemplateForest<forsize, real, treesize, uint, Tree<treesize, uint> >> IO;
    uint index = 1;

    templatePhdf5<TREESIZE, real, PROCSIZE, uint, Tree<PROCSIZE, uint>> IO;

if ( WR == 1 )
    {
        IO.writePolyvertex( forest, index );
    }
else if( WR == 2 )
    {
        IO.writeMultiBlock( forest, index );
    }


#if(1)
    MPI_Barrier( MPI_COMM_WORLD );
    double t9 = MPI_Wtime();

    if ( my_rank == 0 )
    {
        ofstream myfile;
        myfile.open( "timing" ,ios_base::app);
        myfile << "\t Topology Level: \t" << proclevel << endl;
        myfile << "\t Tree Level: \t \t" << meshlevel << endl;
        myfile << "\t Topology Timing: \t " << t2 - t1 << endl;
        myfile << "\t Zoltan Timing: \t" << t4 - t3 << endl;

        myfile << "\t GeomEncoding Timing: \t" << t6 - t5 << endl;
        myfile << "\t Refine Timing: \t" << t7 - t6 << endl;
        myfile << "\t Global Mesh size: \t" << t8 - t7 << endl;
        myfile << "\t I/O Timing: \t \t" << t9 - t8 << endl;

        myfile<<"boolean for debug" <<(is_same< Tree<PROCSIZE,uint>, Tree<TREESIZE,uint> > ::value)<<endl;
        myfile.close();
    }
#endif

//     conditional<!is_same< Tree<PROCSIZE,uint>, Tree<TREESIZE,uint> > ::value,Tree<10,uint>,Tree<PROCSIZE,uint>>::type type ;
    // !is_same< PROCSIZE, TREESIZE> ::value;

}
//  treeProcessorTopology( argcs,pArgs);
//    fTreeProcessorTopology(argcs,pArgs);

/*
 *
 *
 *
 *
 * */

if ( WEAK == 1 )
{
#if ( 1 )
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
    // need to specify WSIZE in defenitions file

    MPI_Barrier( MPI_COMM_WORLD );

    double t4 = MPI_Wtime();

    uint refinelevel;

    const uint proclevel = WSIZE / 3;

    if ( my_rank == 0 )
    {
        //        cout << GREEN << "enter the refinement level for processor topology" << RESET << endl;
        cout << GREEN << "enter the refinement level for each tree" << RESET << endl;

        //      cin >> proclevel;

        cin >> refinelevel;
    }

    MPI_Bcast( &refinelevel, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    MPI_Barrier( MPI_COMM_WORLD );

    double t5 = MPI_Wtime();

    // This is working version

    const uint treesize = TREESIZE;
    const uint weaksize = TREESIZE;
    const uint forsize =  TREESIZE;

//    uint npx = 2, npy = 2, npz = 2;

    for ( uint i = 0; i < 3; i++ )
    {
        ancestorlength[i] = pow( 2., proclevel + 1 );
        xyz1[2 * i] = -ancestorlength[i];
        xyz1[2 * i + 1] = ancestorlength[i];
    }


   // FullTree<PROCSIZE, uint> proc( ancestorlength, ancestorcoords );
    FullTree<WSIZE, uint> proc( ancestorlength, ancestorcoords );

    proc.setLevel( proclevel );

    morton<proclevel * 3> seedKey( my_rank );

    cout << seedKey << endl;

//    morton<PROCSIZE> seed;

    morton<WSIZE> seed;


    for ( i = 0; i < 3 * proclevel; i++ )
    {
       // seed[TREESIZE - 1 - i] = seedKey[3 * proclevel - i - 1];
        seed[WSIZE - 1 - i] = seedKey[3 * proclevel - i - 1];
    }
    cout << seed << endl;

    proc.insertKey( seed );

    /*
       auto it=proc.begin();
       it->second=new uint[1];
       it->second[0]=(uint)my_rank;

       cout<<GREEN<<my_rank<<" "<<it->first<<" "<<it->second[0]<<RESET<<endl;
    */

    vector<uint> Nbrs;

    proc.nbrsConstrcut( Nbrs, my_rank );

    //  std::sort(Nbrs.begin(),Nbrs.end());

    proc.assignProcs( Nbrs, my_rank );

    if ( comsize != pow( 8, proclevel ) )
    {
        throw std::runtime_error( RED "inconsistent full tree toplogy" RESET );
    }

    readSTLGeom( argcs, pArgs, &geom_xyz, &geom_nn, xyz1 );

    double xc[3];
    proc.centroid( seed, xc );

    // cout<<" "<<xc[0]<<" "<<xc[1]<<" "<<xc[2]<<endl;

    double xmax = 0.0, xmin = 0.0;
    double ymax = 0.0, ymin = 0.0;
    double zmax = 0.0, zmin = 0.0;

    for ( int i = 0; i < geom_nn; i++ )
    {
        geom_xyz[3 * i] = geom_xyz[3 * i] + xc[0];
        geom_xyz[3 * i + 1] = geom_xyz[3 * i + 1] + xc[1];
        geom_xyz[3 * i + 2] = geom_xyz[3 * i + 2] + xc[2];

        if ( xmin > geom_xyz[3 * i] )
        {
            xmin = geom_xyz[3 * i];
        }
        if ( xmax < geom_xyz[3 * i] )
        {
            xmax = geom_xyz[3 * i];
        }

        if ( ymin > geom_xyz[3 * i + 1] )
        {
            ymin = geom_xyz[3 * i + 1];
        }
        if ( ymax < geom_xyz[3 * i + 1] )
        {
            ymax = geom_xyz[3 * i + 1];
        }
        if ( zmin > geom_xyz[3 * i + 2] )
        {
            zmin = geom_xyz[3 * i + 2];
        }
        if ( zmax < geom_xyz[3 * i + 2] )
        {
            zmax = geom_xyz[3 * i + 2];
        }
    }

    // cout<<my_rank<<" x "<<xmax<<" "<<xmin<< " "<<ymin <<" "<<ymax<<" "<<zmin<<" "<<zmax<<endl;

  //  uint disc = 16;

    //     cout<<GREEN"proc.size "<<proc.size()<<RESET<<endl;

  //  TemplateForest<PROCSIZE, real, treesize, uint, FullTree<treesize, uint>> forest( proc, ancestorlength, ancestorcoords, npx,npy, npz );
   // TemplateForest<TREESIZE, real, PROCSIZE, uint, FullTree<PROCSIZE, uint>> forest( proc, ancestorlength, ancestorcoords, npx,npy, npz );
    TemplateForest<TREESIZE, real, WSIZE, uint, FullTree<WSIZE, uint>> forest( proc, ancestorlength, ancestorcoords, npx,npy, npz );

    forest.getMaxSeedsLevel( proc );

    forest.setMaxProcLevel( proclevel );

    forest.comPatternConstruct( proc, Nbrs );

    real xx[3] = {0.0, 0.0, 0.0};

    forest.moveGeom( proc, geom_xyz, geom_nn, xx );

    MPI_Barrier( MPI_COMM_WORLD );

    double t6 = MPI_Wtime();

    /*   cout << proc.size() << endl;
     *   fstream myfile;

       for ( uint i = 0; i < Nbrs.size(); i++ )
       {
          cout<<RED<<"my_rank "<<my_rank<<" "<<Nbrs.at(i)<<RESET<<endl;
      }
    */
    /*
    uint direction = 2;
    if ( Nbrs.size() == 6 )
    {
        cout << RED << "my_rank " << my_rank << " " << Nbrs.size() << RESET << endl;
    }
    //      cout<<"my_rank "<<my_rank<<"        "<<proc.size()<<"  "<< Nbrs.size()<<endl;
    */

    forest.refineForestBalanced( refinelevel, proc );

    MPI_Barrier( MPI_COMM_WORLD );

    double t7 = MPI_Wtime();
    /*
    if(my_rank==0)
    {
    cout<<"time refinement "<<t7-t6<<endl;
    cout<<"time preparation "<<t6-t5<<endl;
    cout<<" time to solution "<<t7-t4<<endl;
    }
    */

    forest.getTotalMeshSize();

   // templatePhdf5<WSIZE, real, treesize, uint, FullTree<treesize, uint>> IO;
    templatePhdf5<TREESIZE, real, WSIZE, uint, FullTree<WSIZE, uint>> IO;
   // templatePhdf5<TREESIZE, real, PROCSIZE, uint, FullTree<PROCSIZE, uint>> IO;


    uint index = 2;

if ( WR == 1 )
    {
        IO.writePolyvertex( forest, index);
    }
else if( WR == 2 )
    {
        IO.writeMultiBlock( forest, index);
    }
#endif

    MPI_Barrier( MPI_COMM_WORLD );

    double t8 = MPI_Wtime();
 
    if ( my_rank == 0 )
    {
        ofstream myfile;
        myfile.open( "timing" );
        myfile << "\t preparation-time = " << t6 - t5 << endl;
        myfile << "\t refinement-time= " << t7 - t6 << endl;
        myfile << "\t total run-time without I/O= " << t7 - t4 << endl;
        myfile << "\t I/O= " << t8 - t7 << endl;
        myfile.close();
    }

}

    MPI_Finalize();
}

//    readSTLGeom( argcs, pArgs, &geom_xyz, &geom_nn, xyz1 );
// Jerry need to pass like this delete_ship(geom_xyz,&geom_nn)
void delete_ship(real *geom_xyz, int *geom_nn)
{
int newsize=0;
real *geom_tmp=nullptr;

for(int i=0;i<(*geom_nn);i++)
{
if((geom_xyz)[3*i+2]>1.e-6)
{
newsize++;
}
}

geom_tmp=new real[3*newsize];

int count=0;
for(int i=0;i<(*geom_nn);i++)
{
if(geom_xyz[3*i+2]>1.e-6)
{
geom_tmp[3*count]=geom_xyz[3*i];
geom_tmp[3*count+1]=geom_xyz[3*i+1];
geom_tmp[3*count+2]=geom_xyz[3*i+2];
count++;
}

}



*geom_nn=count;

//not necessary
//**geom_xyz=(double*)realloc((**geom_xyz),count*sizeof(double));


for(int i=0;i<3*count;i++)
{
geom_xyz[i]=geom_tmp[i];
}


delete[] geom_tmp;

}






