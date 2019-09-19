#include "typedefs.h"
#include "communicate.h"
#include "datatype.h"
#include "forest.h"
#include "phdf5.h"
#include "tree.h"
#include "scale.h"
#include "definitions.h"

#if(0)
#define PROCSIZE 64
#define TREESIZE 64

/* ZOLTAN =1 turns partitioning on set to zero to turn it off*/

#define ZOLTAN_ON 1
#define WEIGHT 1
#define ZOLTAN_GEOMETRIC_PARTITION 1
#endif

void treeProcessorTopology( int argcs, char *pArgs[] )
{
#if(0)
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

    const uint forsize = PROCSIZE;
    const uint treesize = TREESIZE;

    readSTLGeom( argcs, pArgs, &geom_xyz, &geom_nn, xyz1 );
    // cout << "geometrysize" << geom_nn << endl;
    uint npx = 2, npy = 2, npz = 2;

    //Tree<64, real> myTree( ancestorlength, ancestorcoords, npx, npy, npz );
    //Tree<64, uint> myTree1( ancestorlength, ancestorcoords, npx, npy, npz );
    //Tree<32, real> myTree2( ancestorlength, ancestorcoords, npx, npy, npz );
    Tree<treesize, uint> proc( ancestorlength, ancestorcoords, npx, npy, npz );

    // cout << "size of my tree " << proc.size() << endl;

    uint mylevel;

    // remove unnecessary files
    if ( my_rank == 0 )
    {
        system( "exec rm -r /sol/*" );
    }

    uint proclevel;
    uint meshlevel;
    uint levels[2];

    if ( my_rank == 0 )
    {
        cout << GREEN << "enter the refinement level for processor topology" << RESET << endl;
        cout << GREEN << "enter the refinement level for each tree" << RESET << endl;

        cin >> proclevel;

        cin >> meshlevel;
        levels[0] = proclevel;
        levels[1] = meshlevel;
    }
    // will be in config file
    //
    //
    MPI_Bcast( levels, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    proclevel = levels[0];
    meshlevel = levels[1];

    if ( proclevel > ( treesize / 3 ) || meshlevel > ( forsize / 3 ) )
    {
        cout << proclevel << " " << treesize / 3 << endl;
        throw std::runtime_error( RED "template size is not big to fit all the levels" RESET );
    }

    Center_coords XYZ;
    double t1;
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
        //   cout<<j<<" : "<<proc.size()<<endl;
    }
    //  cout << "proc done \n" << endl;
    double t2;
    t2 = MPI_Wtime();

    //****************************************************************
    // clock_t t1 = clock();
    // cout<<"rank "<<my_rank<<" topology creation time "<<float(t1-t0)/CLOCKS_PER_SEC<<endl;

    // find the weights for each cube

    real xyz[6];
    uint count;
    morton<treesize> key;
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

     ======================================================================*/

    double t3;
    double t4;

    real xyzc[3];
    uint co = 0;
    real *weight = nullptr;
    weight = new real[proc.size()];

    t3 = MPI_Wtime();

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

    if ( ZOLTAN_ON )
    {
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
    }

    //
    //
    //  this is only to check to see if the distribution was OK
    // note that Zoltan might give zero elements to one processes
    // leading to idling or possibly dead-locks
    //
    //
    //
    //

    t4 = MPI_Wtime();
    //
    uint disc = 2;

#if ( 1 )
    Forest<forsize, real, treesize, uint> forest( ancestorlength, ancestorcoords, proc, 0, disc, disc, disc );

    //    forest.checkZoltanPartConsistency(proc);

    // notice that the we first add the list to derefeine and then move
    // this saves on memory since I do not have to save the mortonSTl from previous step

    // proc.printMesh();
    double t5;
    double t6;

    forest.getMaxSeedsLevel( proc );

    forest.comPatternConstruct( proc );

    uint refinelevel = meshlevel;

    //   forest.createNbrsOfNbrs();
    //   forest.checkNbrsOfNbrsConsistency();
    // forest.createCommGraph( 0 );
    //   forest.checkGraphConsistency();
    //   forest.rcvrMessageSize( );
    // initial location
    /*
       uint incount,outcount;
         uint * weigh;
        MPI_Dist_graph_neighbors_count(graphComm, &incount,&outcount,weigh);
    */
    real xx[3] = {0.0, 0, 0};

    t5 = MPI_Wtime();
    // bottle-neck in move geometry, grows rapidly with meshlevel

    forest.moveGeom( proc, 0, geom_xyz, geom_nn, xx );

    t6 = MPI_Wtime();
    forest.refineForestBalanced( refinelevel, proc );

    double t7 = MPI_Wtime();
//    forest.zoltanGeomrepart(proc, 1 );

// cout<<"??????????????????????????????????????????"<<endl;
// forest.pushToDerefineEachTree( refinelevel,proc);

#if ( 0 )
    FullOctreeTop<6> myscale;
    myscale.convertRank2Bits( my_rank );

    myscale.constructNbrProcs();
    myscale.checkGraphConsistency( my_rank );

#endif
/*

#if ( 0 )
    xx[1] = -0.5;

    forest.moveGeom( proc, geom_xyz, geom_nn, xx );
    //    refinelevel=refinelevel-2;
    forest.refineForestBalanced( refinelevel, proc );

#endif
    //
    //
    // Debug derefine
    //
    //

    // forest.debugDerefine(proc );

    /*
         xx[0]=1.0;
        forest.moveGeom( proc, geom_xyz, geom_nn,xx);
        forest.refineForestBalanced( refinelevel ,proc);
    /*
    for(uint i=0;i<proc.size();i++)
    {
    cout<<<<endl;
    }
    */
//
//
// refine by hand to debug
//
//

//  forest.debug(proc);
//     cout<<"total_size "<<forest.getTotalSize()<<endl;
#if ( 1 )
    forest.getTotalMeshSize();
    Phdf5<forsize, real, treesize, uint> IO;
    uint index = 0;

    /*

       if(refinelevel>5)
      {   IO.writePolyvertex(forest,index);
     }
     else
     {
         IO.writeMultiBlock( forest, 0 );
     }
       */

    IO.writeMultiBlock( forest, 1 );
#endif
    //      IO.writePolyvertex(forest,index);

    //     IO.writeMultiBlock( forest, 0 );

    delete[] geom_xyz;
    delete[] weight;
    int id;

    if ( my_rank == 0 )
    {
        cout << "Topology Level: " << proclevel << endl;
        cout << "Tree Level: " << meshlevel << endl;
        cout << "Topology Timing: " << t2 - t1 << endl;
        cout << "Zoltan Timing: " << t4 - t3 << endl;

        cout << "GeomEncoding Timing: " << t6 - t5 << endl;
        cout << "Refine Timing: " << t7 - t6 << endl;
    }

#endif
#endif

}
