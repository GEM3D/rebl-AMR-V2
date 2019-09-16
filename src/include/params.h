#ifndef _PARAMS_H_
#define _PARAMS_H_

/*!
  * \brief  parameters to vconfigure the code for strong and weak scalings
 */

#define METHOD 2
#define MYSCALE 0.4    /* Input geometry scale unit*/
#define PPDTL 0.0  /* Box tolerance for object derefinement */

enum ReblAmrParameters: uint
{
    PROCSIZE = 32, /*!< size of container in bits for each seed, i.e. topology*/
    TREESIZE = 32, /*!< size of the container in bits for each tree */  
    ZOLTAN_ON =1 , /*!< turn zoltan on and off*/
    ZOLTAN_GEOMETRIC_PARTITION = 1, /*!< partition using geometric methods HSFC  */
    WEIGHT = 1, /*!< use the eighted version of HSFC */
    WR  = 2, /*!<controls the output, set 0 to turn it off, set to 1 for polyvertex and set 2 for multiblock*/
    WEAK = 0, /*!< set to 0 for weak scaling to 1 for strong, default is 0  */
    WSIZE = 3, /*!< size of the container for full tree  */
    REORDER = 0, /*!< redoredr processors in MPI_Ineighbir.... */
    OVERLAP = 0, /*!< overlap communication and computation using non-blocking collectives */
    PART_METHOD= 1, /*!< Specift the partitoning method*/     
  
//    ZOLTAN_DEBUG=1, /*!< defines the debug level of zoltan*/  
    /*!< discretization in x,y, and zdirections */
    npx=5,
    npy=3,
    npz=3,
    CHECK_MESH=0,
    MPI_ERROR_DISABLE  = 0,        /*! if set to 1 it will reset MPI_ERROR_FATAl to MPI_ERROR_RETURN */
};


typedef enum ReblAmrErrorCodes {
    SUCCESS                        = 0,
    NUM_INPUT_ARGS                 = 1,
    MPI_INIT_CHECK_FAIL            = 2,
    MPI_INIT_FAIL                  = 3,
    MPI_DUP_FAIL                   = 4,
    COMSIZE_FAIL                   = 5,
    PROC_LEVEL                     = 6,
    MESH_LEVEL                     = 7,
    MPI_GET_RANK_FAIL              = 8,
    MPI_COMSIZE_FAIL               = 9,
    COMBINED_SIZE                  = 10,
    NO_SEED                        = 11,
    GRAPH_CREATE_FAIL              = 12,
    MPI_INEIGHBOR_FAIL_ZX          = 13,
    MPI_INEIGHBOR_FAIL_XY          = 14,
    BLOCK_NUMBER_FAIL              = 15,
    ALLOCATION_FAIL                = 16,
    MPI_FINALIZE_FAIL              = 17,
    MPI_ERROR_HANDLE_FAIL          = 18,
    CONNECTIVITY_CONSTRUCTION_FAIL = 19,
    THOMAS_FAIL                    = 20,
} ReblAmrResult;


const char *   ReblAmrGetErrorEnum( ReblAmrResult error ); /*!<Information in tex form for Exit Codes */


// #define PROCSIZE 32
// #define TREESIZE 32

/* ZOLTAN =1 turns partitioning on set to zero to turn it off*/

//#define ZOLTAN_ON 1
//#define WEIGHT 1

/* if you set the value to 1 to use geometirc based partitioning */
// #define ZOLTAN_GEOMETRIC_PARTITION 1

// #define WEAK 6
// #define WSIZE 6
// #define WR 0 /*! set WR to zero it will not write out, set to 1 writes out polyVertex and 2 will write multiBloc  */

//#define SENDBOOL 0
// #define REORDER 0
// #define OVERLAP 0

//#if ( METHOD == 3 || METHOD == 4 )
// #undef OVERLAP
// #define OVERLAP 1
//#endif

#endif
