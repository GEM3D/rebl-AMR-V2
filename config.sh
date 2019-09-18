#!/bin/bash

GREEN='\033[22;32m'
RESET='\033[0m'
RED='\033[01;31m'
#module purge


string=$HOSTNAME
if [[ $string == *"comet"* ]]; then
    export cluster=comet
elif [[ $string == *"bridges"* ]]; then
    export cluster=bridges
elif [[ $string == *"stampede"* ]]; then
    export cluster=stampede
elif [[ $string == *"crc"* ]]; then
    export cluster=crc
elif [[ $string == *"titan"* ]]; then
    export cluster=titan
fi


printf "${RED}\n \n       System is =  $cluster ${RESET}\n \n"


read -p "To Configure for GPU press 1 else enter 0:  " GPU
#echo $GPU

printf "${GREEN} \n Configuring ... ... ${RESET}\n \n"


if [ $cluster == "bridges" ]; then

    export HDF5_ROOT=/pylon5/ac561fp/hasbesta/packages/hdf5-1.8.17/
    module load psc_path/1.1
    module load slurm/default
    module load xdusage/2.0-4
    module load cmake/3.7.2
    module load parmetis
    module load zoltan/3.83


    if [ $GPU -lt 1 ]; then
        module load /opt/packages/openmpi/openmpi-2.0.1/modulefiles/gcc_openmpi-2.0.1
        export CXX=/opt/packages/openmpi/openmpi-2.0.1/bin/mpiCC;
        export CC=/opt/packages/openmpi/openmpi-2.0.1/bin/mpicc;
        export GPU1=0;
    else
        module load mpi/pgi_openmpi.18.1
        export MPI_ROOT=/opt/packages/pgi/linux86-64/2018/mpi/openmpi-2.1.2/
        export CC=$MPI_ROOT/bin/mpicc;
        export CXX=$MPI_ROOT/bin/mpiCC;
        export GPU1=1;
    fi

fi


if [ $cluster == "comet" ]; then

    export MODULEPATH=/opt/modulefiles/applications/.gnu:/opt/modulefiles/mpi/.gnu:/opt/modulefiles/mpi:/opt/modulefiles/compilers:/opt/modulefiles/applications:/usr/share/Modules/modulefiles:/etc/modulefiles:/share/apps/compute/modulefiles/mpi:


    export MODULEPATH=/opt/modulefiles/applications/.gnu:/opt/modulefiles/mpi/.gnu:/opt/modulefiles/mpi:/opt/modulefiles/compilers:/opt/modulefiles/applications:/usr/share/Modules/modulefiles:/etc/modulefiles:/share/apps/compute/modulefiles/mpi:

    #system needs
    module load cmake/3.9.1
    module load hdf5
    module load intel/2013_sp1.2.144
    module load gnutools/2.69
    module load parmetis
    module load trilinos

    export HDF5HOME=/opt/hdf5/intel/openmpi_ib
    ##module load .intel/mvapich2_ib/2.1

    if [ $GPU -lt 1 ]; then
        module load openmpi_ib/2.0.4_gnu
        #module load fftw/3.3.4
        export CC=$MPIHOME/bin/mpicc
        export CXX=$MPIHOME/bin/mpiCC
        export GPU1=0;
    else
        module load mpi/pgi_openmpi.18.1
        export MPI_ROOT=/opt/packages/pgi/linux86-64/2018/mpi/openmpi-2.1.2/
        export CC=$MPI_ROOT/bin/mpicc;
        export CXX=$MPI_ROOT/bin/mpiCC;
        export GPU1=1;
    fi

fi



if [ $cluster == "stampede" ]; then

    #module load intel/18.0.2  impi/18.0.2
    module load phdf5/1.8.16

    module load gcc/7.1.0  impi/17.0.3
    module load intel/17.0.4  impi/17.0.3
    module load zoltan/3.83


    #export HDF5HOME=/opt/hdf5/intel/openmpi_ib
    ##module load .intel/mvapich2_ib/2.1

    if [ $GPU -lt 1 ]; then
        #module load openmpi_ib/2.0.4_gnu
        #module load fftw3/3.3.8
        export TACC_PARMETIS_DIR=/home1/apps/intel17/impi17_0/petsc/3.7/knightslanding/;
        # stand alone zoltan gives funky link errors
        #export TACC_ZOLTAN_DIR=/home1/apps/intel17/impi17_0/petsc/3.9/skylake/;
        export TACC_ZOLTAN_DIR=/home1/apps/intel17/impi17_0/trilinos/12.10.1/;
        export CC=$TACC_IMPI_DIR/intel64/bin/mpicc
        export CXX=$TACC_IMPI_DIR/intel64/bin/mpicxx
        export GPU1=0;
    else
        echo "No GPU is available on Stampede"
    fi

fi


if [ $cluster == "crc" ]; then

    module purge
   # module load intel/2017.1.132
   # module load intel-mpi/2017.1.132
    module load  intel/2018.2.199
    module load intel-mpi/2018.2.199
    module load hdf5/1.10.0
    module load  cmake/3.7.1
    module load parmetis/4.0.3
    module load zoltan/3.83.0



    if [ $GPU -lt 1 ]; then

        export CC=$I_MPI_ROOT/intel64/bin/mpicc
        export CXX=$I_MPI_ROOT/intel64/bin/mpiicpc
        export GPU1=0;
        export MPI_ROOT=$I_MPI_ROOT/intel64;
    else
        echo "No GPU is not coded yet for CRC"
    fi

fi


if [ $cluster == "titan" ]; then

    #system needs
    module load cmake3
    module unload PrgEnv-pgi
    module unload pgi
    module load PrgEnv-gnu/5.2.82
    module load cray-mpich/7.6.3
    module load cray-hdf5-parallel
    module load cray-trilinos
    module load cray-tpsl

    export CRAYPE_LINK_TYPE=dynamic

    #if [ $GPU -lt 1 ]; then
    #module load openmpi_ib/2.0.4_gnu
    #module load fftw/3.3.4
    #export CC=$MPIHOME/bin/mpicc
    #export CXX=$MPIHOME/bin/mpiCC
    #export GPU1=0;
    #else
    #module load mpi/pgi_openmpi.18.1
    #export MPI_ROOT=/opt/packages/pgi/linux86-64/2018/mpi/openmpi-2.1.2/
    #export CC=$MPI_ROOT/bin/mpicc;
    #export CXX=$MPI_ROOT/bin/mpiCC;
    #export GPU1=1;
    #fi

fi

if [ ! -d "./soln" ] || [ ! -d "./build" ]; then 
    printf "${RED} Creating Required Directories ${RESET}\n \n"
fi

if [ ! -d "./soln" ]; then
    mkdir soln
fi

if [ ! -d "./build" ]; then
    mkdir build
fi

printf "${GREEN} Required modules have been loaded for your system! ${RESET}\n \n"
printf "${GREEN} Change directory to buid/ and type cmake .. and make sure it finds all packages ${RESET}\n \n"
printf "${GREEN} After successful cmake, type make to compile the code. make -j n to compile in parallel with n processes ${RESET}\n \n"

################################################################

#                  FORMATING the SHELL SCRIPT

################################################################

#1 Load the file into Emacs "emacs config.sh"
#2 Press Ctrl-space at the top of the file
#3 Move the cursor to the bottom of the file
#4 Press Alt-X and type untabify then return
#4 Press Alt-X and type indent-region then return


