#!/bin/csh -f
unalias *
set echo
#--------------------------------------------------------------------------------------------------------
set platform  = taiyuan.intel                        # A unique identifier for your platform
set template  = $cwd/../bin/mkmf.template.$platform     # path to template for your platform
set mkmf      = $cwd/../bin/mkmf                        # path to executable mkmf
set sourcedir = $cwd/../src                             # path to directory containing model source code
set pathnames = $cwd/path_names                         # path to file containing list of source paths
set ppdir     = $cwd/../postprocessing                  # path to directory containing the tool for combining distributed diagnostic output files
#--------------------------------------------------------------------------------------------------------
set execdir = $cwd/exec.$platform  # where code is compiled and executable is created
set executable = $execdir/full_dynamics.x

#source $MODULESHOME/init/csh
#module rm netcdf-3.6.1 netcdf-3.6.2 netcdf-3.6.3 netcdf-4.0.1 netcdfAll netcdf
#module rm hdf5
#module load ifort/11.1.073
#module load intel_compilers
#module use /home/sdu/publicmodules
#module load netcdf/4.1.2
#module load mpich2/1.5b1
#module list

#--------------------------------------------------------------------------------------------------------
# compile combine tool
cd $ppdir
icc -O -c -I/public/software/mathlib/netcdf/4.4.1/intel_parallel/include mppnccombine.c
#icc -O -c -I/usr/include/linux mppnccombine.c
if ( $status != 0 ) exit 1
icc -O -o mppnccombine.x -L/public/software/mathlib/netcdf/4.4.1/intel_parallel/lib -lnetcdf -lnetcdff mppnccombine.o
if ( $status != 0 ) exit 1
#--------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------
# setup directory structure
if ( ! -d $execdir ) mkdir -p $execdir
cd $execdir
#--------------------------------------------------------------------------------------------------------

# execute mkmf to create makefile
set cppDefs = "-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8"
$mkmf -a $sourcedir -t $template -p $executable:t -c "$cppDefs" $pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include
if ( $status != 0 ) then
   unset echo
   echo "ERROR: mkmf failed for full_dynamics" 
   exit 1
endif

# --- execute make ---
make $executable:t
if ( $status != 0 ) then
   unset echo
   echo "ERROR: make failed for full_dynamics" 
   exit 1
endif

unset echo
echo "NOTE: make successful for full_dynamics"
