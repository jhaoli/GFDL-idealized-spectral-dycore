#!/bin/csh -f
#Minimal runscript
set echo 
#--------------------------------------------------------------------------------------------------------
#set experiment = rossby_haurwitz_wave
set experiment = HSt42
# The available experiments are:
# HSt42
# t42_polvani_2004
# t42_polvani_2007_LC1
# t42_jablonowski_2006
# rossby_haurwitz_wave
# mountain_wave
#--------------------------------------------------------------------------------------------------------
# define variables
set platform  = taiyuan.intel                            # A unique identifier for your platform
set npes = 16                                                # Number of processors
set num_executions = 1                                      # Number of times the model is run. Each run restarts from previous run.
set time_stamp = $cwd/../bin/time_stamp.csh                 # Path to timestamp.csh
set model_executable = $cwd/exec.$platform/full_dynamics.x  # Path to model executable
set mpiexec = $cwd/../bin/mpiexec.hydra                     # Enables mpi on workstations
set mppnccombine = $cwd/../postprocessing/mppnccombine.x    # The tool for combining distributed diagnostic output files
set workdir = /public/home/wab/lijh/GFDL/test/$experiment  # Where model is run and model output is produced
#--------------------------------------------------------------------------------------------------------
#source $MODULESHOME/init/csh
#module use -a /home/fms/local/modulefiles
#module use /home/Seth.Underwood/publicmodules
#module purge
#module rm netcdf hdf5
#module load ifort/11.1.073
#module load icc/11.1.073
#module load hdf5/1.8.6
#module load netcdf/4.1.2
#module load mpich2/1.5b1
#module list
limit stacksize unlimited

set namelist   = $cwd/input/$experiment/input.nml    # path to namelist file (contains all namelists)
set diagtable  = $cwd/input/$experiment/diag_table   # path to diagnositics table (specifies fields and files for diagnostic output)
set fieldtable = $cwd/input/$experiment/field_table  # path to field table (specifies tracers)
#--------------------------------------------------------------------------------------------------------

# setup directory structure
if ( -d $workdir ) then
  /bin/rm -rf $workdir/*
else
  mkdir -p $workdir
endif
cd $workdir
mkdir RESTART
#--------------------------------------------------------------------------------------------------------
# get input data and executable
cp $namelist   input.nml
cp $diagtable  diag_table
cp $fieldtable field_table
cp $model_executable .
#touch data_table # Is this still necessary?
if ( $experiment == t42_polvani_2004 ) then
  cp $namelist:h/lmp_gqntc_t42l20_ess_pert.nc INPUT
endif

set irun = 1
while ( $irun <= $num_executions )
#--------------------------------------------------------------------------------------------------------
# run the model
$mpiexec -np $npes ./$model_executable:t
if ($status != 0) then
  echo "Error in execution of $cwd/$model_executable:t"
  exit 1
endif
#--------------------------------------------------------------------------------------------------------
set date_name = `$time_stamp -bf digital`
foreach outfile ( *.out )
  mv $outfile $date_name.$outfile
end
#--------------------------------------------------------------------------------------------------------
# combine diagnostic files, then remove the uncombined files.
if ( $npes > 1 ) then
  foreach ncfile (`/bin/ls *.nc.0000`)
    $mppnccombine $ncfile:r
    if ($status == 0) then
      rm -f $ncfile:r.[0-9][0-9][0-9][0-9]
      mv $ncfile:r $date_name.$ncfile:r
    else
      echo "Error in execution of $mppnccombine while working on $ncfile:r"
      exit 1
    endif
  end
endif
#--------------------------------------------------------------------------------------------------------
# Prepare to run the model again
cd $workdir
/bin/rm INPUT/*.res   INPUT/*.res.nc   INPUT/*.res.nc.0???   INPUT/*.res.tile?.nc   INPUT/*.res.tile?.nc.0???
mv    RESTART/*.res RESTART/*.res.nc RESTART/*.res.nc.0??? RESTART/*.res.tile?.nc RESTART/*.res.tile?.nc.0??? INPUT
#--------------------------------------------------------------------------------------------------------
@ irun ++
end
echo "NOTE: $experiment completed successfully"
exit 0
