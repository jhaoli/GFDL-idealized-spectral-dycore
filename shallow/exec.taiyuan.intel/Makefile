# Makefile created by mkmf $Id: mkmf,v 16.1 2010/05/19 18:49:19 fms Exp $ 

SRCROOT = /public/home/wab/lijh/GFDL/shallow/../src/

CPPDEFS = -Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DOVERLOAD_C8


include /public/home/wab/lijh/GFDL/shallow/../bin/mkmf.template.taiyuan.intel


.DEFAULT:
	-echo $@ does not exist.
all: shallow.x
affinity.o: $(SRCROOT)shared/mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/affinity.c
atmos_model.o: $(SRCROOT)atmos_solo/atmos_model.F90 mpp.o fms.o atmosphere.o time_manager.o fms_io.o mpp_domains.o mpp_io.o diag_manager.o field_manager.o tracer_manager.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_solo/atmos_model.F90
atmosphere.o: $(SRCROOT)atmos_spectral_shallow/atmosphere.F90 fms.o mpp.o constants.o transforms.o time_manager.o shallow_dynamics.o shallow_dynamics_type.o shallow_physics.o shallow_diagnostics.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral_shallow/atmosphere.F90
constants.o: $(SRCROOT)shared/constants/constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/constants/constants.F90
create_xgrid.o: $(SRCROOT)shared/mosaic/create_xgrid.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/create_xgrid.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic	$(SRCROOT)shared/mosaic/create_xgrid.c
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.F90 mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)shared/diag_manager/diag_grid.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/diag_manager/diag_grid.F90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h time_manager.o mpp_io.o mpp.o fms.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o mpp_domains.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.F90 mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)shared/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_util.F90
fft.o: $(SRCROOT)shared/fft/fft.F90 platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fft/fft.F90
fft99.o: $(SRCROOT)shared/fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fft/fft99.F90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)shared/field_manager/fm_util.F90 field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/field_manager/fm_util.F90
fms.o: $(SRCROOT)shared/fms/fms.F90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fms/fms.F90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_2d.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_3d.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_4d.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include -I/public/home/wab/lijh/GFDL/shallow/../src/shared/fms	$(SRCROOT)shared/fms/fms_io.F90
fv_advection.o: $(SRCROOT)atmos_spectral/model/fv_advection.F90 fms.o constants.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/fv_advection.F90
gauss_and_legendre.o: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.F90 fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/get_cal_time.F90
gradient.o: $(SRCROOT)shared/mosaic/gradient.F90 mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)shared/mosaic/gradient_c2l.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic	$(SRCROOT)shared/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)shared/mosaic/grid.F90 mpp.o constants.o fms.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/grid.F90
grid_fourier.o: $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 fms.o constants.o fft.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/grid_fourier.F90
interp.o: $(SRCROOT)shared/mosaic/interp.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/interp.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic	$(SRCROOT)shared/mosaic/interp.c
leapfrog.o: $(SRCROOT)atmos_spectral/model/leapfrog.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/leapfrog.F90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
mosaic.o: $(SRCROOT)shared/mosaic/mosaic.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/mosaic.F90
mosaic_util.o: $(SRCROOT)shared/mosaic/mosaic_util.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic	$(SRCROOT)shared/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/system_clock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_mpi.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_nocomm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_a.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_s.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_a.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_s.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_sma.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_sma.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_sma.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_mpi.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_mpi.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_mpi.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_mpi.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_nocomm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_int.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_scalar.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_gather.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_mpi.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_define_nest_domains.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_comm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_define.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_misc.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_check.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_checkV.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_nest_domains.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nest.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_redistribute.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_get_boundary.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_get_boundary.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_reduce.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_reduce.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum_tl.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_field.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_global_field.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_misc.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_connect.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_read.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_read_2Ddecomp.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_write.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write_2Ddecomp.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
mpp_memutils.o: $(SRCROOT)shared/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)shared/mpp/mpp_utilities.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_utilities.F90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
platform.o: $(SRCROOT)shared/platform/platform.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/platform/platform.F90
read_mosaic.o: $(SRCROOT)shared/mosaic/read_mosaic.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/read_mosaic.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic	$(SRCROOT)shared/mosaic/read_mosaic.c
shallow_diagnostics.o: $(SRCROOT)atmos_spectral_shallow/shallow_diagnostics.F90 fms.o transforms.o diag_manager.o time_manager.o shallow_physics.o shallow_dynamics_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral_shallow/shallow_diagnostics.F90
shallow_dynamics.o: $(SRCROOT)atmos_spectral_shallow/shallow_dynamics.F90 fms.o time_manager.o constants.o transforms.o spectral_damping.o leapfrog.o fv_advection.o shallow_dynamics_type.o williamson92_test_cases.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral_shallow/shallow_dynamics.F90
shallow_dynamics_type.o: $(SRCROOT)atmos_spectral_shallow/shallow_dynamics_type.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral_shallow/shallow_dynamics_type.F90
shallow_physics.o: $(SRCROOT)atmos_spectral_shallow/shallow_physics.F90 fms.o transforms.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral_shallow/shallow_physics.F90
spec_mpp.o: $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 fms.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spec_mpp.F90
spectral_damping.o: $(SRCROOT)atmos_spectral/model/spectral_damping.F90 fms.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/spectral_damping.F90
spherical.o: $(SRCROOT)atmos_spectral/tools/spherical.F90 fms.o spec_mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical.F90
spherical_fourier.o: $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 fms.o mpp.o mpp_domains.o constants.o spec_mpp.o spherical.o gauss_and_legendre.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical_fourier.F90
test_fms_io.o: $(SRCROOT)shared/fms/test_fms_io.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/fms/test_fms_io.F90
test_mpp.o: $(SRCROOT)shared/mpp/test_mpp.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp.F90
test_mpp_domains.o: $(SRCROOT)shared/mpp/test_mpp_domains.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_domains.F90
test_mpp_io.o: $(SRCROOT)shared/mpp/test_mpp_io.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_io.F90
test_mpp_pset.o: $(SRCROOT)shared/mpp/test_mpp_pset.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/shallow/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_pset.F90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/time_manager.F90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
transforms.o: $(SRCROOT)atmos_spectral/tools/transforms.F90 mpp.o fms.o mpp_domains.o spec_mpp.o constants.o spherical_fourier.o grid_fourier.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/transforms.F90
williamson92_test_cases.o: $(SRCROOT)atmos_spectral_shallow/init/williamson92_test_cases.F90 constants.o shallow_dynamics_type.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral_shallow/init/williamson92_test_cases.F90
./fms_platform.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h .
./create_xgrid.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/create_xgrid.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/create_xgrid.h .
./mpp_do_check.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_check.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_check.h .
./mpp_do_update_nonblock.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nonblock.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nonblock.h .
./mosaic_util.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./mpp_write.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write.h .
./test_mpp_pset.F90: $(SRCROOT)shared/mpp/test_mpp_pset.F90
	cp $(SRCROOT)shared/mpp/test_mpp_pset.F90 .
./mpp_reduce_sma.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_sma.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_sma.h .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./mpp_do_updateV.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV.h .
./mpp_do_updateV_nonblock.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV_nonblock.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV_nonblock.h .
./mpp_transmit_mpi.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_mpi.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_mpi.h .
./mpp_do_get_boundary.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_get_boundary.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_get_boundary.h .
./mpp_io_read.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_read.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_read.inc .
./mpp_do_redistribute.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_redistribute.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_redistribute.h .
./mpp_chksum_int.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_int.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_int.h .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./interp.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/interp.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/interp.h .
./mosaic.F90: $(SRCROOT)shared/mosaic/mosaic.F90
	cp $(SRCROOT)shared/mosaic/mosaic.F90 .
./gradient_c2l.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/gradient_c2l.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/gradient_c2l.h .
./mpp_data_nocomm.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_nocomm.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_nocomm.inc .
./mpp_io_misc.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_misc.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_misc.inc .
./mpp_do_update.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update.h .
./fft.F90: $(SRCROOT)shared/fft/fft.F90
	cp $(SRCROOT)shared/fft/fft.F90 .
./mpp_util_sma.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_sma.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_sma.inc .
./mpp_data_mpi.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_mpi.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_mpi.inc .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./mpp_transmit.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit.inc .
./mpp_sum.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum.inc .
./mpp_global_field.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_field.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_field.h .
./gradient.F90: $(SRCROOT)shared/mosaic/gradient.F90
	cp $(SRCROOT)shared/mosaic/gradient.F90 .
./mpp_error_s_s.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_s.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_s.h .
./transforms.F90: $(SRCROOT)atmos_spectral/tools/transforms.F90
	cp $(SRCROOT)atmos_spectral/tools/transforms.F90 .
./mpp_comm_nocomm.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_nocomm.inc .
./mpp_utilities.F90: $(SRCROOT)shared/mpp/mpp_utilities.F90
	cp $(SRCROOT)shared/mpp/mpp_utilities.F90 .
./mpp_update_domains2D.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D.h .
./read_mosaic.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/read_mosaic.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/read_mosaic.h .
./mpp_do_global_field.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_global_field.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_global_field.h .
./mpp_get_boundary.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_get_boundary.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_get_boundary.h .
./mpp_io_util.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_util.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_util.inc .
./mpp_define_nest_domains.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_define_nest_domains.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_define_nest_domains.inc .
./system_clock.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/system_clock.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/system_clock.h .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./mpp_transmit_sma.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_sma.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_sma.h .
./shallow_dynamics_type.F90: $(SRCROOT)atmos_spectral_shallow/shallow_dynamics_type.F90
	cp $(SRCROOT)atmos_spectral_shallow/shallow_dynamics_type.F90 .
./test_mpp_domains.F90: $(SRCROOT)shared/mpp/test_mpp_domains.F90
	cp $(SRCROOT)shared/mpp/test_mpp_domains.F90 .
./spherical_fourier.F90: $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90
	cp $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 .
./mpp_domains_misc.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_misc.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_misc.inc .
./mpp_comm_mpi.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_mpi.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_mpi.inc .
./diag_axis.F90: $(SRCROOT)shared/diag_manager/diag_axis.F90
	cp $(SRCROOT)shared/diag_manager/diag_axis.F90 .
./mpp_comm_sma.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_sma.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_sma.inc .
./mpp_read_2Ddecomp.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_read_2Ddecomp.h .
./parse.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/field_manager/parse.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/field_manager/parse.inc .
./mpp_util_nocomm.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_nocomm.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_nocomm.inc .
./fft99.F90: $(SRCROOT)shared/fft/fft99.F90
	cp $(SRCROOT)shared/fft/fft99.F90 .
./mpp_reduce_mpi.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_mpi.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_mpi.h .
./williamson92_test_cases.F90: $(SRCROOT)atmos_spectral_shallow/init/williamson92_test_cases.F90
	cp $(SRCROOT)atmos_spectral_shallow/init/williamson92_test_cases.F90 .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./diag_manager.F90: $(SRCROOT)shared/diag_manager/diag_manager.F90
	cp $(SRCROOT)shared/diag_manager/diag_manager.F90 .
./fms.F90: $(SRCROOT)shared/fms/fms.F90
	cp $(SRCROOT)shared/fms/fms.F90 .
./read_data_3d.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_3d.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_3d.inc .
./mpp_update_domains2D_nonblock.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h .
./leapfrog.F90: $(SRCROOT)atmos_spectral/model/leapfrog.F90
	cp $(SRCROOT)atmos_spectral/model/leapfrog.F90 .
./mpp_global_sum_tl.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum_tl.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum_tl.h .
./mpp_io_connect.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_connect.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_connect.inc .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./fm_util.F90: $(SRCROOT)shared/field_manager/fm_util.F90
	cp $(SRCROOT)shared/field_manager/fm_util.F90 .
./mpp_sum_nocomm.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_nocomm.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_nocomm.h .
./write_data.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/write_data.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/write_data.inc .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./mpp_error_s_a.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_a.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_a.h .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./shallow_physics.F90: $(SRCROOT)atmos_spectral_shallow/shallow_physics.F90
	cp $(SRCROOT)atmos_spectral_shallow/shallow_physics.F90 .
./mpp_util.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util.inc .
./mpp_domains_reduce.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_reduce.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_reduce.inc .
./mpp_comm.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm.inc .
./mpp_update_nest_domains.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_nest_domains.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_nest_domains.h .
./diag_table.F90: $(SRCROOT)shared/diag_manager/diag_table.F90
	cp $(SRCROOT)shared/diag_manager/diag_table.F90 .
./mpp_reduce_nocomm.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_nocomm.h .
./mpp_error_a_a.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_a.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_a.h .
./read_data_2d.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_2d.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_2d.inc .
./spherical.F90: $(SRCROOT)atmos_spectral/tools/spherical.F90
	cp $(SRCROOT)atmos_spectral/tools/spherical.F90 .
./mpp_write_2Ddecomp.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write_2Ddecomp.h .
./gradient_c2l.c: $(SRCROOT)shared/mosaic/gradient_c2l.c
	cp $(SRCROOT)shared/mosaic/gradient_c2l.c .
./test_mpp_io.F90: $(SRCROOT)shared/mpp/test_mpp_io.F90
	cp $(SRCROOT)shared/mpp/test_mpp_io.F90 .
./diag_grid.F90: $(SRCROOT)shared/diag_manager/diag_grid.F90
	cp $(SRCROOT)shared/diag_manager/diag_grid.F90 .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./shallow_diagnostics.F90: $(SRCROOT)atmos_spectral_shallow/shallow_diagnostics.F90
	cp $(SRCROOT)atmos_spectral_shallow/shallow_diagnostics.F90 .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./mpp_io_write.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_write.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_write.inc .
./spec_mpp.F90: $(SRCROOT)atmos_spectral/tools/spec_mpp.F90
	cp $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 .
./mpp_memutils.F90: $(SRCROOT)shared/mpp/mpp_memutils.F90
	cp $(SRCROOT)shared/mpp/mpp_memutils.F90 .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./read_data_4d.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_4d.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_4d.inc .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./test_fms_io.F90: $(SRCROOT)shared/fms/test_fms_io.F90
	cp $(SRCROOT)shared/fms/test_fms_io.F90 .
./mpp_chksum_scalar.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_scalar.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_scalar.h .
./gauss_and_legendre.F90: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90
	cp $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 .
./constants.F90: $(SRCROOT)shared/constants/constants.F90
	cp $(SRCROOT)shared/constants/constants.F90 .
./mpp_do_checkV.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_checkV.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_checkV.h .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./test_mpp.F90: $(SRCROOT)shared/mpp/test_mpp.F90
	cp $(SRCROOT)shared/mpp/test_mpp.F90 .
./mpp_sum_sma.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_sma.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_sma.h .
./mpp_do_update_nest.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nest.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nest.h .
./mpp_error_a_s.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_s.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_s.h .
./affinity.c: $(SRCROOT)shared/mpp/affinity.c
	cp $(SRCROOT)shared/mpp/affinity.c .
./mosaic_util.c: $(SRCROOT)shared/mosaic/mosaic_util.c
	cp $(SRCROOT)shared/mosaic/mosaic_util.c .
./shallow_dynamics.F90: $(SRCROOT)atmos_spectral_shallow/shallow_dynamics.F90
	cp $(SRCROOT)atmos_spectral_shallow/shallow_dynamics.F90 .
./read_mosaic.c: $(SRCROOT)shared/mosaic/read_mosaic.c
	cp $(SRCROOT)shared/mosaic/read_mosaic.c .
./interp.c: $(SRCROOT)shared/mosaic/interp.c
	cp $(SRCROOT)shared/mosaic/interp.c .
./mpp_global_reduce.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_reduce.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_reduce.h .
./grid_fourier.F90: $(SRCROOT)atmos_spectral/tools/grid_fourier.F90
	cp $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 .
./mpp_gather.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_gather.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_gather.h .
./spectral_damping.F90: $(SRCROOT)atmos_spectral/model/spectral_damping.F90
	cp $(SRCROOT)atmos_spectral/model/spectral_damping.F90 .
./mpp_chksum.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum.h .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./mpp_domains_define.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_define.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_define.inc .
./mpp_domains_util.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_util.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_util.inc .
./mpp_data_sma.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_sma.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_sma.inc .
./mpp_util_mpi.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_mpi.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_mpi.inc .
./get_cal_time.F90: $(SRCROOT)shared/time_manager/get_cal_time.F90
	cp $(SRCROOT)shared/time_manager/get_cal_time.F90 .
./mpp_transmit_nocomm.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_nocomm.h .
./diag_output.F90: $(SRCROOT)shared/diag_manager/diag_output.F90
	cp $(SRCROOT)shared/diag_manager/diag_output.F90 .
./atmosphere.F90: $(SRCROOT)atmos_spectral_shallow/atmosphere.F90
	cp $(SRCROOT)atmos_spectral_shallow/atmosphere.F90 .
./mpp_domains_comm.inc: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_comm.inc
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_comm.inc .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./mpp_sum_mpi.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_mpi.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_mpi.h .
./grid.F90: $(SRCROOT)shared/mosaic/grid.F90
	cp $(SRCROOT)shared/mosaic/grid.F90 .
./constant.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h .
./fv_advection.F90: $(SRCROOT)atmos_spectral/model/fv_advection.F90
	cp $(SRCROOT)atmos_spectral/model/fv_advection.F90 .
./atmos_model.F90: $(SRCROOT)atmos_solo/atmos_model.F90
	cp $(SRCROOT)atmos_solo/atmos_model.F90 .
./mpp_global_sum.h: /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum.h
	cp /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum.h .
./create_xgrid.c: $(SRCROOT)shared/mosaic/create_xgrid.c
	cp $(SRCROOT)shared/mosaic/create_xgrid.c .
SRC = $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)atmos_spectral_shallow/shallow_diagnostics.F90 $(SRCROOT)atmos_solo/atmos_model.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)atmos_spectral_shallow/atmosphere.F90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)atmos_spectral_shallow/init/williamson92_test_cases.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)atmos_spectral/model/spectral_damping.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)atmos_spectral_shallow/shallow_dynamics_type.F90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)atmos_spectral_shallow/shallow_physics.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)atmos_spectral/tools/transforms.F90 $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)atmos_spectral/tools/spherical.F90 $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)atmos_spectral/model/leapfrog.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)atmos_spectral/model/fv_advection.F90 $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)atmos_spectral_shallow/shallow_dynamics.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/diag_manager/diag_table.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/create_xgrid.h /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_3d.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum_tl.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_check.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_connect.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/write_data.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_a.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_sma.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_reduce.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_mpi.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_get_boundary.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_read.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_redistribute.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_nest_domains.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_int.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_a.h /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_2d.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/interp.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write_2Ddecomp.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/gradient_c2l.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_nocomm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_misc.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_write.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_4d.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_scalar.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_checkV.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_mpi.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_field.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_sma.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nest.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_s.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_s.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_nocomm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_reduce.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_gather.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/read_mosaic.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_global_field.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_get_boundary.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_define.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/system_clock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_define_nest_domains.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_sma.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_mpi.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_misc.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_mpi.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_read_2Ddecomp.h /public/home/wab/lijh/GFDL/shallow/../src/shared/field_manager/parse.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_nocomm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_comm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_mpi.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_mpi.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum.h
OBJ = time_manager.o spherical_fourier.o mpp_data.o affinity.o gradient_c2l.o shallow_diagnostics.o atmos_model.o mpp_io.o atmosphere.o fms_io.o tracer_manager.o mpp_memutils.o williamson92_test_cases.o diag_util.o mosaic_util.o spectral_damping.o gradient.o shallow_dynamics_type.o fm_util.o fft.o shallow_physics.o mosaic.o memuse.o memutils.o threadloc.o diag_data.o grid.o diag_axis.o get_cal_time.o read_mosaic.o diag_manager.o transforms.o mpp_utilities.o spherical.o constants.o mpp_pset.o fms.o test_mpp_pset.o test_mpp_io.o spec_mpp.o test_fms_io.o create_xgrid.o grid_fourier.o field_manager.o interp.o diag_output.o gauss_and_legendre.o nsclock.o leapfrog.o mpp_parameter.o platform.o mpp.o fv_advection.o fft99.o mpp_domains.o test_mpp_domains.o test_mpp.o shallow_dynamics.o diag_grid.o diag_table.o
OFF = /public/home/wab/lijh/GFDL/shallow/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/create_xgrid.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_check.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/mosaic_util.h $(SRCROOT)shared/mpp/mpp_pset.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write.h $(SRCROOT)shared/mpp/test_mpp_pset.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_sma.h $(SRCROOT)shared/fms/fms_io.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_mpi.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_get_boundary.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_read.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_redistribute.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_int.h $(SRCROOT)shared/platform/platform.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/interp.h $(SRCROOT)shared/mosaic/mosaic.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/gradient_c2l.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_nocomm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_misc.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update.h $(SRCROOT)shared/fft/fft.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_mpi.inc $(SRCROOT)shared/mpp/mpp_data.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_field.h $(SRCROOT)shared/mosaic/gradient.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_s.h $(SRCROOT)atmos_spectral/tools/transforms.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)shared/mpp/mpp_utilities.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/read_mosaic.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_global_field.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_get_boundary.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_define_nest_domains.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/system_clock.h $(SRCROOT)shared/time_manager/time_manager.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_sma.h $(SRCROOT)atmos_spectral_shallow/shallow_dynamics_type.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_misc.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_mpi.inc $(SRCROOT)shared/diag_manager/diag_axis.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_read_2Ddecomp.h /public/home/wab/lijh/GFDL/shallow/../src/shared/field_manager/parse.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_nocomm.inc $(SRCROOT)shared/fft/fft99.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)atmos_spectral_shallow/init/williamson92_test_cases.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)shared/fms/fms.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_3d.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)atmos_spectral/model/leapfrog.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum_tl.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_connect.inc $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/field_manager/fm_util.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/write_data.inc $(SRCROOT)shared/mpp/mpp_io.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_s_a.h $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)atmos_spectral_shallow/shallow_physics.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_reduce.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_comm.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_update_nest_domains.h $(SRCROOT)shared/diag_manager/diag_table.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_reduce_nocomm.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_a.h /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_2d.inc $(SRCROOT)atmos_spectral/tools/spherical.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)atmos_spectral_shallow/shallow_diagnostics.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_io_write.inc $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/fms/read_data_4d.inc $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/fms/test_fms_io.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum_scalar.h $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 $(SRCROOT)shared/constants/constants.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_checkV.h $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/mpp/test_mpp.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_sma.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_do_update_nest.h /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_error_a_s.h $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)atmos_spectral_shallow/shallow_dynamics.F90 $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/mosaic/interp.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_reduce.h $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_gather.h $(SRCROOT)atmos_spectral/model/spectral_damping.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_chksum.h $(SRCROOT)shared/mpp/threadloc.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_define.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_util.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_data_sma.inc /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_util_mpi.inc $(SRCROOT)shared/time_manager/get_cal_time.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)atmos_spectral_shallow/atmosphere.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_domains_comm.inc $(SRCROOT)shared/mpp/nsclock.c /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)shared/mosaic/grid.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mosaic/constant.h $(SRCROOT)atmos_spectral/model/fv_advection.F90 $(SRCROOT)atmos_solo/atmos_model.F90 /public/home/wab/lijh/GFDL/shallow/../src/shared/mpp/include/mpp_global_sum.h $(SRCROOT)shared/mosaic/create_xgrid.c
clean: neat
	-rm -f .shallow.x.cppdefs $(OBJ) shallow.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
shallow.x: $(OBJ) 
	$(LD) $(OBJ) -o shallow.x  $(LDFLAGS)
