# Makefile created by mkmf $Id: mkmf,v 16.1 2010/05/19 18:49:19 fms Exp $ 

SRCROOT = /public/home/wab/lijh/GFDL/full_dynamics/../src/

CPPDEFS = -Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DINTERNAL_FILE_NML -DOVERLOAD_C8


include /public/home/wab/lijh/GFDL/full_dynamics/../bin/mkmf.template.taiyuan.intel


.DEFAULT:
	-echo $@ does not exist.
all: full_dynamics.x
affinity.o: $(SRCROOT)shared/mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/affinity.c
atmos_model.o: $(SRCROOT)atmos_solo/atmos_model.F90 mpp.o fms.o atmosphere.o time_manager.o fms_io.o mpp_domains.o mpp_io.o diag_manager.o field_manager.o tracer_manager.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_solo/atmos_model.F90
atmosphere.o: $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 fms.o constants.o transforms.o time_manager.o press_and_geopot.o spectral_dynamics.o tracer_type.o hs_forcing.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90
constants.o: $(SRCROOT)shared/constants/constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/constants/constants.F90
create_xgrid.o: $(SRCROOT)shared/mosaic/create_xgrid.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/create_xgrid.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic	$(SRCROOT)shared/mosaic/create_xgrid.c
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.F90 mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)shared/diag_manager/diag_grid.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/diag_manager/diag_grid.F90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h time_manager.o mpp_io.o mpp.o fms.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o mpp_domains.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.F90 mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)shared/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_util.F90
every_step_diagnostics.o: $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 fms.o time_manager.o field_manager.o tracer_manager.o diag_manager.o press_and_geopot.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90
fft.o: $(SRCROOT)shared/fft/fft.F90 platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fft/fft.F90
fft99.o: $(SRCROOT)shared/fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fft/fft99.F90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)shared/field_manager/fm_util.F90 field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/field_manager/fm_util.F90
fms.o: $(SRCROOT)shared/fms/fms.F90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fms/fms.F90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_2d.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_3d.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_4d.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms	$(SRCROOT)shared/fms/fms_io.F90
fv_advection.o: $(SRCROOT)atmos_spectral/model/fv_advection.F90 fms.o constants.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/fv_advection.F90
gauss_and_legendre.o: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90
gaussian_topog.o: $(SRCROOT)shared/topography/gaussian_topog.F90 fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/topography/gaussian_topog.F90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.F90 fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/get_cal_time.F90
global_integral.o: $(SRCROOT)atmos_spectral/model/global_integral.F90 fms.o press_and_geopot.o transforms.o constants.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/global_integral.F90
gradient.o: $(SRCROOT)shared/mosaic/gradient.F90 mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)shared/mosaic/gradient_c2l.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic	$(SRCROOT)shared/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)shared/mosaic/grid.F90 mpp.o constants.o fms.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/grid.F90
grid_fourier.o: $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 fms.o constants.o fft.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/grid_fourier.F90
horiz_interp.o: $(SRCROOT)shared/horiz_interp/horiz_interp.F90 fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
hs_forcing.o: $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 mpp.o fms.o constants.o time_manager.o diag_manager.o field_manager.o tracer_manager.o interpolator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90
ic_from_external_file.o: $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 mpp.o fms.o transforms.o tracer_type.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/ic_from_external_file.F90
implicit.o: $(SRCROOT)atmos_spectral/model/implicit.F90 fms.o constants.o press_and_geopot.o matrix_invert.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/implicit.F90
interp.o: $(SRCROOT)shared/mosaic/interp.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/interp.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic	$(SRCROOT)shared/mosaic/interp.c
interpolator.o: $(SRCROOT)atmos_shared/interpolator/interpolator.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h mpp.o mpp_io.o mpp_domains.o diag_manager.o fms.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)atmos_shared/interpolator/interpolator.F90
jablonowski_2006.o: $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 mpp.o fms.o constants.o vert_coordinate.o transforms.o press_and_geopot.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/jablonowski_2006.F90
leapfrog.o: $(SRCROOT)atmos_spectral/model/leapfrog.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/leapfrog.F90
matrix_invert.o: $(SRCROOT)atmos_spectral/model/matrix_invert.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/matrix_invert.F90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
mosaic.o: $(SRCROOT)shared/mosaic/mosaic.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/mosaic.F90
mosaic_util.o: $(SRCROOT)shared/mosaic/mosaic_util.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic	$(SRCROOT)shared/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/system_clock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_nocomm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_a.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_s.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_a.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_s.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_sma.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_sma.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_sma.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_mpi.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_mpi.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_mpi.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_nocomm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_nocomm.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_nocomm.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_nocomm.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_int.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_scalar.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_gather.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_define_nest_domains.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_comm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_define.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_misc.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nonblock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_check.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_checkV.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_nest_domains.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nest.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_redistribute.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_get_boundary.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_get_boundary.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_reduce.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_reduce.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum_tl.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_field.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_global_field.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_misc.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_connect.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_read.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_read_2Ddecomp.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_write.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write_2Ddecomp.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
mpp_memutils.o: $(SRCROOT)shared/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)shared/mpp/mpp_utilities.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_utilities.F90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
platform.o: $(SRCROOT)shared/platform/platform.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/platform/platform.F90
polvani_2004.o: $(SRCROOT)atmos_spectral/init/polvani_2004.F90 mpp.o fms.o constants.o vert_coordinate.o transforms.o press_and_geopot.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/polvani_2004.F90
polvani_2007.o: $(SRCROOT)atmos_spectral/init/polvani_2007.F90 mpp.o fms.o constants.o vert_coordinate.o transforms.o press_and_geopot.o diag_manager.o tracer_type.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/polvani_2007.F90
press_and_geopot.o: $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/press_and_geopot.F90
read_mosaic.o: $(SRCROOT)shared/mosaic/read_mosaic.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/read_mosaic.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic	$(SRCROOT)shared/mosaic/read_mosaic.c
spec_mpp.o: $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 fms.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spec_mpp.F90
spectral_damping.o: $(SRCROOT)atmos_spectral/model/spectral_damping.F90 fms.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/spectral_damping.F90
spectral_dynamics.o: $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 mpp.o fms.o constants.o time_manager.o field_manager.o tracer_manager.o diag_manager.o transforms.o vert_advection.o implicit.o press_and_geopot.o spectral_damping.o leapfrog.o fv_advection.o water_borrowing.o global_integral.o spectral_init_cond.o tracer_type.o every_step_diagnostics.o mpp_domains.o polvani_2004.o polvani_2007.o jablonowski_2006.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/spectral_dynamics.F90
spectral_init_cond.o: $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 mpp.o fms.o mpp_domains.o constants.o vert_coordinate.o transforms.o press_and_geopot.o spectral_initialize_fields.o ic_from_external_file.o topog_regularization.o topography.o tracer_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/spectral_init_cond.F90
spectral_initialize_fields.o: $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 fms.o constants.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90
spherical.o: $(SRCROOT)atmos_spectral/tools/spherical.F90 fms.o spec_mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical.F90
spherical_fourier.o: $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 fms.o mpp.o mpp_domains.o constants.o spec_mpp.o spherical.o gauss_and_legendre.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/spherical_fourier.F90
test_fms_io.o: $(SRCROOT)shared/fms/test_fms_io.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/fms/test_fms_io.F90
test_mpp.o: $(SRCROOT)shared/mpp/test_mpp.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp.F90
test_mpp_domains.o: $(SRCROOT)shared/mpp/test_mpp_domains.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_domains.F90
test_mpp_io.o: $(SRCROOT)shared/mpp/test_mpp_io.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_io.F90
test_mpp_pset.o: $(SRCROOT)shared/mpp/test_mpp_pset.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_pset.F90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
time_interp.o: $(SRCROOT)shared/time_interp/time_interp.F90 time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_interp/time_interp.F90
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/time_manager.F90
topog_regularization.o: $(SRCROOT)atmos_spectral/init/topog_regularization.F90 fms.o mpp.o constants.o transforms.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/topog_regularization.F90
topography.o: $(SRCROOT)shared/topography/topography.F90 gaussian_topog.o horiz_interp.o fms.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/topography/topography.F90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
tracer_type.o: $(SRCROOT)atmos_spectral/model/tracer_type.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/tracer_type.F90
transforms.o: $(SRCROOT)atmos_spectral/tools/transforms.F90 mpp.o fms.o mpp_domains.o spec_mpp.o constants.o spherical_fourier.o grid_fourier.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/tools/transforms.F90
vert_advection.o: $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_shared/vert_advection/vert_advection.F90
vert_coordinate.o: $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/init/vert_coordinate.F90
water_borrowing.o: $(SRCROOT)atmos_spectral/model/water_borrowing.F90 fms.o transforms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_spectral/model/water_borrowing.F90
./mpp_update_domains2D_nonblock.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h .
./spectral_initialize_fields.F90: $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90
	cp $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 .
./interpolator.F90: $(SRCROOT)atmos_shared/interpolator/interpolator.F90
	cp $(SRCROOT)atmos_shared/interpolator/interpolator.F90 .
./parse.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/field_manager/parse.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/field_manager/parse.inc .
./mpp_util_mpi.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_mpi.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_mpi.inc .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./test_mpp_pset.F90: $(SRCROOT)shared/mpp/test_mpp_pset.F90
	cp $(SRCROOT)shared/mpp/test_mpp_pset.F90 .
./polvani_2004.F90: $(SRCROOT)atmos_spectral/init/polvani_2004.F90
	cp $(SRCROOT)atmos_spectral/init/polvani_2004.F90 .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./mpp_do_update.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update.h .
./vert_coordinate.F90: $(SRCROOT)atmos_spectral/init/vert_coordinate.F90
	cp $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 .
./mpp_reduce_sma.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_sma.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_sma.h .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./mosaic.F90: $(SRCROOT)shared/mosaic/mosaic.F90
	cp $(SRCROOT)shared/mosaic/mosaic.F90 .
./mpp_error_a_s.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_s.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_s.h .
./gaussian_topog.F90: $(SRCROOT)shared/topography/gaussian_topog.F90
	cp $(SRCROOT)shared/topography/gaussian_topog.F90 .
./fft.F90: $(SRCROOT)shared/fft/fft.F90
	cp $(SRCROOT)shared/fft/fft.F90 .
./read_data_4d.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_4d.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_4d.inc .
./mpp_io_misc.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_misc.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_misc.inc .
./global_integral.F90: $(SRCROOT)atmos_spectral/model/global_integral.F90
	cp $(SRCROOT)atmos_spectral/model/global_integral.F90 .
./mpp_chksum_scalar.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_scalar.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_scalar.h .
./atmosphere.F90: $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90
	cp $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 .
./spectral_dynamics.F90: $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90
	cp $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 .
./mpp_sum.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum.inc .
./read_data_2d.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_2d.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_2d.inc .
./polvani_2007.F90: $(SRCROOT)atmos_spectral/init/polvani_2007.F90
	cp $(SRCROOT)atmos_spectral/init/polvani_2007.F90 .
./mpp_domains_define.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_define.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_define.inc .
./matrix_invert.F90: $(SRCROOT)atmos_spectral/model/matrix_invert.F90
	cp $(SRCROOT)atmos_spectral/model/matrix_invert.F90 .
./mpp_io_write.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_write.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_write.inc .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./mpp_sum_nocomm.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_nocomm.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_nocomm.h .
./gradient.F90: $(SRCROOT)shared/mosaic/gradient.F90
	cp $(SRCROOT)shared/mosaic/gradient.F90 .
./horiz_interp_type.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 .
./mpp_comm_sma.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_sma.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_sma.inc .
./transforms.F90: $(SRCROOT)atmos_spectral/tools/transforms.F90
	cp $(SRCROOT)atmos_spectral/tools/transforms.F90 .
./read_mosaic.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/read_mosaic.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/read_mosaic.h .
./mosaic_util.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h .
./mpp_utilities.F90: $(SRCROOT)shared/mpp/mpp_utilities.F90
	cp $(SRCROOT)shared/mpp/mpp_utilities.F90 .
./fms_platform.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h .
./mpp_do_update_nest.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nest.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nest.h .
./mpp_do_checkV.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_checkV.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_checkV.h .
./create_xgrid.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/create_xgrid.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/create_xgrid.h .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./mpp_write_2Ddecomp.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write_2Ddecomp.h .
./test_mpp_domains.F90: $(SRCROOT)shared/mpp/test_mpp_domains.F90
	cp $(SRCROOT)shared/mpp/test_mpp_domains.F90 .
./spherical_fourier.F90: $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90
	cp $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 .
./constant.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h .
./mpp_do_updateV_nonblock.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV_nonblock.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV_nonblock.h .
./mpp_chksum.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum.h .
./mpp_error_s_s.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_s.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_s.h .
./mpp_domains_comm.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_comm.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_comm.inc .
./diag_axis.F90: $(SRCROOT)shared/diag_manager/diag_axis.F90
	cp $(SRCROOT)shared/diag_manager/diag_axis.F90 .
./mpp_do_get_boundary.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_get_boundary.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_get_boundary.h .
./horiz_interp_conserve.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 .
./mpp_io_util.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_util.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_util.inc .
./fft99.F90: $(SRCROOT)shared/fft/fft99.F90
	cp $(SRCROOT)shared/fft/fft99.F90 .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./diag_manager.F90: $(SRCROOT)shared/diag_manager/diag_manager.F90
	cp $(SRCROOT)shared/diag_manager/diag_manager.F90 .
./mpp_sum_mpi.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_mpi.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_mpi.h .
./horiz_interp_bicubic.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 .
./water_borrowing.F90: $(SRCROOT)atmos_spectral/model/water_borrowing.F90
	cp $(SRCROOT)atmos_spectral/model/water_borrowing.F90 .
./fms.F90: $(SRCROOT)shared/fms/fms.F90
	cp $(SRCROOT)shared/fms/fms.F90 .
./mpp_transmit.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit.inc .
./leapfrog.F90: $(SRCROOT)atmos_spectral/model/leapfrog.F90
	cp $(SRCROOT)atmos_spectral/model/leapfrog.F90 .
./mpp_gather.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_gather.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_gather.h .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./hs_forcing.F90: $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90
	cp $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 .
./interp.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/interp.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/interp.h .
./system_clock.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/system_clock.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/system_clock.h .
./fm_util.F90: $(SRCROOT)shared/field_manager/fm_util.F90
	cp $(SRCROOT)shared/field_manager/fm_util.F90 .
./mpp_util_sma.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_sma.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_sma.inc .
./mpp_global_reduce.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_reduce.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_reduce.h .
./read_data_3d.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_3d.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_3d.inc .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./gradient_c2l.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/gradient_c2l.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/gradient_c2l.h .
./mpp_comm.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm.inc .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./mpp_transmit_nocomm.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_nocomm.h .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./mpp_io_read.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_read.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_read.inc .
./implicit.F90: $(SRCROOT)atmos_spectral/model/implicit.F90
	cp $(SRCROOT)atmos_spectral/model/implicit.F90 .
./mpp_do_check.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_check.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_check.h .
./diag_table.F90: $(SRCROOT)shared/diag_manager/diag_table.F90
	cp $(SRCROOT)shared/diag_manager/diag_table.F90 .
./spherical.F90: $(SRCROOT)atmos_spectral/tools/spherical.F90
	cp $(SRCROOT)atmos_spectral/tools/spherical.F90 .
./mpp_reduce_mpi.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_mpi.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_mpi.h .
./mpp_do_redistribute.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_redistribute.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_redistribute.h .
./tracer_type.F90: $(SRCROOT)atmos_spectral/model/tracer_type.F90
	cp $(SRCROOT)atmos_spectral/model/tracer_type.F90 .
./gradient_c2l.c: $(SRCROOT)shared/mosaic/gradient_c2l.c
	cp $(SRCROOT)shared/mosaic/gradient_c2l.c .
./mpp_do_update_nonblock.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nonblock.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nonblock.h .
./test_mpp_io.F90: $(SRCROOT)shared/mpp/test_mpp_io.F90
	cp $(SRCROOT)shared/mpp/test_mpp_io.F90 .
./diag_grid.F90: $(SRCROOT)shared/diag_manager/diag_grid.F90
	cp $(SRCROOT)shared/diag_manager/diag_grid.F90 .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./mpp_domains_reduce.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_reduce.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_reduce.inc .
./spec_mpp.F90: $(SRCROOT)atmos_spectral/tools/spec_mpp.F90
	cp $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./mpp_comm_nocomm.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_nocomm.inc .
./mpp_memutils.F90: $(SRCROOT)shared/mpp/mpp_memutils.F90
	cp $(SRCROOT)shared/mpp/mpp_memutils.F90 .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./mpp_write.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write.h .
./write_data.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/write_data.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/write_data.inc .
./ic_from_external_file.F90: $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90
	cp $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 .
./mpp_transmit_sma.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_sma.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_sma.h .
./test_fms_io.F90: $(SRCROOT)shared/fms/test_fms_io.F90
	cp $(SRCROOT)shared/fms/test_fms_io.F90 .
./horiz_interp_spherical.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./gauss_and_legendre.F90: $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90
	cp $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 .
./constants.F90: $(SRCROOT)shared/constants/constants.F90
	cp $(SRCROOT)shared/constants/constants.F90 .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./every_step_diagnostics.F90: $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90
	cp $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 .
./mpp_io_connect.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_connect.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_connect.inc .
./vert_advection.F90: $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90
	cp $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 .
./topography.F90: $(SRCROOT)shared/topography/topography.F90
	cp $(SRCROOT)shared/topography/topography.F90 .
./test_mpp.F90: $(SRCROOT)shared/mpp/test_mpp.F90
	cp $(SRCROOT)shared/mpp/test_mpp.F90 .
./mpp_sum_sma.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_sma.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_sma.h .
./mpp_domains_util.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_util.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_util.inc .
./mpp_do_updateV.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV.h .
./mpp_error_a_a.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_a.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_a.h .
./horiz_interp_bilinear.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 .
./affinity.c: $(SRCROOT)shared/mpp/affinity.c
	cp $(SRCROOT)shared/mpp/affinity.c .
./mosaic_util.c: $(SRCROOT)shared/mosaic/mosaic_util.c
	cp $(SRCROOT)shared/mosaic/mosaic_util.c .
./read_mosaic.c: $(SRCROOT)shared/mosaic/read_mosaic.c
	cp $(SRCROOT)shared/mosaic/read_mosaic.c .
./interp.c: $(SRCROOT)shared/mosaic/interp.c
	cp $(SRCROOT)shared/mosaic/interp.c .
./spectral_init_cond.F90: $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90
	cp $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 .
./grid_fourier.F90: $(SRCROOT)atmos_spectral/tools/grid_fourier.F90
	cp $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 .
./mpp_data_mpi.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_mpi.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_mpi.inc .
./mpp_update_domains2D.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D.h .
./mpp_global_field.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_field.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_field.h .
./topog_regularization.F90: $(SRCROOT)atmos_spectral/init/topog_regularization.F90
	cp $(SRCROOT)atmos_spectral/init/topog_regularization.F90 .
./spectral_damping.F90: $(SRCROOT)atmos_spectral/model/spectral_damping.F90
	cp $(SRCROOT)atmos_spectral/model/spectral_damping.F90 .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./mpp_data_nocomm.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_nocomm.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_nocomm.inc .
./mpp_transmit_mpi.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_mpi.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_mpi.h .
./get_cal_time.F90: $(SRCROOT)shared/time_manager/get_cal_time.F90
	cp $(SRCROOT)shared/time_manager/get_cal_time.F90 .
./mpp_util_nocomm.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_nocomm.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_nocomm.inc .
./mpp_error_s_a.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_a.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_a.h .
./mpp_comm_mpi.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_mpi.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_mpi.inc .
./mpp_update_nest_domains.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_nest_domains.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_nest_domains.h .
./mpp_global_sum.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum.h .
./mpp_read_2Ddecomp.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_read_2Ddecomp.h .
./time_interp.F90: $(SRCROOT)shared/time_interp/time_interp.F90
	cp $(SRCROOT)shared/time_interp/time_interp.F90 .
./press_and_geopot.F90: $(SRCROOT)atmos_spectral/model/press_and_geopot.F90
	cp $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 .
./horiz_interp.F90: $(SRCROOT)shared/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp.F90 .
./mpp_data_sma.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_sma.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_sma.inc .
./mpp_reduce_nocomm.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_nocomm.h .
./mpp_domains_misc.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_misc.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_misc.inc .
./diag_output.F90: $(SRCROOT)shared/diag_manager/diag_output.F90
	cp $(SRCROOT)shared/diag_manager/diag_output.F90 .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./jablonowski_2006.F90: $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90
	cp $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 .
./grid.F90: $(SRCROOT)shared/mosaic/grid.F90
	cp $(SRCROOT)shared/mosaic/grid.F90 .
./mpp_util.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util.inc .
./mpp_chksum_int.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_int.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_int.h .
./fv_advection.F90: $(SRCROOT)atmos_spectral/model/fv_advection.F90
	cp $(SRCROOT)atmos_spectral/model/fv_advection.F90 .
./mpp_define_nest_domains.inc: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_define_nest_domains.inc
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_define_nest_domains.inc .
./atmos_model.F90: $(SRCROOT)atmos_solo/atmos_model.F90
	cp $(SRCROOT)atmos_solo/atmos_model.F90 .
./mpp_do_global_field.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_global_field.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_global_field.h .
./mpp_get_boundary.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_get_boundary.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_get_boundary.h .
./create_xgrid.c: $(SRCROOT)shared/mosaic/create_xgrid.c
	cp $(SRCROOT)shared/mosaic/create_xgrid.c .
./mpp_global_sum_tl.h: /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum_tl.h
	cp /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum_tl.h .
SRC = $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)atmos_spectral/init/polvani_2004.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)shared/topography/gaussian_topog.F90 $(SRCROOT)atmos_spectral/model/implicit.F90 $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 $(SRCROOT)atmos_spectral/tools/transforms.F90 $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)atmos_spectral/tools/spherical.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)atmos_spectral/model/fv_advection.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)shared/diag_manager/diag_table.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)atmos_solo/atmos_model.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)atmos_spectral/init/polvani_2007.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)atmos_spectral/model/spectral_damping.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)atmos_spectral/model/matrix_invert.F90 $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)atmos_spectral/model/tracer_type.F90 $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)atmos_shared/interpolator/interpolator.F90 $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)atmos_spectral/model/leapfrog.F90 $(SRCROOT)atmos_spectral/model/global_integral.F90 $(SRCROOT)atmos_spectral/init/topog_regularization.F90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)atmos_spectral/model/water_borrowing.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_mpi.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_gather.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/interp.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/system_clock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_reduce.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/field_manager/parse.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_3d.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/gradient_c2l.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_nocomm.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_read.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_check.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_redistribute.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_sma.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_mpi.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nonblock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_s.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_reduce.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_4d.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_misc.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_nocomm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_scalar.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/write_data.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_sma.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_2d.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_define.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_connect.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_write.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_nocomm.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_sma.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_a.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/read_mosaic.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_field.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nest.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_checkV.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/create_xgrid.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write_2Ddecomp.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_nocomm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_mpi.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_a.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_nocomm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_s.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_comm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_read_2Ddecomp.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_nest_domains.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_get_boundary.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_nocomm.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_misc.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_int.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_define_nest_domains.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_global_field.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_get_boundary.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum_tl.h
OBJ = spectral_dynamics.o mpp_data.o gradient_c2l.o affinity.o atmosphere.o ic_from_external_file.o fms_io.o every_step_diagnostics.o tracer_manager.o mpp_memutils.o diag_util.o mosaic_util.o time_interp.o polvani_2004.o mosaic.o memuse.o jablonowski_2006.o memutils.o threadloc.o horiz_interp_bicubic.o gaussian_topog.o implicit.o diag_axis.o diag_manager.o vert_coordinate.o transforms.o mpp_utilities.o spherical.o fms.o mpp_pset.o test_mpp_io.o create_xgrid.o interp.o diag_output.o horiz_interp_conserve.o mpp_parameter.o fv_advection.o platform.o mpp_domains.o diag_table.o diag_grid.o time_manager.o spherical_fourier.o vert_advection.o atmos_model.o mpp_io.o polvani_2007.o horiz_interp_spherical.o topography.o spectral_damping.o horiz_interp.o gradient.o spectral_initialize_fields.o fm_util.o fft.o horiz_interp_bilinear.o diag_data.o grid.o matrix_invert.o hs_forcing.o get_cal_time.o read_mosaic.o constants.o tracer_type.o spec_mpp.o test_mpp_pset.o test_fms_io.o interpolator.o grid_fourier.o field_manager.o gauss_and_legendre.o nsclock.o leapfrog.o global_integral.o topog_regularization.o mpp.o fft99.o press_and_geopot.o horiz_interp_type.o water_borrowing.o test_mpp_domains.o test_mpp.o spectral_init_cond.o
OFF = /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)atmos_spectral/init/spectral_initialize_fields.F90 $(SRCROOT)atmos_shared/interpolator/interpolator.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/field_manager/parse.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_mpi.inc $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)atmos_spectral/init/polvani_2004.F90 $(SRCROOT)shared/fms/fms_io.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update.h $(SRCROOT)atmos_spectral/init/vert_coordinate.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_sma.h $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/mosaic/mosaic.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_s.h $(SRCROOT)shared/topography/gaussian_topog.F90 $(SRCROOT)shared/fft/fft.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_4d.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_misc.inc $(SRCROOT)atmos_spectral/model/global_integral.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_scalar.h $(SRCROOT)atmos_spectral/driver/solo/atmosphere.F90 $(SRCROOT)atmos_spectral/model/spectral_dynamics.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_2d.inc $(SRCROOT)atmos_spectral/init/polvani_2007.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_define.inc $(SRCROOT)atmos_spectral/model/matrix_invert.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_write.inc $(SRCROOT)shared/mpp/mpp_data.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_sma.inc $(SRCROOT)atmos_spectral/tools/transforms.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/read_mosaic.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/mosaic_util.h $(SRCROOT)shared/mpp/mpp_utilities.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/include/fms_platform.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nest.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_checkV.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/create_xgrid.h $(SRCROOT)shared/time_manager/time_manager.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)atmos_spectral/tools/spherical_fourier.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/constant.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_s.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_comm.inc $(SRCROOT)shared/diag_manager/diag_axis.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_get_boundary.h $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_util.inc $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)atmos_spectral/model/water_borrowing.F90 $(SRCROOT)shared/fms/fms.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit.inc $(SRCROOT)atmos_spectral/model/leapfrog.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_gather.h $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)atmos_param/hs_forcing/hs_forcing.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/interp.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/system_clock.h $(SRCROOT)shared/field_manager/fm_util.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_reduce.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/read_data_3d.inc $(SRCROOT)shared/mpp/mpp_io.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mosaic/gradient_c2l.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm.inc $(SRCROOT)shared/mpp/mpp.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)shared/memutils/memuse.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_read.inc $(SRCROOT)atmos_spectral/model/implicit.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_check.h $(SRCROOT)shared/diag_manager/diag_table.F90 $(SRCROOT)atmos_spectral/tools/spherical.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_mpi.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_redistribute.h $(SRCROOT)atmos_spectral/model/tracer_type.F90 $(SRCROOT)shared/mosaic/gradient_c2l.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_reduce.inc $(SRCROOT)atmos_spectral/tools/spec_mpp.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_write.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/fms/write_data.inc $(SRCROOT)atmos_spectral/init/ic_from_external_file.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_sma.h $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)atmos_spectral/tools/gauss_and_legendre.F90 $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)atmos_spectral/model/every_step_diagnostics.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_io_connect.inc $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)shared/mpp/test_mpp.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_sum_sma.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_updateV.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_a_a.h $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)atmos_spectral/init/spectral_init_cond.F90 $(SRCROOT)atmos_spectral/tools/grid_fourier.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_domains2D.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_field.h $(SRCROOT)atmos_spectral/init/topog_regularization.F90 $(SRCROOT)atmos_spectral/model/spectral_damping.F90 $(SRCROOT)shared/mpp/threadloc.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_nocomm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_transmit_mpi.h $(SRCROOT)shared/time_manager/get_cal_time.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util_nocomm.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_error_s_a.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_comm_mpi.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_update_nest_domains.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)atmos_spectral/model/press_and_geopot.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_data_sma.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_reduce_nocomm.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)atmos_spectral/init/jablonowski_2006.F90 $(SRCROOT)shared/mosaic/grid.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_util.inc /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_chksum_int.h $(SRCROOT)atmos_spectral/model/fv_advection.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)atmos_solo/atmos_model.F90 /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_do_global_field.h /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_get_boundary.h $(SRCROOT)shared/mosaic/create_xgrid.c /public/home/wab/lijh/GFDL/full_dynamics/../src/shared/mpp/include/mpp_global_sum_tl.h
clean: neat
	-rm -f .full_dynamics.x.cppdefs $(OBJ) full_dynamics.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
full_dynamics.x: $(OBJ) 
	$(LD) $(OBJ) -o full_dynamics.x  $(LDFLAGS)
