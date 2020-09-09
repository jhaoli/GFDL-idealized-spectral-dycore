module mountain_wave_mod

!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! This program is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------

! Jablonowski, C., P. H. Lauritzen, R. Nair and M. Taylor, 2008:
! Idealized test cases for the dynamical cores of Atmospheric General Circulation Models :
! A proposal for the NCAR ASP 2008 summer colloquium Idealized test cases for 3D dynamical cores

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use               fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, stdlog, &
                                 write_version_number, close_file, check_nml_error

use         constants_mod, only: PI,    &
                                 GRAV,  & ! GRAV=9.800, not 9.806 as specified by Polvani
                                 RDGAS, & ! RDGAS=287.04, not 287.00 as specified by Polvani
                                 OMEGA, & ! same as specified by Polvani
                                 RADIUS,& ! same as specified by Polvani
                                 KAPPA    ! same as specified by Polvani

use   vert_coordinate_mod, only: compute_vert_coord

use        transforms_mod, only: get_grid_boundaries, get_sin_lat, get_cos_lat, trans_grid_to_spherical, trans_spherical_to_grid,&
                                 get_deg_lon, get_deg_lat, &
                                 vor_div_from_uv_grid, uv_grid_from_vor_div, get_grid_domain

use  press_and_geopot_mod, only: press_and_geopot_init, pressure_variables

use      diag_manager_mod, only: diag_axis_init, register_static_field, send_data

implicit none
private

character(len=128), parameter :: version = &
'$Id: mountain_wave.F90,v 1.0 2020/08/15 add by lijh$'

character(len=128), parameter :: tagname = &
'$Name: siena_201207 $'

public :: mountain_wave

real :: T0   = 288.d0      ! K
real :: h0   = 2000.d0     ! m
real :: d    = 1.5e6 
real :: u0   = 20.d0       ! m s-1
real :: lonc = 0.5 * PI
real :: latc = PI / 6.0
real :: psp  = 93000.d0    ! Pa
real :: N    = 0.0182      ! s-1

namelist / mountain_wave_nml / T0, h0, d, u0, lonc, latc, psp, N

contains
!=========================================================================================================================
subroutine mountain_wave(sea_level_press, triang_trunc, &
                   vert_coord_option, vert_difference_option, scale_heights, surf_res, &
                   p_press, p_sigma, exponent, pk, bk, &
                   vors, divs, ts, ln_ps, ug, vg, tg, psg, vorg, divg, surf_geopotential)

  real,    intent(in) :: sea_level_press
  logical, intent(in) :: triang_trunc 
  character(len=*), intent(in) :: vert_coord_option, vert_difference_option
  real,    intent(in) :: scale_heights, surf_res, p_press, p_sigma, exponent
  real,    intent(out), dimension(:) :: pk, bk
  complex, intent(out), dimension(:,:,:) :: vors, divs, ts
  complex, intent(out), dimension(:,:  ) :: ln_ps
  real,    intent(out), dimension(:,:,:) :: ug, vg, tg
  real,    intent(out), dimension(:,:  ) :: psg
  real,    intent(out), dimension(:,:,:) :: vorg, divg
  real,    intent(out), dimension(:,:  ) :: surf_geopotential
  
  integer :: unit, ierr, io, i, j, k, num_lon, num_lat, num_levels
  real, dimension(size(ug,1), size(ug,2)) :: ln_psg
  real, dimension(size(ug,1), size(ug,2), size(ug,3)) :: p_full, ln_p_full
  real, dimension(size(ug,1), size(ug,2), size(ug,3)+1) :: p_half, ln_p_half

  ! The allocatable arrays below must be allocated for the global domain
  real, allocatable, dimension(:)   :: deg_lon, deg_lat, rad_lon, rad_lat, sin_lat, cos_lat
  real r
!------------------------------------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=mountain_wave_nml, iostat=io)
  ierr = check_nml_error(io, 'mountain_wave_nml')
#else
  unit = open_namelist_file()
  ierr=1
  do while (ierr /= 0)
    read(unit, nml=mountain_wave_nml, iostat=io, end=20)
    ierr = check_nml_error (io, 'mountain_wave_nml')
  enddo
20  call close_file (unit)
#endif
  call write_version_number(version, tagname)
  if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=mountain_wave_nml)

  call compute_vert_coord(vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, sea_level_press, pk, bk)

  call press_and_geopot_init(pk, bk, .false., vert_difference_option)

  num_lon    = size(ug,1)
  num_lat    = size(ug,2)
  num_levels = size(ug,3)
  allocate(rad_lon(num_lon), deg_lon(num_lon))
  allocate(rad_lat(num_lat), deg_lat(num_lat))
  allocate(sin_lat(num_lat), cos_lat(num_lat))


  call get_deg_lon(deg_lon)
  rad_lon = PI*deg_lon/180.
  call get_deg_lat(deg_lat)
  rad_lat = PI*deg_lat/180.
  call get_sin_lat(sin_lat)
  call get_cos_lat(cos_lat)

!! lijh begin
  do j = 1, num_lat
    do i = 1, num_lon
      r = RADIUS * acos(sin(latc) * sin_lat(j) + cos(latc) * cos_lat(j) * cos(rad_lon(i) - lonc))
      surf_geopotential(i,j) = GRAV * h0 * exp(-(r / d)**2)
    end do
  end do

  do k = 1, num_levels
    do j = 1, num_lat
      ug(:,j,k) = u0 * cos_lat(j)
    end do
  end do
  
  vg = 0.0

  do j = 1, num_lat
    do i = 1, num_lon
      psg(i,j) = psp * exp(-0.5 * RADIUS * N**2 * u0 / GRAV**2 / kappa * (u0 / RADIUS + 2.0 * OMEGA) * &
                 (sin_lat(j)**2 - 1.0) - N**2 / GRAV**2 / kappa * surf_geopotential(i,j))
    end do
  end do
  
  call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, psg)

  tg = 288.0

!! lijh end
! Transform grid fields to spherical waves then transform
! back to ensure that spectral and grid representations match.

  call trans_grid_to_spherical(tg, ts)
  call trans_spherical_to_grid(ts, tg)
  call vor_div_from_uv_grid(ug, vg, vors, divs, triang=triang_trunc)
  call uv_grid_from_vor_div(vors, divs, ug, vg)
  call trans_spherical_to_grid(vors, vorg)
  call trans_spherical_to_grid(divs, divg)
  
  ln_psg = alog(psg)
  call trans_grid_to_spherical(ln_psg, ln_ps)
  call trans_spherical_to_grid(ln_ps,  ln_psg)
  psg = exp(ln_psg)
  
  deallocate(deg_lon, deg_lat, rad_lon, rad_lat, sin_lat, cos_lat)
  
  return
end subroutine mountain_wave
!================================================================================

end module mountain_wave_mod
