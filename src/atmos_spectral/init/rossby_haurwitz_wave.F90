module rossby_haurwitz_wave_mod

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
'$Id: rossby_haurwitz_wave.F90,v 1.0 2020/07/21 add by lijh$'

character(len=128), parameter :: tagname = &
'$Name: siena_201207 $'

public :: rossby_haurwitz_wave

real, parameter :: P00    = 1.e5 ! Used to compute potential temperature. sea_level_press may be the same, but one should not assume that.
real, parameter :: halfpi = .5*PI

real :: R     = 4.       ! wave number 
real :: omg   = 1.962e-6 ! rotation parameter
real :: pref  = 9.55e4   ! reference pressure
real :: gamma = 0.0065   ! lapse rate below tropopause (deg K/m)
real :: T0    = 288.     ! global mean surface temperature (deg K)

namelist / rossby_haurwitz_wave_nml / R, omg, pref, gamma, T0

contains
!=========================================================================================================================
subroutine rossby_haurwitz_wave(sea_level_press, triang_trunc, &
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
  ! integer :: id_lon, id_lat, id_phalf, id_pfull, id_basic_flow, id_basic_temp, id_pot_temp, id_perturbation
  ! logical :: used
  real, dimension(size(ug,1), size(ug,2)) :: ln_psg
  real, dimension(size(ug,1), size(ug,2), size(ug,3)) :: p_full, ln_p_full
  real, dimension(size(ug,1), size(ug,2), size(ug,3)+1) :: p_half, ln_p_half

  ! The allocatable arrays below must be allocated for the global domain
  real, allocatable, dimension(:)   :: deg_lon, deg_lat, rad_lon, rad_lat, sin_lat, cos_lat, coriolis
  ! real, allocatable, dimension(:,:) :: basic_flow, basic_temp, pot_temp
  real a, b, c, gzp
!------------------------------------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=rossby_haurwitz_wave_nml, iostat=io)
  ierr = check_nml_error(io, 'rossby_haurwitz_wave_nml')
#else
  unit = open_namelist_file()
  ierr=1
  do while (ierr /= 0)
    read(unit, nml=rossby_haurwitz_wave_nml, iostat=io, end=20)
    ierr = check_nml_error (io, 'rossby_haurwitz_wave_nml')
  enddo
20  call close_file (unit)
#endif
  call write_version_number(version, tagname)
  if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=rossby_haurwitz_wave_nml)

  call compute_vert_coord(vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, sea_level_press, pk, bk)

  call press_and_geopot_init(pk, bk, .false., vert_difference_option)

  num_lon    = size(ug,1)
  num_lat    = size(ug,2)
  num_levels = size(ug,3)
  allocate(rad_lon(num_lon), deg_lon(num_lon))
  allocate(rad_lat(num_lat), deg_lat(num_lat))
  allocate(sin_lat(num_lat), cos_lat(num_lat), coriolis(num_lat))
! allocate(basic_flow(num_lat,num_levels))
! allocate(basic_temp(num_lat,num_levels), pot_temp(num_lat,num_levels))

  call get_deg_lon(deg_lon)
  rad_lon = PI*deg_lon/180.
  call get_deg_lat(deg_lat)
  rad_lat = PI*deg_lat/180.
  call get_sin_lat(sin_lat)
  call get_cos_lat(cos_lat)
  coriolis = 2*OMEGA*sin_lat

  surf_geopotential = 0

! id_phalf = diag_axis_init('phalf',.01*p_half,'hPa','z','approx half pressure level',direction=-1)
! id_pfull = diag_axis_init('pfull',.01*p_full,'hPa','z','approx full pressure level',direction=-1,edges=id_phalf)
! id_lon   = diag_axis_init('lon',deg_lon,'degrees E','x','longitude')
! id_lat   = diag_axis_init('lat',deg_lat,'degrees N','y','latitude')
! id_basic_flow   = register_static_field('jablonowski_2006', 'basic_flow',  (/id_lat,id_pfull/),'initial zonal wind', 'm/sec')
! id_basic_temp   = register_static_field('jablonowski_2006', 'basic_temp',  (/id_lat,id_pfull/),'initial basic_temp', 'K')
! id_pot_temp     = register_static_field('jablonowski_2006', 'pot_temp',    (/id_lat,id_pfull/),'initial potential temp', 'K')
! if(id_basic_flow > 0) used = send_data(id_basic_flow, basic_flow)


! if(id_basic_temp > 0) used = send_data(id_basic_temp, basic_temp)
! if(id_pot_temp   > 0) used = send_data(id_pot_temp,     pot_temp)

!! lijh begin
  do k = 1, num_levels
    do j = 1, num_lat
      do i = 1, num_lon 
        a = cos_lat(j)
        b = R * cos_lat(j)**(R - 1) * sin_lat(j)**2 * cos(R * rad_lon(i))
        c = cos_lat(j)**(R + 1) * cos(R * rad_lon(i))
        ug(i,j,k) = RADIUS * omg * (a + b - c)
      end do
    end do
  end do
  
  do k = 1, num_levels
    do j = 1, num_lat
      do i = 1, num_lon
        a = R * cos_lat(j)**(R - 1) * sin_lat(j) * sin(R * rad_lon(i))
        vg(i,j,k) = - RADIUS * omg * a
      end do
    end do
  end do

  do j = 1, num_lat
    a = 0.5 * omg * (2 * OMEGA + omg) * cos_lat(j)**2 + &
        0.25 * omg**2 * cos_lat(j)**(2 * R) * ((R + 1) * cos_lat(j)**2 + (2 * R**2 - R - 2)) -&
        0.5 * R**2 * omg**2 * cos_lat(j)**(2 * R - 2)
    b = 2 * (OMEGA + omg) * omg * cos_lat(j)**R * &
        (R**2 + 2 * R + 2 - (R + 1)**2 * cos_lat(j)**2) / (R + 1) / (R + 2)
    c = 0.25 * omg**2 * cos_lat(j)**(2 + R) * ((R + 1) * cos_lat(j)**2 - R - 2)
    do i = 1, num_lon
      gzp = RADIUS**2 * (a + b * cos(R * rad_lon(i)) + c * cos(2 * R * rad_lon(i)))
      psg(i,j) = pref * (1 + gamma / GRAV / T0 * gzp)**(GRAV / gamma / RDGAS)
    end do
  end do
  
  call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, psg)

  do k = 1, num_levels
    do j = 1, num_lat
      do i = 1, num_lon
       tg(i,j,k) = T0 * (p_full(i,j,k) / pref) ** (gamma * RDGAS / GRAV)
       ! pot_temp(i,j,k) = tg(i,j,k) * (P00 / p_full(i,j,k))**KAPPA
      end do
    enddo
  enddo

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
  
  deallocate(deg_lon, deg_lat, rad_lon, rad_lat, sin_lat, cos_lat, coriolis)
  ! deallocate(basic_flow, basic_temp, pot_temp)
  
  return
end subroutine rossby_haurwitz_wave
!================================================================================

end module rossby_haurwitz_wave_mod
