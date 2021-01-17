module jet_zonal_flow_test_case_mod
  
  use fms_mod, only: error_mesg,                     &
                     FATAL
  use constants_mod, only : RADIUS, OMEGA, GRAV, PI
  use  shallow_dynamics_type_mod, only : dynamics_type
  use transforms_mod, only: get_sin_lat, get_cos_lat, &
                            get_deg_lat, get_deg_lon, &
                            vor_div_from_uv_grid,     &
                            trans_grid_to_spherical,  &
                            trans_spherical_to_grid,  &
                            uv_grid_from_vor_div
  use quadpack_mod, only: qags

  implicit none

  private
  public jet_zonal_flow

  real, parameter :: u_max = 80.d0 ! m s-1
  real, parameter :: lat0 = PI / 7.d0
  real, parameter :: lat1 = PI / 2.d0 - lat0
  real, parameter :: en = exp(-4.0 / (lat1 - lat0)**2.0)
  real, parameter :: gh0 = GRAV * 1.0e4 ! m2 s-2
  real, parameter :: ghd = GRAV * 120.0
  real, parameter :: lat2 = PI / 4.d0
  real, parameter :: alpha = 1.d0 / 3.d0
  real, parameter :: beta = 1.d0 / 15.d0

contains

  subroutine jet_zonal_flow(Dyn, cos_lat, sin_lat, is, ie, js, je)

    type(dynamics_type)   , intent(inout) :: Dyn
    integer               , intent(in   ) :: is, ie, js, je
    real, dimension(js:je), intent(in   ) :: sin_lat, cos_lat
    real, dimension(is:ie) :: deg_lon, rad_lon
    real, dimension(js:je) :: deg_lat, rad_lat
    integer i, j, neval, ierr
    real abserr
    
    call get_deg_lon(deg_lon)
    rad_lon = PI * deg_lon / 180.
    call get_deg_lat(deg_lat)
    rad_lat = PI * deg_lat / 180. 

    Dyn%Grid%hs(:,:)  = 0.d0

    do j = js, je
      Dyn%Grid%u(:,j,1) = u_function(rad_lat(j))
    end do
   
    Dyn%Grid%v(:,:,1) = 0.d0

    do j = js, je
      i = is
      if (j == js) then
        Dyn%Grid%h(i,j,1) = gh0
      else
        call qags(gh_integrand, -0.5*PI, rad_lat(j), 1.0e-12, 1.0d-3, Dyn%Grid%h(i,j,1), abserr, neval, ierr)
        if (ierr /= 0) then
          call error_mesg('atmosphere', &
              'Failed to calculate integration in jet_zonal_flow_test_case_mod', FATAL)
        end if
        Dyn%Grid%h(i,j,1) = gh0 - Dyn%Grid%h(i,j,1)
      end if
      do i = is, ie
        Dyn%Grid%h(i,j,1) = Dyn%Grid%h(is,j,1)
        ! Add perturbation
        Dyn%Grid%h(i,j,1) = Dyn%Grid%h(i,j,1) + ghd * &
          cos(rad_lat(j)) * &
          exp(-(merge(rad_lon(i) - 2*PI, rad_lon(i), rad_lon(i) > PI) / alpha)**2) * &
          exp(-((lat2 - rad_lat(j)) / beta)**2)
      end do
    end do

    call vor_div_from_uv_grid(Dyn%Grid%u  (:,:,1), Dyn%Grid%v  (:,:,1), &
                              Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1))
    call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
    call trans_spherical_to_grid(Dyn%Spec%vor(:,:,1), Dyn%Grid%vor(:,:,1))
    call trans_spherical_to_grid(Dyn%Spec%div(:,:,1), Dyn%Grid%div(:,:,1))
    
  end subroutine jet_zonal_flow
  
  real function gh_integrand(lat) result(res)
    
    real, intent(in) :: lat
    
    real u, f

    u = u_function(lat)
    f = 2 * OMEGA * sin(lat)
    res = RADIUS * u * (f + tan(lat) / RADIUS * u)

  end function gh_integrand

  real function u_function(lat) result(res)
    
    real, intent(in) :: lat

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0.0
    else
      res = u_max / en * exp(1 / (lat - lat0) / (lat - lat1))
    end if
  end function u_function

end module jet_zonal_flow_test_case_mod
