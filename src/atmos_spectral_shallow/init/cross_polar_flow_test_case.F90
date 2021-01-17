module cross_polar_flow_test_case_mod
  
  use constants_mod, only : RADIUS, OMEGA, GRAV, PI
  use  shallow_dynamics_type_mod, only : dynamics_type
  use transforms_mod, only: get_sin_lat, get_cos_lat, &
                            get_deg_lat, get_deg_lon, &
                            vor_div_from_uv_grid,     &
                            trans_grid_to_spherical,  &
                            trans_spherical_to_grid,  &
                            uv_grid_from_vor_div
  implicit none

  private
  public :: cross_polar_flow

contains
  subroutine cross_polar_flow(Dyn, cos_lat, sin_lat, is, ie, js, je)

    type(dynamics_type)   , intent(inout) :: Dyn
    integer               , intent(in   ) :: is, ie, js, je
    real, dimension(js:je), intent(in   ) :: sin_lat, cos_lat
    real, parameter :: u0 = 20 ! m s-1
    real, parameter :: gz0 = 5.7684e4 ! m2 s-2
    real, dimension(is:ie) :: deg_lon, rad_lon
    real, dimension(js:je) :: deg_lat, rad_lat
    integer i, j
    
    call get_deg_lon(deg_lon)
    rad_lon = PI * deg_lon / 180.
    call get_deg_lat(deg_lat)
    rad_lat = PI * deg_lat / 180. 

    do j = js, je
      do i = is, ie
        Dyn%Grid%u(i,j,1) = -u0 * sin(rad_lon(i)) * sin_lat(j) * (4.0 * cos_lat(j)**2 - 1.0)
        Dyn%Grid%v(i,j,1) = u0 * sin_lat(j)**2 * cos(rad_lon(i))
        Dyn%Grid%h(i,j,1) = gz0 + 2.0 * RADIUS * OMEGA * u0  * sin_lat(j)**3 * cos_lat(j) * sin(rad_lon(i))
      end do
    end do

    Dyn%Grid%hs(:,:)  = 0.d0
  
    call vor_div_from_uv_grid(Dyn%Grid%u  (:,:,1), Dyn%Grid%v  (:,:,1), &
                              Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1))
    call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
    call trans_spherical_to_grid(Dyn%Spec%vor(:,:,1), Dyn%Grid%vor(:,:,1))
    call trans_spherical_to_grid(Dyn%Spec%div(:,:,1), Dyn%Grid%div(:,:,1))
    
  end subroutine cross_polar_flow

end module cross_polar_flow_test_case_mod
