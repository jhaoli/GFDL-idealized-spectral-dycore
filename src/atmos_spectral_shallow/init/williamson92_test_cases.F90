module williamson92_test_cases_mod
  
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
  public steady_geostrophic_flow
  public mountain_zonal_flow
  public rossby_haurwitz_wave

contains
  subroutine steady_geostrophic_flow(Dyn, cos_lat, sin_lat, is, ie, js, je)

    type(dynamics_type)   , intent(inout) :: Dyn
    integer               , intent(in   ) :: is, ie, js, je
    real, dimension(js:je), intent(in   ) :: sin_lat, cos_lat
    real, parameter :: u0 = 2 * PI * RADIUS / (12.0 * 86400.0)
    real, parameter :: gz0 = 2.94e4
    integer i, j
   
   ! do j = js, je
   !   do i = is, ie
   !     Dyn%Grid%u(i,j,1) = u0 * cos_lat(j)
   !     Dyn%Grid%h(i,j,1) = gz0 - (RADIUS * OMEGA * u0 + u0**2 * 0.5) * sin_lat(j)**2
   !   end do
   ! end do

   ! Dyn%Grid%v(:,:,1) = 0.d0
   ! Dyn%Grid%hs(:,:)  = 0.d0
  
   ! call vor_div_from_uv_grid(Dyn%Grid%u  (:,:,1), Dyn%Grid%v  (:,:,1), &
   !                           Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1))
   ! call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
   ! call trans_spherical_to_grid(Dyn%Spec%vor(:,:,1), Dyn%Grid%vor(:,:,1))
   ! call trans_spherical_to_grid(Dyn%Spec%div(:,:,1), Dyn%Grid%div(:,:,1))
    
    Dyn%Grid%div(:,:,1) = 0.0
    do j = js, je
      Dyn%Grid%vor(:,j,1) = 2.0 * u0 / RADIUS * sin_lat(j)
      Dyn%Grid%h(:,j,1) = gz0 - (RADIUS * OMEGA * u0 + u0**2 * 0.5) * sin_lat(j)**2
    end do
    Dyn%Grid%hs(:,:) = 0.0
    call trans_grid_to_spherical(Dyn%Grid%vor(:,:,1), Dyn%Spec%vor(:,:,1))
    call trans_grid_to_spherical(Dyn%Grid%div(:,:,1), Dyn%Spec%div(:,:,1))
    call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
    call uv_grid_from_vor_div(Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1), &
                              Dyn%Grid%u(:,:,1), Dyn%Grid%v(:,:,1))

  end subroutine steady_geostrophic_flow

  subroutine mountain_zonal_flow(Dyn, cos_lat, sin_lat, is, ie, js, je)

    type(dynamics_type), intent(inout) :: Dyn
    integer, intent(in) :: is, ie, js, je
    real, dimension(js:je), intent(in) :: sin_lat, cos_lat
    real, parameter :: u0 = 20.0
    real, parameter :: gh0 = 5960.0 * GRAV
    real, parameter :: lon0 = PI * 1.5
    real, parameter :: lat0 = PI / 6.0
    real, parameter :: hs0  = 2000.0
    real, parameter :: R = PI / 9.0
    
    real, dimension(is:ie) :: deg_lon, rad_lon
    real, dimension(js:je) :: deg_lat, rad_lat
    integer num_lon, num_lat
    real dlon, d
    integer i, j
  
    call get_deg_lon(deg_lon)
    rad_lon = PI * deg_lon / 180.
    call get_deg_lat(deg_lat)
    rad_lat = PI * deg_lat / 180.

    do j = js, je
      do i = is, ie
        dlon = abs(rad_lon(i) - lon0)
        dlon = min(dlon, 2.0 * PI - dlon)
        d = min(R, sqrt(dlon**2 + (rad_lat(j) - lat0)**2))
        Dyn%Grid%hs(i,j) = hs0 * (1.0 - d / R) * GRAV
      end do
    end do

    do j = js, je
      Dyn%Grid%u(:,j,1) = u0 * cos_lat(j)
    end do
    Dyn%Grid%v(:,:,1) = 0.0
  
    do j = js, je
      do i = is, ie
        Dyn%Grid%h(i,j,1) = gh0 - (RADIUS * OMEGA * u0 + u0**2 * 0.5) * sin_lat(j)**2 - Dyn%Grid%hs(i,j)
      end do
    end do
    Dyn%Grid%div(:,:,1) = 0.0
    do j = js, je
      Dyn%Grid%vor(:,j,1) = 2.0 * u0 / RADIUS * sin_lat(j)
    end do
    call vor_div_from_uv_grid(Dyn%Grid%u  (:,:,1), Dyn%Grid%v  (:,:,1), &
                              Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1))
    call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
!or  
    ! call trans_grid_to_spherical(Dyn%Grid%vor(:,:,1), Dyn%Spec%vor(:,:,1))
    ! call trans_grid_to_spherical(Dyn%Grid%div(:,:,1), Dyn%Spec%div(:,:,1))
    ! call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
    ! call uv_grid_from_vor_div(Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1), &
                              ! Dyn%Grid%u(:,:,1), Dyn%Grid%v(:,:,1))

  end subroutine mountain_zonal_flow
  
  subroutine rossby_haurwitz_wave(Dyn, cos_lat, sin_lat, is, ie, js, je)
    
    type(dynamics_type), intent(inout) :: Dyn
    integer, intent(in) :: is, ie, js, je
    real, dimension(js:je), intent(in) :: sin_lat, cos_lat
    real, dimension(is:ie) :: deg_lon, rad_lon
    real, dimension(js:je) :: deg_lat, rad_lat
    real, parameter :: R = 4.0
    real, parameter :: omg = 7.848e-6
    real, parameter :: gh0 = 8.0e3 * GRAV
    real a, b, c
    integer i, j

    call get_deg_lon(deg_lon)
    rad_lon = PI * deg_lon / 180.
    call get_deg_lat(deg_lat)
    rad_lat = PI * deg_lat / 180.

    do j = js, je
      do i = is, ie
        a = cos_lat(j)
        b = R * cos_lat(j)**(R - 1) * sin_lat(j)**2 * cos(R * rad_lon(i))
        c = cos_lat(j)**(R + 1) * cos(R * rad_lon(i))
        Dyn%Grid%u(i,j,1) = RADIUS * omg * (a + b - c)
      end do
    end do

    do j = js, je
      do i = is, ie
        a = R * cos_lat(j)**(R - 1) * sin_lat(j) * sin(R * rad_lon(i))
        Dyn%Grid%v(i,j,1) = -RADIUS * omg * a
      end do
    end do

    do j = js, je 
      a = 0.5 * omg * (2 * OMEGA + omg) * cos_lat(j)**2 + &
         0.25 * omg**2 * ((R + 1) * cos_lat(j)**(2 * R + 2) + (2 * R**2 - R - 2) * cos_lat(j)**(2 * R) - 2 * R**2 * cos_lat(j)**(2 * R - 2))
      b = 2 * (OMEGA + omg) * omg * cos_lat(j)**R * &
          (R**2 + 2 * R + 2 - (R + 1)**2 * cos_lat(j)**2) / (R + 1) / (R + 2)
      c = 0.25 * omg**2 * cos_lat(j)**(2 * R) * ((R + 1) * cos_lat(j)**2 - R - 2)
      do i = is, ie
        Dyn%Grid%h(i,j,1) = gh0 + RADIUS**2 * (a + b * cos(R * rad_lon(i)) + c * cos(2 * R * rad_lon(i)))
      end do
    end do

    Dyn%Grid%hs(:,:) = 0.0 
    Dyn%Grid%div(:,:,1) = 0.0
    do j = js, je
      do i = is, ie
        Dyn%Grid%vor(i,j,1) = 2.0 * omg * sin_lat(j) - omg * sin_lat(j) * cos_lat(j)**R * (R**2 + 3.0 * R + 2.0) * cos(R * rad_lon(i))
      end do
    end do
    call vor_div_from_uv_grid(Dyn%Grid%u  (:,:,1), Dyn%Grid%v  (:,:,1), &
                              Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1))
    call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
!or
    ! call trans_grid_to_spherical(Dyn%Grid%vor(:,:,1), Dyn%Spec%vor(:,:,1))
    ! call trans_grid_to_spherical(Dyn%Grid%div(:,:,1), Dyn%Spec%div(:,:,1))
    ! call trans_grid_to_spherical(Dyn%Grid%h(:,:,1), Dyn%Spec%h(:,:,1))
    ! call uv_grid_from_vor_div(Dyn%Spec%vor(:,:,1), Dyn%Spec%div(:,:,1), &
    !                           Dyn%Grid%u(:,:,1), Dyn%Grid%v(:,:,1))

  end subroutine rossby_haurwitz_wave

end module williamson92_test_cases_mod
