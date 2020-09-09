module shallow_dynamics_type_mod

implicit none

type grid_type
   real, pointer, dimension(:,:,:) :: u=>NULL(), v=>NULL(), vor=>NULL(), div=>NULL(), h=>NULL(), trs=>NULL(), tr=>NULL()
   real, pointer, dimension(:,:)   :: stream=>NULL(), pv=>NULL(), hs=>NULL()
end type
type spectral_type
   complex, pointer, dimension(:,:,:) :: vor=>NULL(), div=>NULL(), h=>NULL(), trs=>NULL()
end type
type tendency_type
   real, pointer, dimension(:,:) :: u=>NULL(), v=>NULL(), h=>NULL(), trs=>NULL(), tr=>NULL()
end type
type dynamics_type
   type(grid_type)     :: grid
   type(spectral_type) :: spec
   type(tendency_type) :: tend
   integer             :: num_lon, num_lat
   logical             :: grid_tracer, spec_tracer
end type

end module shallow_dynamics_type_mod