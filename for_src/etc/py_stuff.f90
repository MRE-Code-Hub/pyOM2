
! dummies for new subroutine and modules

subroutine panic_snap
end subroutine panic_snap

module rossmix_module
logical :: enable_rossmix =.false., enable_rossmix_mean_flow_interaction=.false.
end module rossmix_module

subroutine rossmix_friction
end subroutine rossmix_friction

subroutine rossmix_eddy_advect
end subroutine rossmix_eddy_advect

module rossmix2_module
logical :: enable_rossmix2 =.false.
real*8, allocatable :: E_r(:,:,:,:),e_back(:,:,:)
end module rossmix2_module

subroutine  rossmix2_friction
end subroutine  rossmix2_friction

subroutine  rossmix2_eddy_advect
end subroutine  rossmix2_eddy_advect

module biharmonic_thickness_module
logical :: enable_biharmonic_thickness_backscatter_integrate_energy =.false.
end module biharmonic_thickness_module

subroutine biharmonic_thickness_friction
end subroutine biharmonic_thickness_friction

subroutine biharmonic_thickness_mixing
end subroutine biharmonic_thickness_mixing


