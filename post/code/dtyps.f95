module dtyps

implicit none

type tag
  character(len=12) wrtdst
end type

type varspace
  real, dimension(:,:,:), allocatable :: var
end type

end module dtyps
