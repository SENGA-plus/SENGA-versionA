module indata

implicit none

!===========================================#
!       Original domain configuration
!===========================================#
double precision, parameter :: Lx=1.0D-2

integer, parameter :: nx=1000
integer, parameter :: ny=20
integer, parameter :: nz=20

integer, parameter :: nxproc=40
integer, parameter :: nyproc=1
integer, parameter :: nzproc=1

!-> NOTE: The number of species (nspec) is defined in the tags module
!integer, parameter ::nspec=36

!-> number of variables in binary file
integer, parameter :: nvar=4   
!-> number of variables in HDF5 file
integer, parameter :: nvars=6

integer, parameter ::nxsize=nx/nxproc
integer, parameter ::nysize=ny/nyproc
integer, parameter ::nzsize=nz/nzproc

double precision ,parameter :: Yfu=0.055d0
double precision ,parameter :: Yfp=3.1d-5

double precision, parameter :: deltax=Lx/(dble(nx-1))
!===========================================#
!          Read/Write Parameters 
!===========================================#
integer, parameter :: snap_start=0
integer, parameter :: snap_step=1
integer, parameter :: snap_end=20
!---0=off 1=on----
integer, parameter :: inplane=0
integer, parameter :: midplane=0
integer, parameter :: outplane=0
integer, parameter :: full_domain=0
integer, parameter :: decomp=0
integer, parameter :: reacting=1
integer, parameter :: full_domain_Binary=0
integer, parameter :: full_domain_HDF5=1
integer, parameter :: Q_criterion=0

end module indata
