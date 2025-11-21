subroutine XDMF_write(snap)

use indata
use chmech
use arrinit

implicit none
integer,intent(in) :: snap
integer i

character(len=4 ) snx
character(len=4 ) sny
character(len=4 ) snz
character(len=4 ) pnsnap
character(len=6 ) fname
character(len=60) fdir

!Domain dimensions
write(snx,"(I4)") nx
write(sny,"(I4)") ny
write(snz,"(I4)") nz

!------------------------------------------------------------------------
write (pnsnap,'(I4.4)') snap
fname='3B'//pnsnap
fdir ='../3D/3B'//pnsnap//'.xmf'

print*, 'writing the XDMF file to: ', fname//'.xmf'

open(unit=10,file=trim(fdir),status="unknown",form="formatted")


! Write header
write(10,"(a)") '<?xml version="1.0" ?>'
write(10,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
write(10,"(a)") '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'


!Definition of the domain and topology
write(10,"(2x,a)") '<Domain>'
write(10,"(4x,a,a,a)") '<Grid Name="',fname,'" GridType="Uniform">'
!In hdf5 files, x/y/z is written as z/y/x
write(10,"(6x,2a,2(x,a),a)")                                         & 
          '<Topology TopologyType="3DRectMesh" Dimensions="',        & 
          trim(adjustl(snz)), trim(adjustl(sny)),trim(adjustl(snx)), & 
          '"/>'

!==============================================================================
! Writing the grid
!------------------------------------------------------------------------------
write(10,"(6x,a)") '<Geometry GeometryType="VxVyVz">'
! x-coordinate
write(10,"(8x,a,a,a)") '<DataItem Dimensions="',trim(adjustl(snx)),           &
                       '" NumberType="Float" Precision="8" Format="HDF">'
write(10,"(10x,a)") trim(fname)//'.h5:/Grid/x-grid'
write(10,"(8x,a)") '</DataItem>'

! y-coordinate
write(10,"(8x,a,a,a)") '<DataItem Dimensions="',trim(adjustl(sny)),           &
                       '" NumberType="Float" Precision="8" Format="HDF">'
write(10,"(10x,a)") trim(fname)//'.h5:/Grid/y-grid'
write(10,"(8x,a)") '</DataItem>'

! z-coordinate
write(10,"(8x,a,a,a)") '<DataItem Dimensions="',trim(adjustl(snz)),           &
                       '" NumberType="Float" Precision="8" Format="HDF">'
write(10,"(10x,a)") trim(fname)//'.h5:/Grid/z-grid'
write(10,"(8x,a)") '</DataItem>'
!------------------------------------------------------------------------------
write(10,"(6x,a)") '</Geometry>'
!==============================================================================
! Writing data locations
!------------------------------------------------------------------------------
do i=1,nvars
  write(10,"(6x,a,a,a)") '<Attribute Name="',trim(wrtvar(i)%wrtdst),          &
                         '" AttributeType="Scalar" Center="Node">'
  ! In hdf5 files, x/y/z is written as z/y/x
  write(10,"(8x,2a,2(x,a),a)") '<DataItem Dimensions="',                      &
                  trim(adjustl(snz)),trim(adjustl(sny)), trim(adjustl(snx)),  &
                          '" NumberType="Float" Precision="8" Format="HDF">'  
  write(10,"(10x,a)") trim(fname)//'.h5:/Field/'//trim(wrtvar(i)%wrtdst)
  write(10,"(8x,a)") '</DataItem>'
  write(10,"(6x,a)") '</Attribute>'
enddo

do i=1,nspec
  write(10,"(6x,a,a,a)") '<Attribute Name="',trim(wrtspc(i)%wrtdst),          &
                         '" AttributeType="Scalar" Center="Node">'
  ! In hdf5 files, x/y/z is written as z/y/x
  write(10,"(8x,2a,2(x,a),a)") '<DataItem Dimensions="',                      &
                  trim(adjustl(snz)),trim(adjustl(sny)), trim(adjustl(snx)),  &
                          '" NumberType="Float" Precision="8" Format="HDF">'  
  write(10,"(10x,a)") trim(fname)//'.h5:/Field/'//trim(wrtspc(i)%wrtdst)
  write(10,"(8x,a)") '</DataItem>'
  write(10,"(6x,a)") '</Attribute>'
enddo

do i=1,nspec
  write(10,"(6x,a,a,a)") '<Attribute Name="',trim(wrtrrt(i)%wrtdst),          &
                         '" AttributeType="Scalar" Center="Node">'
  ! In hdf5 files, x/y/z is written as z/y/x
  write(10,"(8x,2a,2(x,a),a)") '<DataItem Dimensions="',                      &
                  trim(adjustl(snz)),trim(adjustl(sny)), trim(adjustl(snx)),  &
                          '" NumberType="Float" Precision="8" Format="HDF">'  
  write(10,"(10x,a)") trim(fname)//'.h5:/Field/'//trim(wrtrrt(i)%wrtdst)
  write(10,"(8x,a)") '</DataItem>'
  write(10,"(6x,a)") '</Attribute>'
enddo
!==============================================================================

write(10,"(4x,a)") '</Grid>'
write(10,"(2x,a)") '</Domain>'
write(10,"(a)") '</Xdmf>'

close(10)

return

end subroutine
