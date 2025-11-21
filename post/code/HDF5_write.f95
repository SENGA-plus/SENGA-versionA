subroutine HDF5_write(snap)
use hdf5
use arrinit
use chmech
use indata
implicit none


integer i,j,k
integer, intent(in) :: snap


double precision dist_x(nx)
double precision dist_y(ny)
double precision dist_z(nz)
double precision dummy(nx,ny,nz)

!===========decleration for HDF5 write===========================
integer :: hdferr
integer(hid_t):: file_id  ! file identifer 
integer(hid_t):: d_id     ! data set identifier
integer(hid_t):: f_id     ! field group identifer
integer(hid_t):: g_id     ! grid group identifer
integer(hid_t):: dspace   ! Dataspace identifier
!-mod-----------------
integer, parameter :: ndim = 3, nxdim=1, nydim=1, nzdim=1
integer(hsize_t) gdim(ndim) 
integer(hsize_t) xdim(1) 
integer(hsize_t) ydim(1) 
integer(hsize_t) zdim(1) 
!=================================================================

character(len=60) fname
character(len=4 ) pnsnap
!-----------------------------------------------
! HDF5 tags for nvars variable defined in indata
!-----------------------------------------------
wrtvar(1)%wrtdst  = "URUN"
wrtvar(2)%wrtdst  = "VRUN"
wrtvar(3)%wrtdst  = "WRUN"
wrtvar(4)%wrtdst  = "PRUN"
wrtvar(5)%wrtdst  = "DRUN"
wrtvar(6)%wrtdst  = "TRUN"

work(1)%var=real(urun)
work(2)%var=real(vrun)
work(3)%var=real(wrun)
work(4)%var=real(prun)
work(5)%var=real(drun)
work(6)%var=real(trun)

!-----------------------------------------------
! HDF5 tags for the species of the chemical mech
!-----------------------------------------------
call set_mech_tags

gdim(1)=nx
gdim(2)=ny
gdim(3)=nz

xdim=nx
ydim=ny
zdim=nz

!-----
do i=1, nx
  dist_x(i)=dble(i-1)*deltax
enddo

do j=1, ny
  dist_y(j)=dble(j-1)*deltax
enddo

do k=1, nz
  dist_z(k)=dble(k-1)*deltax
enddo

call XDMF_write(snap)

!============================= HDF5 write ================================
write (pnsnap,'(I4.4)') snap
fname=trim('../3D/3B'//pnsnap//'.h5')

print*, 'writing the HDF5 file to: ', '3B'//pnsnap//'.h5'
print*, '==================================='

! initialize fortran interface.
call h5open_f(hdferr)

! create a new file using default properties.
call h5fcreate_f(fname, h5f_acc_trunc_f, file_id, hdferr)

!=========================================================================
! Write grid data
!=========================================================================
!create group at the root using absolute name containing the x/y/z coordinates datasets
call h5gcreate_f(file_id,'/Grid',g_id,hdferr)
!*************************************************************************
!---x coordinate---
! create the dataspace for the grid dataset in x
call h5screate_simple_f(nxdim,xdim,dspace,hdferr)
! create the dataset in group "/Grid" with default properties
call h5dcreate_f(g_id,'x-grid',h5t_native_real,dspace,d_id,hdferr)
! write in the dataset
call h5dwrite_f(d_id,h5t_native_real,real(dist_x),xdim,hdferr)
! close the dataspace for the first dataset.
call h5sclose_f(dspace, hdferr)
! close the first dataset.
call h5dclose_f(d_id, hdferr)

!---y coordinate---
! create the dataspace for the grid dataset in y
call h5screate_simple_f(nydim,ydim,dspace,hdferr)
! create the dataset in group "/Grid" with default properties
call h5dcreate_f(g_id,'y-grid',h5t_native_real,dspace,d_id,hdferr)
! write in the dataset
call h5dwrite_f(d_id,h5t_native_real,real(dist_y),ydim,hdferr)
! close the dataspace for the first dataset.
call h5sclose_f(dspace, hdferr)
! close the first dataset.
call h5dclose_f(d_id, hdferr)

!---z coordinate---
! create the dataspace for the grid dataset in z
call h5screate_simple_f(nzdim,zdim,dspace,hdferr)
! create the dataset in group "/Grid" with default properties
call h5dcreate_f(g_id,'z-grid',h5t_native_real,dspace,d_id,hdferr)
! write in the dataset
call h5dwrite_f(d_id,h5t_native_real,real(dist_z),zdim,hdferr)
! close the dataspace for the grid dataset.
call h5sclose_f(dspace, hdferr)
! close the grid dataset.
call h5dclose_f(d_id, hdferr)

! close group /Grid
call h5gclose_f(g_id,hdferr)

!=========================================================================
! Start writing datasets
!=========================================================================
! create group at the root using absolute name containing the x/y/z coordinates datasets
call h5gcreate_f(file_id,'/Field',f_id,hdferr)
!*************************************************************************
do i=1,nvars
  ! create the dataspace for the dataset.
  call h5screate_simple_f(ndim, gdim, dspace, hdferr)
  ! create a dataset in group "mygroup" with default properties.
  call h5dcreate_f(f_id, trim(wrtvar(i)%wrtdst), h5t_native_real, dspace,d_id, hdferr)
  ! write the dataset.
  call h5dwrite_f(d_id, h5t_native_real, work(i)%var, gdim, hdferr)
  ! close the dataspace for the dataset.
  call h5sclose_f(dspace, hdferr)
  ! close the dataset.
  call h5dclose_f(d_id, hdferr)
enddo

do i=1,nspec
  dummy = 0.0d0
  dummy(:,:,:) = yrun(:,:,:,i)
  ! create the dataspace for the dataset.
  call h5screate_simple_f(ndim, gdim, dspace, hdferr)
  ! create a dataset in group "mygroup" with default properties.
  call h5dcreate_f(f_id, trim(wrtspc(i)%wrtdst), h5t_native_real, dspace, d_id, hdferr)
  ! write the dataset.
  call h5dwrite_f(d_id, h5t_native_real, real(dummy), gdim, hdferr)
  ! close the dataspace for the dataset.
  call h5sclose_f(dspace, hdferr)
  ! close the dataset.
  call h5dclose_f(d_id, hdferr)
enddo

do i=1,nspec
  dummy = 0.0d0
  dummy(:,:,:) = rrte(:,:,:,i)
  ! create the dataspace for the dataset.
  call h5screate_simple_f(ndim, gdim, dspace, hdferr)
  ! create a dataset in group "mygroup" with default properties.
  call h5dcreate_f(f_id, trim(wrtrrt(i)%wrtdst), h5t_native_real, dspace, d_id, hdferr)
  ! write the dataset.
  call h5dwrite_f(d_id, h5t_native_real, real(dummy), gdim, hdferr)
  ! close the dataspace for the dataset.
  call h5sclose_f(dspace, hdferr)
  ! close the dataset.
  call h5dclose_f(d_id, hdferr)
enddo

! close group /Field
call h5gclose_f(f_id,hdferr)
! close the file.
call h5fclose_f(file_id, hdferr)
! close fortran interface.
call h5close_f(hdferr)

return
end
