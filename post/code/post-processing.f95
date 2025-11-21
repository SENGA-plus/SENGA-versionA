program post

use indata
use arrinit

implicit none

real temp_write(nx,ny,nz,nvar)

integer xoffset, yoffset, zoffset
integer ixproc,iyproc,izproc
integer i,j,k
integer ip,jp,kp
integer istart,ifinish,jstart,jfinish,kstart,kfinish
integer iproc
integer snap
integer v,n

character*4 pnproc
character*5 pnsnap
character*4 snapname
character*60 fname

!=========================================================================
! Start of time snapshot loop
!=========================================================================
do snap=snap_start,snap_end,snap_step
!=========================================================================
  write (snapname,'(I4.4)') snap
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Read Snapshot (snap)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  print*,'Reading files at snapshot: ',snapname
  !------------------------------------
  do iproc=0,(nxproc*nyproc*nzproc)-1
    ixproc = mod(iproc,nxproc)
    iyproc = iproc/nxproc
    izproc = iyproc/nyproc
    iyproc = mod(iyproc,nyproc)
    xoffset = ixproc*nx/nxproc
    yoffset = iyproc*ny/nyproc
    zoffset = izproc*nz/nzproc
  
    write(pnproc,'(I4.4)') iproc
    write(pnsnap,'(I5.5)') snap
    fname = '../../output/out'//pnsnap//pnproc//'.res'

    if (reacting.eq.1) then
      open(unit=13,file=trim(fname),status='old',form='unformatted')
              !drun, urun, vrun, wrun, erun, trun, prun, yrun, RRTE, time
      read(13) drun4,drun1,drun2,drun3,drun6,drun5,drun7,drun8,drun9,drun10
      close(13)
    else
      print*,'does not exist'
    endif

    !------------------
    ! Shifting the data
    !------------------
    do k = 1,(nz/nzproc)
      do j = 1,(ny/nyproc)
        do i = 1,(nx/nxproc)
    
          drun(xoffset+i,yoffset+j,zoffset+k)= drun4(i,j,k)    
          urun(xoffset+i,yoffset+j,zoffset+k)= drun1(i,j,k) 
          vrun(xoffset+i,yoffset+j,zoffset+k)= drun2(i,j,k)            
          wrun(xoffset+i,yoffset+j,zoffset+k)= drun3(i,j,k)
          trun(xoffset+i,yoffset+j,zoffset+k)= drun5(i,j,k)
          erun(xoffset+i,yoffset+j,zoffset+k)= drun6(i,j,k)
          prun(xoffset+i,yoffset+j,zoffset+k)= drun7(i,j,k)
           
          if (reacting.eq.1) then
            do n= 1,nspec
              yrun(xoffset+i,yoffset+j,zoffset+k,n)= drun8(i,j,k,n)
              rrte(xoffset+i,yoffset+j,zoffset+k,n)= drun9(i,j,k,n)
            enddo
          endif
    
        enddo
      enddo
    enddo
  
  enddo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! End of Snapshot Read
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  write(*,'("   ETIME = ", 1E12.4)') drun10

  do k = 1,nz
    do j = 1,ny 
      do i = 1,nx
         rrte_YF(i,j,k)=rrte(i,j,k,1)
      enddo
    enddo
  enddo

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Write Snapshot (snap)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !===================
  !Write Tecplot files
  !===================
  ! Midplane
  !---------
  if(midplane.eq.1) then  
    print*,'writing the midplane at snapshot: ',snapname                        
    k=nz/2

    fname=trim('../midplane/midplane'//snapname//'.dat')
    open(unit=10, file=fname,status='unknown')
  
    if (reacting.eq.1) then
      write (10,*) 'variables="x","y","z","rho","u","v","w","Yfuel","T","p","E"'
      write (10,*) 'zone i=',nx,'j=',ny ,'k=',nz/nz,'f=point'
      do j=1,ny
        do i=1,nx
          write(10,'(12e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                                urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                                yrun(i,j,k,1),trun(i,j,k),              &
                                prun(i,j,k),erun(i,j,k)
                        
        enddo
      enddo
    else
      write (10,*) 'variables="x","y","z","rho","u","v","w","T","p","E"'
      write (10,*) 'zone i=',nx,'j=',ny ,'k=',nz/nz,'f=point'
      do j=1,ny
        do i=1,nx          
          write(10,'(10e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                               urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                               trun(i,j,k),prun(i,j,k),erun(i,j,k)
                       
        enddo
      enddo
    endif
    close(10)
  endif

  !---------
  ! Inplane
  !---------
  if(inplane.eq.1)then
    print*,'writing the inlet plane at snapshot: ',snapname                        
    i=nx/nx

    fname=trim('../inplane/inplane'//snapname//'.dat')
    open(unit=10, file=fname,status='unknown')
    
    if (reacting.eq.1) then
      write (10,*) 'variables="x","y","z","rho","u","v","w","Yfuel","T","p","E"'   
      write (10,*) 'zone i=',nx/nx,'j=',ny ,'k=',nz,'f=point'
      do k=1,nz
        do j=1,ny
          write(10,'(12e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                                urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                                yrun(i,j,k,1),trun(i,j,k),&
                                prun(i,j,k),erun(i,j,k)
                         
        enddo
      enddo
    else
      write (10,*) 'variables="x","y","z","rho","u","v","w","T","p","E"'
      write (10,*) 'zone i=',nx/nx,'j=',ny ,'k=',nz,'f=point'
      do k=1,nz
        do j=1,ny
          write(10,'(10e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                                urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                                trun(i,j,k),prun(i,j,k),erun(i,j,k)
                         
        enddo
      enddo
    endif 
    close(10)
  endif

  !---------
  ! Outplane
  !---------
  if(outplane.eq.1) then  
    print*,'writing the outlet plane at snapshot: ',snapname                        
    i=nx

    fname=trim('../outplane/outplane'//snapname//'.dat')
    open(unit=10, file=fname,status='unknown')

    if (reacting.eq.1) then
      write (10,*) 'variables="x","y","z","rho","u","v","w","Yfuel","t","p","e"'
      write (10,*) 'zone i=',nx/nx,'j=',ny ,'k=',nz,'f=point'
      do k=1,nz
        do j=1,ny
          write(10,'(12e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                                urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                                yrun(i,j,k,1),trun(i,j,k),&
                                prun(i,j,k),erun(i,j,k)
                         
        enddo
      enddo
    else
      write (10,*) 'variables="x","y","z","rho","u","v","w","T","p","E"'
      write (10,*) 'zone i=',nx/nx,'j=',ny ,'k=',nz,'f=point'
      do k=1,nz
        do j=1,ny
          write(10,'(10e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                                urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                                trun(i,j,k),prun(i,j,k),erun(i,j,k)
                         
        enddo
      enddo
    endif 
    close(10)
  endif

  !-----------
  ! 3D Domain
  !-----------
  if(full_domain.eq.1) then  
    print*,'writing the full 3D domain at snapshot: ',snapname                        
     
    fname=trim('../3D/3D'//snapname//'.dat')
    open(unit=10, file=fname,status='unknown')

    if (reacting.eq.1) then
      write (10,*) 'variables="x","y","z","rho","u","v","w","Yfuel","T","p","E"'
      write (10,*) 'zone i=',nx,'j=',ny ,'k=',nz,'f=point'
      do k=1,nz
        do j=1,ny
          do i=1,nx
            write(10,'(12e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                                  urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                                  yrun(i,j,k,1),trun(i,j,k),              &
                                  prun(i,j,k),erun(i,j,k)
                      
          enddo
        enddo
      enddo
    else
      write (10,*) 'variables="x","y","z","d","u","v","w","t","p","e"'
      write (10,*) 'zone i=',nx,'j=',ny ,'k=',nz,'f=point'
      do k=1,nz
        do j=1,ny
          do i=1,nx
            write(10,'(10e17.9)') dble(i),dble(j),dble(k),drun(i,j,k),    &
                                  urun(i,j,k),vrun(i,j,k),wrun(i,j,k),    & 
                                  trun(i,j,k),prun(i,j,k),erun(i,j,k)
                      
          enddo
        enddo
      enddo
    endif
    close(10)
  endif

  !---------------------
  ! Full 3D Binary Write
  !---------------------
  if(full_domain_Binary.eq.1) then
    if (reacting.eq.1) then
    !++++++++++++++++++++++
      print*,'writing the full 3D domain (binary) at snapshot: ',snapname
      !---------temp array for writing the file-------
      do k = 1,nz
        do j = 1,ny 
          do i = 1,nx
            crun_YF(i,j,k)=(Yfu-yrun(i,j,k,1))/(Yfu-Yfp)
          enddo
        enddo
      enddo

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            temp_write(i,j,k,1)=real(crun_YF(i,j,k))
            temp_write(i,j,k,2)=real(urun(i,j,k))
            temp_write(i,j,k,3)=real(vrun(i,j,k))
            temp_write(i,j,k,4)=real(wrun(i,j,k))
          enddo
        enddo
      enddo
      !-----------------------
      fname = trim('../3D/grid'//snapname//'.g')
      open(unit=10,file=fname,form='unformatted',status='new')
      ! Write the grid data file
      write(10) 1
      write(10) nx,ny,nz
      write(10) (((real(i-1)*deltax,i=1,nx),j=1,ny),k=1,nz), &
                (((real(j-1)*deltax,i=1,nx),j=1,ny),k=1,nz), &
                (((real(k-1)*deltax,i=1,nx),j=1,ny),k=1,nz)
      close(10)
      !---------------------------------------------- 

      fname = trim('../3D/3B'//snapname//'.all.f')
      open(unit=10, file=fname,form='unformatted',status='new')

      ! Write data
      write(10) 1
      write(10) nx,ny,nz,nvar
      write(10) ((((temp_write(i,j,k,v),i=1,nx),j=1,ny),k=1,nz),v=1,nvar)
      close(10)
    !++++++++++++++++++++++
    else
    !++++++++++++++++++++++
      print*,'NON-REACTING case not developed in this code'
    !++++++++++++++++++++++
    endif
    !++++++++++++++++++++++
  endif

  !===================
  ! End Tecplot Write
  !===================

  !---------------------------------------------- 
  !                HDF5 WRITE
  !---------------------------------------------- 
  if(full_domain_HDF5.eq.1) then
    if (reacting.eq.1) then
      call HDF5_write((snap))
    else
      print*,'NON-REACTING case not developed in this code'
    endif
  endif

  !---------------------------------------------- 
  !             Full 3D Decomp
  !---------------------------------------------- 
  if(decomp.eq.1) then  
    print*,'writing file full 3D Decomp for snapshot',snap                        
    open(unit=10, file='../3D/inflow.dat',form='unformatted',status='new')
    write(10) urun,vrun,wrun
    close(10)
  endif

!=========================================================================
! End of Time Snapshot Loop
!=========================================================================
enddo   
!=========================================================================
end
