module arrinit

use dtyps
use indata
use chmech
implicit none

double precision drun1(nxsize,nysize,nzsize)
double precision drun2(nxsize,nysize,nzsize)
double precision drun3(nxsize,nysize,nzsize)
double precision drun4(nxsize,nysize,nzsize)
double precision drun5(nxsize,nysize,nzsize)
double precision drun6(nxsize,nysize,nzsize)
double precision drun7(nxsize,nysize,nzsize)
double precision drun8(nxsize,nysize,nzsize,nspec)
double precision drun9(nxsize,nysize,nzsize,nspec)
double precision drun10
      
double precision urun(nx,ny,nz)
double precision vrun(nx,ny,nz)
double precision wrun(nx,ny,nz)
double precision trun(nx,ny,nz)
double precision drun(nx,ny,nz)
double precision erun(nx,ny,nz)
double precision prun(nx,ny,nz)
double precision crun_YF(nx,ny,nz)
double precision rrte_YF(nx,ny,nx)
double precision yrun(nx,ny,nz,nspec)
double precision rrte(nx,ny,nz,nspec)

type (tag) :: wrtvar(nvars)
type (varspace) :: work(nvars)

end module 
