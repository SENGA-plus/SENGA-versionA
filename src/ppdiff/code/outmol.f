      SUBROUTINE OUTMOL
 
C     *************************************************************************
C
C     OUTMOL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     06-OCT-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     WRITES OUT THE MOLECULAR TRANSPORT DATA FOR PPDIFF
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_diffin.h'
      INCLUDE 'com_ppdcom.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER ISPEC,JSPEC
      INTEGER ICOEFF
 
      CHARACTER*50 ERRSTR


C     BEGIN
C     =====

C     =========================================================================

      OPEN(UNIT=NCOUTP,FILE=FNOUTP,STATUS='UNKNOWN',FORM='FORMATTED',
     +     ERR=1010)

C     =========================================================================

C     HEADER
      WRITE(NCOUTP,*)'*****************************'
      WRITE(NCOUTP,*)'*                           *'
      WRITE(NCOUTP,*)'*  Output file from ppdiff  *'
      WRITE(NCOUTP,*)'*                           *'
      WRITE(NCOUTP,*)'*****************************'

C     =========================================================================

C     SPECIES LIST
      WRITE(NCOUTP,*)'Species list:'
      WRITE(NCOUTP,'(I5)')NSPEC
      DO ISPEC = 1, NSPEC
        WRITE(NCOUTP,'(I5,3X,A)')ISPEC,SPCSYM(ISPEC)
      ENDDO

C     SPECIES DATA
      WRITE(NCOUTP,*)'Molecular transport reference data:'
      WRITE(NCOUTP,'(2(1PE12.4))')PREFGB,TREFGB

      WRITE(NCOUTP,*)'Viscosity'
      DO ISPEC = 1,NSPEC
        WRITE(NCOUTP,'(2I5)')ISPEC,NCOVIS
        WRITE(NCOUTP,'(5(1PE15.7))')
     +       (VISCCO(ICOEFF,ISPEC),ICOEFF=1,NCOVIS)
      ENDDO

      WRITE(NCOUTP,*)'Thermal conductivity'
      DO ISPEC = 1,NSPEC
        WRITE(NCOUTP,'(2I5)')ISPEC,NCOCON
        WRITE(NCOUTP,'(5(1PE15.7))')
     +       (CONDCO(ICOEFF,ISPEC),ICOEFF=1,NCOCON)
      ENDDO

      WRITE(NCOUTP,*)'Binary diffusion coefficient'
      DO ISPEC = 1,NSPEC
        WRITE(NCOUTP,'(2I5)')ISPEC,NCODIF
        DO JSPEC = 1,ISPEC
          WRITE(NCOUTP,'(I5,5(1PE15.7))')JSPEC,
     +         (DIFFCO(ICOEFF,JSPEC,ISPEC),ICOEFF=1,NCODIF)
        ENDDO
      ENDDO

      WRITE(NCOUTP,*)'Thermal diffusion ratio'
      DO ISPEC = 1,NSPEC
        WRITE(NCOUTP,'(2I5)')ISPEC,NCOTDR
        DO JSPEC = 1,ISPEC-1
          WRITE(NCOUTP,'(I5,5(1PE15.7))')JSPEC,
     +         (TDRCCO(ICOEFF,JSPEC,ISPEC),ICOEFF=1,NCOTDR)
        ENDDO
      ENDDO

C     =========================================================================

      WRITE(NCOUTP,*)'End of file'
      CLOSE(NCOUTP)

C     =========================================================================

C     FILENAME ERROR HANDLER
      GOTO 1020
1010  ERRSTR = 'opening output file:'
      CALL ERHAND(ERRSTR,IEFATL)
1020  CONTINUE
C     END OF FILENAME ERROR HANDLER

C     =========================================================================


      RETURN
      END
