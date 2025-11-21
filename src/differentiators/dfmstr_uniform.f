      SUBROUTINE DFMSTR
 
C     *************************************************************************
C
C     DFMSTR
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     10-NOV-2013:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES MESH STRETCHING
C     IN ONE DIMENSION ONLY
C     REGULAR SPACING
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION DELTAG
      INTEGER IGSTAL,IGSTOL,IGOFST
      INTEGER IC,ICPROC


C     BEGIN
C     =====
      
C     =========================================================================

C     SPECIFY THE DIMENSION FOR COORDINATE STRETCHING
C     Y-DIRECTION
      IGOFST = 0
      DO ICPROC = 0, IYPROC-1
        IGOFST = IGOFST + NPMAPY(ICPROC)
      ENDDO
      IGSTAL = IGOFST + 1
      IGSTOL = IGOFST + NPMAPY(IYPROC)
      DELTAG = DELTAY

C     REGULAR COORDINATES
      DO IC = IGSTAL,IGSTOL
        GCMREG(IC) = REAL(IC-1)*DELTAG
        GCMSTR(IC) = GCMREG(IC)
        DGDHAT(IC) = ONE
        DGDHSQ(IC) = ONE
        D2GDH2(IC) = ZERO
      ENDDO     

C     DIAGNOSTICS
      WRITE(6,'(2I5)')IGSTAL,IGSTOL
C      DO IC = IGSTAL,IGSTOL
C        WRITE(6,'(6(1PE12.4))')GCMREG(IC),GCMSTR(IC),
C     +                         DHATDG(IC),D2HDG2(IC),
C     +                         DGDHAT(IC),D2GDH2(IC)
C      ENDDO     

      OPEN(UNIT=9,FILE='meshpoints.res',FORM='FORMATTED')
      DO IC = IGSTAL,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMREG(IC),GCMSTR(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

      OPEN(UNIT=9,FILE='dgdhat.res',FORM='FORMATTED')
      DO IC = IGSTAL+1,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMREG(IC),DGDHAT(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

      OPEN(UNIT=9,FILE='d2gdh2.res',FORM='FORMATTED')
      DO IC = IGSTAL+1,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMREG(IC),D2GDH2(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

C     =========================================================================


      RETURN
      END
