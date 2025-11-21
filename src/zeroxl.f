      SUBROUTINE ZEROXL(FARRAY)
 
C     *************************************************************************
C
C     ZEROXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     06-JUL-2003:  RSC MODIFIED FOR SENGA2
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     ZEROS THE X-WISE LEFT END ELEMENTS OF THE SPECIFIED ARRAY
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION FARRAY(NXSIZE,NYSIZE,NZSIZE)


C     LOCAL DATA
C     ==========
      INTEGER JC,KC


C     BEGIN
C     =====

C     =========================================================================

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL

          FARRAY(ISTAL,JC,KC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
