      SUBROUTINE CXHALO(BIGARR,JNDEXL,JNDEXR,KNDEXL,KNDEXR)

C     *************************************************************************
C
C     CXHALO
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT
C
C     CHANGE RECORD
C     -------------
C     14-MAY-2003:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     CARRIES OUT HALO EXCHANGE FOR PERIODIC BCS IN X DIRECTION
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION BIGARR(NXBIGL:NXBIGR,NYBIGL:NYBIGR,NZBIGL:NZBIGR)
      INTEGER JNDEXL,JNDEXR,KNDEXL,KNDEXR
      

C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC
      INTEGER IS


C     BEGIN
C     =====

C     =========================================================================

C     RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
      DO KC = KNDEXL,KNDEXR
        DO JC = JNDEXL,JNDEXR

          IS = ISTALI - 1
          DO IC = ISTARO,ISTORO

            IS = IS + 1
            BIGARR(IC,JC,KC) = BIGARR(IS,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================

C     LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
      DO KC = KNDEXL,KNDEXR
        DO JC = JNDEXL,JNDEXR

          IS = ISTARI - 1
          DO IC = ISTALO,ISTOLO

            IS = IS + 1
            BIGARR(IC,JC,KC) = BIGARR(IS,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
