      SUBROUTINE CXHALS(BIGARR,JNDEXL,JNDEXR,KNDEXL,KNDEXR)

C     *************************************************************************
C
C     CXHALS
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT
C
C     CHANGE RECORD
C     -------------
C     16-MAY-2003:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     CARRIES OUT HALO EXCHANGE FOR PERIODIC BCS IN X DIRECTION
C     FOR SPECIES MASS FRACTIONS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION BIGARR(NXBIGL:NXBIGR,NYBIGL:NYBIGR,NZBIGL:NZBIGR,
     +                        NSPCMX)
      INTEGER JNDEXL,JNDEXR,KNDEXL,KNDEXR
      

C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC
      INTEGER IS,ISPEC


C     BEGIN
C     =====

C     =========================================================================

C     RUN THROUGH ALL SPECIES
C     -----------------------
      DO ISPEC = 1,NSPEC

C       =======================================================================

C       RIGHT OUTER HALO SET EQUAL TO LEFT INNER HALO
        DO KC = KNDEXL,KNDEXR
          DO JC = JNDEXL,JNDEXR

            IS = ISTALI - 1
            DO IC = ISTARO,ISTORO

              IS = IS + 1
              BIGARR(IC,JC,KC,ISPEC) = BIGARR(IS,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO

C       =======================================================================

C       LEFT OUTER HALO SET EQUAL TO RIGHT INNER HALO
        DO KC = KNDEXL,KNDEXR
          DO JC = JNDEXL,JNDEXR

            IS = ISTARI - 1
            DO IC = ISTALO,ISTOLO

              IS = IS + 1
              BIGARR(IC,JC,KC,ISPEC) = BIGARR(IS,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO

C       =======================================================================

      ENDDO
C     END OF SPECIES LOOP

C     =========================================================================


      RETURN
      END
