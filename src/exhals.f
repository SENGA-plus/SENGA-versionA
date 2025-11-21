      SUBROUTINE EXHALS(BIGARR,ISPEC,BUFFER,INDEXL,INDEXR,
     +                                      JNDEXL,JNDEXR,
     +                                      KNDEXL,KNDEXR)

C     *************************************************************************
C
C     EXHALS
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
C     EXTRACTS HALO DATA FOR PARALLEL TRANSFER BUFFER
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
      DOUBLE PRECISION BUFFER(NPARAY)
      INTEGER INDEXL,INDEXR,JNDEXL,JNDEXR,KNDEXL,KNDEXR
      INTEGER ISPEC


C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC
      INTEGER ICOUNT


C     BEGIN
C     =====

C     =========================================================================

      ICOUNT = 0
      DO KC = KNDEXL,KNDEXR
        DO JC = JNDEXL,JNDEXR
          DO IC = INDEXL,INDEXR

            ICOUNT = ICOUNT + 1
            BUFFER(ICOUNT) = BIGARR(IC,JC,KC,ISPEC)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
