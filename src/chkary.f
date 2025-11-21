      SUBROUTINE CHKARY(ARRAY,
     +                 NCHKXL,NCHKXR,NCHKYL,NCHKYR,NCHKZL,NCHKZR,NCPCMX,
     +                  ISTAC,ISTOC,JSTAC,JSTOC,KSTAC,KSTOC,ICPEC)
 
C     *************************************************************************
C
C     CHKARY
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT - CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     04-JUL-2004:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     DIAGNOSTIC ROUTINE
C     CHECKS UNIFORMITY OF THE SPECIFIED ARRAY
C     VARIANT FOR SPECIES ARRAYS
C
C     *************************************************************************


C     ARGUMENTS
C     =========
      INTEGER NCHKXL,NCHKXR,NCHKYL,NCHKYR,NCHKZL,NCHKZR,NCPCMX
      DOUBLE PRECISION
     +       ARRAY(NCHKXL:NCHKXR,NCHKYL:NCHKYR,NCHKZL:NCHKZR,NCPCMX)
      INTEGER ISTAC,ISTOC,JSTAC,JSTOC,KSTAC,KSTOC,ICPEC


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ARRMIN,ARRMAX
      INTEGER IC,JC,KC
      INTEGER ICMIN,JCMIN,KCMIN
      INTEGER ICMAX,JCMAX,KCMAX


C     BEGIN
C     =====

      WRITE(6,*)NCHKXL,NCHKXR,NCHKYL,NCHKYR,NCHKZL,NCHKZR,NCPCMX
      WRITE(6,*)ISTAC,ISTOC,JSTAC,JSTOC,KSTAC,KSTOC,ICPEC

C     =========================================================================

      ARRMIN = ARRAY(ISTAC,JSTAC,KSTAC,ICPEC)
      ICMIN = ISTAC
      JCMIN = JSTAC
      KCMIN = KSTAC
      ARRMAX = ARRAY(ISTAC,JSTAC,KSTAC,ICPEC)
      ICMAX = ISTAC
      JCMAX = JSTAC
      KCMAX = KSTAC

      DO KC = KSTAC,KSTOC
        DO JC = JSTAC,JSTOC
          DO IC = ISTAC,ISTOC

            IF(ARRAY(IC,JC,KC,ICPEC).LT.ARRMIN)THEN
              ARRMIN = ARRAY(IC,JC,KC,ICPEC)
              ICMIN = IC
              JCMIN = JC
              KCMIN = KC
            ENDIF
            IF(ARRAY(IC,JC,KC,ICPEC).GT.ARRMAX)THEN
              ARRMAX = ARRAY(IC,JC,KC,ICPEC)
              ICMAX = IC
              JCMAX = JC
              KCMAX = KC
            ENDIF

          ENDDO
        ENDDO
      ENDDO

      WRITE(6,'(2(1PE12.4,3I5))')ARRMIN,ICMIN,JCMIN,KCMIN,
     +                           ARRMAX,ICMAX,JCMAX,KCMAX

C     =========================================================================


      RETURN
      END
