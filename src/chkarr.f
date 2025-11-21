      SUBROUTINE CHKARR(ARRAY,
     +                  NCHKXL,NCHKXR,NCHKYL,NCHKYR,NCHKZL,NCHKZR,
     +                  ISTAC,ISTOC,JSTAC,JSTOC,KSTAC,KSTOC)
 
C     *************************************************************************
C
C     CHKARR
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
C
C     *************************************************************************


C     ARGUMENTS
C     =========
      INTEGER NCHKXL,NCHKXR,NCHKYL,NCHKYR,NCHKZL,NCHKZR
      DOUBLE PRECISION ARRAY(NCHKXL:NCHKXR,NCHKYL:NCHKYR,NCHKZL:NCHKZR)
      INTEGER ISTAC,ISTOC,JSTAC,JSTOC,KSTAC,KSTOC


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ARRMIN,ARRMAX
      INTEGER IC,JC,KC
      INTEGER ICMIN,JCMIN,KCMIN
      INTEGER ICMAX,JCMAX,KCMAX


C     BEGIN
C     =====

      WRITE(6,*)NCHKXL,NCHKXR,NCHKYL,NCHKYR,NCHKZL,NCHKZR
      WRITE(6,*)ISTAC,ISTOC,JSTAC,JSTOC,KSTAC,KSTOC

C     =========================================================================

      ARRMIN = ARRAY(ISTAC,JSTAC,KSTAC)
      ICMIN = ISTAC
      JCMIN = JSTAC
      KCMIN = KSTAC
      ARRMAX = ARRAY(ISTAC,JSTAC,KSTAC)
      ICMAX = ISTAC
      JCMAX = JSTAC
      KCMAX = KSTAC

      DO KC = KSTAC,KSTOC
        DO JC = JSTAC,JSTOC
          DO IC = ISTAC,ISTOC

            IF(ARRAY(IC,JC,KC).LT.ARRMIN)THEN
              ARRMIN = ARRAY(IC,JC,KC)
              ICMIN = IC
              JCMIN = JC
              KCMIN = KC
            ENDIF
            IF(ARRAY(IC,JC,KC).GT.ARRMAX)THEN
              ARRMAX = ARRAY(IC,JC,KC)
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
