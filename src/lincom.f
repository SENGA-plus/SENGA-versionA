      SUBROUTINE LINCOM
 
C     *************************************************************************
C
C     LINCOM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     15-JAN-2003:  CREATED
C     08-AUG-2012:  RSC EVALUATE ALL SPECIES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES INTERMEDIATE SOLUTION VALUES IN ERK SCHEME
C     BY DOING LINEAR COMBINATIONS OF LEFT- AND RIGHT-HAND SIDES
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FORNOW,TEMP1,TEMP2,TEMP3,TEMP4
      INTEGER IC,JC,KC,ISPEC


C     BEGIN
C     =====

C     =========================================================================

C     ERK SUBSTEP
C     ===========

C     -------------------------------------------------------------------------
C     NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
C     -------------------------------------------------------------------------

C     DENSITY
C     -------
      DO KC = KSTALD,KSTOLD
        DO JC = JSTALD,JSTOLD
          DO IC = ISTALD,ISTOLD

            DERR(IC,JC,KC) = DERR(IC,JC,KC)
     +                              + RKERR(IRKSTP)*DRHS(IC,JC,KC)

            FORNOW = DRUN(IC,JC,KC)
            DRUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*DRHS(IC,JC,KC)
            DRHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*DRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
                       
C     -------------------------------------------------------------------------

C     U-VELOCITY
C     ----------
      DO KC = KSTALU,KSTOLU
        DO JC = JSTALU,JSTOLU
          DO IC = ISTALU,ISTOLU

            UERR(IC,JC,KC) = UERR(IC,JC,KC)
     +                              + RKERR(IRKSTP)*URHS(IC,JC,KC)

            FORNOW = URUN(IC,JC,KC)
            URUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*URHS(IC,JC,KC)
            URHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*URHS(IC,JC,KC)
            
          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------
            
C     V-VELOCITY
C     ----------
      DO KC = KSTALV,KSTOLV
        DO JC = JSTALV,JSTOLV
          DO IC = ISTALV,ISTOLV

            VERR(IC,JC,KC) = VERR(IC,JC,KC)
     +                              + RKERR(IRKSTP)*VRHS(IC,JC,KC)

            FORNOW = VRUN(IC,JC,KC)
            VRUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*VRHS(IC,JC,KC)
            VRHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*VRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
      
C     -------------------------------------------------------------------------

C     W-VELOCITY
C     ----------
      DO KC = KSTALW,KSTOLW
        DO JC = JSTALW,JSTOLW
          DO IC = ISTALW,ISTOLW

            WERR(IC,JC,KC) = WERR(IC,JC,KC)
     +                              + RKERR(IRKSTP)*WRHS(IC,JC,KC)

            FORNOW = WRUN(IC,JC,KC)
            WRUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*WRHS(IC,JC,KC)
            WRHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*WRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
      
C     -------------------------------------------------------------------------

C     STAGNATION INTERNAL ENERGY
C     --------------------------
      DO KC = KSTALE,KSTOLE
        DO JC = JSTALE,JSTOLE
          DO IC = ISTALE,ISTOLE

            EERR(IC,JC,KC) = EERR(IC,JC,KC)
     +                              + RKERR(IRKSTP)*ERHS(IC,JC,KC)

            FORNOW = ERUN(IC,JC,KC)
            ERUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*ERHS(IC,JC,KC)
            ERHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*ERHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
      
C     -------------------------------------------------------------------------

C     SPECIES MASS FRACTIONS
C     ----------------------
C     RSC 08-AUG-2012 EVALUATE ALL SPECIES
C      DO ISPEC = 1,NSPM1
      DO ISPEC = 1,NSPEC

        DO KC = KSTALY,KSTOLY
          DO JC = JSTALY,JSTOLY
            DO IC = ISTALY,ISTOLY

              YERR(IC,JC,KC,ISPEC) = YERR(IC,JC,KC,ISPEC)
     +                             + RKERR(IRKSTP)*YRHS(IC,JC,KC,ISPEC)

              FORNOW = YRUN(IC,JC,KC,ISPEC)
              YRUN(IC,JC,KC,ISPEC) = FORNOW
     +                             + RKLHS(IRKSTP)*YRHS(IC,JC,KC,ISPEC)
              YRHS(IC,JC,KC,ISPEC) = FORNOW
     +                             + RKRHS(IRKSTP)*YRHS(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     VM & NC: GRADIENT OF SPECIES AT WALL EQUAL TO ZERO
      IF((NSBCXL.EQ.NSBCW2).OR.(NSBCXL.EQ.NSBCW1))THEN
        DO ISPEC=1,NSPEC
          DO KC=KSTAL,KSTOL
            DO JC=JSTAL,JSTOL
              TEMP1=YRHS(ISTAL+1,JC,KC,ISPEC)/DRHS(ISTAL+1,JC,KC)
              TEMP2=YRHS(ISTAL+2,JC,KC,ISPEC)/DRHS(ISTAL+2,JC,KC)
              TEMP3=YRHS(ISTAL+3,JC,KC,ISPEC)/DRHS(ISTAL+3,JC,KC)
              TEMP4=YRHS(ISTAL+4,JC,KC,ISPEC)/DRHS(ISTAL+4,JC,KC)
              YRUN(ISTAL,JC,KC,ISPEC)=(12.0/25.0)*(4.0*TEMP1-3.0*TEMP2+
     +                           (4.0/3.0)*TEMP3-(1.0/4.0*TEMP4))
              YRUN(ISTAL,JC,KC,ISPEC)=MIN(1.0,YRUN(ISTAL,JC,KC,ISPEC))
              YRUN(ISTAL,JC,KC,ISPEC)=MAX(0.0,YRUN(ISTAL,JC,KC,ISPEC))
              YRUN(ISTAL,JC,KC,ISPEC)=DRHS(ISTAL,JC,KC)
     +            *YRUN(ISTAL,JC,KC,ISPEC)
              YRHS(ISTAL,JC,KC,ISPEC)=YRUN(ISTAL,JC,KC,ISPEC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF((NSBCXR.EQ.NSBCW2).OR.(NSBCXR.EQ.NSBCW1))THEN
        DO ISPEC=1,NSPEC
          DO KC=KSTAL,KSTOL
            DO JC=JSTAL,JSTOL
              TEMP1=YRHS(ISTOL-1,JC,KC,ISPEC)/DRHS(ISTOL-1,JC,KC)
              TEMP2=YRHS(ISTOL-2,JC,KC,ISPEC)/DRHS(ISTOL-2,JC,KC)
              TEMP3=YRHS(ISTOL-3,JC,KC,ISPEC)/DRHS(ISTOL-3,JC,KC)
              TEMP4=YRHS(ISTOL-4,JC,KC,ISPEC)/DRHS(ISTOL-4,JC,KC)
              YRUN(ISTOL,JC,KC,ISPEC)=(12.0/25.0)*(4.0*TEMP1-3.0*TEMP2+
     +                           (4.0/3.0)*TEMP3-(1.0/4.0*TEMP4))
              YRUN(ISTOL,JC,KC,ISPEC)=MIN(1.0,YRUN(ISTOL,JC,KC,ISPEC))
              YRUN(ISTOL,JC,KC,ISPEC)=MAX(0.0,YRUN(ISTOL,JC,KC,ISPEC))
              YRUN(ISTOL,JC,KC,ISPEC)=DRHS(ISTOL,JC,KC)
     +                         *YRUN(ISTOL,JC,KC,ISPEC)
              YRHS(ISTOL,JC,KC,ISPEC)=YRUN(ISTOL,JC,KC,ISPEC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF((NSBCYL.EQ.NSBCW2).OR.(NSBCYL.EQ.NSBCW1))THEN
        DO ISPEC=1,NSPEC
          DO KC=KSTAL,KSTOL
            DO IC=ISTAL,ISTOL
              TEMP1=YRHS(IC,JSTAL+1,KC,ISPEC)/DRHS(IC,JSTAL+1,KC)
              TEMP2=YRHS(IC,JSTAL+2,KC,ISPEC)/DRHS(IC,JSTAL+2,KC)
              TEMP3=YRHS(IC,JSTAL+3,KC,ISPEC)/DRHS(IC,JSTAL+3,KC)
              TEMP4=YRHS(IC,JSTAL+4,KC,ISPEC)/DRHS(IC,JSTAL+4,KC)
              YRUN(IC,JSTAL,KC,ISPEC)=(12.0/25.0)*(4.0*TEMP1-3.0*TEMP2+
     +                           (4.0/3.0)*TEMP3-(1.0/4.0*TEMP4))
              YRUN(IC,JSTAL,KC,ISPEC)=MIN(1.0,YRUN(IC,JSTAL,KC,ISPEC))
              YRUN(IC,JSTAL,KC,ISPEC)=MAX(0.0,YRUN(IC,JSTAL,KC,ISPEC))
              YRUN(IC,JSTAL,KC,ISPEC)=DRHS(IC,JSTAL,KC)
     +                               *YRUN(IC,JSTAL,KC,ISPEC)
              YRHS(IC,JSTAL,KC,ISPEC)=YRUN(IC,JSTAL,KC,ISPEC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF((NSBCYR.EQ.NSBCW2).OR.(NSBCYR.EQ.NSBCW1))THEN
        DO ISPEC=1,NSPEC
          DO KC=KSTAL,KSTOL
            DO IC=ISTAL,ISTOL
              TEMP1=YRHS(IC,JSTOL-1,KC,ISPEC)/DRHS(IC,JSTOL-1,KC)
              TEMP2=YRHS(IC,JSTOL-2,KC,ISPEC)/DRHS(IC,JSTOL-2,KC)
              TEMP3=YRHS(IC,JSTOL-3,KC,ISPEC)/DRHS(IC,JSTOL-3,KC)
              TEMP4=YRHS(IC,JSTOL-4,KC,ISPEC)/DRHS(IC,JSTOL-4,KC)
              YRUN(IC,JSTOL,KC,ISPEC)=(12.0/25.0)*(4.0*TEMP1-3.0*TEMP2+
     +                           (4.0/3.0)*TEMP3-(1.0/4.0*TEMP4))
              YRUN(IC,JSTOL,KC,ISPEC)=MIN(1.0,YRUN(IC,JSTOL,KC,ISPEC))
              YRUN(IC,JSTOL,KC,ISPEC)=MAX(0.0,YRUN(IC,JSTOL,KC,ISPEC))
              YRUN(IC,JSTOL,KC,ISPEC)=DRHS(IC,JSTOL,KC)
     +                               *YRUN(IC,JSTOL,KC,ISPEC)
              YRHS(IC,JSTOL,KC,ISPEC)=YRUN(IC,JSTOL,KC,ISPEC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF((NSBCZL.EQ.NSBCW2).OR.(NSBCZL.EQ.NSBCW1))THEN
        DO ISPEC=1,NSPEC
          DO JC=JSTAL,JSTOL
            DO IC=ISTAL,ISTOL
              TEMP1=YRHS(IC,JC,KSTAL+1,ISPEC)/DRHS(IC,JC,KSTAL+1)
              TEMP2=YRHS(IC,JC,KSTAL+2,ISPEC)/DRHS(IC,JC,KSTAL+2)
              TEMP3=YRHS(IC,JC,KSTAL+3,ISPEC)/DRHS(IC,JC,KSTAL+3)
              TEMP4=YRHS(IC,JC,KSTAL+4,ISPEC)/DRHS(IC,JC,KSTAL+4)
              YRUN(IC,JC,KSTAL,ISPEC)=(12.0/25.0)*(4.0*TEMP1-3.0*TEMP2+
     +                           (4.0/3.0)*TEMP3-(1.0/4.0*TEMP4))
              YRUN(IC,JC,KSTAL,ISPEC)=MIN(1.0,YRUN(IC,JC,KSTAL,ISPEC))
              YRUN(IC,JC,KSTAL,ISPEC)=MAX(0.0,YRUN(IC,JC,KSTAL,ISPEC))
              YRUN(IC,JC,KSTAL,ISPEC)=DRHS(IC,JC,KSTAL)
     +                               *YRUN(IC,JC,KSTAL,ISPEC)
              YRHS(IC,JC,KSTAL,ISPEC)=YRUN(IC,JC,KSTAL,ISPEC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF((NSBCZR.EQ.NSBCW2).OR.(NSBCZR.EQ.NSBCW1))THEN
        DO ISPEC=1,NSPEC
          DO JC=JSTAL,JSTOL
            DO IC=ISTAL,ISTOL
              TEMP1=YRHS(IC,JC,KSTOL-1,ISPEC)/DRHS(IC,JC,KSTOL-1)
              TEMP2=YRHS(IC,JC,KSTOL-2,ISPEC)/DRHS(IC,JC,KSTOL-2)
              TEMP3=YRHS(IC,JC,KSTOL-3,ISPEC)/DRHS(IC,JC,KSTOL-3)
              TEMP4=YRHS(IC,JC,KSTOL-4,ISPEC)/DRHS(IC,JC,KSTOL-4)
              YRUN(IC,JC,KSTOL,ISPEC)=(12.0/25.0)*(4.0*TEMP1-3.0*TEMP2+
     +                           (4.0/3.0)*TEMP3-(1.0/4.0*TEMP4))
              YRUN(IC,JC,KSTOL,ISPEC)=MIN(1.0,YRUN(IC,JC,KSTOL,ISPEC))
              YRUN(IC,JC,KSTOL,ISPEC)=MAX(0.0,YRUN(IC,JC,KSTOL,ISPEC))
              YRUN(IC,JC,KSTOL,ISPEC)=DRHS(IC,JC,KSTOL)
     +                               *YRUN(IC,JC,KSTOL,ISPEC)
              YRHS(IC,JC,KSTOL,ISPEC)=YRUN(IC,JC,KSTOL,ISPEC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C     -------------------------------------------------------------------------

CC     NTH SPECIES
C      DO KC = KSTALY,KSTOLY
C        DO JC = JSTALY,JSTOLY
C          DO IC = ISTALY,ISTOLY
C
C            YRUN(IC,JC,KC,NSPEC) = ZERO
C            YRHS(IC,JC,KC,NSPEC) = ZERO
C
C          ENDDO
C        ENDDO
C      ENDDO
C
C      DO ISPEC = 1,NSPM1
C        DO KC = KSTALY,KSTOLY
C          DO JC = JSTALY,JSTOLY
C            DO IC = ISTALY,ISTOLY
C
C              YRUN(IC,JC,KC,NSPEC) = YRUN(IC,JC,KC,NSPEC)
C     +                             + YRUN(IC,JC,KC,ISPEC)
C              YRHS(IC,JC,KC,NSPEC) = YRHS(IC,JC,KC,NSPEC)
C     +                             + YRHS(IC,JC,KC,ISPEC)
C
C            ENDDO
C          ENDDO
C        ENDDO
C      ENDDO
C
C      DO KC = KSTALY,KSTOLY
C        DO JC = JSTALY,JSTOLY
C          DO IC = ISTALY,ISTOLY
C
C            YRUN(IC,JC,KC,NSPEC)
C     +        = DRUN(IC,JC,KC)*(ONE-YRUN(IC,JC,KC,NSPEC)/DRUN(IC,JC,KC))
C
C            YRHS(IC,JC,KC,NSPEC)
C     +        = DRHS(IC,JC,KC)*(ONE-YRHS(IC,JC,KC,NSPEC)/DRHS(IC,JC,KC))
C
C          ENDDO
C        ENDDO
C      ENDDO

C     =========================================================================


      RETURN
      END
