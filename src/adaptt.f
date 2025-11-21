      SUBROUTINE ADAPTT
 
C     *************************************************************************
C
C     ADAPTT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     19-JAN-2003:  CREATED
C     23-AUG-2009:  RSC REVISE ERROR NORM EVALUATION
C     08-AUG-2012:  RSC EVALUATE ALL SPECIES
C     09-AUG-2012   RSC/RACG USE GLOBAL ERROR
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES NEW TIMESTEP FOR ERK SCHEME
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ERYTOT(NSPCMX)
      DOUBLE PRECISION ERDTOT,ERUTOT,ERVTOT,ERWTOT,ERETOT
      DOUBLE PRECISION ERRMAX,TRATIO,TSTOLD
C     RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
C      DOUBLE PRECISION TSTLOC
      DOUBLE PRECISION ERRLOC
      DOUBLE PRECISION FORNOW
      INTEGER IC,JC,KC,ISPEC


C     BEGIN
C     =====

C     =========================================================================

C     CHECK ADAPTION FLAG
C     -------------------
      IF(FLADPT)THEN

C       =======================================================================

C       INITIALISE THE ERROR NORM TOTALS
C       --------------------------------
        ERDTOT = ZERO
        ERUTOT = ZERO
        ERVTOT = ZERO
        ERWTOT = ZERO
        ERETOT = ZERO
C       RSC 08-AUG-2012 EVALUATE ALL SPECIES
C        DO ISPEC = 1,NSPM1
        DO ISPEC = 1,NSPEC
          ERYTOT(ISPEC) = ZERO
        ENDDO

C       =======================================================================

C       ERK ERROR EVALUATION
C       --------------------
C       USING ERK ERROR ARRAYS
C       RSC 23 AUG-2009 REVISE ERROR NORM EVALUATION

C       EVALUATE ERROR NORMS
        DO KC = KSTALD,KSTOLD
          DO JC = JSTALD,JSTOLD
            DO IC = ISTALD,ISTOLD

C              FORNOW = ABS(DERR(IC,JC,KC))
              FORNOW = ABS(DERR(IC,JC,KC))/(ABS(DRUN(IC,JC,KC))+ERDNRM)
              IF(FORNOW.GT.ERDTOT)ERDTOT = FORNOW

            ENDDO
          ENDDO
        ENDDO

        DO KC = KSTALU,KSTOLU
          DO JC = JSTALU,JSTOLU
            DO IC = ISTALU,ISTOLU

C              FORNOW = ABS(UERR(IC,JC,KC))
              FORNOW = ABS(UERR(IC,JC,KC))/(ABS(URUN(IC,JC,KC))+ERUNRM)
              IF(FORNOW.GT.ERUTOT)ERUTOT = FORNOW

            ENDDO
          ENDDO
        ENDDO

        DO KC = KSTALV,KSTOLV
          DO JC = JSTALV,JSTOLV
            DO IC = ISTALV,ISTOLV

C              FORNOW = ABS(VERR(IC,JC,KC))
              FORNOW = ABS(VERR(IC,JC,KC))/(ABS(VRUN(IC,JC,KC))+ERVNRM)
              IF(FORNOW.GT.ERVTOT)ERVTOT = FORNOW

            ENDDO
          ENDDO
        ENDDO

        DO KC = KSTALW,KSTOLW
          DO JC = JSTALW,JSTOLW
            DO IC = ISTALW,ISTOLW

C              FORNOW = ABS(WERR(IC,JC,KC))
              FORNOW = ABS(WERR(IC,JC,KC))/(ABS(WRUN(IC,JC,KC))+ERWNRM)
              IF(FORNOW.GT.ERWTOT)ERWTOT = FORNOW

            ENDDO
          ENDDO
        ENDDO

        DO KC = KSTALE,KSTOLE
          DO JC = JSTALE,JSTOLE
            DO IC = ISTALE,ISTOLE

C              FORNOW = ABS(EERR(IC,JC,KC))
              FORNOW = ABS(EERR(IC,JC,KC))/(ABS(ERUN(IC,JC,KC))+ERENRM)
              IF(FORNOW.GT.ERETOT)ERETOT = FORNOW

            ENDDO
          ENDDO
        ENDDO

C       RSC 08-AUG-2012 EVALUATE ALL SPECIES
C        DO ISPEC = 1,NSPM1
        DO ISPEC = 1,NSPEC

          DO KC = KSTALY,KSTOLY
            DO JC = JSTALY,JSTOLY
              DO IC = ISTALY,ISTOLY

C                FORNOW = ABS(YERR(IC,JC,KC,ISPEC))
                FORNOW = ABS(YERR(IC,JC,KC,ISPEC))
     +                 /(ABS(YRUN(IC,JC,KC,ISPEC))+ERYNRM(ISPEC))
                IF(FORNOW.GT.ERYTOT(ISPEC))ERYTOT(ISPEC) = FORNOW

              ENDDO
            ENDDO
          ENDDO

        ENDDO

C       =======================================================================

CC       NORMALISE THE ERROR NORMS
CC       -------------------------
C        ERDTOT = ABS(ERDTOT)*ERDNRM
C        ERUTOT = ABS(ERUTOT)*ERUNRM
C        ERVTOT = ABS(ERVTOT)*ERVNRM
C        ERWTOT = ABS(ERWTOT)*ERWNRM
C        ERETOT = ABS(ERETOT)*ERENRM
C        DO ISPEC = 1,NSPM1
C          ERYTOT(ISPEC) = ABS(ERYTOT(ISPEC))*ERYNRM(ISPEC)
C        ENDDO

C       =======================================================================

C       FIND THE MAXIMUM 
C       ----------------
        ERRMAX = ZERO
        IF(ERDTOT.GT.ERRMAX)THEN
          ERRMAX = ERDTOT
          INDERR = -4
        ENDIF
        IF(ERUTOT.GT.ERRMAX)THEN
          ERRMAX = ERUTOT
          INDERR = -1
        ENDIF
        IF(ERVTOT.GT.ERRMAX)THEN
          ERRMAX = ERVTOT
          INDERR = -2
        ENDIF
        IF(ERWTOT.GT.ERRMAX)THEN
          ERRMAX = ERWTOT
          INDERR = -3
        ENDIF
        IF(ERETOT.GT.ERRMAX)THEN
          ERRMAX = ERETOT
          INDERR = 0
        ENDIF
C       RSC 08-AUG-2012 EVALUATE ALL SPECIES
C        DO ISPEC = 1,NSPM1
        DO ISPEC = 1,NSPEC
          IF(ERYTOT(ISPEC).GT.ERRMAX)THEN
            ERRMAX = ERYTOT(ISPEC)
            INDERR = ISPEC
          ENDIF
        ENDDO

C       =======================================================================

C       FIND THE LARGEST GLOBAL ERROR
C       -----------------------------
C       RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
        ERRLOC = ERRMAX
        CALL P_GMAX(ERRLOC,ERRMAX)

C       =======================================================================

C       EVALUATE THE NEW TIME STEP
C       --------------------------
C       ZERO CHECK
        IF(ERRMAX.LT.ERRLOW)ERRMAX = ERRLOW

C       -----------------------------------------------------------------------

CC       I-CONTROLLER
C        TRATIO = CTMULT*EXP(CTALPH*LOG(ERRTOL/ERRMAX))

C       -----------------------------------------------------------------------

CC       PI-CONTROLLER
C        TRATIO = CTMULT*EXP(CTALPH*LOG(ERRTOL/ERRMAX)
C     +                     +CTBETA*LOG(ERROLD/ERRTOL))
C        ERROLD = ERRMAX

C       -----------------------------------------------------------------------

C       PID-CONTROLLER
        TRATIO = CTMULT*EXP(CTALPH*LOG(ERRTOL/ERRMAX)
     +                     +CTBETA*LOG(ERROLD/ERRTOL)
     +                     +CTGAMA*LOG(ERRTOL/ERRLDR))
        ERRLDR = ERROLD
        ERROLD = ERRMAX

C       -----------------------------------------------------------------------

C       LIMIT CHANGES TO TIME STEP
        IF(TRATIO.GT.TRMAX)TRATIO = TRMAX
        IF(TRATIO.LT.TRMIN)TRATIO = TRMIN

C       -----------------------------------------------------------------------

C       SAVE THE OLD TIME STEP
        TSTOLD = TSTEP

C       SET THE NEW TIME STEP
        TSTEP = TSTEP*TRATIO

C       LIMIT THE TIME STEP
        IF(TSTEP.GT.TSMAX)TSTEP = TSMAX
        IF(TSTEP.LT.TSMIN)TSTEP = TSMIN

C       =======================================================================

C       PARALLEL TRANSFER TO SET NEW GLOBAL TIME STEP
C       ---------------------------------------------
C       NEW TIME STEP IS THE GLOBAL MINIMUM OVER ALL PROCESSORS
C       RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
C        TSTLOC = TSTEP
C        CALL P_GMIN(TSTLOC,TSTEP)
        
C       =======================================================================

C       UPDATE THE TIME ADVANCEMENT COEFFICIENTS
C       AND THE RK SUBSTEP TIME LEVELS
        TRATIO = TSTEP/TSTOLD
        DO IRKSTP = 1, NRKSTP
          RKLHS(IRKSTP) = RKLHS(IRKSTP)*TRATIO
          RKRHS(IRKSTP) = RKRHS(IRKSTP)*TRATIO
          RKERR(IRKSTP) = RKERR(IRKSTP)*TRATIO
          RKTIM(IRKSTP) = RKTIM(IRKSTP)*TRATIO
        ENDDO

C       =======================================================================

C       (RE)INITIALISE ERK ERROR ARRAYS
C       -------------------------------
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DERR(IC,JC,KC) = ZERO
              UERR(IC,JC,KC) = ZERO
              VERR(IC,JC,KC) = ZERO
              WERR(IC,JC,KC) = ZERO
              EERR(IC,JC,KC) = ZERO

            ENDDO
          ENDDO
        ENDDO
C       RSC 08-AUG-2012 EVALUATE ALL SPECIES
C        DO ISPEC = 1,NSPM1
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YERR(IC,JC,KC,ISPEC) = ZERO

              ENDDO
            ENDDO
          ENDDO

        ENDDO

C       =======================================================================

C       (RE)INITIALISE ERK SUBSTEP ERROR NORMS
C       --------------------------------------
        DO IRKSTP = 1,NRKSTP

          ERDRHS(IRKSTP) = ZERO
          ERURHS(IRKSTP) = ZERO
          ERVRHS(IRKSTP) = ZERO
          ERWRHS(IRKSTP) = ZERO
          ERERHS(IRKSTP) = ZERO
C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC
            ERYRHS(ISPEC,IRKSTP) = ZERO
          ENDDO

        ENDDO

C       =======================================================================

      ENDIF
C     ADAPTION FLAG

C     =========================================================================


      RETURN
      END
