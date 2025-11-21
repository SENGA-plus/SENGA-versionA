      SUBROUTINE DTINIT
 
C     *************************************************************************
C
C     DTINIT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT -- CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     18-MAY-2003:  CREATED
C     07-JUL-2009:  RSC BUG FIX ERROR NORMS
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES TIME-STEPPING FOR ERK SCHEME
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ACOFRK(NRKMAX),BCOFRK(NRKMAX),BHATRK(NRKMAX)
      INTEGER IC,JC,KC
      INTEGER ISPEC


C     BEGIN
C     =====

C     =========================================================================

C     VALUES DEFINED IN THIS BLOCK FOR

C     ERK COEFFICIENTS
C     ERK ERROR ESTIMATION NORMALISING VALUES
C     ERK TIME STEP CONTROLLER COEFFICIENTS
C     ERK TOLERANCES AND LIMITS

CC     -------------------------------------------------------------------------
CC
CC     ERK SCHEME RK3(2)4[2R+]C
CC     ------------------------
C      NRKSTP = 4
C
C      ACOFRK(1) = 11847461282814D0/36547543011857D0
C      ACOFRK(2) = 3943225443063D0/7078155732230D0
C      ACOFRK(3) = -346793006927D0/4029903576067D0
C      ACOFRK(4) = ZERO
C
C      BCOFRK(1) = 1017324711453D0/9774461848756D0
C      BCOFRK(2) = 8237718856693D0/13685301971492D0
C      BCOFRK(3) = 57731312506979D0/19404895981398D0
C      BCOFRK(4) = -101169746363290D0/37734290219643D0
C
C      BHATRK(1) = 15763415370699D0/46270243929542D0
C      BHATRK(2) = 514528521746D0/5659431552419D0
C      BHATRK(3) = 27030193851939D0/9429696342944D0
C      BHATRK(4) = -69544964788955D0/30262026368149D0
C
C     -------------------------------------------------------------------------

C     ERK SCHEME RK4(3)5[2R+]C
C     ------------------------
      NRKSTP = 5

      ACOFRK(1) = 970286171893D0/4311952581923D0
      ACOFRK(2) = 6584761158862D0/12103376702013D0
      ACOFRK(3) = 2251764453980D0/15575788980749D0
      ACOFRK(4) = 26877169314380D0/34165994151039D0
      ACOFRK(5) = ZERO

      BCOFRK(1) = 1153189308089D0/22510343858157D0
      BCOFRK(2) = 1772645290293D0/4653164025191D0
      BCOFRK(3) = -1672844663538D0/4480602732383D0
      BCOFRK(4) = 2114624349019D0/3568978502595D0
      BCOFRK(5) = 5198255086312D0/14908931495163D0

      BHATRK(1) = 1016888040809D0/7410784769900D0
      BHATRK(2) = 11231460423587D0/58533540763752D0
      BHATRK(3) = -1563879915014D0/6823010717585D0
      BHATRK(4) = 606302364029D0/971179775848D0
      BHATRK(5) = 1097981568119D0/3980877426909D0

C     -------------------------------------------------------------------------

C     SET TOLERANCES AND LIMITS FOR THE RK TIME STEP CONTROLLER
C     ---------------------------------------------------------
C     RSC 07-JUL-2009 BUG FIX ERROR NORMS
C      ERRTOL = 1.0D-3
      ERRTOL = 1.0D-4
      ERRLOW = 1.0D-30
      IF(NCDMPI.EQ.0)THEN
        ERROLD = ERRTOL
        ERRLDR = ERRTOL
      ENDIF
      TRMAX = 1.01D0
      TRMIN = 1.0D-2
      TSMAX = 1.0D0
      TSMIN = 1.0D-15

C     -------------------------------------------------------------------------

C     INITIALISE THE NORMALISING VALUES FOR RK ERROR ESTIMATION
C     ---------------------------------------------------------
      ERDNRM = DRIN
      ERUNRM = ONE
      IF(ABS(URIN).GT.ERRTOL)ERUNRM = DRIN*ABS(URIN)
      ERVNRM = ONE
      IF(ABS(VRIN).GT.ERRTOL)ERVNRM = DRIN*ABS(VRIN)
      ERWNRM = ONE
      IF(ABS(WRIN).GT.ERRTOL)ERWNRM = DRIN*ABS(WRIN)
      ERENRM = ONE
      IF(ABS(ERIN).GT.ERRTOL)ERENRM = DRIN*ABS(ERIN)
      DO ISPEC = 1,NSPEC
        ERYNRM(ISPEC) = ONE
      ENDDO

C     RSC 07-JUL-2009 BUG FIX ERROR NORMS
      ERDNRM = ONE/ERDNRM
      ERUNRM = ONE/ERUNRM
      ERVNRM = ONE/ERVNRM
      ERWNRM = ONE/ERWNRM
      ERENRM = ONE/ERENRM

C     RSC 23-AUG-2009 REFORMULATE ERROR NORMS
      ERDNRM = ZERO
      ERUNRM = 1.0D-6
      ERVNRM = 1.0D-6
      ERWNRM = 1.0D-6
      ERENRM = 1.0D-2
      DO ISPEC = 1,NSPEC
        ERYNRM(ISPEC) = 1.0D-10
      ENDDO

C     -------------------------------------------------------------------------

C     SET THE COEFFICIENTS FOR THE RK TIME STEP CONTROLLER
C     ----------------------------------------------------
C     I CONTROLLER
C     CTALPH = 1.0/ORDER OF EMBEDDED SCHEME+1
      CTMULT = 9.0D-1
      CTALPH = ONE/FOUR
C     PI CONTROLLER
C     CTALPH = 0.7/ORDER OF EMBEDDED SCHEME
C     CTBETA = 0.4/ORDER OF EMBEDDED SCHEME
      CTMULT = 9.0D-1
      CTALPH = 7.0D-1/THREE
      CTBETA = 4.0D-1/THREE
C     PID CONTROLLER
C     CTALPH = 0.49/ORDER OF EMBEDDED SCHEME
C     CTBETA = 0.34/ORDER OF EMBEDDED SCHEME
C     CTGAMA = 0.10/ORDER OF EMBEDDED SCHEME
      CTMULT = 9.0D-1
      CTALPH = 4.9D-1/THREE
      CTBETA = 3.4D-1/THREE
      CTGAMA = 1.0D-1/THREE

C     -------------------------------------------------------------------------

C     SET TIME STEP ADAPTION FLAG
C     ---------------------------
      FLADPT = .TRUE.
      IF(NSTPSW.EQ.0)FLADPT = .FALSE.

C     -------------------------------------------------------------------------

C     END OF BLOCK FOR

C     ERK COEFFICIENTS
C     ERK ERROR ESTIMATION NORMALISING VALUES
C     ERK TIME STEP CONTROLLER COEFFICIENTS
C     ERK TOLERANCES AND LIMITS

C     =========================================================================

C     NO OF SUBSTEPS MINUS ONE
C     ------------------------
      NRKSM1 = NRKSTP - 1

C     =========================================================================

C     RK COEFFICIENTS
C     ---------------
      DO IRKSTP = 1,NRKSTP
        RKLHS(IRKSTP) = BCOFRK(IRKSTP)
        RKRHS(IRKSTP) = ACOFRK(IRKSTP)
        RKERR(IRKSTP) = BCOFRK(IRKSTP) - BHATRK(IRKSTP)
      ENDDO

C     =========================================================================

C     TIME-ADVANCEMENT COEFFICIENTS
C     -----------------------------
      DO IRKSTP = 1,NRKSTP
        RKLHS(IRKSTP) = RKLHS(IRKSTP)*TSTEP
        RKRHS(IRKSTP) = RKRHS(IRKSTP)*TSTEP
        RKERR(IRKSTP) = RKERR(IRKSTP)*TSTEP
      ENDDO

C     =========================================================================

C     TIME LEVELS FOR EACH RK SUBSTEP
C     -------------------------------
      RKTIM(1) = ZERO
      DO IRKSTP = 2,NRKSTP
        RKTIM(IRKSTP) = ACOFRK(IRKSTP-1)
        DO JC = 1,IRKSTP-2
          RKTIM(IRKSTP) = RKTIM(IRKSTP) + BCOFRK(JC)
        ENDDO
        RKTIM(IRKSTP) = RKTIM(IRKSTP)*TSTEP
      ENDDO

C     =========================================================================

C     INITIALISE ERK ERROR ARRAYS
C     ---------------------------
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
      DO ISPEC = 1,NSPEC

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              YERR(IC,JC,KC,ISPEC) = ZERO

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     =========================================================================

C     INITIALISE ERK SUBSTEP ERROR NORMS
C     ----------------------------------
      DO IRKSTP = 1,NRKSTP

        ERDRHS(IRKSTP) = ZERO
        ERURHS(IRKSTP) = ZERO
        ERVRHS(IRKSTP) = ZERO
        ERWRHS(IRKSTP) = ZERO
        ERERHS(IRKSTP) = ZERO
        DO ISPEC = 1,NSPEC
          ERYRHS(ISPEC,IRKSTP) = ZERO
        ENDDO

      ENDDO

C     =========================================================================


      RETURN
      END
