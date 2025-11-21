      SUBROUTINE INTEGF(FUNCTN,ARGMIN,ARGMAX,ANSWER)

C     *************************************************************************
C
C     INTEGF
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     20-APR-1997:  CREATED
C     24-MAY-2003:  UPDATED FOR SENGA2
C     15-SEP-2012:  RSC IMPROVE INITIAL SAMPLING OF FUNCTION
C
C     DESCRIPTION
C     -----------
C     INTEGRATES THE SPECIFIED FUNCTION BETWEEN THE SPECIFIED LIMITS
C     USING SUCCESSIVELY-REFINED SIMPSON'S RULE
C
C     REFERENCE
C     ---------
C     NUMERICAL RECIPES 1st Ed, CUP (1986), pp110-113
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION TOLINT,VSMALL,ZERO,HALF,THREE,FOUR
      INTEGER INTMAX
      PARAMETER(TOLINT=1.0D-6, INTMAX=20, VSMALL = 1.0D-30)
      PARAMETER(ZERO=0.0D0, HALF=0.5D0, THREE=3.0D0, FOUR=4.0D0)


C     ARGUMENTS
C     =========
      DOUBLE PRECISION FUNCTN,ARGMIN,ARGMAX,ANSWER
      EXTERNAL FUNCTN


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION OLDTMP,OLDANS,TMPANS,ARGTMP,DELTRG,ADDPTS
      INTEGER ICOUNT,IC,INTERV


C     BEGIN
C     =====

C     =========================================================================

C     RSC 15-SEP-2012 IMPROVE INITIAL SAMPLING OF FUNCTION
C      INTERV = 1
      INTERV = 4
      ICOUNT = 1
      ANSWER = VSMALL
      TMPANS = HALF*(ARGMAX-ARGMIN)*(FUNCTN(ARGMIN)+FUNCTN(ARGMAX))

C     =========================================================================

C     MAIN LOOP TO REFINE THE INTEGRAL
C     --------------------------------

1000  CONTINUE

        OLDANS = ANSWER
        OLDTMP = TMPANS

        DELTRG = (ARGMAX-ARGMIN)/REAL(INTERV)
        ARGTMP = ARGMIN-HALF*DELTRG

        ADDPTS = ZERO
        DO IC = 1, INTERV
          ARGTMP = ARGTMP + DELTRG
          ADDPTS = ADDPTS + FUNCTN(ARGTMP)
        ENDDO

        TMPANS = HALF*(TMPANS+ADDPTS*DELTRG)
        ANSWER = (FOUR*TMPANS-OLDTMP)/THREE

        IF(ABS(ANSWER-OLDANS).GT.(TOLINT*ABS(OLDANS)))THEN
          ICOUNT = ICOUNT + 1
          IF(ICOUNT.LE.INTMAX)THEN
            INTERV = 2*INTERV
            GO TO 1000
          ELSE
C           FINISH EVEN IF INTEGRAL IS NOT CONVERGED
C           RSC 15-SEP-2012 MINOR BUG FIX
            WRITE(6,*)'Warning: INTEGF: integral not converged'
            WRITE(6,*)'at iteration:',INTMAX
            WRITE(6,*)'with values',ANSWER,OLDANS
          ENDIF
        ENDIF

C     END OF MAIN LOOP
C     ----------------

C     =========================================================================


      RETURN
      END
