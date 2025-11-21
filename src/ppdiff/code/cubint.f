      FUNCTION CUBINT(XARRAY,YARRAY,NARRAY,XREQRD)

C     *************************************************************************
C
C     CUBINT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT, AFTER NUMERICAL RECIPES
C
C     CHANGE RECORD
C     -------------
C     23-JUL-1990: CREATED
C     22-JUL-1994: RSC MODIFIED FOR NEW RUNS AT CTR 1994
C     17-OCT-2012: RSC MODIFIED FOR MOLECULAR TRANSPORT EVALUATION
C
C     DESCRIPTION
C     -----------
C     EVALUATES THE KNOTS FOR A SPLINE FIT TO MAX. 261 DATA POINTS
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ZERO,ONE,TWO,SIX
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, SIX = 6.0D0)

      INTEGER NWRKMX
      PARAMETER(NWRKMX = 1000)


C     RETURN VALUE
C     ============
      DOUBLE PRECISION CUBINT


C     ARGUMENTS
C     =========
      DOUBLE PRECISION XARRAY(NARRAY),YARRAY(NARRAY)
      DOUBLE PRECISION XREQRD
      INTEGER NARRAY


C     LOCAL DATA
      DOUBLE PRECISION SKNOTS(NWRKMX),WORKSP(NWRKMX)
      DOUBLE PRECISION SUPPER,WUPPER,SIG,PVAL
      DOUBLE PRECISION DELTAX,ACOEFF,BCOEFF
      INTEGER ICOUNT
      INTEGER KCOUNT,KLO,KHI


C     BEGIN
C     =====

C     =========================================================================

C     CALCULATE THE KNOTS
C     -------------------

C     SET NATURAL BCS
      SKNOTS(1) = ZERO
      WORKSP(1) = ZERO
      SUPPER = ZERO
      WUPPER = ZERO

C     DECOMPOSITION LOOP OF TRID ALGORITHM
      DO 1000 ICOUNT = 2,NARRAY-1
        SIG = (XARRAY(ICOUNT)-XARRAY(ICOUNT-1))
     +       /(XARRAY(ICOUNT+1)-XARRAY(ICOUNT-1))
        PVAL = SIG*SKNOTS(ICOUNT-1)+TWO
        SKNOTS(ICOUNT) = (SIG-ONE)/PVAL
        WORKSP(ICOUNT) = (SIX*((YARRAY(ICOUNT+1)-YARRAY(ICOUNT))
     +                        /(XARRAY(ICOUNT+1)-XARRAY(ICOUNT))
     +                        -(YARRAY(ICOUNT)-YARRAY(ICOUNT-1))
     +                        /(XARRAY(ICOUNT)-XARRAY(ICOUNT-1)))
     +                        /(XARRAY(ICOUNT+1)-XARRAY(ICOUNT-1))
     +                        -SIG*WORKSP(ICOUNT-1))/PVAL
1000  CONTINUE

      SKNOTS(NARRAY) = (WUPPER-SUPPER*WORKSP(NARRAY-1))
     +                /(SUPPER*SKNOTS(NARRAY-1)+ONE)

C     BACKSUBSTITUTION
      DO 2000 ICOUNT = NARRAY-1,1,-1
        SKNOTS(ICOUNT) = SKNOTS(ICOUNT)*SKNOTS(ICOUNT+1)+WORKSP(ICOUNT)
2000  CONTINUE

C     =========================================================================

C     EVALUATE THE SPLINE
C     -------------------

C     BISECTION
      KLO = 1
      KHI = NARRAY

3000  CONTINUE
        IF(KHI-KLO.GT.1)THEN
          KCOUNT = (KHI+KLO)/2
          IF(XARRAY(KCOUNT).GT.XREQRD)THEN
            KHI = KCOUNT
          ELSE
            KLO = KCOUNT
          ENDIF
          GOTO 3000
        ENDIF
C     END OF BISECTION LOOP 3000

      DELTAX = XARRAY(KHI) - XARRAY(KLO)
      ACOEFF = (XARRAY(KHI)-XREQRD)/DELTAX
      BCOEFF = (XREQRD-XARRAY(KLO))/DELTAX

      CUBINT = ACOEFF*YARRAY(KLO) + BCOEFF*YARRAY(KHI)
     1       +((ACOEFF**3-ACOEFF)*SKNOTS(KLO)
     2       +(BCOEFF**3-BCOEFF)*SKNOTS(KHI))*(DELTAX**2)/SIX

C     =========================================================================


      RETURN
      END    
