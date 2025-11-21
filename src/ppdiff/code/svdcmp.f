      SUBROUTINE SVDCMP(AMATRX,MSIZE,NSIZE,MPHYS,NPHYS,WVECTR,VMATRX)

C     *************************************************************************
C
C     SVDCMP
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     28-OCT-2012:  CREATED
C
C     DESCRIPTION
C     -----------
C     COMPUTES THE SINGULAR VALUE DECOMPOSITION OF A MATRIX
C
C     IN THE FORM A = U.W.V^T
C     
C     INPUT:
C       MATRIX AMATRX(MSIZE,NSIZE) OF PHYSICAL DIMENSIONS MPHYS,NPHYS
C
C     OUTPUT:
C       MATRIX U(MSIZE,NSIZE) STORED IN AMATRX
C       MATRIX W = DIAG(NSIZE) STORED AS WVECTR OF PHYSICAL DIMENSION NPHYS
C       MATRIX V(NSIZE,NSIZE) (NOT THE TRANSPOSE) STORED IN VMATRX(NPHYS,NPHYS)
C
C     BASED ON NUMERICAL RECIPES SUBROUTINE SVDCMP
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      INTEGER NMAX
      PARAMETER(NMAX=500)


C     ARGUMENTS
C     =========
      DOUBLE PRECISION AMATRX(MPHYS,NPHYS)
      DOUBLE PRECISION VMATRX(NPHYS,NPHYS)
      DOUBLE PRECISION WVECTR(NPHYS)
      INTEGER MSIZE,MPHYS,NSIZE,NPHYS


C     EXTERNAL FUNCTION
C     =================
      DOUBLE PRECISION PYTHAG
      EXTERNAL PYTHAG


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION RV1(NMAX)
      DOUBLE PRECISION ANORM,CVAR,FVAR,GVAR,HVAR,SVAR,SSCALE
      DOUBLE PRECISION XVAR,YVAR,ZVAR
      INTEGER IC,ITS,JC,JJ,KC,LC,NM


C     BEGIN
C     =====
C     HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM
      GVAR = ZERO
      SSCALE = ZERO
      ANORM = ZERO
      DO IC = 1,NSIZE
        LC = IC + 1
        RV1(IC) = SSCALE*GVAR
        GVAR = ZERO
        SVAR = ZERO
        SSCALE = ZERO
        IF(IC.LE.MSIZE)THEN
          DO KC = IC,MSIZE
            SSCALE = SSCALE + ABS(AMATRX(KC,IC))
          ENDDO
          IF(SSCALE.NE.ZERO)THEN
            DO KC = IC,MSIZE
              AMATRX(KC,IC) = AMATRX(KC,IC)/SSCALE
              SVAR = SVAR + AMATRX(KC,IC)*AMATRX(KC,IC)
            ENDDO
            FVAR = AMATRX(IC,IC)
            GVAR = -SIGN(SQRT(SVAR),FVAR)
            HVAR = FVAR*GVAR - SVAR
            AMATRX(IC,IC) = FVAR - GVAR
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = IC,MSIZE
                SVAR = SVAR + AMATRX(KC,IC)*AMATRX(KC,JC)
              ENDDO
              FVAR = SVAR/HVAR
              DO KC = IC,MSIZE
                AMATRX(KC,JC) = AMATRX(KC,JC) + FVAR*AMATRX(KC,IC)
              ENDDO
            ENDDO
            DO KC = IC,MSIZE
              AMATRX(KC,IC) = SSCALE*AMATRX(KC,IC)
            ENDDO
          ENDIF
        ENDIF
        WVECTR(IC) = SSCALE*GVAR
        GVAR = ZERO
        SVAR = ZERO
        SSCALE = ZERO
        IF((IC.LE.MSIZE).AND.(IC.NE.NSIZE))THEN
          DO KC = LC,NSIZE
            SSCALE = SSCALE + ABS(AMATRX(IC,KC))
          ENDDO
          IF(SSCALE.NE.ZERO)THEN
            DO KC = LC,NSIZE
              AMATRX(IC,KC) = AMATRX(IC,KC)/SSCALE
              SVAR = SVAR + AMATRX(IC,KC)*AMATRX(IC,KC)
            ENDDO
            FVAR = AMATRX(IC,LC)
            GVAR = -SIGN(SQRT(SVAR),FVAR)
            HVAR = FVAR*GVAR - SVAR
            AMATRX(IC,LC) = FVAR - GVAR
            DO KC = LC,NSIZE
              RV1(KC)=AMATRX(IC,KC)/HVAR
            ENDDO
            DO JC = LC,MSIZE
              SVAR = ZERO
              DO KC = LC,NSIZE
                SVAR = SVAR + AMATRX(JC,KC)*AMATRX(IC,KC)
              ENDDO
              DO KC = LC,NSIZE
                AMATRX(JC,KC) = AMATRX(JC,KC) + SVAR*RV1(KC)
              ENDDO
            ENDDO
            DO KC = LC,NSIZE
              AMATRX(IC,KC)=SSCALE*AMATRX(IC,KC)
            ENDDO
          ENDIF
        ENDIF
        ANORM = MAX(ANORM,(ABS(WVECTR(IC))+ABS(RV1(IC))))
      ENDDO

C     ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS
      DO IC = NSIZE,1,-1
        IF(IC.LT.NSIZE)THEN
          IF(GVAR.NE.ZERO)THEN
C           DOUBLE DIVISION TO AVOID POSSIBLE UNDERFLOW
            DO JC = LC,NSIZE
              VMATRX(JC,IC) = (AMATRX(IC,JC)/AMATRX(IC,LC))/GVAR
            ENDDO
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = LC,NSIZE
                SVAR = SVAR + AMATRX(IC,KC)*VMATRX(KC,JC)
              ENDDO
              DO KC = LC,NSIZE
                VMATRX(KC,JC) = VMATRX(KC,JC) + SVAR*VMATRX(KC,IC)
              ENDDO
            ENDDO
          ENDIF
          DO JC = LC,NSIZE
            VMATRX(IC,JC) = ZERO
            VMATRX(JC,IC) = ZERO
          ENDDO
        ENDIF
        VMATRX(IC,IC) = ONE
        GVAR = RV1(IC)
        LC = IC
      ENDDO

C     ACCUMULATION OF LEFT-HAND TRANSFORMATIONS
      DO IC = MIN(MSIZE,NSIZE),1,-1
        LC = IC + 1
        GVAR = WVECTR(IC)
        DO JC = LC,NSIZE
          AMATRX(IC,JC) = ZERO
        ENDDO
        IF(GVAR.NE.ZERO)THEN
          GVAR = ONE/GVAR
          DO JC=LC,NSIZE
            SVAR = ZERO
            DO KC = LC,MSIZE
              SVAR = SVAR + AMATRX(KC,IC)*AMATRX(KC,JC)
            ENDDO
            FVAR = (SVAR/AMATRX(IC,IC))*GVAR
            DO KC = IC,MSIZE
              AMATRX(KC,JC) = AMATRX(KC,JC) + FVAR*AMATRX(KC,IC)
            ENDDO
          ENDDO
          DO JC = IC,MSIZE
            AMATRX(JC,IC) = AMATRX(JC,IC)*GVAR
          ENDDO
        ELSE
          DO JC= IC,MSIZE
            AMATRX(JC,IC) = ZERO
          ENDDO
        ENDIF
        AMATRX(IC,IC) = AMATRX(IC,IC) + ONE
      ENDDO

C     DIAGONALIZATION OF THE BIDIAGONAL FORM
C     LOOP OVER SINGULAR VALUES
      NM = 1
      DO KC = NSIZE,1,-1
C       LOOP OVER ALLOWED ITERATIONS
        DO ITS = 1,30
          DO LC = KC,1,-1
C           TEST FOR SPLITTING
            NM = LC-1
C           NOTE THAT RV1(1) IS ALWAYS ZERO
            IF((ABS(RV1(LC))+ANORM).EQ.ANORM) GOTO 2000
            IF((ABS(WVECTR(NM))+ANORM).EQ.ANORM) GOTO 1000
          ENDDO

C         CANCELLATION OF RV1(LC) IF LC > 1
1000      CVAR = ZERO
          SVAR = ONE
          DO IC = LC,KC
            FVAR = SVAR*RV1(IC)
            RV1(IC) = CVAR*RV1(IC)
            IF((ABS(FVAR)+ANORM).EQ.ANORM) GOTO 2000
            GVAR = WVECTR(IC)
            HVAR = PYTHAG(FVAR,GVAR)
            WVECTR(IC) = HVAR
            HVAR = ONE/HVAR
            CVAR =  (GVAR*HVAR)
            SVAR = -(FVAR*HVAR)
            DO JC = 1,MSIZE
              YVAR = AMATRX(JC,NM)
              ZVAR = AMATRX(JC,IC)
              AMATRX(JC,NM) =  (YVAR*CVAR)+(ZVAR*SVAR)
              AMATRX(JC,IC) = -(YVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
          ENDDO

2000      ZVAR = WVECTR(KC)
          IF(LC.EQ.KC)THEN
C           CONVERGENCE
            IF(ZVAR.LT.ZERO)THEN
C             SINGULAR VALUE IS MADE NONNEGATIVE
              WVECTR(KC) = -ZVAR
              DO JC = 1,NSIZE
                VMATRX(JC,KC) = -VMATRX(JC,KC)
              ENDDO
            ENDIF
            GOTO 3000
          ENDIF
          IF(ITS.EQ.30)THEN
            WRITE(6,*)'SVDCMP: no convergence'
            STOP
          ENDIF

C         SHIFT FROM BOTTOM 2-BY-2 MINOR
          XVAR = WVECTR(LC)
          NM = KC-1
          YVAR = WVECTR(NM)
          GVAR = RV1(NM)
          HVAR = RV1(KC)
          FVAR = ((YVAR-ZVAR)*(YVAR+ZVAR)
     +         +  (GVAR-HVAR)*(GVAR+HVAR))/(TWO*HVAR*YVAR)
          GVAR = PYTHAG(FVAR,ONE)
          FVAR = ((XVAR-ZVAR)*(XVAR+ZVAR)
     +         +  HVAR*((YVAR/(FVAR+SIGN(GVAR,FVAR)))-HVAR))/XVAR
C         NEXT QR TRANSFORMATION
          CVAR = ONE
          SVAR = ONE
          DO JC = LC,NM
            IC = JC+1
            GVAR = RV1(IC)
            YVAR = WVECTR(IC)
            HVAR = SVAR*GVAR
            GVAR = CVAR*GVAR
            ZVAR = PYTHAG(FVAR,HVAR)
            RV1(JC) = ZVAR
            CVAR = FVAR/ZVAR
            SVAR = HVAR/ZVAR
            FVAR =  (XVAR*CVAR)+(GVAR*SVAR)
            GVAR = -(XVAR*SVAR)+(GVAR*CVAR)
            HVAR = YVAR*SVAR
            YVAR = YVAR*CVAR
            DO JJ = 1,NSIZE
              XVAR = VMATRX(JJ,JC)
              ZVAR = VMATRX(JJ,IC)
              VMATRX(JJ,JC) =  (XVAR*CVAR)+(ZVAR*SVAR)
              VMATRX(JJ,IC) = -(XVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
            ZVAR = PYTHAG(FVAR,HVAR)
            WVECTR(JC) = ZVAR
C           ROTATION CAN BE ARBITRARY IF ZVAR = 0
            IF(ZVAR.NE.ZERO)THEN
              ZVAR = ONE/ZVAR
              CVAR = FVAR*ZVAR
              SVAR = HVAR*ZVAR
            ENDIF
            FVAR =  (CVAR*GVAR)+(SVAR*YVAR)
            XVAR = -(SVAR*GVAR)+(CVAR*YVAR)
            DO JJ = 1,MSIZE
              YVAR = AMATRX(JJ,JC)
              ZVAR = AMATRX(JJ,IC)
              AMATRX(JJ,JC) =  (YVAR*CVAR)+(ZVAR*SVAR)
              AMATRX(JJ,IC) = -(YVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
          ENDDO
          RV1(LC) = ZERO
          RV1(KC) = FVAR
          WVECTR(KC) = XVAR
        ENDDO

3000    CONTINUE

      ENDDO


      RETURN
      END

