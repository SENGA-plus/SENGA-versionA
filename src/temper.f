      SUBROUTINE TEMPER
 
C     *************************************************************************
C
C     TEMPER
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     16-NOV-2002:  CREATED
C     
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES TEMPERATURE AND PRESSURE
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
      DOUBLE PRECISION TOLTMP
      PARAMETER(TOLTMP = 1.0E-10)
      INTEGER NTITRS
      PARAMETER(NTITRS = 100)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TCOEFF(0:NCTMAX),TDERIV(1:NCTMAX)
      DOUBLE PRECISION UKUK
      DOUBLE PRECISION TEMPOR,TFPOLY,TDPOLY,DELTMP,CPFORY
      INTEGER IC,JC,KC,ISPEC,ITINT,ICP,ITITRS
      INTEGER IINDEX,IPOWER,ICOEF1,ICOEF2
 

C     BEGIN
C     =====

C     =========================================================================

C     TEMPERATURE AND PRESSURE
C     ------------------------

C     TEMPERATURE AND PRESSURE ARE PARALLEL

      DO KC = KSTALT,KSTOLT
        DO JC = JSTALT,JSTOLT
          DO IC = ISTALT,ISTOLT

C           ===================================================================

C           KINETIC ENERGY TERM
            UKUK = (URHS(IC,JC,KC)*URHS(IC,JC,KC)
     +           +  VRHS(IC,JC,KC)*VRHS(IC,JC,KC)
     +           +  WRHS(IC,JC,KC)*WRHS(IC,JC,KC))/DRHS(IC,JC,KC)

C           ===================================================================

C           INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
C           AND ITS DERIVATIVE
            TCOEFF(0) = HALF*UKUK - ERHS(IC,JC,KC)
            DO ICP = 1, NCTMAX
              TCOEFF(ICP) = ZERO
              TDERIV(ICP) = ZERO
            ENDDO

C           ===================================================================

C           USE STORE7 TO ACCUMULATE MIXTURE SPECIFIC GAS CONSTANT
C           INITIALISE STORE7
            STORE7(IC,JC,KC) = ZERO

C           ===================================================================

C           RUN THROUGH ALL SPECIES
            DO ISPEC = 1,NSPEC

C             =================================================================

C             LOCATE TEMPERATURE IN AN INTERVAL
              IINDEX = 1 + (ISPEC-1)/NSPIMX
              IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
              ICOEF2 = NTBASE**IPOWER
              ICOEF1 = ICOEF2*NTBASE
              ITINT = 1 + MOD(ITNDEX(IC,JC,KC,IINDEX),ICOEF1)/ICOEF2

C             =================================================================

C             CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
              TCOEFF(0) = TCOEFF(0)
     +                  + YRHS(IC,JC,KC,ISPEC)*
     +                    AMASCP(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
              TCOEFF(1) = TCOEFF(1)
     +                  + YRHS(IC,JC,KC,ISPEC)*AMASCT(1,ITINT,ISPEC)
              TDERIV(1) = TCOEFF(1)
              DO ICP = 2, NCPOLY(ITINT,ISPEC)
                TCOEFF(ICP) = TCOEFF(ICP)
     +                      + YRHS(IC,JC,KC,ISPEC)*
     +                        AMASCT(ICP,ITINT,ISPEC)
                TDERIV(ICP) = TDERIV(ICP)
     +                      + YRHS(IC,JC,KC,ISPEC)*
     +                        AMASCP(ICP,ITINT,ISPEC)
              ENDDO

C             =================================================================

C             USE STORE7
C             TO ACCUMULATE (DENSITY TIMES) MIXTURE SPECIFIC GAS CONSTANT
              STORE7(IC,JC,KC) = STORE7(IC,JC,KC)
     +                         + YRHS(IC,JC,KC,ISPEC)*RGSPEC(ISPEC)

C             =================================================================

            ENDDO
C           END OF RUN THROUGH ALL SPECIES

C           ===================================================================

C           SOLVE FOR TEMPERATURE
C           USING NEWTON-RAPHSON
            TEMPOR = TRUN(IC,JC,KC)
            ITITRS = 1
1000        CONTINUE

C             EVALUATE TEMPERATURE POLYNOMIAL AND ITS DERIVATIVE
              TFPOLY = TCOEFF(NCTMAX)
              TDPOLY = TDERIV(NCTMAX)
              DO ICP = NCTMM1,1,-1
                TFPOLY = TCOEFF(ICP) + TFPOLY*TEMPOR
                TDPOLY = TDERIV(ICP) + TDPOLY*TEMPOR
              ENDDO
              TFPOLY = TCOEFF(0) + TFPOLY*TEMPOR

C             EVALUATE TEMPERATURE CORRECTION
              DELTMP = -TFPOLY/TDPOLY

C             CHECK FOR CONVERGENCE
              IF(ABS(DELTMP).GT.TOLTMP)THEN
                IF(ITITRS.LT.NTITRS)THEN
                  TEMPOR = TEMPOR + DELTMP
                  ITITRS = ITITRS + 1
                  GOTO 1000
                ELSE
                  WRITE(6,*)
     +            'Fatal: TEMPER: T iteration failed to converge'
                  WRITE(6,*)'processor:',IPROC
                  WRITE(6,*)'at point:',IC,JC,KC
                  WRITE(6,*)'with values:',TEMPOR,DELTMP
                  WRITE(6,*)DRHS(IC,JC,KC)
                  WRITE(6,*)URHS(IC,JC,KC)
                  WRITE(6,*)VRHS(IC,JC,KC)
                  WRITE(6,*)WRHS(IC,JC,KC)
                  WRITE(6,*)ERHS(IC,JC,KC)
                  STOP
                ENDIF
              ENDIF

C           END OF LOOP 1000

C           ===================================================================

C           SET THE NEW TEMPERATURE
            TRUN(IC,JC,KC) = TEMPOR

C           ===================================================================

C           FOR ALL SPECIES RELOCATE TEMPERATURE IN AN INTERVAL
C           EVALUATE MIXTURE SPECIFIC HEAT CP
            DO IINDEX = 1,NINTMX
              ITNDEX(IC,JC,KC,IINDEX) = 0
            ENDDO
            TRANSP(IC,JC,KC) = ZERO
            DO ISPEC = 1,NSPEC

              ITINT = 1
1100          CONTINUE
                IF(TRUN(IC,JC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                  IF(ITINT.LT.NTINT(ISPEC))THEN
                    ITINT = ITINT + 1
                    GOTO 1100
                  ENDIF
                ENDIF
C             END OF LOOP 1100

C             SET THE TEMPERATURE INTERVAL INDEX
              IINDEX = 1 + (ISPEC-1)/NSPIMX
              IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
              ITNDEX(IC,JC,KC,IINDEX) = ITNDEX(IC,JC,KC,IINDEX)
     +                                +(ITINT-1)*NTBASE**IPOWER

C             =================================================================

C             EVALUATE MIXTURE SPECIFIC HEAT CP
              CPFORY = AMASCP(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
              DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                CPFORY = CPFORY*TRUN(IC,JC,KC) + AMASCP(ICP,ITINT,ISPEC)
              ENDDO
              TRANSP(IC,JC,KC) = TRANSP(IC,JC,KC)
     +                         + YRHS(IC,JC,KC,ISPEC)*CPFORY

            ENDDO
            TRANSP(IC,JC,KC) = TRANSP(IC,JC,KC)/DRHS(IC,JC,KC)

C           ===================================================================

C           EVALUATE MIXTURE PRESSURE
            PRUN(IC,JC,KC) = TRUN(IC,JC,KC)*STORE7(IC,JC,KC)

C           ===================================================================

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
