      SUBROUTINE NONDIM
 
C     *************************************************************************
C
C     NONDIM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     16-DEC-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     EVALUATES NON-DIMENSIONAL NUMBERS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_diffin.h'
      INCLUDE 'com_ppdcom.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
      INTEGER NRTMAX
      PARAMETER(NRTMAX = 50)

      INTEGER NRTEMP
      PARAMETER(NRTEMP = 21)
      DOUBLE PRECISION TINC
      PARAMETER(TINC = 1.0D2)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TVALUE(NRTMAX),TRATIO(NRTMAX)
      DOUBLE PRECISION VISCOS(NRTMAX),DIFFUS(NRTMAX),CONDUC(NRTMAX)
      DOUBLE PRECISION DVALUE(NRTMAX),CPTOTL(NRTMAX)
      DOUBLE PRECISION PRANNO(NRTMAX),SHMDNO(NRTMAX),CLEWNO(NRTMAX)
      DOUBLE PRECISION FORNOW
      INTEGER ISPEC,JSPEC
      INTEGER ITINT
      INTEGER IRTEMP
      INTEGER ICOEFF,NCOEFF
 

C     BEGIN
C     =====

C     =========================================================================

C     TEMPERATURE
      TVALUE(1) = 2.73D2
      TVALUE(2) = TREFGB
      DO IRTEMP = 3, NRTEMP
        TVALUE(IRTEMP) = TVALUE(IRTEMP-1) + TINC
      ENDDO

C     NON-DIMENSIONAL NUMBERS
      DO ISPEC = 1,NSPEC

        DO IRTEMP = 1, NRTEMP

C         TEMPERATURE RATIO
          TRATIO(IRTEMP) = LOG(TVALUE(IRTEMP)/TREFGB)

C         VISCOSITY
          FORNOW = VISCCO(NCOVIS,ISPEC)
          DO ICOEFF = NCOVIS-1,1,-1
            FORNOW = FORNOW*TRATIO(IRTEMP) + VISCCO(ICOEFF,ISPEC)
          ENDDO
          VISCOS(IRTEMP) = EXP(FORNOW)

C         DIFFUSION COEFFICIENT
          JSPEC = ISPEC
          FORNOW = DIFFCO(NCODIF,ISPEC,JSPEC)
          DO ICOEFF = NCODIF-1,1,-1
            FORNOW = FORNOW*TRATIO(IRTEMP) + DIFFCO(ICOEFF,ISPEC,JSPEC)
          ENDDO
          DIFFUS(IRTEMP) = EXP(FORNOW)

C         THERMAL CONDUCTIVITY
          FORNOW = CONDCO(NCOCON,ISPEC)
          DO ICOEFF = NCOCON-1,1,-1
            FORNOW = FORNOW*TRATIO(IRTEMP) + CONDCO(ICOEFF,ISPEC)
          ENDDO
          CONDUC(IRTEMP) = EXP(FORNOW)

C         DENSITY
          FORNOW = RGUNIV*TVALUE(IRTEMP)/WMOLAR(ISPEC)
          DVALUE(IRTEMP) = PREFGB/FORNOW

C         MASS-BASED SPECIFIC HEAT CAPACITY AT CONSTANT PRESSURE

C         LOCATE TEMPERATURE WITHIN AN INTERVAL
          ITINT = 1
1000      CONTINUE
            IF(TVALUE(IRTEMP).GT.TINTHI(ITINT,ISPEC))THEN
              IF(ITINT.LT.NTINT(ISPEC))THEN
                ITINT = ITINT + 1
                GOTO 1000
              ENDIF
            ENDIF
C         END OF LOOP 1000

C         EVALUATE CP
          NCOEFF = NCOFCP(ITINT,ISPEC)-2
          FORNOW = ACOFCP(NCOEFF,ITINT,ISPEC)
          DO ICOEFF = NCOEFF-1,1,-1
            FORNOW = FORNOW*TVALUE(IRTEMP) + ACOFCP(ICOEFF,ITINT,ISPEC)
          ENDDO
          CPTOTL(IRTEMP) = FORNOW*RGUNIV/WMOLAR(ISPEC)

C         NON-DIMENSIONAL NUMBERS
          PRANNO(IRTEMP) = VISCOS(IRTEMP)*CPTOTL(IRTEMP)/CONDUC(IRTEMP)

          FORNOW = DVALUE(IRTEMP)*DIFFUS(IRTEMP)
          SHMDNO(IRTEMP) = VISCOS(IRTEMP)/FORNOW

          FORNOW = FORNOW*CPTOTL(IRTEMP)
          CLEWNO(IRTEMP) = CONDUC(IRTEMP)/FORNOW

        ENDDO

C       WRITE OUT THE RESULTS
        WRITE(6,'(I5)')ISPEC
        DO IRTEMP = 1, NRTEMP

          WRITE(6,'(5(1PE15.7))')TVALUE(IRTEMP),
     +       VISCOS(IRTEMP),DIFFUS(IRTEMP),CONDUC(IRTEMP),CPTOTL(IRTEMP)

C          WRITE(6,'(5(1PE15.7))')TVALUE(IRTEMP),
C     +       PRANNO(IRTEMP),SHMDNO(IRTEMP),CLEWNO(IRTEMP)

        ENDDO

      ENDDO

C     =========================================================================


      RETURN
      END
