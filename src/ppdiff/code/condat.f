      SUBROUTINE CONDAT
 
C     *************************************************************************
C
C     CONDAT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     06-OCT-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     CREATES THE THERMAL CONDUCTIVITY DATA FOR PPDIFF
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_diffin.h'
      INCLUDE 'com_ppdcom.h'
C     -------------------------------------------------------------------------


C     EXTERNAL FUNCTION
C     =================
      DOUBLE PRECISION CUBINT
      EXTERNAL CUBINT


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TMVALS(NTRDMX),OMVALS(NTRDMX)
      DOUBLE PRECISION TM2FIT(NTRDMX),OM2FIT(NTRDMX)
      DOUBLE PRECISION OMDELT(NDELMX)
      DOUBLE PRECISION ACOEFS(NFCOMX)
      DOUBLE PRECISION FCTRAN(NTRDMX),FCVROT(NTRDMX),FCVVIB(NTRDMX)
      DOUBLE PRECISION CVTOTL(NTRDMX),CVRVIB(NTRDMX)
      DOUBLE PRECISION CVTRAN,CVRROT,CVROVT
      DOUBLE PRECISION FMTRAN
      DOUBLE PRECISION DELRED
      DOUBLE PRECISION TREDZD,ZFNREF,ZFNRED,ZRATIO
      DOUBLE PRECISION RRATIO,AFACTR,BFACTR,AOVERB
      DOUBLE PRECISION CHISQR
      DOUBLE PRECISION ZROTCF(4)
      DOUBLE PRECISION FORNOW,TMLOCL,VISFIT
      DOUBLE PRECISION THRHLF,FIVHLF,FIVTHD,TWOVPI
      DOUBLE PRECISION DVFACT
      INTEGER ICOEFF,NCOEFF,ITINT,IRTEMP
      INTEGER ISPEC
      INTEGER IFITLO,IFITHI,IR2FIT,NR2FIT
      CHARACTER*50 ERRSTR

CC     DIAGNOSTICS
C      DOUBLE PRECISION TINCRM
C      INTEGER NCPLOF
C      CHARACTER*50 FNPLOF
C      CHARACTER*2 CSPEC


C     BEGIN
C     =====

C     =========================================================================

C     PRELIMINARIES
C     -------------

C     NUMBERS
      THRHLF = 3.0D0/2.0D0
      FIVHLF = 5.0D0/2.0D0
      FIVTHD = 5.0D0/3.0D0
      TWOVPI = TWO/PI
      DVFACT = DIFACT/VSFACT

C     ROTATIONAL RELAXATION COLLISION NUMBER TEMP FUNCTION COEFFICIENTS
      ZROTCF(4) = ROOTPI*ROOTPI*ROOTPI
      ZROTCF(3) = HALF*HALF*PI*PI + TWO
      ZROTCF(2) = HALF*ZROTCF(4)
      ZROTCF(1) = ONE

C     UNIVERSAL GAS CONSTANT
      RGUNIV = AVOGNO*BOLTZC

C     CV TRANSLATIONAL
      CVTRAN = THRHLF

C     CV ROTATIONAL
      CVRROT = ONE

C     NUMBER OF COEFFICIENTS FOR CONDUCTIVITY
      NCOCON = 4

C     =========================================================================

C     RUN THROUGH ALL SPECIES
C     -----------------------
      DO ISPEC = 1, NSPEC

C       =======================================================================

CC       DIAGNOSTICS
C        WRITE(6,*)'CONDAT',ISPEC

C       INITIALISE THE CONDUCTIVITY COEFFICIENTS
C       ----------------------------------------
        DO ICOEFF = 1,NCOCON
          CONDCO(ICOEFF,ISPEC) = ZERO
        ENDDO

C       =======================================================================

C       CHECK MOLECULAR GEOMETRY INDICATOR
C       ----------------------------------
        IF(IMGEOM(ISPEC).EQ.0)THEN

C         =====================================================================

C         MOLECULE IS MONATOMIC
C         ---------------------
          FMTRAN = FIVHLF
          FORNOW = FMTRAN*CVTRAN*RGUNIV/WMOLAR(ISPEC)
          CONDCO(1,ISPEC) = VISCCO(1,ISPEC) + LOG(FORNOW)
          DO ICOEFF = 2, MIN(NCOCON,NCOVIS)
            CONDCO(ICOEFF,ISPEC) = VISCCO(ICOEFF,ISPEC)
          ENDDO

CC         DIAGNOSTICS
C          WRITE(6,*)'CONDAT:',ISPEC,'  Monatomic'
C        WRITE(6,'(4(1PE12.4))')(CONDCO(ICOEFF,ISPEC),ICOEFF = 1, NCOCON)
 
CC         DIAGNOSTICS
C          NCPLOF = 7
C          WRITE(CSPEC,'(I2.2)')ISPEC
C          FNPLOF = 'confit'//CSPEC//'.res'
C          OPEN(UNIT=NCPLOF,FILE=FNPLOF,FORM='FORMATTED')
C          TINCRM = 1.0D2
C          NR2FIT = 1 + NINT((TFITHI-TFITLO)/TINCRM)
C          TM2FIT(1) = TFITLO
C          DO IR2FIT = 2, NR2FIT
C            TM2FIT(IR2FIT) = TM2FIT(IR2FIT-1) + TINCRM
C            FORNOW = CONDCO(NCOCON,ISPEC)
C            DO ICOEFF = NCOCON-1,1,-1
C              FORNOW = FORNOW*LOG(TM2FIT(IR2FIT)/TREFGB)
C     +               + CONDCO(ICOEFF,ISPEC)
C            ENDDO
C            WRITE(NCPLOF,'(2(1PE12.4))')TM2FIT(IR2FIT),FORNOW
C          ENDDO
C          CLOSE(NCPLOF)

CC         DIAGNOSTICS
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = CONDCO(NCOCON,ISPEC)
C            DO ICOEFF = NCOCON-1,1,-1
C              FORNOW = FORNOW*LOG(TM2FIT(IR2FIT)/TREFGB)
C     +               + CONDCO(ICOEFF,ISPEC)
C            ENDDO
C            WRITE(6,'(2(1PE12.4))')TM2FIT(IR2FIT),FORNOW
C          ENDDO
C          WRITE(6,*)
    
C         =====================================================================

        ELSE

C         =====================================================================

C         MOLECULE IS NOT MONATOMIC
C         -------------------------

C         ---------------------------------------------------------------------

C         ROTATIONAL RELAXATION COLLISION NUMBER
C         TEMPERATURE FUNCTION AT REFERENCE TEMPERATURE
          TREDZD = SQRT(EPSOKB(ISPEC)/TREFZD)
          ZFNREF = ZROTCF(4)
          DO ICOEFF = 3,1,-1
            ZFNREF = ZFNREF*TREDZD + ZROTCF(ICOEFF)
          ENDDO

CC         DIAGNOSTICS
C          WRITE(6,'(3(1PE12.4))')TREFZD,TREDZD,ZFNREF

C         ---------------------------------------------------------------------

C         TEMPERATURE-DEPENDENT FACTORS
C         -----------------------------

C         CHECK FOR POLAR MOLECULE
          IF(FPOLAR(ISPEC))THEN

C           -------------------------------------------------------------------

C           POLAR MOLECULE
C           --------------

CC           DIAGNOSTICS
C            WRITE(6,*)'CONDAT:',ISPEC,'  Polar'
 
C           REDUCED DIPOLE MOMENT
            FORNOW = SIGMAD(ISPEC)
            FORNOW = EPSOKB(ISPEC)*BOLTZC*FORNOW*FORNOW*FORNOW
            DELRED = DIMOMU(ISPEC)
            DELRED = HALF*DELRED*DELRED/FORNOW

C           RATIO OF COLLISION INTEGRALS: OMEGA(2,2)*/OMEGA(1,1)* (= ASTAR)
C           COLLECT THE VALUES OF ASTAR AT NON-ZERO DIPOLE MOMENT
C           EVALUATE THE NON-REDUCED TEMPERATURE
            DO IRTEMP = 1, NTRDAA

              DO IDELTA = 1, NDELAA

                OMDELT(IDELTA) = ASTAR(IRTEMP,IDELTA)

              ENDDO

              OMVALS(IRTEMP) = CUBINT(DELRAA,OMDELT,NDELAA,DELRED)
              TMVALS(IRTEMP) = TREDAA(IRTEMP)*EPSOKB(ISPEC)

CC             DIAGNOSTICS
C      WRITE(6,'(2I5,2F12.3)')ISPEC,IRTEMP,TMVALS(IRTEMP),OMVALS(IRTEMP)

            ENDDO

C           -------------------------------------------------------------------

          ELSE

C           -------------------------------------------------------------------

C           NON-POLAR MOLECULE
C           ------------------
C           RATIO OF COLLISION INTEGRALS: OMEGA(2,2)*/OMEGA(1,1)* (= ASTAR)
C           COLLECT THE VALUES OF ASTAR AT ZERO DIPOLE MOMENT
C           EVALUATE THE NON-REDUCED TEMPERATURE

CC           DIAGNOSTICS
C            WRITE(6,*)'CONDAT:',ISPEC,'  Non-polar'
 
            IDELTA = 1
            DO IRTEMP = 1, NTRDAA

              OMVALS(IRTEMP) = ASTAR(IRTEMP,IDELTA)
              TMVALS(IRTEMP) = TREDAA(IRTEMP)*EPSOKB(ISPEC)

CC             DIAGNOSTICS
C      WRITE(6,'(2I5,2F12.3)')ISPEC,IRTEMP,TMVALS(IRTEMP),OMVALS(IRTEMP)

            ENDDO

C           -------------------------------------------------------------------

          ENDIF

C         =====================================================================

C         FIT THE TEMPERATURE-DEPENDENT DATA
C         ----------------------------------
C         WITHIN THE REQUIRED RANGE OF (NON-REDUCED) TEMPERATURE

C         LOW END
          IFITLO = 0
          IRTEMP = 0
2000      CONTINUE
            IRTEMP = IRTEMP + 1
            IF(TMVALS(IRTEMP).LT.TFITLO)THEN
              IF(IRTEMP.LT.NTRDAA)GOTO 2000
            ELSE
              IFITLO = IRTEMP - 1
            ENDIF
C         END OF LOOP 2000
          IFITLO = MAX(IFITLO,1)

C         HIGH END
          IFITHI = NTRDAA + 1
          IRTEMP = NTRDAA + 1
2010      CONTINUE
            IRTEMP = IRTEMP - 1
            IF(TMVALS(IRTEMP).GT.TFITHI)THEN
              IF(IRTEMP.GT.1)GOTO 2010
            ELSE
              IFITHI = IRTEMP + 1
            ENDIF
C         END OF LOOP 2010
          IFITHI = MIN(IFITHI,NTRDAA)

CC         DIAGNOSTICS
C          WRITE(6,'(2I5)')IFITLO,IFITHI

C         =====================================================================

C         CV TOTAL
C         --------
          DO IRTEMP = IFITLO, IFITHI

C           LOCATE TEMPERATURE WITHIN AN INTERVAL
            TMLOCL = TMVALS(IRTEMP)
            ITINT = 1
1000        CONTINUE
              IF(TMLOCL.GT.TINTHI(ITINT,ISPEC))THEN
                IF(ITINT.LT.NTINT(ISPEC))THEN
                  ITINT = ITINT + 1
                  GOTO 1000
                ENDIF
              ENDIF
C           END OF LOOP 1000

CC           DIAGNOSTICS
C            WRITE(6,'(I5,F12.4,I5)')IRTEMP,TMLOCL,ITINT

C           EVALUATE CV/R0 TOTAL
            NCOEFF = NCOFCP(ITINT,ISPEC)-2
            FORNOW = ACOFCP(NCOEFF,ITINT,ISPEC)
            DO ICOEFF = NCOEFF-1,1,-1
              FORNOW = FORNOW*TMLOCL + ACOFCP(ICOEFF,ITINT,ISPEC)
            ENDDO
            CVTOTL(IRTEMP) = FORNOW - ONE

CC           DIAGNOSTICS
C          WRITE(6,'(2I5,2(1PE12.4))')IRTEMP,ITINT,TMLOCL,CVTOTL(IRTEMP)

          ENDDO

C         =====================================================================

C         CHECK LINEARITY OF MOLECULE
C         ---------------------------
          IF(IMGEOM(ISPEC).EQ.1)THEN

C           -------------------------------------------------------------------

C           MOLECULE IS LINEAR
C           ------------------
            CVRROT = ONE
            DO IRTEMP = IFITLO, IFITHI
              CVRVIB(IRTEMP) = CVTOTL(IRTEMP) - FIVHLF
            ENDDO

C           DIAGNOSTICS
            WRITE(6,*)'       ',ISPEC,'  Linear'
 
C           -------------------------------------------------------------------

          ELSE IF(IMGEOM(ISPEC).EQ.2)THEN

C           -------------------------------------------------------------------

C           MOLECULE IS NON-LINEAR
C           ----------------------
            CVRROT = THRHLF
            DO IRTEMP = IFITLO, IFITHI
              CVRVIB(IRTEMP) = CVTOTL(IRTEMP) - THREE
            ENDDO

C           DIAGNOSTICS
            WRITE(6,*)'       ',ISPEC,'  Non-linear'
 
C           -------------------------------------------------------------------

          ELSE

C           INCORRECT MOLECULAR GEOMETRY INDICATOR
            ERRSTR = 'incorrect molecular geometry indicator'
            CALL ERHAND(ERRSTR,IEFATL)

          ENDIF

C         RATIO OF ROTATIONAL TO TRANSLATIONAL CV
          CVROVT = CVRROT/CVTRAN

C         =====================================================================

C         WEIGHTING FACTORS
C         -----------------
          DO IRTEMP = IFITLO, IFITHI

C           ROTATIONAL RELAXATION COLLISION NUMBER
C           TEMPERATURE FUNCTION AT REDUCED TEMPERATURE
            TREDZD = SQRT(EPSOKB(ISPEC)/TMVALS(IRTEMP))
            ZFNRED = ZROTCF(4)
            DO ICOEFF = 3,1,-1
              ZFNRED = ZFNRED*TREDZD + ZROTCF(ICOEFF)
            ENDDO
            ZRATIO = ZEDROT(ISPEC)*ZFNREF/ZFNRED

CC           DIAGNOSTICS
C            WRITE(6,'(I5,4(1PE12.4))')IRTEMP,TMVALS(IRTEMP),TREDZD,
C     +                                       ZEDROT(ISPEC),ZRATIO

C           RATIO OF MASS DIFFUSIVITY TO VISCOSITY
            RRATIO = DVFACT*TWO*OMVALS(IRTEMP)

C           FACTORS A AND B
            AFACTR = FIVHLF - RRATIO
            BFACTR = ZRATIO + TWOVPI*(FIVTHD*CVRROT + RRATIO)
            AOVERB = AFACTR/BFACTR

C           WEIGHTING FACTORS
            FCTRAN(IRTEMP) = FIVHLF*(ONE - TWOVPI*CVROVT*AOVERB)
            FCVROT(IRTEMP) = RRATIO*(ONE + TWOVPI*AOVERB)
            FCVVIB(IRTEMP) = RRATIO

          ENDDO

C         =====================================================================

C         COLLECT THE VALUES TO FIT
C         -------------------------
          IR2FIT = 0
          DO IRTEMP = IFITLO, IFITHI

C           LOG OF TEMPERATURE
            IR2FIT = IR2FIT + 1
            TM2FIT(IR2FIT) = LOG(TMVALS(IRTEMP)/TREFGB)

C           VISCOSITY

CC           DIAGNOSTICS
C            WRITE(6,'(I5,4(1PE12.4))')ISPEC,
C     +                 (VISCCO(ICOEFF,ISPEC),ICOEFF=1,NCOVIS)
           
            VISFIT = VISCCO(NCOVIS,ISPEC)
            DO ICOEFF = NCOVIS-1,1,-1
              VISFIT = VISFIT*TM2FIT(IR2FIT) + VISCCO(ICOEFF,ISPEC)
            ENDDO
            VISFIT = EXP(VISFIT)

CC           DIAGNOSTICS
C            WRITE(6,'(I5,2(1PE12.4))')IRTEMP,EXP(TM2FIT(IR2FIT)),VISFIT
           
C           CONDUCTIVITY
            FORNOW = FCTRAN(IRTEMP)*CVTRAN
     +             + FCVROT(IRTEMP)*CVRROT
     +             + FCVVIB(IRTEMP)*CVRVIB(IRTEMP)
            FORNOW = FORNOW*VISFIT*RGUNIV/WMOLAR(ISPEC)
            OM2FIT(IR2FIT) = LOG(FORNOW)

          ENDDO
          NR2FIT = IR2FIT

CC         DIAGNOSTICS
C          DO IRTEMP = 1, NR2FIT
C      WRITE(6,'(2I5,2F10.3)')ISPEC,IRTEMP,TM2FIT(IRTEMP),OM2FIT(IRTEMP)
C          ENDDO

C         FIT THE POLYNOMIAL
          CHISQR = ZERO
          CALL POLFIT(TM2FIT,OM2FIT,NR2FIT,NFITMX,
     +                ACOEFS,NCOCON,NFCOMX,CHISQR)


CC         DIAGNOSTICS
C          WRITE(6,'(I5,1PE12.4)')ISPEC,CHISQR
C        WRITE(6,'(I5,6(1PE12.4))')ISPEC,(ACOEFS(ICOEFF),ICOEFF=1,NCOCON)
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = ACOEFS(NCOCON)
C            DO ICOEFF = NCOCON-1,1,-1
C              FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C            ENDDO
C        WRITE(6,'(I5,5(1PE12.4))')IR2FIT,TM2FIT(IR2FIT),OM2FIT(IR2FIT),
C     +                            FORNOW,OM2FIT(IR2FIT)-FORNOW,
C     +                            ONE-FORNOW/OM2FIT(IR2FIT)
C          ENDDO

CC         DIAGNOSTICS
C          NCPLOF = 7
C          WRITE(CSPEC,'(I2.2)')ISPEC
C          FNPLOF = 'confit'//CSPEC//'.res'
C          OPEN(UNIT=NCPLOF,FILE=FNPLOF,FORM='FORMATTED')
C          DO IR2FIT = 1, NR2FIT
C          WRITE(NCPLOF,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT)),OM2FIT(IR2FIT)
C          ENDDO
C          WRITE(NCPLOF,*)
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = ACOEFS(NCOCON)
C            DO ICOEFF = NCOCON-1,1,-1
C              FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C            ENDDO
C            WRITE(NCPLOF,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT)),FORNOW
C          ENDDO
C          CLOSE(NCPLOF)

C         =====================================================================

C         COMPLETE THE POLYNOMIAL COEFFICIENTS FOR THERMAL CONDUCTIVITY
C         -------------------------------------------------------------
          DO ICOEFF = 1,NCOCON
            CONDCO(ICOEFF,ISPEC) = ACOEFS(ICOEFF)
          ENDDO

C         =====================================================================

CC         DIAGNOSTICS
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = CONDCO(NCOCON,ISPEC)
C            DO ICOEFF = NCOCON-1,1,-1
C              FORNOW = FORNOW*TM2FIT(IR2FIT) + CONDCO(ICOEFF,ISPEC)
C            ENDDO
C          WRITE(6,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT))*TREFGB,EXP(FORNOW)
C          ENDDO
C          WRITE(6,*)
    
C         =====================================================================

        ENDIF
C       CHECK MOLECULAR GEOMETRY INDICATOR

C       =======================================================================

      ENDDO

C     =========================================================================

C     DIAGNOSTICS
      CALL NONDIM

C     =========================================================================


      RETURN
      END
