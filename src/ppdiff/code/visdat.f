      SUBROUTINE VISDAT
 
C     *************************************************************************
C
C     VISDAT
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
C     CREATES THE VISCOSITY DATA FOR PPDIFF
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
      DOUBLE PRECISION DELRED,CROSEC
      DOUBLE PRECISION CHISQR
      DOUBLE PRECISION FORNOW
      INTEGER ISPEC
      INTEGER IFITLO,IFITHI,IR2FIT,NR2FIT

CC     DIAGNOSTICS
C      INTEGER NCPLOF
C      CHARACTER*50 FNPLOF
C      CHARACTER*2 CSPEC


C     BEGIN
C     =====

C     =========================================================================

C     PRELIMINARIES
C     -------------

C     MOLECULAR MASS
      DO ISPEC = 1, NSPEC

        EMOMAS(ISPEC) = WMOLAR(ISPEC)/AVOGNO        

      ENDDO

C     POLAR MOLECULE FLAG
      DO ISPEC = 1, NSPEC

        FPOLAR(ISPEC) = .FALSE.
        IF(DIMOMU(ISPEC).GT.DIMTOL)FPOLAR(ISPEC) = .TRUE.

      ENDDO

C     NUMBER OF COEFFICIENTS FOR VISCOSITY
      NCOVIS = 4

C     =========================================================================

C     RUN THROUGH ALL SPECIES
C     -----------------------
      DO ISPEC = 1, NSPEC

C       INITIALISE THE VISCOSITY COEFFICIENTS
        DO ICOEFF = 1,NCOVIS
          VISCCO(ICOEFF,ISPEC) = ZERO
        ENDDO

CC       DIAGNOSTICS
C        WRITE(6,*)'VISDAT',ISPEC,FPOLAR(ISPEC)
C        WRITE(6,*)'VISDAT:',ISPEC

C       TEMPERATURE-DEPENDENT FACTORS
C       -----------------------------

C       CHECK FOR POLAR MOLECULE
        IF(FPOLAR(ISPEC))THEN

C         ---------------------------------------------------------------------

C         POLAR MOLECULE
C         --------------

C         REDUCED DIPOLE MOMENT
          FORNOW = SIGMAD(ISPEC)
          FORNOW = EPSOKB(ISPEC)*BOLTZC*FORNOW*FORNOW*FORNOW
          DELRED = DIMOMU(ISPEC)
          DELRED = HALF*DELRED*DELRED/FORNOW

CC         DIAGNOSTICS
C          WRITE(6,'(I5,F10.3)')ISPEC,DELRED

C         COLLISION INTEGRAL
C         COLLECT THE VALUES OF OMEGA(2,2)* AT NON-ZERO DIPOLE MOMENT
C         EVALUATE THE NON-REDUCED TEMPERATURE
          DO IRTEMP = 1, NTRD22

            DO IDELTA = 1, NDEL22

              OMDELT(IDELTA) = OMEG22(IRTEMP,IDELTA)

CC             DIAGNOSTICS
C              WRITE(6,'(I5,2F10.3)')IDELTA,DELR22(IDELTA),OMDELT(IDELTA)

            ENDDO            

            OMVALS(IRTEMP) = CUBINT(DELR22,OMDELT,NDEL22,DELRED)
            TMVALS(IRTEMP) = TRED22(IRTEMP)*EPSOKB(ISPEC)

CC         DIAGNOSTICS
C      WRITE(6,'(2I5,2F10.3)')ISPEC,IRTEMP,TMVALS(IRTEMP),OMVALS(IRTEMP)

          ENDDO

C         ---------------------------------------------------------------------

        ELSE

C         ---------------------------------------------------------------------

C         NON-POLAR MOLECULE
C         ------------------

C         COLLISION INTEGRAL
C         COLLECT THE OMEGA VALUES AT ZERO DIPOLE MOMENT
C         EVALUATE THE NON-REDUCED TEMPERATURE
          IDELTA = 1
          DO IRTEMP = 1, NTRD22

            OMVALS(IRTEMP) = OMEG22(IRTEMP,IDELTA)
            TMVALS(IRTEMP) = TRED22(IRTEMP)*EPSOKB(ISPEC)

CC         DIAGNOSTICS
C      WRITE(6,'(2I5,2F10.3)')ISPEC,IRTEMP,TMVALS(IRTEMP),OMVALS(IRTEMP)

          ENDDO

C         ---------------------------------------------------------------------

        ENDIF

C       =======================================================================

C       FIT TEMPERATURE-DEPENDENT DATA
C       ------------------------------
C       WITHIN THE REQUIRED RANGE OF (NON-REDUCED) TEMPERATURE

C       LOW END
        IFITLO = 0
        IRTEMP = 0
1000    CONTINUE
          IRTEMP = IRTEMP + 1
          IF(TMVALS(IRTEMP).LT.TFITLO)THEN
            IF(IRTEMP.LT.NTRD22)GOTO 1000
          ELSE
            IFITLO = IRTEMP - 1
          ENDIF
C       END OF LOOP 1000
        IFITLO = MAX(IFITLO,1)

C       HIGH END
        IFITHI = NTRD22 + 1
        IRTEMP = NTRD22 + 1
1010    CONTINUE
          IRTEMP = IRTEMP - 1
          IF(TMVALS(IRTEMP).GT.TFITHI)THEN
            IF(IRTEMP.GT.1)GOTO 1010
          ELSE
            IFITHI = IRTEMP + 1
          ENDIF
C       END OF LOOP 1010
        IFITHI = MIN(IFITHI,NTRD22)

CC       DIAGNOSTICS
C        WRITE(6,'(2I5)')IFITLO,IFITHI

C       COLLECT THE VALUES TO FIT
        IR2FIT = 0
        DO IRTEMP = IFITLO, IFITHI
          IR2FIT = IR2FIT + 1
          TM2FIT(IR2FIT) = LOG(TMVALS(IRTEMP)/TREFGB)
          OM2FIT(IR2FIT) = HALF*TM2FIT(IR2FIT) - LOG(OMVALS(IRTEMP))
        ENDDO
        NR2FIT = IR2FIT

CC       DIAGNOSTICS
C        DO IRTEMP = 1, NR2FIT
C      WRITE(6,'(2I5,2F10.3)')ISPEC,IRTEMP,TM2FIT(IRTEMP),OM2FIT(IRTEMP)
C        ENDDO

C       FIT THE POLYNOMIAL
        CHISQR = ZERO
        CALL POLFIT(TM2FIT,OM2FIT,NR2FIT,NFITMX,
     +              ACOEFS,NCOVIS,NFCOMX,CHISQR)


CC       DIAGNOSTICS
C        WRITE(6,'(I5,1PE12.4)')ISPEC,CHISQR
C        WRITE(6,'(I5,6(1PE12.4))')ISPEC,(ACOEFS(ICOEFF),ICOEFF=1,NCOVIS)
C        DO IR2FIT = 1, NR2FIT
C          FORNOW = ACOEFS(NCOVIS)
C          DO ICOEFF = NCOVIS-1,1,-1
C            FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C          ENDDO
C        WRITE(6,'(I5,5(1PE12.4))')IR2FIT,TM2FIT(IR2FIT),OM2FIT(IR2FIT),
C     +                            FORNOW,OM2FIT(IR2FIT)-FORNOW,
C     +                            ONE-FORNOW/OM2FIT(IR2FIT)
C        ENDDO

CC       DIAGNOSTICS
C        NCPLOF = 7
C        WRITE(CSPEC,'(I2.2)')ISPEC
C        FNPLOF = 'visfit'//CSPEC//'.res'
C        OPEN(UNIT=NCPLOF,FILE=FNPLOF,FORM='FORMATTED')
C        DO IR2FIT = 1, NR2FIT
C          WRITE(NCPLOF,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT)),OM2FIT(IR2FIT)
C        ENDDO
C        WRITE(NCPLOF,*)
C        DO IR2FIT = 1, NR2FIT
C          FORNOW = ACOEFS(NCOVIS)
C          DO ICOEFF = NCOVIS-1,1,-1
C            FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C          ENDDO
C          WRITE(NCPLOF,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT)),FORNOW
C        ENDDO
C        CLOSE(NCPLOF)

C       =======================================================================

C       TEMPERATURE-INDEPENDENT FACTOR
C       ------------------------------
        FORNOW = SIGMAD(ISPEC)
        CROSEC = PI*FORNOW*FORNOW
        FORNOW = PI*EMOMAS(ISPEC)*BOLTZC*TREFGB
        VISCCO(1,ISPEC) = VSFACT*SQRT(FORNOW)/CROSEC

C       =======================================================================

C       COMPLETE THE POLYNOMIAL COEFFICIENTS FOR VISCOSITY
C       --------------------------------------------------
        VISCCO(1,ISPEC) = LOG(VISCCO(1,ISPEC)) + ACOEFS(1)
        DO ICOEFF = 2,NCOVIS
          VISCCO(ICOEFF,ISPEC) = ACOEFS(ICOEFF)
        ENDDO

CC       DIAGNOSTICS
C        DO IR2FIT = 1, NR2FIT
C          FORNOW = VISCCO(NCOVIS,ISPEC)
C          DO ICOEFF = NCOVIS-1,1,-1
C            FORNOW = FORNOW*TM2FIT(IR2FIT) + VISCCO(ICOEFF,ISPEC)
C          ENDDO
C          WRITE(6,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT))*TREFGB,EXP(FORNOW)
C        ENDDO
C        WRITE(6,*)

C       =======================================================================

      ENDDO

C     =========================================================================


      RETURN
      END
