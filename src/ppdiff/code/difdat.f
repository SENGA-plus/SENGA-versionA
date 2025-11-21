      SUBROUTINE DIFDAT
 
C     *************************************************************************
C
C     DIFDAT
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
C     CREATES THE BINARY DIFFUSION COEFFICIENT DATA FOR PPDIFF
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
      DOUBLE PRECISION SIGRED,EOKRED,DELRED,DIMAMU,DIMBMU
      DOUBLE PRECISION DIMPMU,ALPHAN,XIFACT
      DOUBLE PRECISION REDMAS,CROSEC
      DOUBLE PRECISION CHISQR
      DOUBLE PRECISION FORNOW
      INTEGER ISPEC,JSPEC
      INTEGER IPSPEC,INSPEC
      INTEGER IFITLO,IFITHI,IR2FIT,NR2FIT

CC     DIAGNOSTICS
C      INTEGER NCPLOF
C      CHARACTER*50 FNPLOF
C      CHARACTER*2 CSPEC1,CSPEC2


C     BEGIN
C     =====

C     =========================================================================

C     PRELIMINARIES
C     -------------

C     NUMBER OF COEFFICIENTS FOR DIFFUSIVITY
      NCODIF = 4

C     INITIALISE DIFFUSIVITY COEFFICIENTS
      DO JSPEC = 1, NSPEC
        DO ISPEC = 1, NSPEC

          DO ICOEFF = 1, NCODIF
            DIFFCO(ICOEFF,ISPEC,JSPEC) = ZERO
          ENDDO

        ENDDO
      ENDDO

C     =========================================================================

C     RUN THROUGH ALL PAIRS OF SPECIES
C     --------------------------------
      DO JSPEC = 1, NSPEC
        DO ISPEC = 1, JSPEC

CC         DIAGNOSTICS
C          WRITE(6,*)'DIFDAT:',ISPEC,JSPEC

C         =====================================================================

C         REDUCED MOLECULAR MASS
          FORNOW = EMOMAS(ISPEC) + EMOMAS(JSPEC)
          REDMAS = EMOMAS(ISPEC)*EMOMAS(JSPEC)/FORNOW

C         REDUCED COLLISION CROSS-SECTION
          SIGRED = HALF*(SIGMAD(ISPEC) + SIGMAD(JSPEC))

C         REDUCED LENNARD-JONES POTENTIAL WELL DEPTH
          EOKRED = SQRT(EPSOKB(ISPEC)*EPSOKB(JSPEC))

C         =====================================================================

C         TEMPERATURE-DEPENDENT FACTORS
C         -----------------------------

C         CHECK COLLISION POLARITY
          IF(FPOLAR(ISPEC).EQV.FPOLAR(JSPEC))THEN

C           -------------------------------------------------------------------

C           POLAR:POLAR OR NON-POLAR:NON-POLAR COLLISION
C           --------------------------------------------

C           CHECK FOR POLAR MOLECULE(S)
            IF(FPOLAR(ISPEC))THEN

C             -----------------------------------------------------------------

C             POLAR MOLECULES
C             ---------------

C             EVALUATE REDUCED DIPOLE MOMENTS
              FORNOW = SIGMAD(ISPEC)
              FORNOW = EPSOKB(ISPEC)*BOLTZC*FORNOW*FORNOW*FORNOW
              DIMAMU = DIMOMU(ISPEC)
              DIMAMU = DIMAMU/SQRT(FORNOW)
              FORNOW = SIGMAD(JSPEC)
              FORNOW = EPSOKB(JSPEC)*BOLTZC*FORNOW*FORNOW*FORNOW
              DIMBMU = DIMOMU(JSPEC)
              DIMBMU = DIMBMU/SQRT(FORNOW)
              DELRED = HALF*DIMAMU*DIMBMU

C             COLLISION INTEGRAL
C             COLLECT THE VALUES OF OMEGA(1,1)* AT NON-ZERO DIPOLE MOMENT
C             EVALUATE THE NON-REDUCED TEMPERATURE
              DO IRTEMP = 1, NTRD11
  
                DO IDELTA = 1, NDEL11
  
                  OMDELT(IDELTA) = OMEG11(IRTEMP,IDELTA)
  
                ENDDO
  
                OMVALS(IRTEMP) = CUBINT(DELR11,OMDELT,NDEL11,DELRED)
                TMVALS(IRTEMP) = TRED11(IRTEMP)*EOKRED

              ENDDO
  
C             -----------------------------------------------------------------

            ELSE

C             -----------------------------------------------------------------

C             NON-POLAR MOLECULES
C             -------------------

C             EVALUATE REDUCED DIPOLE MOMENT
              DELRED = ZERO

C             COLLISION INTEGRAL
C             COLLECT THE VALUES OF OMEGA(1,1)* AT ZERO DIPOLE MOMENT
C             EVALUATE THE NON-REDUCED TEMPERATURE
              IDELTA = 1
              DO IRTEMP = 1, NTRD11
  
                OMVALS(IRTEMP) = OMEG11(IRTEMP,IDELTA)
                TMVALS(IRTEMP) = TRED11(IRTEMP)*EOKRED

              ENDDO

C             -----------------------------------------------------------------

            ENDIF

C           -------------------------------------------------------------------

          ELSE

C           -------------------------------------------------------------------

C           POLAR:NON-POLAR COLLISION
C           -------------------------

C           CHECK WHICH MOLECULE IS POLAR
            IPSPEC = ISPEC
            INSPEC = JSPEC
            IF(FPOLAR(JSPEC))THEN
              IPSPEC = JSPEC
              INSPEC = ISPEC
            ENDIF

C           REDUCED POLARISABILITY
            FORNOW = SIGMAD(INSPEC)
            FORNOW = FORNOW*FORNOW*FORNOW
            ALPHAN = POLALF(INSPEC)/FORNOW

C           REDUCED DIPOLE MOMENT
            FORNOW = SIGMAD(IPSPEC)
            FORNOW = EPSOKB(IPSPEC)*BOLTZC*FORNOW*FORNOW*FORNOW
            DIMPMU = DIMOMU(IPSPEC)
            DIMPMU = DIMPMU/SQRT(FORNOW)

C           XI FACTOR
            FORNOW = EPSOKB(IPSPEC)/EPSOKB(INSPEC)
            XIFACT = ONE + HALF*HALF*ALPHAN*DIMPMU*SQRT(FORNOW)

C           MODIFIED REDUCED COLLISION CROSS-SECTION
            SIGRED = SIGRED*EXP(-6.0D0*LOG(XIFACT))

C           MODIFIED REDUCED LENNARD-JONES POTENTIAL WELL DEPTH
            EOKRED = EOKRED*XIFACT*XIFACT

C           EVALUATE REDUCED DIPOLE MOMENT
            DELRED = ZERO

C           COLLISION INTEGRAL
C           COLLECT THE VALUES OF OMEGA(1,1)* AT ZERO DIPOLE MOMENT
C           EVALUATE THE NON-REDUCED TEMPERATURE
            IDELTA = 1
            DO IRTEMP = 1, NTRD11
  
              OMVALS(IRTEMP) = OMEG11(IRTEMP,IDELTA)
              TMVALS(IRTEMP) = TRED11(IRTEMP)*EOKRED

            ENDDO

C           -------------------------------------------------------------------

          ENDIF

C         =====================================================================

C         FIT TEMPERATURE-DEPENDENT DATA
C         ------------------------------
C         WITHIN THE REQUIRED RANGE OF (NON-REDUCED) TEMPERATURE

C         LOW END
          IFITLO = 0
          IRTEMP = 0
1000      CONTINUE
            IRTEMP = IRTEMP + 1
            IF(TMVALS(IRTEMP).LT.TFITLO)THEN
              IF(IRTEMP.LT.NTRD11)GOTO 1000
            ELSE
              IFITLO = IRTEMP - 1
            ENDIF
C         END OF LOOP 1000
          IFITLO = MAX(IFITLO,1)

C         HIGH END
          IFITHI = NTRD11 + 1
          IRTEMP = NTRD11 + 1
1010      CONTINUE
            IRTEMP = IRTEMP - 1
            IF(TMVALS(IRTEMP).GT.TFITHI)THEN
              IF(IRTEMP.GT.1)GOTO 1010
            ELSE
              IFITHI = IRTEMP + 1
            ENDIF
C         END OF LOOP 1010
          IFITHI = MIN(IFITHI,NTRD11)

CC         DIAGNOSTICS
C          WRITE(6,'(2I5)')IFITLO,IFITHI

C         COLLECT THE VALUES TO FIT
          IR2FIT = 0
          DO IRTEMP = IFITLO, IFITHI
            IR2FIT = IR2FIT + 1
            TM2FIT(IR2FIT) = LOG(TMVALS(IRTEMP)/TREFGB)
            OM2FIT(IR2FIT) = THREE*HALF*TM2FIT(IR2FIT)
     +                     - LOG(OMVALS(IRTEMP))
          ENDDO
          NR2FIT = IR2FIT

CC         DIAGNOSTICS
C          DO IRTEMP = 1, NR2FIT
C      WRITE(6,'(2I5,2F10.3)')ISPEC,IRTEMP,TM2FIT(IRTEMP),OM2FIT(IRTEMP)
C          ENDDO

C         FIT THE POLYNOMIAL
          CHISQR = ZERO
          CALL POLFIT(TM2FIT,OM2FIT,NR2FIT,NFITMX,
     +                ACOEFS,NCODIF,NFCOMX,CHISQR)


CC         DIAGNOSTICS
C          WRITE(6,'(I5,1PE12.4)')ISPEC,CHISQR
C        WRITE(6,'(I5,6(1PE12.4))')ISPEC,(ACOEFS(ICOEFF),ICOEFF=1,NCODIF)
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = ACOEFS(NCODIF)
C            DO ICOEFF = NCODIF-1,1,-1
C              FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C            ENDDO
C        WRITE(6,'(I5,5(1PE12.4))')IR2FIT,TM2FIT(IR2FIT),OM2FIT(IR2FIT),
C     +                            FORNOW,OM2FIT(IR2FIT)-FORNOW,
C     +                            ONE-FORNOW/OM2FIT(IR2FIT)
C          ENDDO

CC         DIAGNOSTICS
C          NCPLOF = 7
C          WRITE(CSPEC1,'(I2.2)')ISPEC
C          WRITE(CSPEC2,'(I2.2)')JSPEC
C          FNPLOF = 'diffit'//CSPEC1//CSPEC2//'.res'
C          OPEN(UNIT=NCPLOF,FILE=FNPLOF,FORM='FORMATTED')
C          DO IR2FIT = 1, NR2FIT
C        WRITE(NCPLOF,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT)),OM2FIT(IR2FIT)
C          ENDDO
C          WRITE(NCPLOF,*)
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = ACOEFS(NCODIF)
C            DO ICOEFF = NCODIF-1,1,-1
C              FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C            ENDDO
C            WRITE(NCPLOF,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT)),FORNOW
C          ENDDO
C          CLOSE(NCPLOF)

C         =====================================================================

C         TEMPERATURE-INDEPENDENT FACTOR
C         ------------------------------
          CROSEC = PI*SIGRED*SIGRED
          FORNOW = BOLTZC*TREFGB
          FORNOW = TWOPI*FORNOW*FORNOW*FORNOW/REDMAS
          DIFFCO(1,ISPEC,JSPEC) = DIFACT*SQRT(FORNOW)/(PREFGB*CROSEC)

C         =====================================================================

C         COMPLETE THE COEFFICIENTS FOR THE BINARY DIFFUSION COEFFICIENTS
          DIFFCO(1,ISPEC,JSPEC) = LOG(DIFFCO(1,ISPEC,JSPEC)) + ACOEFS(1)
          DO ICOEFF = 2,NCODIF
            DIFFCO(ICOEFF,ISPEC,JSPEC) = ACOEFS(ICOEFF)
          ENDDO

CC         DIAGNOSTICS
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = DIFFCO(NCODIF,ISPEC,JSPEC)
C            DO ICOEFF = NCODIF-1,1,-1
C          FORNOW = FORNOW*TM2FIT(IR2FIT) + DIFFCO(ICOEFF,ISPEC,JSPEC)
C            ENDDO
C          WRITE(6,'(2(1PE12.4))')EXP(TM2FIT(IR2FIT))*TREFGB,EXP(FORNOW)
C          ENDDO
C          WRITE(6,*)

C         =====================================================================

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
