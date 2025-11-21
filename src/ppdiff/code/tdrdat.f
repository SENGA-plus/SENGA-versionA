      SUBROUTINE TDRDAT
 
C     *************************************************************************
C
C     TDRDAT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     02-DEC-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     CREATES THE THERMAL DIFFUSION RATIO COEFFICIENT DATA FOR PPDIFF
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
      DOUBLE PRECISION FIVE,SIX,TWELVE,FFTEEN,SXTEEN,FFTYFV
      PARAMETER(FIVE = 5.0D0, SIX = 6.0D0, TWELVE = 1.2D1)
      PARAMETER(FFTEEN = 1.5D1, SXTEEN = 1.6D1, FFTYFV = 5.5D1)


C     EXTERNAL FUNCTION
C     =================
      DOUBLE PRECISION CUBINT
      EXTERNAL CUBINT


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TMVALS(NTRDMX)
      DOUBLE PRECISION TMTMAA(NTRDMX),OMTMAA(NTRDMX)
      DOUBLE PRECISION TMTMBB(NTRDMX),OMTMBB(NTRDMX)
      DOUBLE PRECISION TMTMCC(NTRDMX),OMTMCC(NTRDMX)
      DOUBLE PRECISION OMVR22(NTRDMX),OMVR12(NTRDMX),OMVR13(NTRDMX)
      DOUBLE PRECISION AASTAR(NTRDMX),BBSTAR(NTRDMX),CCSTAR(NTRDMX)
      DOUBLE PRECISION TM2FIT(NTRDMX),OM2FIT(NTRDMX)
      DOUBLE PRECISION OMDELT(NDELMX)
      DOUBLE PRECISION ACOEFS(NFCOMX)
      DOUBLE PRECISION EOKRED,DELRED,DIMAMU,DIMBMU
      DOUBLE PRECISION DIMPMU,ALPHAN,XIFACT
      DOUBLE PRECISION REDMAS
      DOUBLE PRECISION CHISQR
      DOUBLE PRECISION FORNOW,THETAR
      DOUBLE PRECISION THRD,QRTR
      INTEGER ISPEC,JSPEC
      INTEGER IPSPEC,INSPEC
      INTEGER IFITLO,IFITHI,IR2FIT,NR2FIT
      INTEGER NTRDCR

CC     DIAGNOSTICS
C      INTEGER NCPLOF
C      CHARACTER*50 FNPLOF
C      CHARACTER*2 CSPEC1,CSPEC2


C     BEGIN
C     =====

C     =========================================================================

C     PRELIMINARIES
C     -------------

C     NUMBERS
      THRD = ONE/THREE
      QRTR = HALF*HALF

C     NUMBER OF COEFFICIENTS FOR THERMAL DIFFUSION RATIO
      NCOTDR = 4

C     INITIALISE THERMAL DIFFUSION RATIO COEFFICIENTS
      DO JSPEC = 1, NSPEC
        DO ISPEC = 1, NSPEC

          DO ICOEFF = 1, NCOTDR
            TDRCCO(ICOEFF,ISPEC,JSPEC) = ZERO
          ENDDO

        ENDDO
      ENDDO

C     =========================================================================

C     RUN THROUGH ALL PAIRS OF SPECIES
C     --------------------------------
      DO JSPEC = 1, NSPEC
        DO ISPEC = 1, JSPEC-1

CC         DIAGNOSTICS
C          WRITE(6,*)'TDRDAT:',ISPEC,JSPEC

C         =====================================================================

C         REDUCED MOLAR MASS
          REDMAS = WMOLAR(JSPEC) - WMOLAR(ISPEC)
          FORNOW = WMOLAR(JSPEC) + WMOLAR(ISPEC)
          REDMAS = REDMAS/FORNOW

C         REDUCED LENNARD-JONES POTENTIAL WELL DEPTH
          EOKRED = SQRT(EPSOKB(ISPEC)*EPSOKB(JSPEC))

C         MODIFIED REDUCED COLLISION CROSS-SECTION
          SIGRED = ZERO

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

C             RATIOS OF COLLISION INTEGRALS
C             OMEGA(2,2)*/OMEGA(1,1)* = ASTAR
C             OMEGA(1,2)*/OMEGA(1,1)* = CSTAR
C             OMEGA(1,3)*/OMEGA(1,1)* = (5*CSTAR-BSTAR)/4
C             COLLECT THE VALUES AT NON-ZERO DIPOLE MOMENT
C             EVALUATE THE NON-REDUCED TEMPERATURES

C             ASTAR
              DO IRTEMP = 1, NTRDAA
  
                DO IDELTA = 1, NDELAA
  
                  OMDELT(IDELTA) = ASTAR(IRTEMP,IDELTA)
  
                ENDDO
  
                OMTMAA(IRTEMP) = CUBINT(DELRAA,OMDELT,NDELAA,DELRED)
                TMTMAA(IRTEMP) = TREDAA(IRTEMP)*EOKRED

CC               DIAGNOSTICS
C                WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDAA(IRTEMP),
C     +                            TMTMAA(IRTEMP),OMTMAA(IRTEMP)

              ENDDO
  
C             BSTAR
              DO IRTEMP = 1, NTRDBB
  
                DO IDELTA = 1, NDELBB
  
                  OMDELT(IDELTA) = BSTAR(IRTEMP,IDELTA)
  
                ENDDO
  
                OMTMBB(IRTEMP) = CUBINT(DELRBB,OMDELT,NDELBB,DELRED)
                TMTMBB(IRTEMP) = TREDBB(IRTEMP)*EOKRED

CC               DIAGNOSTICS
C                WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDBB(IRTEMP),
C     +                            TMTMBB(IRTEMP),OMTMBB(IRTEMP)

              ENDDO

C             CSTAR
              DO IRTEMP = 1, NTRDCC
  
                DO IDELTA = 1, NDELCC
  
                  OMDELT(IDELTA) = CSTAR(IRTEMP,IDELTA)
  
                ENDDO
  
                OMTMCC(IRTEMP) = CUBINT(DELRCC,OMDELT,NDELCC,DELRED)
                TMTMCC(IRTEMP) = TREDCC(IRTEMP)*EOKRED

CC               DIAGNOSTICS
C                WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDCC(IRTEMP),
C     +                            TMTMCC(IRTEMP),OMTMCC(IRTEMP)

              ENDDO
 
C             -----------------------------------------------------------------

            ELSE

C             -----------------------------------------------------------------

C             NON-POLAR MOLECULES
C             -------------------

C             EVALUATE REDUCED DIPOLE MOMENT
              DELRED = ZERO

C             RATIOS OF COLLISION INTEGRALS
C             OMEGA(2,2)*/OMEGA(1,1)* = ASTAR
C             OMEGA(1,2)*/OMEGA(1,1)* = CSTAR
C             OMEGA(1,3)*/OMEGA(1,1)* = (5*CSTAR-BSTAR)/4
C             COLLECT THE VALUES AT ZERO DIPOLE MOMENT
C             EVALUATE THE NON-REDUCED TEMPERATURES
              IDELTA = 1

C             ASTAR
              DO IRTEMP = 1, NTRDAA
  
                OMTMAA(IRTEMP) = ASTAR(IRTEMP,IDELTA)
                TMTMAA(IRTEMP) = TREDAA(IRTEMP)*EOKRED

CC               DIAGNOSTICS
C                WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDAA(IRTEMP),
C     +                            TMTMAA(IRTEMP),OMTMAA(IRTEMP)

              ENDDO
  
C             BSTAR
              DO IRTEMP = 1, NTRDBB
  
                OMTMBB(IRTEMP) = BSTAR(IRTEMP,IDELTA)
                TMTMBB(IRTEMP) = TREDBB(IRTEMP)*EOKRED

CC               DIAGNOSTICS
C                WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDBB(IRTEMP),
C     +                            TMTMBB(IRTEMP),OMTMBB(IRTEMP)

              ENDDO
  
C             CSTAR
              DO IRTEMP = 1, NTRDCC
  
                OMTMCC(IRTEMP) = CSTAR(IRTEMP,IDELTA)
                TMTMCC(IRTEMP) = TREDCC(IRTEMP)*EOKRED

CC               DIAGNOSTICS
C                WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDCC(IRTEMP),
C     +                            TMTMCC(IRTEMP),OMTMCC(IRTEMP)

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

C           RATIOS OF COLLISION INTEGRALS
C           OMEGA(2,2)*/OMEGA(1,1)* = ASTAR
C           OMEGA(1,2)*/OMEGA(1,1)* = CSTAR
C           OMEGA(1,3)*/OMEGA(1,1)* = (5*CSTAR-BSTAR)/4
C           COLLECT THE VALUES AT ZERO DIPOLE MOMENT
C           EVALUATE THE NON-REDUCED TEMPERATURES
            IDELTA = 1

C           ASTAR
            DO IRTEMP = 1, NTRDAA
  
              OMTMAA(IRTEMP) = ASTAR(IRTEMP,IDELTA)
              TMTMAA(IRTEMP) = TREDAA(IRTEMP)*EOKRED

CC             DIAGNOSTICS
C              WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDAA(IRTEMP),
C     +                          TMTMAA(IRTEMP),OMTMAA(IRTEMP)

            ENDDO
  
C           BSTAR
            DO IRTEMP = 1, NTRDBB
  
              OMTMBB(IRTEMP) = BSTAR(IRTEMP,IDELTA)
              TMTMBB(IRTEMP) = TREDBB(IRTEMP)*EOKRED

CC             DIAGNOSTICS
C              WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDBB(IRTEMP),
C     +                          TMTMBB(IRTEMP),OMTMBB(IRTEMP)

            ENDDO
  
C           CSTAR
            DO IRTEMP = 1, NTRDCC
  
              OMTMCC(IRTEMP) = CSTAR(IRTEMP,IDELTA)
              TMTMCC(IRTEMP) = TREDCC(IRTEMP)*EOKRED

CC             DIAGNOSTICS
C              WRITE(6,'(I5,3(1PE12.4))')IRTEMP,TREDCC(IRTEMP),
C     +                          TMTMCC(IRTEMP),OMTMCC(IRTEMP)

            ENDDO
 
C           -------------------------------------------------------------------

          ENDIF

C         =====================================================================

C         STANDARDISE THE TEMPERATURE BASIS FOR ALL THE COLLISION INTEGRALS
C         -----------------------------------------------------------------
C         USING TEMPERATURE VALUES FOR ASTAR AS THE BASIS
          NTRDCR = NTRDAA

C         OMEGA(2,2)*/OMEGA(1,1)* = ASTAR
          DO IRTEMP = 1, NTRDCR
  
            TMVALS(IRTEMP) = TMTMAA(IRTEMP)
            OMVR22(IRTEMP) = OMTMAA(IRTEMP)

          ENDDO
  
C         OMEGA(1,2)*/OMEGA(1,1)* = CSTAR
          DO IRTEMP = 1, NTRDCR
  
            OMVR12(IRTEMP) = CUBINT(TMTMCC,OMTMCC,NTRDCC,TMVALS(IRTEMP))

          ENDDO
 
C         OMEGA(1,3)*/OMEGA(1,1)* = (5*CSTAR-BSTAR)/4
          DO IRTEMP = 1, NTRDCR
  
            OMVR13(IRTEMP) = CUBINT(TMTMBB,OMTMBB,NTRDBB,TMVALS(IRTEMP))
            OMVR13(IRTEMP) = QRTR*(FIVE*OMVR12(IRTEMP) - OMVR13(IRTEMP))

          ENDDO
  
C         LINEAR COMBINATIONS OF COLLISION INTEGRALS
C         AASTAR(ALPHA,BETA),BBSTAR(ALPHA,BETA),CCSTAR(ALPHA,BETA)
C         NOT TO BE CONFUSED WITH ASTAR,BSTAR,CSTAR
          DO IRTEMP = 1, NTRDCR
  
            AASTAR(IRTEMP) = HALF*OMVR22(IRTEMP)
            BBSTAR(IRTEMP) = THRD*(FIVE*OMVR12(IRTEMP) - OMVR13(IRTEMP))
            CCSTAR(IRTEMP) = THRD*OMVR12(IRTEMP)

          ENDDO
  
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
              IF(IRTEMP.LT.NTRDCR)GOTO 1000
            ELSE
              IFITLO = IRTEMP - 1
            ENDIF
C         END OF LOOP 1000
          IFITLO = MAX(IFITLO,1)

C         HIGH END
          IFITHI = NTRDCR + 1
          IRTEMP = NTRDCR + 1
1010      CONTINUE
            IRTEMP = IRTEMP - 1
            IF(TMVALS(IRTEMP).GT.TFITHI)THEN
              IF(IRTEMP.GT.1)GOTO 1010
            ELSE
              IFITHI = IRTEMP + 1
            ENDIF
C         END OF LOOP 1010
          IFITHI = MIN(IFITHI,NTRDCR)

CC         DIAGNOSTICS
C          WRITE(6,'(2I5)')IFITLO,IFITHI

C         COLLECT THE VALUES TO FIT
          IR2FIT = 0
          DO IRTEMP = IFITLO, IFITHI
            IR2FIT = IR2FIT + 1
            TM2FIT(IR2FIT) = TMVALS(IRTEMP)/TREFGB
            THETAR = TWO*AASTAR(IRTEMP) + FIVE
            THETAR = THETAR*(SIX*CCSTAR(IRTEMP) - FIVE)
            FORNOW = SXTEEN*AASTAR(IRTEMP)-TWELVE*BBSTAR(IRTEMP)+FFTYFV
            FORNOW = FORNOW*AASTAR(IRTEMP)
            OM2FIT(IR2FIT) = HALF*FFTEEN*THETAR*REDMAS/FORNOW 
          ENDDO
          NR2FIT = IR2FIT

CC         DIAGNOSTICS
C          DO IRTEMP = 1, NR2FIT
C      WRITE(6,'(2I5,2F10.3)')ISPEC,IRTEMP,TM2FIT(IRTEMP),OM2FIT(IRTEMP)
C          ENDDO

C         FIT THE POLYNOMIAL
          CHISQR = ZERO
          CALL POLFIT(TM2FIT,OM2FIT,NR2FIT,NFITMX,
     +                ACOEFS,NCOTDR,NFCOMX,CHISQR)


CC         DIAGNOSTICS
C          WRITE(6,'(I5,1PE12.4)')ISPEC,CHISQR
C        WRITE(6,'(I5,6(1PE12.4))')ISPEC,(ACOEFS(ICOEFF),ICOEFF=1,NCOTDR)
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = ACOEFS(NCOTDR)
C            DO ICOEFF = NCOTDR-1,1,-1
C              FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C            ENDDO
C        WRITE(6,'(I5,5(1PE12.4))')IR2FIT,TM2FIT(IR2FIT),OM2FIT(IR2FIT),
C     +                            FORNOW,OM2FIT(IR2FIT)-FORNOW,
C     +          (OM2FIT(IR2FIT)-FORNOW)/(OM2FIT(IR2FIT)+1.0D-30)
C          ENDDO

CC         DIAGNOSTICS
CC          IF(WMOLAR(ISPEC).LE.FIVE)THEN
C          NCPLOF = 7
C          WRITE(CSPEC1,'(I2.2)')ISPEC
C          WRITE(CSPEC2,'(I2.2)')JSPEC
C          FNPLOF = 'tdrfit'//CSPEC1//CSPEC2//'.res'
C          OPEN(UNIT=NCPLOF,FILE=FNPLOF,FORM='FORMATTED')
C          DO IR2FIT = 1, NR2FIT
C        WRITE(NCPLOF,'(2(1PE12.4))')TM2FIT(IR2FIT)*TREFGB,OM2FIT(IR2FIT)
C          ENDDO
C          WRITE(NCPLOF,*)
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = ACOEFS(NCOTDR)
C            DO ICOEFF = NCOTDR-1,1,-1
C              FORNOW = FORNOW*TM2FIT(IR2FIT) + ACOEFS(ICOEFF)
C            ENDDO
C            WRITE(NCPLOF,'(2(1PE12.4))')TM2FIT(IR2FIT)*TREFGB,FORNOW
C          ENDDO
C          CLOSE(NCPLOF)
CC          ENDIF

C         =====================================================================

C         COMPLETE THE COEFFICIENTS FOR THE THERMAL DIFFUSION RATIO
          DO ICOEFF = 1,NCOTDR
            TDRCCO(ICOEFF,ISPEC,JSPEC) = ACOEFS(ICOEFF)
          ENDDO

CC         DIAGNOSTICS
C          DO IR2FIT = 1, NR2FIT
C            FORNOW = TDRCCO(NCOTDR,ISPEC,JSPEC)
C            DO ICOEFF = NCOTDR-1,1,-1
C          FORNOW = FORNOW*TM2FIT(IR2FIT) + TDRCCO(ICOEFF,ISPEC,JSPEC)
C            ENDDO
C            WRITE(6,'(2(1PE12.4))')TM2FIT(IR2FIT)*TREFGB,FORNOW
C          ENDDO
C          WRITE(6,*)

C         =====================================================================

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
