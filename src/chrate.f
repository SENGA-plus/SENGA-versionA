      SUBROUTINE CHRATE 
 
C     *************************************************************************
C
C     CHRATE
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     04-JAN-2003:  CREATED
C     22-MAR-2008:  RSC BUG FIX LINDEMANN RATE EXPRESSION
C     04-SEP-2009:  RSC RECODE REACTANT AND PRODUCT COEFFICIENT LISTS
C     13-SEP-2009:  RSC REFORMULATE FOR LINDEMANN RATE EXPRESSIONS
C     04-OCT-2009:  RSC REFORMULATE FOR TROE/SRI RATE EXPRESSIONS
C     08-AUG-2012:  RSC/ZN BUG FIX PRESSURE-DEPENDENT RATES
C      
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES CHEMICAL REACTION RATES
C     USING MULTI-SPECIES MULTI-STEP CHEMISTRY
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
      DOUBLE PRECISION YSMALL,YDENOM
      PARAMETER(YSMALL = 1.0D-30,YDENOM = 1.0D-15)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION RACNST,RNCNST,REOVRR
      DOUBLE PRECISION OVWMAS,GIBBSP,SCOEF,PREDUC
      DOUBLE PRECISION OVTST1,TSTAR2,OVTST3,OMALPH
      DOUBLE PRECISION FBROAD,FTCENT,TRATS1,TRATS2,TRATS3,CFACTR,ENFACT
      DOUBLE PRECISION ACFSRI,BCFSRM,OVCSRM,DCFSRI,ECFSRI
      DOUBLE PRECISION FORNOW
      INTEGER ISPEC,ISSPEC
      INTEGER ISTEP,IBODY
      INTEGER IINDEX,IPOWER,ICOEF1,ICOEF2,ITINT,ICP
      INTEGER IC,JC,KC
      LOGICAL FLTHRD

CVM      IF (IRKSTP==1) THEN
CVM        OPEN(UNIT=113,FILE='chrate.dat',STATUS='NEW',FORM='FORMATTED')
CVM      ENDIF

 

C     BEGIN
C     =====

C     =========================================================================

C     GLOBAL DATA
C     -----------
C     NSPEC total no of species
C     NSTEP total no of steps
C     NBODY total no of third bodies
C     NGIBB total no of Gibbs steps
C     NLIND total no of Lindemann steps
C     NTROE total no of Troe steps
C     NSRIF total no of SRI steps
C     NSPCMX max no of species
C     NSTPMX max no of steps
C     NSSMAX max length of step species-list
C     NRSMAX max length of step reactant- and product-lists
C     NBDYMX max no of third bodies
C     NLLMAX max no of Lindemann steps
C     RPARAM(1,NSTEP) pre-exponential factor ln(A)
C     RPARAM(2,NSTEP) temperature exponent n
C     RPARAM(3,NSTEP) activation energy E/R0
C     NSSPEC(NSSMAX,NSTEP) is the step species-list
C     NRSPEC(NRSMAX,NSTEP) is the step reactant-list
C     NPSPEC(NRSMAX,NSTEP) is the step product-list
C     NRCPEC(NRSMAX,NSTEP) is the step coefficient reactant-list
C     NPCPEC(NRSMAX,NSTEP) is the step coefficient product-list
C     NSSLEN(NSTPMX) contains the length of the species-list for each step
C     NRSLEN(NSTPMX) contains the length of the reactant-list for each step
C     NPSLEN(NSTPMX) contains the length of the product-list for each step
C     NRCLEN(NSTPMX) contains the length of the reactant coefficient-list for
C                    each step
C     NPCLEN(NSTPMX) contains the length of the product coefficient-list for
C                    each step
C     WMOLAR(NSPCMX) molar mass of species
C     OVWMOL(NSPCMX) reciprocal of molar mass of species
C     DIFFMU(NSSMAX,NSTPMX) is the species delta-mu for each step
C     DIFFMW(NSSMAX,NSTPMX) is the species (delta-mu times molar mass)
C                           for each step
C     CRSPEC(NSSMAX,NSTPMX) is the reactant stoichiometric coefficient-list
C                           for each step
C     CPSPEC(NSSMAX,NSTPMX) is the product stoichiometric coefficient-list
C                           for each step
C     MBLIST(NSTPMX) contains 0 if no third body in this step
C                    else contains the index number of the third body
C     MGLIST(NSTPMX) contains 0 if no Gibbs function evaluation in this step
C                    else contains the index number of the Gibbs step
C     MLLIST(NSTPMX) contains 0 if no Lindemann rate in this step
C                    else contains the index number of the Lindemann step
C     MTLIST(NSTPMX) contains 0 if no Troe rate in this step
C                    else contains the index number of the Troe step
C     MSLIST(NSTPMX) contains 0 if no SRI rate in this step
C                    else contains the index number of the SRI step
C     EFFY3B(NSPEC,NBDYMX) contains the third-body efficiencies
C     RCLIND(1,MLLIST(NSTPMX)) Lindemann rate pre-exponential factor ln(A)
C     RCLIND(2,MLLIST(NSTPMX)) Lindemann rate temperature exponent n
C     RCLIND(3,MLLIST(NSTPMX)) Lindemann rate activation energy E/R0
C     RCLIND(4,MLLIST(NSTPMX)) Lindemann rate broadening factor
C     RCTROE(1,MTLIST(NSTPMX)) Troe form rate pre-exponential factor ln(A)
C     RCTROE(2,MTLIST(NSTPMX)) Troe form rate temperature exponent n
C     RCTROE(3,MTLIST(NSTPMX)) Troe form rate activation energy E/R0
C     RCTROE(4,MTLIST(NSTPMX)) Troe form rate temperature parameter alpha
C     RCTROE(5,MTLIST(NSTPMX)) Troe form rate temperature parameter no 1*
C     RCTROE(6,MTLIST(NSTPMX)) Troe form rate temperature parameter no 2*
C     RCTROE(7,MTLIST(NSTPMX)) Troe form rate temperature parameter no 3*
C     RCTROE(8,MTLIST(NSTPMX)) Troe form rate c-const no 1
C     RCTROE(9,MTLIST(NSTPMX)) Troe form rate c-const no 2
C     RCTROE(10,MTLIST(NSTPMX)) Troe form rate n-const no 1
C     RCTROE(11,MTLIST(NSTPMX)) Troe form rate n-const no 2
C     RCTROE(12,MTLIST(NSTPMX)) Troe form rate d-const
C     RCSRIF(1,MSLIST(NSTPMX)) SRI form rate pre-exponential factor ln(A)
C     RCSRIF(2,MSLIST(NSTPMX)) SRI form rate temperature exponent n
C     RCSRIF(3,MSLIST(NSTPMX)) SRI form rate activation energy E/R0
C     RCSRIF(4,MSLIST(NSTPMX)) SRI form rate parameter a
C     RCSRIF(5,MSLIST(NSTPMX)) SRI form rate parameter b
C     RCSRIF(6,MSLIST(NSTPMX)) SRI form rate parameter c
C     RCSRIF(7,MSLIST(NSTPMX)) SRI form rate parameter d
C     RCSRIF(8,MSLIST(NSTPMX)) SRI form rate parameter e

C     =========================================================================

C     ZERO THE REACTION RATE ACCUMULATOR ARRAYS
      DO ISPEC = 1,NSPEC

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              RATE(IC,JC,KC,ISPEC) = ZERO

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     =========================================================================

C     RUN THROUGH ALL STEPS IN THE REACTION MECHANISM
C     -----------------------------------------------

      DO ISTEP = 1, NSTEP

C       =======================================================================

C       INCLUDE THIRD BODY CONCENTRATIONS, IF ANY
C       -----------------------------------------

C       SET THIRD-BODY EVALUATION FLAG
C       RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
        FLTHRD = .TRUE.

C       CHECK THIRD-BODY LIST
        IF(MBLIST(ISTEP).NE.0)THEN

C         GET THE INDEX OF THE THIRD BODY
          IBODY = MBLIST(ISTEP)

C         ZERO THE THIRD-BODY CONCENTRATION ACCUMULATOR
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                STORE3(IC,JC,KC) = ZERO

              ENDDO
            ENDDO
          ENDDO

C         USE THE THIRD-BODY SPECIES-LIST
C         AND THE THIRD-BODY EFFICIENCY LIST
C         TO GET THE THIRD-BODY CONCENTRATION
          DO ISPEC = 1,NSPEC

C           GET THE THIRD-BODY EFFICIENCY FROM THE LIST
C           AND COMBINE WITH THE RECIPROCAL MOLAR MASS
            OVWMAS = OVWMOL(ISPEC)*EFFY3B(ISPEC,IBODY)

C           EVALUATE THE THIRD BODY CONCENTRATION
C           AND STORE IN TEMPORARY ARRAY
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  STORE3(IC,JC,KC) = STORE3(IC,JC,KC)
     +                             + YRHS(IC,JC,KC,ISPEC)*OVWMAS

                ENDDO
              ENDDO
            ENDDO

          ENDDO
C         THIRD-BODY SPECIES-LIST

CC         DIAGNOSTICS
C          WRITE(6,*)'3body',ISTEP,IBODY,FLTHRD,STORE3(250,1,1)

        ENDIF
C                                             STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       EVALUATE THE SPECIFIC REACTION RATE CONSTANT FOR THIS STEP
C       ----------------------------------------------------------
        RACNST = RPARAM(1,ISTEP)
        RNCNST = RPARAM(2,ISTEP)
        REOVRR = RPARAM(3,ISTEP)
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STORE2(IC,JC,KC) = RACNST
     +                         + RNCNST*LOG(TRUN(IC,JC,KC))
     +                         - REOVRR/TRUN(IC,JC,KC)
              STORE1(IC,JC,KC) = EXP(STORE2(IC,JC,KC))


            ENDDO
          ENDDO
        ENDDO

CC       DIAGNOSTICS
C        WRITE(6,'("kf",I5,5(1PE12.4))')
C     +    ISTEP,RACNST,RNCNST,REOVRR,STORE1(250,1,1),STORE2(250,1,1)

C                                            STORE1 = FORWARD RATE CONSTANT
C                                            STORE2 = LN(FORWARD RATE CONSTANT)
C                                            STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       INCLUDE LINDEMANN RATE EVALUATION, IF ANY
C       -----------------------------------------
C       RSC 13-SEP-2009 REFORMULATE FOR LINDEMANN RATE EXPRESSIONS
        IF(MLLIST(ISTEP).NE.0)THEN

          RACNST = RCLIND(1,MLLIST(ISTEP))
          RNCNST = RCLIND(2,MLLIST(ISTEP))
          REOVRR = RCLIND(3,MLLIST(ISTEP))
          FLCNST = RCLIND(4,MLLIST(ISTEP))

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               EVALUATE K0
                STORE4(IC,JC,KC) = RACNST
     +                           + RNCNST*LOG(TRUN(IC,JC,KC))
     +                           - REOVRR/TRUN(IC,JC,KC)
                STORE4(IC,JC,KC) = EXP(STORE4(IC,JC,KC))

C               EVALUATE REDUCED PRESURE
                PREDUC = STORE3(IC,JC,KC)*STORE4(IC,JC,KC)
     +                                   /STORE1(IC,JC,KC)

C               EVALUATE UPDATED FORWARD RATE CONSTANT
                FORNOW = FLCNST*PREDUC/(ONE+PREDUC)
                STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*FORNOW

C               RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
                STORE2(IC,JC,KC) = LOG(STORE1(IC,JC,KC))

CC               DIAGNOSTICS
C                WRITE(6,'("L",2I5,5(1PE12.4))')ISTEP,IC,
C     +          STORE1(IC,JC,KC),STORE2(IC,JC,KC),STORE3(IC,JC,KC),
C     +          STORE4(IC,JC,KC),FORNOW

              ENDDO
            ENDDO
          ENDDO

C         RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
C         RESET THIRD BODY EVALUATION FLAG
          FLTHRD = .FALSE.

CC         DIAGNOSTICS
C          WRITE(6,'("L",I5,6(1PE12.4))')
C     +      ISTEP,FLCNST,RACNST,RNCNST,REOVRR,
C     +            STORE1(250,1,1),STORE2(250,1,1)

        ENDIF
C       END OF LINDEMANN RATE EVALUATION
C                                            STORE1 = FORWARD RATE CONSTANT
C                                            STORE2 = LN(FORWARD RATE CONSTANT)
C                                            STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       INCLUDE TROE FORM RATE EVALUATION, IF ANY
C       -----------------------------------------
C       RSC 04-OCT-2009
        IF(MTLIST(ISTEP).NE.0)THEN

          RACNST = RCTROE(1,MTLIST(ISTEP))
          RNCNST = RCTROE(2,MTLIST(ISTEP))
          REOVRR = RCTROE(3,MTLIST(ISTEP))
          TALPHA = RCTROE(4,MTLIST(ISTEP))
          OVTST1 = RCTROE(5,MTLIST(ISTEP))
          TSTAR2 = RCTROE(6,MTLIST(ISTEP))
          OVTST3 = RCTROE(7,MTLIST(ISTEP))
          CFCST1 = RCTROE(8,MTLIST(ISTEP))
          CFCST2 = RCTROE(9,MTLIST(ISTEP))
          ENCST1 = RCTROE(10,MTLIST(ISTEP))
          ENCST2 = RCTROE(11,MTLIST(ISTEP))
          DTCNST = RCTROE(12,MTLIST(ISTEP))

          OMALPH = ONE - TALPHA

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               EVALUATE K0
                STORE4(IC,JC,KC) = RACNST
     +                           + RNCNST*LOG(TRUN(IC,JC,KC))
     +                           - REOVRR/TRUN(IC,JC,KC)
                STORE4(IC,JC,KC) = EXP(STORE4(IC,JC,KC))

C               EVALUATE REDUCED PRESURE
                PREDUC = STORE3(IC,JC,KC)*STORE4(IC,JC,KC)
     +                                   /STORE1(IC,JC,KC)

C               EVALUATE FCENT
                TRATS1 = TRUN(IC,JC,KC)*OVTST1
                TRATS2 = -TSTAR2/TRUN(IC,JC,KC)
C                TRATS2 = TSTAR2/TRUN(IC,JC,KC)
                TRATS3 = TRUN(IC,JC,KC)*OVTST3
                FTCENT = OMALPH*EXP(TRATS3)
     +                 + TALPHA*EXP(TRATS1)
     +                 + EXP(TRATS2)
                FTCENT = LOG10(FTCENT)

C               EVALUATE BROADENING FACTOR
                CFACTR = CFCST1 + CFCST2*FTCENT
                ENFACT = ENCST1 + ENCST2*FTCENT
                FORNOW = LOG10(PREDUC) + CFACTR
                FORNOW = FORNOW/(ENFACT - DTCNST*FORNOW)
                FORNOW = ONE + FORNOW*FORNOW
                FBROAD = FTCENT/FORNOW
                FBROAD = EXP(FBROAD*CLNTEN)

C               EVALUATE UPDATED FORWARD RATE CONSTANT
                FORNOW = FBROAD*PREDUC/(ONE+PREDUC)
                STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*FORNOW

C               RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
                STORE2(IC,JC,KC) = LOG(STORE1(IC,JC,KC))

              ENDDO
            ENDDO
          ENDDO

C         RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
C         RESET THIRD BODY EVALUATION FLAG
          FLTHRD = .FALSE.

        ENDIF
C       END OF TROE FORM RATE EVALUATION
C                                            STORE1 = FORWARD RATE CONSTANT
C                                            STORE2 = LN(FORWARD RATE CONSTANT)
C                                            STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       INCLUDE SRI FORM RATE EVALUATION, IF ANY
C       ----------------------------------------
C       RSC 04-OCT-2009
        IF(MSLIST(ISTEP).NE.0)THEN

          RACNST = RCSRIF(1,MSLIST(ISTEP))
          RNCNST = RCSRIF(2,MSLIST(ISTEP))
          REOVRR = RCSRIF(3,MSLIST(ISTEP))
          ACFSRI = RCSRIF(4,MSLIST(ISTEP))
          BCFSRM = RCSRIF(5,MSLIST(ISTEP))
          OVCSRM = RCSRIF(6,MSLIST(ISTEP))
          DCFSRI = RCSRIF(7,MSLIST(ISTEP))
          ECFSRI = RCSRIF(8,MSLIST(ISTEP))

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               EVALUATE K0
                STORE4(IC,JC,KC) = RACNST
     +                           + RNCNST*LOG(TRUN(IC,JC,KC))
     +                           - REOVRR/TRUN(IC,JC,KC)
                STORE4(IC,JC,KC) = EXP(STORE4(IC,JC,KC))

C               EVALUATE REDUCED PRESURE
                PREDUC = STORE3(IC,JC,KC)*STORE4(IC,JC,KC)
     +                                   /STORE1(IC,JC,KC)
C               EVALUATE FCENT
                FORNOW = BCFSRM/TRUN(IC,JC,KC)
                FBROAD = ACFSRI*EXP(FORNOW)
                FORNOW = TRUN(IC,JC,KC)*OVCSRM
                FBROAD = FBROAD + EXP(FORNOW)
                FORNOW = LOG10(PREDUC)
                FORNOW = ONE/(ONE+LOG10(FORNOW))
                FBROAD = DCFSRI*EXP(FORNOW*LOG(FBROAD))
                FORNOW = FBROAD*EXP(ECFSRI*LOG(TRUN(IC,JC,KC)))

C               EVALUATE UPDATED FORWARD RATE CONSTANT
                FORNOW = FBROAD*PREDUC/(ONE+PREDUC)
                STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*FORNOW

C               RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
                STORE2(IC,JC,KC) = LOG(STORE1(IC,JC,KC))

              ENDDO
            ENDDO
          ENDDO

C         RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
C         RESET THIRD BODY EVALUATION FLAG
          FLTHRD = .FALSE.

        ENDIF
C       END OF SRI FORM RATE EVALUATION
C                                            STORE1 = FORWARD RATE CONSTANT
C                                            STORE2 = LN(FORWARD RATE CONSTANT)
C                                            STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       USE STEP REACTANT-LIST
C       ----------------------
C       TO FIND REACTANT SPECIES FOR THIS STEP
C       HAVING POSITIVE INTEGER STOICHIOMETRIC COEFFICIENTS

C       NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP REACTANT-LIST
C       AND ARE TREATED SEPARATELY (SEE BELOW)

        DO ISSPEC = 1, NRSLEN(ISTEP)

C         GET THE SPECIES NUMBER FROM THE REACTANT-LIST
          ISPEC = NRSPEC(ISSPEC,ISTEP)

C         EVALUATE REACTANT CONCENTRATIONS AND MULTIPLY UP
C         CONCENTRATION IS RHO*Y/WMOLAR
C         YRHS CONTAINS RHO*Y
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = MAX(YRHS(IC,JC,KC,ISPEC)*OVWMOL(ISPEC),ZERO)
                STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*FORNOW

              ENDDO
            ENDDO
          ENDDO

        ENDDO
C       END OF STEP REACTANT-LIST
C                                            STORE1 = FORWARD REACTION RATE
C                                            STORE2 = LN(FORWARD RATE CONSTANT)
C                                            STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       USE STEP REACTANT COEFFICIENT-LIST
C       ----------------------------------
C       TO FIND REACTANT SPECIES FOR THIS STEP
C       HAVING NON-POSITIVE-INTEGER STOICHIOMETRIC COEFFICIENTS
C       AND OBTAIN THE COEFFICIENTS
C       RSC 04-SEP-2009 RECODE REACTANT AND PRODUCT COEFFICIENT LISTS

C       NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP REACTANT-LIST
C       AND ARE TREATED SEPARATELY (SEE BELOW)

        IF(NRCLEN(ISTEP).GT.0)THEN

          DO ISSPEC = 1, NRCLEN(ISTEP)

C           GET THE SPECIES NUMBER FROM THE STEP COEFFICIENT REACTANT-LIST
            ISPEC = NRCPEC(ISSPEC,ISTEP)

C           GET THE REACTANT STOICHIOMETRIC COEFFICIENT
            SCOEF = CRSPEC(ISSPEC,ISTEP)

C           EVALUATE REACTANT CONCENTRATIONS AND MULTIPLY UP
C           CONCENTRATION IS RHO*Y/WMOLAR
C           YRHS CONTAINS RHO*Y
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  FORNOW = YRHS(IC,JC,KC,ISPEC)*OVWMOL(ISPEC)
                  FORNOW = MAX(FORNOW,YSMALL)
                  FORNOW = EXP(SCOEF*LOG(FORNOW))
                  FORNOW = FORNOW/(ONE+YDENOM*FORNOW)
                  STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*FORNOW

                ENDDO
              ENDDO
            ENDDO

          ENDDO

        ENDIF
C       END OF STEP REACTANT COEFFICIENT-LIST
C                                            STORE1 = FORWARD REACTION RATE
C                                            STORE2 = LN(FORWARD RATE CONSTANT)
C                                            STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       INCLUDE BACKWARD RATE EVALUATION, IF ANY
C       ----------------------------------------
        IF(MGLIST(ISTEP).NE.0)THEN

C         =====================================================================

C         GIBBS FUNCTION EVALUATION
C         -------------------------

C         RUN THROUGH ALL SPECIES FOR THIS STEP
          DO ISSPEC = 1, NSSLEN(ISTEP)

C           GET THE SPECIES NUMBER FROM THE STEP SPECIES-LIST
            ISPEC = NSSPEC(ISSPEC,ISTEP) 

C           EVALUATE GIBBS FUNCTION FOR EACH SPECIES
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

C                 LOCATE THE TEMPERATURE IN AN INTERVAL
                  IINDEX = 1 + (ISPEC-1)/NSPIMX
                  IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                  ICOEF2 = NTBASE**IPOWER
                  ICOEF1 = ICOEF2*NTBASE
                  ITINT = 1 + MOD(ITNDEX(IC,JC,KC,IINDEX),ICOEF1)/ICOEF2

C                 CONSTRUCT GIBBS FUNCTION FROM ITS POLYNOMIAL COEFFICIENTS
                  GIBBSP = AMOLGB(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                  DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                    GIBBSP = AMOLGB(ICP,ITINT,ISPEC)
     +                     + GIBBSP*TRUN(IC,JC,KC)
                  ENDDO
                  GIBBSP = AMOLGB(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                    /TRUN(IC,JC,KC)
     +                   - AMOLGB(NCENPY(ITINT,ISPEC),ITINT,ISPEC)
     +                    *LOG(TRUN(IC,JC,KC))
     +                   - GIBBSP

C                 ADD GIBBS FUNCTION CONTRIBUTION TO RATE COEFFICIENT
C                 USING STEP SPECIES DELTA-LIST
C                 TO GET BACKWARD RATE COEFFICIENT
                  STORE2(IC,JC,KC) = STORE2(IC,JC,KC)
     +                             + DIFFMU(ISSPEC,ISTEP)*GIBBSP

                ENDDO
              ENDDO
            ENDDO

          ENDDO
C         STEP SPECIES-LIST
C                                           STORE1 = FORWARD REACTION RATE
C                                           STORE2 = LN(BACKWARD RATE CONSTANT)
C                                           STORE3 = THIRD BODY CONCENTRATION
C         =====================================================================

C         FINALISE BACKWARD RATE COEFFICIENT
C         ----------------------------------
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                STORE2(IC,JC,KC) = EXP(STORE2(IC,JC,KC))

              ENDDO
            ENDDO
          ENDDO
C                                             STORE1 = FORWARD REACTION RATE
C                                             STORE2 = BACKWARD RATE CONSTANT
C                                             STORE3 = THIRD BODY CONCENTRATION
C         =====================================================================

C         USE STEP PRODUCT-LIST
C         ---------------------
C         TO FIND PRODUCT SPECIES FOR THIS STEP
C         HAVING POSITIVE INTEGER STOICHIOMETRIC COEFFICIENTS

C         NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP PRODUCT-LIST
C         AND ARE TREATED SEPARATELY (SEE BELOW)

          DO ISSPEC = 1, NPSLEN(ISTEP)

C           GET THE SPECIES NUMBER FROM THE PRODUCT-LIST
            ISPEC = NPSPEC(ISSPEC,ISTEP)

C           EVALUATE PRODUCT CONCENTRATIONS AND MULTIPLY UP
C           CONCENTRATION IS RHO*Y/WMOLAR
C           YRHS CONTAINS RHO*Y

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  FORNOW = MAX(YRHS(IC,JC,KC,ISPEC)*OVWMOL(ISPEC),ZERO)
                  STORE2(IC,JC,KC) = STORE2(IC,JC,KC)*FORNOW

                ENDDO
              ENDDO
            ENDDO

          ENDDO
C         END OF STEP PRODUCT-LIST
C                                             STORE1 = FORWARD REACTION RATE
C                                             STORE2 = BACKWARD REACTION RATE
C                                             STORE3 = THIRD BODY CONCENTRATION
C         =====================================================================

C         USE STEP PRODUCT COEFFICIENT-LIST
C         ---------------------------------
C         TO FIND PRODUCT SPECIES FOR THIS STEP
C         HAVING NON-POSITIVE-INTEGER STOICHIOMETRIC COEFFICIENTS
C         AND OBTAIN THE COEFFICIENTS
C         RSC 04-SEP-2009 RECODE REACTANT AND PRODUCT COEFFICIENT LISTS

C         NOTE: THIRD BODIES ARE NOT INCLUDED IN THE STEP PRODUCT-LIST
C         AND ARE TREATED SEPARATELY (SEE BELOW)

          IF(NPCLEN(ISTEP).GT.0)THEN

            DO ISSPEC = 1, NPCLEN(ISTEP)

C             GET THE SPECIES NUMBER FROM THE STEP COEFFICIENT PRODUCT-LIST
              ISPEC = NPCPEC(ISSPEC,ISTEP)

C             GET THE PRODUCT STOICHIOMETRIC COEFFICIENT
              SCOEF = CPSPEC(ISSPEC,ISTEP)

C             EVALUATE PRODUCT CONCENTRATIONS AND MULTIPLY UP
C             CONCENTRATION IS RHO*Y/WMOLAR
C             YRHS CONTAINS RHO*Y
              DO KC = KSTAL,KSTOL
                DO JC = JSTAL,JSTOL
                  DO IC = ISTAL,ISTOL

                    FORNOW = YRHS(IC,JC,KC,ISPEC)*OVWMOL(ISPEC)
                    FORNOW = MAX(FORNOW,YSMALL)
                    FORNOW = EXP(SCOEF*LOG(FORNOW))
                    FORNOW = FORNOW/(ONE+YDENOM*FORNOW)
                    STORE2(IC,JC,KC) = STORE2(IC,JC,KC)*FORNOW

                  ENDDO
                ENDDO
              ENDDO

            ENDDO

          ENDIF
C         END OF STEP PRODUCT COEFFICIENT-LIST
C                                             STORE1 = FORWARD REACTION RATE
C                                             STORE2 = BACKWARD REACTION RATE
C                                             STORE3 = THIRD BODY CONCENTRATION
C         =====================================================================

C         EVALUATE NET FORWARD REACTION RATE
C         ----------------------------------
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                STORE1(IC,JC,KC) = STORE1(IC,JC,KC) - STORE2(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO

C         =====================================================================

        ENDIF
C       END OF BACKWARD RATE EVALUATION
C                                            STORE1 = NET FORWARD REACTION RATE
C                                            STORE3 = THIRD BODY CONCENTRATION
C       =======================================================================

C       INCLUDE THIRD BODY CONCENTRATIONS, IF ANY
C       -----------------------------------------

C       CHECK THIRD-BODY LIST
        IF(MBLIST(ISTEP).NE.0)THEN

C         RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
C         CHECK THIRD-BODY EVALUATION FLAG
          IF(FLTHRD)THEN

C           INCLUDE THIRD BODY CONCENTRATION IN REACTION RATE FOR THE STEP
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*STORE3(IC,JC,KC)

                ENDDO
              ENDDO
            ENDDO

          ENDIF
C         END OF THIRD BODY EVALUATION FLAG

        ENDIF
C       END OF THIRD BODY TREATMENT
C                                            STORE1 = NET FORWARD REACTION RATE
C       =======================================================================

C       USE STEP SPECIES-LIST AND STEP SPECIES DELTA-LIST
C       -------------------------------------------------
C       TO GET DIFFERENCE OF STOICHIOMETRIC COEFFICIENTS FOR EACH SPECIES 
C       (ACTUALLY DELTA-MU TIMES WMOLAR)
C       ADD STEP CONTRIBUTION TO TOTAL RATE FOR EACH SPECIES INVOLVED

        DO ISSPEC = 1, NSSLEN(ISTEP)

C         GET THE SPECIES NUMBER FROM THE STEP SPECIES-LIST
          ISPEC = NSSPEC(ISSPEC,ISTEP) 
 
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC)
     +              + STORE1(IC,JC,KC)*DIFFMW(ISSPEC,ISTEP)

              ENDDO
            ENDDO
          ENDDO

        ENDDO
C       STEP SPECIES-LIST
CVM        IF (IRKSTP==1) THEN
CVM          DO IC=ISTAL,ISTOL
CVM            DO ISPEC=1,NSPEC
CVM              WRITE(113,*)IC,ISPEC,ISTEP,RATE(IC,1,1,ISPEC)
CVM            ENDDO
CVM          ENDDO
CVM        ENDIF


C     =========================================================================

      ENDDO
CVM      IF (IRKSTP==1) THEN
CVM        CLOSE(113)
CVM      ENDIF


C     END OF RUN THROUGH ALL STEPS IN THE REACTION MECHANISM
C     ------------------------------------------------------

C     =========================================================================


      RETURN
      END
