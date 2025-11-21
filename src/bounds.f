      SUBROUTINE BOUNDS 
 
C     *************************************************************************
C
C     BOUNDS
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     13-JUL-2003:  RSC MODIFIED FOR SENGA2
C     08-AUG-2012:  RSC EVALUATE ALL SPECIES
C     26-OCT-2013:  RSC ACTIVATE ALL BCS ON ALL SIDES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES CHARACTERISTIC BOUNDARY CONDITIONS FOR ALL PDES
C
C     *************************************************************************

! VISHNU MOHAN (VM), KHALIL ABO AMSHA (KAA) AND NILANJAN CHAKRABORTY (NC)
! IMPLEMENTED THE TRANSVERSE DERIVATIVE TERMS FOR BOUNDARY CONDITION AS PER METHOD
! OF YOO AND IM


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FORNOW
      INTEGER IC,JC,KC
      INTEGER ISPEC
C     THESE ARE INTRODUCED FOR LODATO'S BC EQN. 3.80- NC
      DOUBLE PRECISION TT1XL(NYSIZE,NZSIZE),
     +                 TT2XL(NYSIZE,NZSIZE),
     +                 TT3XL(NYSIZE,NZSIZE),
     +                 TT4XL(NYSIZE,NZSIZE),
     +                 TT5XL(NYSIZE,NZSIZE),
     +                 TT6XL(NYSIZE,NZSIZE,NSPEC)
      DOUBLE PRECISION TT1XR(NYSIZE,NZSIZE),
     +                 TT2XR(NYSIZE,NZSIZE),
     +                 TT3XR(NYSIZE,NZSIZE),
     +                 TT4XR(NYSIZE,NZSIZE),
     +                 TT5XR(NYSIZE,NZSIZE),
     +                 TT6XR(NYSIZE,NZSIZE,NSPEC)
      !UMAIRS CORRECTION HERE
      DOUBLE PRECISION TT1YL(NXSIZE,NZSIZE),
     +                 TT2YL(NXSIZE,NZSIZE),
     +                 TT3YL(NXSIZE,NZSIZE),
     +                 TT4YL(NXSIZE,NZSIZE),
     +                 TT5YL(NXSIZE,NZSIZE),
     +                 TT6YL(NXSIZE,NZSIZE,NSPEC)
      DOUBLE PRECISION TT1YR(NXSIZE,NZSIZE),
     +                 TT2YR(NXSIZE,NZSIZE),
     +                 TT3YR(NXSIZE,NZSIZE),
     +                 TT4YR(NXSIZE,NZSIZE),
     +                 TT5YR(NXSIZE,NZSIZE),
     +                 TT6YR(NXSIZE,NZSIZE,NSPEC)
      DOUBLE PRECISION TT1ZL(NXSIZE,NYSIZE),
     +                 TT2ZL(NXSIZE,NYSIZE),
     +                 TT3ZL(NXSIZE,NYSIZE),
     +                 TT4ZL(NXSIZE,NYSIZE),
     +                 TT5ZL(NXSIZE,NYSIZE),
     +                 TT6ZL(NXSIZE,NYSIZE,NSPEC)
      DOUBLE PRECISION TT1ZR(NXSIZE,NYSIZE),
     +                 TT2ZR(NXSIZE,NYSIZE),
     +                 TT3ZR(NXSIZE,NYSIZE),
     +                 TT4ZR(NXSIZE,NYSIZE),
     +                 TT5ZR(NXSIZE,NYSIZE),
     +                 TT6ZR(NXSIZE,NYSIZE,NSPEC)
      DOUBLE PRECISION BET
      INTEGER FLAG_BET_XL,FLAG_BET_XR,FLAG_BET_YL,FLAG_BET_YR,
     +        FLAG_BET_ZL,FLAG_BET_ZR
      !VM: SWITCH TO SET BET=LOCAL_MACH
      !FLAG_BET=0 POINTS TO GLOBAL VALUE OF BET
      !FLAG_BET=1 POINTS TO LOCAL VALUE OF BET

      INTEGER FLAG_PIO_XL,FLAG_PIO_XR,FLAG_PIO_YL,FLAG_PIO_YR,
     +        FLAG_PIO_ZL,FLAG_PIO_ZR
      !VM: SWITCH TO TURN POINTWISE INFLOW-OUTFLOW ON/OFF
      !1 IS ON, 0 IS OFF

C     FY - FOR NR INFLOW
      DOUBLE PRECISION M2MAX! = 0.001D0 ! SET VALUE OF MAXIMUM MACH NUMBER SQUARED  
      DOUBLE PRECISION NRIETA1,NRIETA2,NRIETA3,NRIETA4,NRIETA5,NRIETA6
      DOUBLE PRECISION STLUXL(NYSIZE,NZSIZE),STLVXL(NYSIZE,NZSIZE),
     +                 STLWXL(NYSIZE,NZSIZE),STLTXL(NYSIZE,NZSIZE)
      DOUBLE PRECISION STLYXL(NYSIZE,NZSIZE,NSPCMX)

C     SET VALUE OF MAXIMUM MACH NUMBER SQUARED  
      M2MAX = 0.001D0

      BET=0.050D0

C     BEGIN
C     =====

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================


C     X-DIRECTION LEFT-HAND END
C     -------------------------
      IF(FXLCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUXL = PRIMITIVE U-VELOCITY COMPONENT
C       STRVXL = PRIMITIVE V-VELOCITY COMPONENT
C       STRWXL = PRIMITIVE W-VELOCITY COMPONENT
C       STRPXL = PRESSURE
C       STRDXL = DENSITY
C       STRTXL = TEMPERATURE
C       STREXL = INTERNAL ENERGY
C       STRGXL = MIXTURE CP
C       STRRXL = MIXTURE SPECIFIC GAS CONSTANT
C       STRYXL(ISPEC) = SPECIES MASS FRACTION
C       RATEXL(ISPEC) = SPECIES REACTION RATE
C       STRHXL(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1XL = DUDX
C       BCL2XL = DRHODX
C       BCL3XL = DVDX
C       BCL4XL = DWDX
C       BCL5XL = DPDX
C       BCLYXL(ISPEC) = DYDX

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRHXL(JC,KC,ISPEC) = STRHXL(JC,KC,ISPEC)
     +         - STRGXL(JC,KC)*STRTXL(JC,KC)*RGSPEC(ISPEC)/STRRXL(JC,KC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            GAM1XL(JC,KC) = STRGXL(JC,KC) - STRRXL(JC,KC)
            STREXL(JC,KC) = STREXL(JC,KC) - GAM1XL(JC,KC)*STRTXL(JC,KC)

            GAM1XL(JC,KC) = STRRXL(JC,KC)/GAM1XL(JC,KC)
            OVGMXL(JC,KC) = ONE/GAM1XL(JC,KC)

          

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            FORNOW = STRGXL(JC,KC)*GAM1XL(JC,KC)*STRTXL(JC,KC)
            ACOUXL(JC,KC) = SQRT(FORNOW)
            OVA2XL(JC,KC) = ONE/FORNOW

C     THESE ARE INTRODUCED FOR LODATO'S BC- NC

            TT1XL(JC,KC)=T51BXL(JC,KC)+
     +                   T52BXL(JC,KC)*(GAM1XL(JC,KC)+1.0)-
     +        STRDXL(JC,KC)* ACOUXL(JC,KC)*T2BXL(JC,KC)         

            TT2XL(JC,KC)=ACOUXL(JC,KC)*ACOUXL(JC,KC)*
     +      T1BXL(JC,KC)-T51BXL(JC,KC)-(GAM1XL(JC,KC)+1.0)*
     +      T52BXL(JC,KC)

            TT3XL(JC,KC)=T3BXL(JC,KC)
            TT4XL(JC,KC)=T4BXL(JC,KC)
            
            TT5XL(JC,KC)=T51BXL(JC,KC)+
     +                   T52BXL(JC,KC)*(GAM1XL(JC,KC)+1.0)+
     +        STRDXL(JC,KC)* ACOUXL(JC,KC)*T2BXL(JC,KC)   

            DO IS=1,NSPEC
            TT6XL(JC,KC,IS)=T6BXL(JC,KC,IS)
            END DO

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCXL.EQ.NSBCO1)THEN
          FLAG_PIO_XL=NXLPRM(1)!0
          FLAG_BET_XL=NXLPRM(2)!1

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXL(JC,KC) = SORPXL(JC,KC)
     +        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L5X AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                      *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
C             OLD VALUE OF L5X
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +        *(BCL5XL(JC,KC)+STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC))

C             SUBTRACT FROM NEW VALUE OF L5X
C             TERM INTRODUCED FOR LODATO'S BC- NC
              IF(FLAG_BET_XL==1) THEN
                BET=STRUXL(JC,KC)*STRUXL(JC,KC)+
     +              STRVXL(JC,KC)*STRVXL(JC,KC)+
     +              STRWXL(JC,KC)*STRWXL(JC,KC)
                BET=SQRT(BET)/ACOUXL(JC,KC)
              END IF
              IF((STRUXL(JC,KC).GT.ZERO).AND.(FLAG_PIO_XL==1))THEN

C             SUBTRACT FROM NEW VALUE OF L5X
              BCL2XL(JC,KC) =
     +                       -BCL2XL(JC,KC)
     +                       -OVA2XL(JC,KC)*SORPXL(JC,KC)

              BCL3XL(JC,KC) =
     +                        0.1*(STRVXL(JC,KC)-0.0D0)
     +                       -BCL3XL(JC,KC)

              BCL4XL(JC,KC) =
     +                        0.1*(STRWXL(JC,KC)-0.0D0)
     +                       -BCL4XL(JC,KC)
              BCL5XL(JC,KC)= HALF*SORPXL(JC,KC)
     +                     + COBCXL*ACOUXL(JC,KC)*(STRPXL(JC,KC)-PINFXL)
     +                     +0.5*(1.0-BET)*TT5XL(JC,KC)- BCL5XL(JC,KC)
     +                     + 100.0*STRDXL(JC,KC)*(ONE/XGDLEN)
     +               *(ACOUXL(JC,KC)**2.0-STRUXL(JC,KC)**2.0)
     +               *(STRUXL(JC,KC)-0.0)
              ELSE

              BCL5XL(JC,KC)= HALF*SORPXL(JC,KC)
     +                     + COBCXL*ACOUXL(JC,KC)*(STRPXL(JC,KC)-PINFXL)
     +                     +0.5*(1.0-BET)*TT5XL(JC,KC)- BCL5XL(JC,KC)
              ENDIF

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              IF((STRUXL(JC,KC).GT.ZERO).AND.(FLAG_PIO_XL==1))THEN

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

              URHS(ISTAL,JC,KC) = URHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRUXL(JC,KC)
     +      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*(STRUXL(JC,KC)+ACOUXL(JC,KC))

              VRHS(ISTAL,JC,KC) = VRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRVXL(JC,KC)
     +                          - BCL3XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRVXL(JC,KC)

              WRHS(ISTAL,JC,KC) = WRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRWXL(JC,KC)
     +                          - BCL4XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRWXL(JC,KC)

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STREXL(JC,KC)
     +                      - BCL3XL(JC,KC)*STRDXL(JC,KC)*STRVXL(JC,KC)
     +                      - BCL4XL(JC,KC)*STRDXL(JC,KC)*STRWXL(JC,KC)
     +                      - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                                     + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                                     + OVGMXL(JC,KC))

            ELSE
              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

              URHS(ISTAL,JC,KC) = URHS(ISTAL,JC,KC)
     +   - BCL5XL(JC,KC)*OVA2XL(JC,KC)*(STRUXL(JC,KC)+ACOUXL(JC,KC))

              VRHS(ISTAL,JC,KC) = VRHS(ISTAL,JC,KC)
     +                   - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRVXL(JC,KC)

              WRHS(ISTAL,JC,KC) = WRHS(ISTAL,JC,KC)
     +                   - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRWXL(JC,KC)

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +               - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                              + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                              + OVGMXL(JC,KC))
              ENDIF

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

              IF((STRUXL(JC,KC).GT.ZERO).AND.(FLAG_PIO_XL==1))THEN
                FORNOW = BCLYXL(JC,KC,ISPEC)*STRDXL(JC,KC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                            - FORNOW*STRHXL(JC,KC,ISPEC)

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     + - (BCL2XL(JC,KC)+BCL5XL(JC,KC)*OVA2XL(JC,KC))*STRYXL(JC,KC,ISPEC)
     +                                  - FORNOW

                ELSE
                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     +                 - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRYXL(JC,KC,ISPEC)
                ENDIF

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCXL.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXL(JC,KC) = SORPXL(JC,KC)
     +        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L2X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                      *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
C             L1X UNCHANGED
              BCL2XL(JC,KC) = -BCL2XL(JC,KC)
              BCL3XL(JC,KC) = -BCL3XL(JC,KC)
              BCL4XL(JC,KC) = -BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*SORPXL(JC,KC)
     +                     + COBCXL*ACOUXL(JC,KC)*(STRPXL(JC,KC)-PINFXL)
     +                      - BCL5XL(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF L's
                BCLYXL(JC,KC,ISPEC) = STRUXL(JC,KC)*BCLYXL(JC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
                BCLYXL(JC,KC,ISPEC) = RATEXL(JC,KC,ISPEC)/STRDXL(JC,KC)
     +                              - BCLYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

              URHS(ISTAL,JC,KC) = URHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRUXL(JC,KC)
     +      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*(STRUXL(JC,KC)+ACOUXL(JC,KC))

              VRHS(ISTAL,JC,KC) = VRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRVXL(JC,KC)
     +                          - BCL3XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRVXL(JC,KC)

              WRHS(ISTAL,JC,KC) = WRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRWXL(JC,KC)
     +                          - BCL4XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRWXL(JC,KC)

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STREXL(JC,KC)
     +                      - BCL3XL(JC,KC)*STRDXL(JC,KC)*STRVXL(JC,KC)
     +                      - BCL4XL(JC,KC)*STRDXL(JC,KC)*STRWXL(JC,KC)
     +                      - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                                     + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                                     + OVGMXL(JC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = BCLYXL(JC,KC,ISPEC)*STRDXL(JC,KC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                            - FORNOW*STRHXL(JC,KC,ISPEC)

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     + - (BCL2XL(JC,KC)+BCL5XL(JC,KC)*OVA2XL(JC,KC))*STRYXL(JC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXL.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXL(JC,KC) = ZERO
              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SYDTXL(JC,KC) = SYDTXL(JC,KC)
     +                        + DYDTXL(JC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPXL(JC,KC) = SORPXL(JC,KC)
     +                        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXL(JC,KC) = SYDTXL(JC,KC)/STRRXL(JC,KC)
              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X,L2X,L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                      *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X UNCHANGED
              BCL5XL(JC,KC) = BCL1XL(JC,KC)
     +                      - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
     +                      - BCL5XL(JC,KC)
              BCL2XL(JC,KC) = GAM1XL(JC,KC)*OVA2XL(JC,KC)
     +                       *(BCL1XL(JC,KC)+BCL5XL(JC,KC))
     +                      + STRDXL(JC,KC)*(DTDTXL(JC,KC)/STRTXL(JC,KC)
     +                      - SORPXL(JC,KC)/STRPXL(JC,KC)
     +                      + SYDTXL(JC,KC))
     +                      - BCL2XL(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXL.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                    *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X UNCHANGED
              FORNOW = BCL1XL(JC,KC)
     +               - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
              BCL2XL(JC,KC) = -DDDTXL(JC,KC)
     +                      - OVA2XL(JC,KC)*(BCL1XL(JC,KC)+FORNOW)
     +                      - BCL2XL(JC,KC)
              BCL3XL(JC,KC) = -DVDTXL(JC,KC) - BCL3XL(JC,KC)
              BCL4XL(JC,KC) = -DWDTXL(JC,KC) - BCL4XL(JC,KC)
              BCL5XL(JC,KC) = FORNOW - BCL5XL(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STREXL(JC,KC)
     +                  - BCL3XL(JC,KC)*STRDXL(JC,KC)*STRVXL(JC,KC)
     +                  - BCL4XL(JC,KC)*STRDXL(JC,KC)*STRWXL(JC,KC)
     +                  - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                                 + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                                 + OVGMXL(JC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                BCLYXL(JC,KC,ISPEC) = RATEXL(JC,KC,ISPEC)/STRDXL(JC,KC)
     +                              - DYDTXL(JC,KC,ISPEC)
     +                              - STRUXL(JC,KC)*BCLYXL(JC,KC,ISPEC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +          - BCLYXL(JC,KC,ISPEC)*STRDXL(JC,KC)*STRHXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  
 
C       =======================================================================

        IF(NSBCXL.EQ.NSBCI4)THEN

C         INFLOW BC No 4 - IMPLEMENTED BY FY
C         SUBSONIC NON-REFLECTING INFLOW AS PER LODATO WITH TRANSVERSE TERMS INCLUDED

C         ASSIGN NON-REFLECTING INFLOW COEFFICIENTS READ FROM CONTROL FILE
          NRIETA1 = RXLPRC(1)
          NRIETA2 = RXLPRC(2)
          NRIETA3 = RXLPRC(3)
          NRIETA4 = RXLPRC(4)
          NRIETA5 = RXLPRC(5)
          NRIETA6 = RXLPRC(6)

C         STORE VELOCITY, TEMPERATURE MASS FRACTION VALUES AT THE BOUNDARY
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STLUXL(JC,KC) = STRUXL(JC,KC)
              STLVXL(JC,KC) = STRVXL(JC,KC)
              STLWXL(JC,KC) = STRWXL(JC,KC)
              STLTXL(JC,KC) = STRTXL(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL
  
                STLYXL(JC,KC,ISPEC) = STRYXL(JC,KC,ISPEC)
  
              ENDDO
            ENDDO

          ENDDO
          
C         SET VELOCITY, TEMPERATURE, MASS FRACTION TARGET VALUES
          CALL BCUTXL
          CALL BCTTXL
          CALL BCYTXL

C         SPECIFY L's AS REQUIRED
C         L2X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
             
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL2XL(JC,KC) = STLUXL(JC,KC)
     +                      *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STLUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STLUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STLUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             STORE VALUE OF BCL2XL AFTER CALC              
              IF ((JC.EQ.1).AND.(KC.EQ.1)) TEMP3=BCL2XL(JC,KC)

C             SUBTRACT FROM NEW VALUE OF L's 
C             L1X UNCHANGED
              BCL2XL(JC,KC) = NRIETA2*STRDXL(JC,KC)
     +                       *STRRXL(JC,KC)/(XGDLEN*ACOUXL(JC,KC))
     +                       *(STLTXL(JC,KC)-STRTXL(JC,KC))
     +                       +TT2XL(JC,KC)*OVA2XL(JC,KC)-BCL2XL(JC,KC)
              BCL3XL(JC,KC) = NRIETA3*ACOUXL(JC,KC)
     +                       /XGDLEN*(STLVXL(JC,KC)-STRVXL(JC,KC))
     +                       +TT3XL(JC,KC)-BCL3XL(JC,KC)
              BCL4XL(JC,KC) = NRIETA4*ACOUXL(JC,KC)
     +                       /XGDLEN*(STLWXL(JC,KC)-STRWXL(JC,KC))
     +                       +TT4XL(JC,KC)-BCL4XL(JC,KC)
              BCL5XL(JC,KC) = NRIETA5
     +                       *STRDXL(JC,KC)*ACOUXL(JC,KC)
     +                       *ACOUXL(JC,KC)*(ONE-M2MAX)/(TWO*XGDLEN)
     +                       *(STLUXL(JC,KC)-STRUXL(JC,KC))
     +                       +HALF*TT5XL(JC,KC)
     +                       -BCL5XL(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF L's
                BCLYXL(JC,KC,ISPEC) = STLUXL(JC,KC)
     +                                *BCLYXL(JC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's
                BCLYXL(JC,KC,ISPEC) = NRIETA6
     +                               *ACOUXL(JC,KC)/XGDLEN
     +                               *(STLYXL(JC,KC,ISPEC)
     +                               -STRYXL(JC,KC,ISPEC))
     +                               +TT6XL(JC,KC,ISPEC)
     +                               -BCLYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
C         APPLYING WAVE AMPLITUDE VARIATIONS TO CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

              URHS(ISTAL,JC,KC) = URHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STLUXL(JC,KC)
     +      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*(STLUXL(JC,KC)+ACOUXL(JC,KC))

              VRHS(ISTAL,JC,KC) = VRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STLVXL(JC,KC)
     +                          - BCL3XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STLVXL(JC,KC)

              WRHS(ISTAL,JC,KC) = WRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STLWXL(JC,KC)
     +                          - BCL4XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STLWXL(JC,KC)

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STREXL(JC,KC)
     +                      - BCL3XL(JC,KC)*STRDXL(JC,KC)*STLVXL(JC,KC)
     +                      - BCL4XL(JC,KC)*STRDXL(JC,KC)*STLWXL(JC,KC)
     +                      - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                                     + STLUXL(JC,KC)/ACOUXL(JC,KC)
     +                                     + OVGMXL(JC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = BCLYXL(JC,KC,ISPEC)*STRDXL(JC,KC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                            - FORNOW*STRHXL(JC,KC,ISPEC)

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     + - (BCL2XL(JC,KC)+BCL5XL(JC,KC)*OVA2XL(JC,KC))*STLYXL(JC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

C        IF(NSBCXL.EQ.NSBCW1)THEN
C
CC         WALL BOUNDARY CONDITION No 1
CC         NO-SLIP WALL - ADIABATIC
C
CC         ALL VELOCITY COMPONENTS IMPOSED
CC         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
CC         SET IN SUBROUTINE BOUNDT
C
CC         SPECIFY L's AS REQUIRED
CC         L1X,L3X-L5X
C          DO KC = KSTAL,KSTOL
C            DO JC = JSTAL,JSTOL
C           
CC             OLD VALUE OF L's
C              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
C              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
C     +                      *(BCL5XL(JC,KC)-FORNOW)
C              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
C              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
C              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
C     +                      *(BCL5XL(JC,KC)+FORNOW)
C
CC             SUBTRACT FROM NEW VALUE OF L's
CC             L1X,L2X UNCHANGED
C              BCL3XL(JC,KC) = -DVDTXL(JC,KC) - BCL3XL(JC,KC)
C              BCL4XL(JC,KC) = -DWDTXL(JC,KC) - BCL4XL(JC,KC)
C              BCL5XL(JC,KC) = BCL1XL(JC,KC)
C     +                      - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
C     +                      - BCL5XL(JC,KC)
C
C            ENDDO
C          ENDDO
C
CC         ADD TO CONSERVATIVE SOURCE TERMS
C          DO KC = KSTAL,KSTOL
C            DO JC = JSTAL,JSTOL
C
C              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
C     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)
C
C              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
C     +                  - BCL3XL(JC,KC)*STRDXL(JC,KC)*STRVXL(JC,KC)
C     +                  - BCL4XL(JC,KC)*STRDXL(JC,KC)*STRWXL(JC,KC)
C     +                  - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
C     +                                 + STRUXL(JC,KC)/ACOUXL(JC,KC)
C     +                                 + OVGMXL(JC,KC))
C
C            ENDDO
C          ENDDO
C
C          DO ISPEC = 1,NSPEC
C
C            DO KC = KSTAL,KSTOL
C              DO JC = JSTAL,JSTOL
C
C                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
C     +                 - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRYXL(JC,KC,ISPEC)
C
C              ENDDO
C            ENDDO
C
C          ENDDO
C
C        ENDIF  

C       VM & NC CHANGING ADIABATIC WALL TREATMENT
        IF(NSBCXL.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

          DO KC=KSTAL,KSTOL
            DO JC=JSTAL,JSTOL
              STRTXL(JC,KC)=(48.0*TRUN(ISTAL+1,JC,KC)
     +                      -36.0*TRUN(ISTAL+2,JC,KC)
     +                      +16.0*TRUN(ISTAL+3,JC,KC)
     +                      -3.0*TRUN(5,JC,KC))/25.0
              TRUN(ISTAL,JC,KC)=STRTXL(JC,KC)
            ENDDO
          ENDDO

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXL(JC,KC) = SORPXL(JC,KC)
     +                        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                    *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X UNCHANGED
              BCL3XL(JC,KC) = ZERO
              BCL4XL(JC,KC) = ZERO
C              BCL3XL(JC,KC) = -DVDTXL(JC,KC) - BCL3XL(JC,KC)
C              BCL4XL(JC,KC) = -DWDTXL(JC,KC) - BCL4XL(JC,KC)
C             NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
              BCL5XL(JC,KC) = BCL1XL(JC,KC)
     +                      - BCL5XL(JC,KC)
C              BCL5XL(JC,KC) = BCL1XL(JC,KC)
C     +                      - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
C     +                      - BCL5XL(JC,KC)
              BCL2XL(JC,KC) = ZERO
C              BCL2XL(JC,KC) = GAM1XL(JC,KC)*OVA2XL(JC,KC)
C     +                       *(BCL1XL(JC,KC)+BCL5XL(JC,KC))
C     +                      + STRDXL(JC,KC)*(DTDTXL(JC,KC)/STRTXL(JC,KC)
C     +                      - SORPXL(JC,KC)/STRPXL(JC,KC))
C     +                      - BCL2XL(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF LYX
                BCLYXL(JC,KC,ISPEC) = ZERO
C                BCLYXL(JC,KC,ISPEC) = STRUXL(JC,KC)*BCLYXL(JC,KC,ISPEC)

C               UPDATE L2X
C                BCL2XL(JC,KC) = BCL2XL(JC,KC)
C     +                        + (RATEXL(JC,KC,ISPEC)
C     +                         - STRDXL(JC,KC)*BCLYXL(JC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRXL(JC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)
              URHS(ISTAL,JC,KC) = ZERO
              VRHS(ISTAL,JC,KC) = ZERO
              WRHS(ISTAL,JC,KC) = ZERO
              ERHS(ISTAL,JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     + - (BCL2XL(JC,KC)+BCL5XL(JC,KC)*OVA2XL(JC,KC))*STRYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXL.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXL(JC,KC) = SORPXL(JC,KC)
     +                        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                    *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X UNCHANGED
              BCL3XL(JC,KC) = ZERO
              BCL4XL(JC,KC) = ZERO
C              BCL3XL(JC,KC) = -DVDTXL(JC,KC) - BCL3XL(JC,KC)
C              BCL4XL(JC,KC) = -DWDTXL(JC,KC) - BCL4XL(JC,KC)
C             NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
              BCL5XL(JC,KC) = BCL1XL(JC,KC)
     +                      - BCL5XL(JC,KC)
C              BCL5XL(JC,KC) = BCL1XL(JC,KC)
C     +                      - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
C     +                      - BCL5XL(JC,KC)
              BCL2XL(JC,KC) = ZERO
C              BCL2XL(JC,KC) = GAM1XL(JC,KC)*OVA2XL(JC,KC)
C     +                       *(BCL1XL(JC,KC)+BCL5XL(JC,KC))
C     +                      + STRDXL(JC,KC)*(DTDTXL(JC,KC)/STRTXL(JC,KC)
C     +                      - SORPXL(JC,KC)/STRPXL(JC,KC))
C     +                      - BCL2XL(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF LYX
                BCLYXL(JC,KC,ISPEC) = ZERO
C                BCLYXL(JC,KC,ISPEC) = STRUXL(JC,KC)*BCLYXL(JC,KC,ISPEC)

C               UPDATE L2X
C                BCL2XL(JC,KC) = BCL2XL(JC,KC)
C     +                        + (RATEXL(JC,KC,ISPEC)
C     +                         - STRDXL(JC,KC)*BCLYXL(JC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRXL(JC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)
              URHS(ISTAL,JC,KC) = ZERO
              VRHS(ISTAL,JC,KC) = ZERO
              WRHS(ISTAL,JC,KC) = ZERO
              ERHS(ISTAL,JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     + - (BCL2XL(JC,KC)+BCL5XL(JC,KC)*OVA2XL(JC,KC))*STRYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF  
C     X-DIRECTION LEFT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     X-DIRECTION RIGHT-HAND END
C     --------------------------
      IF(FXRCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUXR = PRIMITIVE U-VELOCITY COMPONENT
C       STRVXR = PRIMITIVE V-VELOCITY COMPONENT
C       STRWXR = PRIMITIVE W-VELOCITY COMPONENT
C       STRPXR = PRESSURE
C       STRDXR = DENSITY
C       STRTXR = TEMPERATURE
C       STREXR = INTERNAL ENERGY
C       STRGXR = MIXTURE CP
C       STRRXR = MIXTURE SPECIFIC GAS CONSTANT
C       STRYXR(ISPEC) = SPECIES MASS FRACTION
C       RATEXR(ISPEC) = SPECIES REACTION RATE
C       STRHXR(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1XR = DUDX
C       BCL2XR = DRHODX
C       BCL3XR = DVDX
C       BCL4XR = DWDX
C       BCL5XR = DPDX
C       BCLYXR(ISPEC) = DYDX

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRHXR(JC,KC,ISPEC) = STRHXR(JC,KC,ISPEC)
     +         - STRGXR(JC,KC)*STRTXR(JC,KC)*RGSPEC(ISPEC)/STRRXR(JC,KC)

            ENDDO

          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            GAM1XR(JC,KC) = STRGXR(JC,KC) - STRRXR(JC,KC)
            STREXR(JC,KC) = STREXR(JC,KC) - GAM1XR(JC,KC)*STRTXR(JC,KC)

            GAM1XR(JC,KC) = STRRXR(JC,KC)/GAM1XR(JC,KC)
            OVGMXR(JC,KC) = ONE/GAM1XR(JC,KC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            FORNOW = STRGXR(JC,KC)*GAM1XR(JC,KC)*STRTXR(JC,KC)
            ACOUXR(JC,KC) = SQRT(FORNOW)
            OVA2XR(JC,KC) = ONE/FORNOW

C     THESE ARE INTRODUCED FOR LODATO'S BC- NC

            TT1XR(JC,KC)=T51BXR(JC,KC)+
     +                   T52BXR(JC,KC)*(GAM1XR(JC,KC)+1.0)-
     +        STRDXR(JC,KC)* ACOUXR(JC,KC)*T2BXR(JC,KC)         

            TT2XR(JC,KC)=ACOUXR(JC,KC)*ACOUXR(JC,KC)*
     +      T1BXR(JC,KC)-T51BXR(JC,KC)-(GAM1XR(JC,KC)+1.0)*
     +      T52BXR(JC,KC)

            TT3XR(JC,KC)=T3BXR(JC,KC)
            TT4XR(JC,KC)=T4BXR(JC,KC)
            
            TT5XR(JC,KC)=T51BXR(JC,KC)+
     +                   T52BXR(JC,KC)*(GAM1XR(JC,KC)+1.0)+
     +        STRDXR(JC,KC)* ACOUXR(JC,KC)*T2BXR(JC,KC)   

            DO IS=1,NSPEC
            TT6XR(JC,KC,IS)=T6BXR(JC,KC,IS)
            END DO


          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCXR.EQ.NSBCO1)THEN
          FLAG_PIO_XR=NXRPRM(1)!0
          FLAG_BET_XR=NXRPRM(2)!1

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXR(JC,KC) = SORPXR(JC,KC)
     +        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L1X AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              BCL2XR(JC,KC) = STRUXR(JC,KC)*(BCL2XR(JC,KC)-OVA2XR(JC,KC)
     +                       *BCL5XR(JC,KC))

              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)

C             OLD VALUE OF L1X
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +        *(BCL5XR(JC,KC)-STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC))

C             SUBTRACT FROM NEW VALUE OF L1X
C     TERMS INTRODUCED FOR LODATO'S BC- NC
              IF(FLAG_BET_XR==1) THEN
                BET=STRUXR(JC,KC)*STRUXR(JC,KC)+
     +              STRVXR(JC,KC)*STRVXR(JC,KC)+
     +              STRWXR(JC,KC)*STRWXR(JC,KC)
                BET=SQRT(BET)/ACOUXR(JC,KC)         
              END IF
              IF((STRUXR(JC,KC).LT.ZERO).AND.(FLAG_PIO_XR==1))THEN
              BCL2XR(JC,KC)=-OVA2XR(JC,KC)*SORPXR(JC,KC)
     +               -BCL2XR(JC,KC)

              BCL3XR(JC,KC)=0.1*(STRVXR(JC,KC)-0.0D0)-BCL3XR(JC,KC)

              BCL4XR(JC,KC)=0.1*(STRWXR(JC,KC)-0.0D0)-BCL4XR(JC,KC)
              BCL1XR(JC,KC)= HALF*SORPXR(JC,KC)
     +                     + COBCXR*ACOUXR(JC,KC)*(STRPXR(JC,KC)-PINFXR)
     +                     +0.5*(1.0-BET)*TT1XR(JC,KC)- BCL1XR(JC,KC)
     +                     -100.0*STRDXR(JC,KC)*(ONE/XGDLEN)
     +                 *(ACOUXR(JC,KC)**2.0-STRUXR(JC,KC)**2.0)
     +                 *(STRUXR(JC,KC)-0.0D0)
              ELSE
              BCL1XR(JC,KC)= HALF*SORPXR(JC,KC)
     +                     + COBCXR*ACOUXR(JC,KC)*(STRPXR(JC,KC)-PINFXR)
     +                     +0.5*(1.0-BET)*TT1XR(JC,KC)- BCL1XR(JC,KC)
              ENDIF
            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              IF((STRUXR(JC,KC).LT.ZERO).AND.(FLAG_PIO_XR==1))THEN
              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)
     +                          - BCL2XR(JC,KC)

              URHS(ISTOL,JC,KC) = URHS(ISTOL,JC,KC)
     +      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                          - BCL2XR(JC,KC)*STRUXR(JC,KC)

              VRHS(ISTOL,JC,KC) = VRHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRVXR(JC,KC)
     +                          - BCL2XR(JC,KC)*STRVXR(JC,KC)
     +                          - BCL3XR(JC,KC)*STRDXR(JC,KC)

              WRHS(ISTOL,JC,KC) = WRHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRWXR(JC,KC)
     +                          - BCL2XR(JC,KC)*STRWXR(JC,KC)
     +                          - BCL4XR(JC,KC)*STRDXR(JC,KC)

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                                     - STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                                     + OVGMXR(JC,KC))
     +                      - BCL2XR(JC,KC)*STREXR(JC,KC)
     +                      - BCL3XR(JC,KC)*STRDXR(JC,KC)*STRVXR(JC,KC)
     +                      - BCL4XR(JC,KC)*STRDXR(JC,KC)*STRWXR(JC,KC)
              

              ELSE

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)

              URHS(ISTOL,JC,KC) = URHS(ISTOL,JC,KC)
     +   - BCL1XR(JC,KC)*OVA2XR(JC,KC)*(STRUXR(JC,KC)-ACOUXR(JC,KC))

              VRHS(ISTOL,JC,KC) = VRHS(ISTOL,JC,KC)
     +                   - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRVXR(JC,KC)

              WRHS(ISTOL,JC,KC) = WRHS(ISTOL,JC,KC)
     +                   - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRWXR(JC,KC)

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +               - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                              - STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                              + OVGMXR(JC,KC))
              ENDIF

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

              IF((STRUXR(JC,KC).LT.ZERO).AND.(FLAG_PIO_XR==1))THEN

                FORNOW = BCLYXR(JC,KC,ISPEC)*STRDXR(JC,KC)

                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                            - FORNOW*STRHXR(JC,KC,ISPEC)

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     + - (BCL2XR(JC,KC)+BCL1XR(JC,KC)*OVA2XR(JC,KC))*STRYXR(JC,KC,ISPEC)
     +                                  - FORNOW

              ELSE
                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     +                 - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRYXR(JC,KC,ISPEC)
              ENDIF

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCXR.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXR(JC,KC) = SORPXR(JC,KC)
     +        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L4X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                      *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
C             L5X UNCHANGED
              BCL1XR(JC,KC) = HALF*SORPXR(JC,KC)
     +                     + COBCXR*ACOUXR(JC,KC)*(STRPXR(JC,KC)-PINFXR)
     +                      - BCL1XR(JC,KC)
              BCL2XR(JC,KC) = -BCL2XR(JC,KC)
              BCL3XR(JC,KC) = -BCL3XR(JC,KC)
              BCL4XR(JC,KC) = -BCL4XR(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF L's
                BCLYXR(JC,KC,ISPEC) = STRUXR(JC,KC)*BCLYXR(JC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
                BCLYXR(JC,KC,ISPEC) = RATEXR(JC,KC,ISPEC)/STRDXR(JC,KC)
     +                              - BCLYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)
     +                          - BCL2XR(JC,KC)

              URHS(ISTOL,JC,KC) = URHS(ISTOL,JC,KC)
     +      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                          - BCL2XR(JC,KC)*STRUXR(JC,KC)

              VRHS(ISTOL,JC,KC) = VRHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRVXR(JC,KC)
     +                          - BCL2XR(JC,KC)*STRVXR(JC,KC)
     +                          - BCL3XR(JC,KC)*STRDXR(JC,KC)

              WRHS(ISTOL,JC,KC) = WRHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRWXR(JC,KC)
     +                          - BCL2XR(JC,KC)*STRWXR(JC,KC)
     +                          - BCL4XR(JC,KC)*STRDXR(JC,KC)

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                                     - STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                                     + OVGMXR(JC,KC))
     +                      - BCL2XR(JC,KC)*STREXR(JC,KC)
     +                      - BCL3XR(JC,KC)*STRDXR(JC,KC)*STRVXR(JC,KC)
     +                      - BCL4XR(JC,KC)*STRDXR(JC,KC)*STRWXR(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = BCLYXR(JC,KC,ISPEC)*STRDXR(JC,KC)

                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                            - FORNOW*STRHXR(JC,KC,ISPEC)

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     + - (BCL2XR(JC,KC)+BCL1XR(JC,KC)*OVA2XR(JC,KC))*STRYXR(JC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXR.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXR(JC,KC) = ZERO
              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SYDTXR(JC,KC) = SYDTXR(JC,KC)
     +                        + DYDTXR(JC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPXR(JC,KC) = SORPXR(JC,KC)
     +                        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXR(JC,KC) = SYDTXR(JC,KC)/STRRXR(JC,KC)
              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X,L2X,L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                      *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5X UNCHANGED
              BCL1XR(JC,KC) = BCL5XR(JC,KC)
     +                      + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
     +                      - BCL1XR(JC,KC)
              BCL2XR(JC,KC) = GAM1XR(JC,KC)*OVA2XR(JC,KC)
     +                       *(BCL1XR(JC,KC)+BCL5XR(JC,KC))
     +                      + STRDXR(JC,KC)*(DTDTXR(JC,KC)/STRTXR(JC,KC)
     +                      - SORPXR(JC,KC)/STRPXR(JC,KC)
     +                      + SYDTXR(JC,KC))
     +                      - BCL2XR(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL2XR(JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXR.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                    *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5X UNCHANGED
              FORNOW = BCL5XR(JC,KC)
     +               + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
              BCL1XR(JC,KC) = FORNOW - BCL1XR(JC,KC)
              BCL2XR(JC,KC) = -DDDTXR(JC,KC)
     +                      - OVA2XR(JC,KC)*(BCL5XR(JC,KC)+FORNOW)
     +                      - BCL2XR(JC,KC)
              BCL3XR(JC,KC) = -DVDTXR(JC,KC) - BCL3XR(JC,KC)
              BCL4XR(JC,KC) = -DWDTXR(JC,KC) - BCL4XR(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                  - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                                 - STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                                 + OVGMXR(JC,KC))
     +                  - BCL2XR(JC,KC)*STREXR(JC,KC)
     +                  - BCL3XR(JC,KC)*STRDXR(JC,KC)*STRVXR(JC,KC)
     +                  - BCL4XR(JC,KC)*STRDXR(JC,KC)*STRWXR(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                BCLYXR(JC,KC,ISPEC) = RATEXR(JC,KC,ISPEC)/STRDXR(JC,KC)
     +                              - DYDTXR(JC,KC,ISPEC)
     +                              - STRUXR(JC,KC)*BCLYXR(JC,KC,ISPEC)

                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +          - BCLYXR(JC,KC,ISPEC)*STRDXR(JC,KC)*STRHXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

C        IF(NSBCXR.EQ.NSBCW1)THEN
C
CC         WALL BOUNDARY CONDITION No 1
CC         NO-SLIP WALL - ADIABATIC
C
CC         ALL VELOCITY COMPONENTS IMPOSED
CC         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
CC         SET IN SUBROUTINE BOUNDT
C
CC         SPECIFY L's AS REQUIRED
CC         L1X,L3X-L5X
C          DO KC = KSTAL,KSTOL
C            DO JC = JSTAL,JSTOL
C           
CC             OLD VALUE OF L's
C              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
C              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
C     +                      *(BCL5XR(JC,KC)-FORNOW)
C              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
C              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)
C              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
C     +                      *(BCL5XR(JC,KC)+FORNOW)
C
CC             SUBTRACT FROM NEW VALUE OF L's
CC             L2X,L5X UNCHANGED
C              BCL1XR(JC,KC) = BCL5XR(JC,KC)
C     +                      + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
C     +                      - BCL1XR(JC,KC)
C              BCL3XR(JC,KC) = -DVDTXR(JC,KC) - BCL3XR(JC,KC)
C              BCL4XR(JC,KC) = -DWDTXR(JC,KC) - BCL4XR(JC,KC)
C
C            ENDDO
C          ENDDO
C
CC         ADD TO CONSERVATIVE SOURCE TERMS
C          DO KC = KSTAL,KSTOL
C            DO JC = JSTAL,JSTOL
C
C              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
C     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)
C
C              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
C     +                  - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
C     +                                 + STRUXR(JC,KC)/ACOUXR(JC,KC)
C     +                                 + OVGMXR(JC,KC))
C     +                  - BCL3XR(JC,KC)*STRDXR(JC,KC)*STRVXR(JC,KC)
C     +                  - BCL4XR(JC,KC)*STRDXR(JC,KC)*STRWXR(JC,KC)
C
C            ENDDO
C          ENDDO
C
C          DO ISPEC = 1,NSPEC
C
C            DO KC = KSTAL,KSTOL
C              DO JC = JSTAL,JSTOL
C
C                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
C     +                 - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRYXR(JC,KC,ISPEC)
C
C              ENDDO
C            ENDDO
C
C          ENDDO
C
C        ENDIF  
C       VM & NC CHANGING ADIABATIC WALL TREATMENT
        IF(NSBCXR.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ADIABATIC

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

          DO KC=KSTAL,KSTOL
            DO JC=JSTAL,JSTOL
              STRTXR(JC,KC)=(48.0*TRUN(ISTOL-1,JC,KC)
     +                      -36.0*TRUN(ISTOL-2,JC,KC)
     +                      +16.0*TRUN(ISTOL-3,JC,KC)
     +                      -3.0*TRUN(ISTOL-4,JC,KC))/25.0
              TRUN(ISTOL,JC,KC)=STRTXR(JC,KC)
            ENDDO
          ENDDO

C         PRECOMPUTE CHEMISTRY TERMS

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXR(JC,KC) = SORPXR(JC,KC)
     +                        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                    *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5X UNCHANGED
C             NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
              BCL1XR(JC,KC) = BCL5XR(JC,KC)
     +                      - BCL1XR(JC,KC)
C              BCL1XR(JC,KC) = BCL5XR(JC,KC)
C     +                      + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
C     +                      - BCL1XR(JC,KC)
              BCL3XR(JC,KC) = ZERO
              BCL4XR(JC,KC) = ZERO
              BCL2XR(JC,KC) = ZERO
C              BCL3XR(JC,KC) = -DVDTXR(JC,KC) - BCL3XR(JC,KC)
C              BCL4XR(JC,KC) = -DWDTXR(JC,KC) - BCL4XR(JC,KC)
C              BCL2XR(JC,KC) = GAM1XR(JC,KC)*OVA2XR(JC,KC)
C     +                       *(BCL1XR(JC,KC)+BCL5XR(JC,KC))
C     +                      + STRDXR(JC,KC)*(DTDTXR(JC,KC)/STRTXR(JC,KC)
C     +                      - SORPXR(JC,KC)/STRPXR(JC,KC))
C     +                      - BCL2XR(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF LYX
                BCLYXR(JC,KC,ISPEC) = ZERO
C                BCLYXR(JC,KC,ISPEC) = STRUXR(JC,KC)*BCLYXR(JC,KC,ISPEC)

C               UPDATE L2X
C                BCL2XR(JC,KC) = BCL2XR(JC,KC)
C     +                        + (RATEXR(JC,KC,ISPEC)
C     +                         - STRDXR(JC,KC)*BCLYXR(JC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRXR(JC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)
     +                          - BCL2XR(JC,KC)
              URHS(ISTOL,JC,KC) = ZERO
              VRHS(ISTOL,JC,KC) = ZERO
              WRHS(ISTOL,JC,KC) = ZERO
              ERHS(ISTOL,JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     + - (BCL2XR(JC,KC)+BCL1XR(JC,KC)*OVA2XR(JC,KC))*STRYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXR.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXR(JC,KC) = SORPXR(JC,KC)
     +                        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                    *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5X UNCHANGED
C             NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
              BCL1XR(JC,KC) = BCL5XR(JC,KC)
     +                      - BCL1XR(JC,KC)
C              BCL1XR(JC,KC) = BCL5XR(JC,KC)
C     +                      + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
C     +                      - BCL1XR(JC,KC)
              BCL3XR(JC,KC) = ZERO
              BCL4XR(JC,KC) = ZERO
              BCL2XR(JC,KC) = ZERO
C              BCL3XR(JC,KC) = -DVDTXR(JC,KC) - BCL3XR(JC,KC)
C              BCL4XR(JC,KC) = -DWDTXR(JC,KC) - BCL4XR(JC,KC)
C              BCL2XR(JC,KC) = GAM1XR(JC,KC)*OVA2XR(JC,KC)
C     +                       *(BCL1XR(JC,KC)+BCL5XR(JC,KC))
C     +                      + STRDXR(JC,KC)*(DTDTXR(JC,KC)/STRTXR(JC,KC)
C     +                      - SORPXR(JC,KC)/STRPXR(JC,KC))
C     +                      - BCL2XR(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF LYX
                BCLYXR(JC,KC,ISPEC) = ZERO
C                BCLYXR(JC,KC,ISPEC) = STRUXR(JC,KC)*BCLYXR(JC,KC,ISPEC)

C               UPDATE L2X
C                BCL2XR(JC,KC) = BCL2XR(JC,KC)
C     +                        + (RATEXR(JC,KC,ISPEC)
C     +                         - STRDXR(JC,KC)*BCLYXR(JC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRXR(JC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)
     +                          - BCL2XR(JC,KC)
              URHS(ISTOL,JC,KC) = ZERO
              VRHS(ISTOL,JC,KC) = ZERO
              WRHS(ISTOL,JC,KC) = ZERO
              ERHS(ISTOL,JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     + - (BCL2XR(JC,KC)+BCL1XR(JC,KC)*OVA2XR(JC,KC))*STRYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     X-DIRECTION RIGHT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Y-DIRECTION LEFT-HAND END
C     -------------------------
      IF(FYLCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUYL = PRIMITIVE U-VELOCITY COMPONENT
C       STRVYL = PRIMITIVE V-VELOCITY COMPONENT
C       STRWYL = PRIMITIVE W-VELOCITY COMPONENT
C       STRPYL = PRESSURE
C       STRDYL = DENSITY
C       STRTYL = TEMPERATURE
C       STREYL = INTERNAL ENERGY
C       STRGYL = MIXTURE CP
C       STRRYL = MIXTURE SPECIFIC GAS CONSTANT
C       STRYYL(ISPEC) = SPECIES MASS FRACTION
C       RATEYL(ISPEC) = SPECIES REACTION RATE
C       STRHYL(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1YL = DVDY
C       BCL2YL = DRHODY
C       BCL3YL = DUDY
C       BCL4YL = DWDY
C       BCL5YL = DPDY
C       BCLYYL(ISPEC) = DYDY

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRHYL(IC,KC,ISPEC) = STRHYL(IC,KC,ISPEC)
     +         - STRGYL(IC,KC)*STRTYL(IC,KC)*RGSPEC(ISPEC)/STRRYL(IC,KC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            GAM1YL(IC,KC) = STRGYL(IC,KC) - STRRYL(IC,KC)
            STREYL(IC,KC) = STREYL(IC,KC) - GAM1YL(IC,KC)*STRTYL(IC,KC)

            GAM1YL(IC,KC) = STRRYL(IC,KC)/GAM1YL(IC,KC)
            OVGMYL(IC,KC) = ONE/GAM1YL(IC,KC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGYL(IC,KC)*GAM1YL(IC,KC)*STRTYL(IC,KC)
            ACOUYL(IC,KC) = SQRT(FORNOW)
            OVA2YL(IC,KC) = ONE/FORNOW

C     THESE ARE INTRODUCED FOR LODATO'S BC- NC

            TT1YL(IC,KC)=T51BYL(IC,KC)+
     +                   T52BYL(IC,KC)*(GAM1YL(IC,KC)+1.0)-
     +        STRDYL(IC,KC)* ACOUYL(IC,KC)*T2BYL(IC,KC)         

            TT2YL(IC,KC)=ACOUYL(IC,KC)*ACOUYL(IC,KC)*
     +      T1BYL(IC,KC)-T51BYL(IC,KC)-(GAM1YL(IC,KC)+1.0)*
     +      T52BYL(IC,KC)

            TT3YL(IC,KC)=T3BYL(IC,KC)
            TT4YL(IC,KC)=T4BYL(IC,KC)
            
            TT5YL(IC,KC)=T51BYL(IC,KC)+
     +                   T52BYL(IC,KC)*(GAM1YL(IC,KC)+1.0)+
     +        STRDYL(IC,KC)* ACOUYL(IC,KC)*T2BYL(IC,KC)   

            DO IS=1,NSPEC
            TT6YL(IC,KC,IS)=T6BYL(IC,KC,IS)
            END DO


          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCYL.EQ.NSBCO1)THEN
          FLAG_PIO_YL=NYLPRM(1)!0
          FLAG_BET_YL=NYLPRM(2)!1

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1, NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYL(IC,KC) = SORPYL(IC,KC)
     +        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L5Y AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                      *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
C             OLD VALUE OF L5Y
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +        *(BCL5YL(IC,KC)+STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC))

C             SUBTRACT FROM NEW VALUE OF L5Y
C             TERMS INTRODUCED FOR LODATO'S BC- NC
              IF(FLAG_BET_YL==1) THEN
                BET=STRUYL(IC,KC)*STRUYL(IC,KC)+
     +              STRVYL(IC,KC)*STRVYL(IC,KC)+
     +              STRWYL(IC,KC)*STRWYL(IC,KC)
                BET=SQRT(BET)/ACOUYL(IC,KC)
              END IF
              IF((STRVYL(IC,KC).GT.ZERO).AND.(FLAG_PIO_YL==1))THEN
              BCL2YL(IC,KC) =
     +                       -BCL2YL(IC,KC)
     +                       -OVA2YL(IC,KC)*SORPYL(IC,KC)

              BCL3YL(IC,KC) =
     +                        0.1*(STRUYL(IC,KC)-0.0D0)
     +                       -BCL3YL(IC,KC)

              BCL4YL(IC,KC) =
     +                        0.1*(STRWYL(IC,KC)-0.0D0)
     +                       -BCL4YL(IC,KC)
              BCL5YL(IC,KC)= HALF*SORPYL(IC,KC)
     +                     + COBCYL*ACOUYL(IC,KC)*(STRPYL(IC,KC)-PINFYL)
     + +0.5*(1.0-BET)*TT5YL(IC,KC)- BCL5YL(IC,KC)
     +                     + 100.0*STRDYL(IC,KC)*(ONE/YGDLEN)!ASK NC: CHANGED TO YGDLEN
     +               *(ACOUYL(IC,KC)**2.0-STRVYL(IC,KC)**2.0)
     +               *(STRVYL(IC,KC)-0.0)
                
              ELSE
              BCL5YL(IC,KC)= HALF*SORPYL(IC,KC)
     +                     + COBCYL*ACOUYL(IC,KC)*(STRPYL(IC,KC)-PINFYL)
     + +0.5*(1.0-BET)*TT5YL(IC,KC)- BCL5YL(IC,KC)
              ENDIF
            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
              IF((STRVYL(IC,KC).GT.ZERO).AND.(FLAG_PIO_YL==1))THEN

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

              URHS(IC,JSTAL,KC) = URHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRUYL(IC,KC)
     +                          - BCL3YL(IC,KC)*STRDYL(IC,KC)
     +                      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRUYL(IC,KC)

              VRHS(IC,JSTAL,KC) = VRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRVYL(IC,KC)
     +      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*(STRVYL(IC,KC)+ACOUYL(IC,KC))

              WRHS(IC,JSTAL,KC) = WRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRWYL(IC,KC)
     +                          - BCL4YL(IC,KC)*STRDYL(IC,KC)
     +                      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRWYL(IC,KC)

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STREYL(IC,KC)
     +                      - BCL3YL(IC,KC)*STRDYL(IC,KC)*STRUYL(IC,KC)
     +                      - BCL4YL(IC,KC)*STRDYL(IC,KC)*STRWYL(IC,KC)
     +                      - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                                     + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                                     + OVGMYL(IC,KC))

              ELSE
              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

              URHS(IC,JSTAL,KC) = URHS(IC,JSTAL,KC)
     +                   - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRUYL(IC,KC)

              VRHS(IC,JSTAL,KC) = VRHS(IC,JSTAL,KC)
     +   - BCL5YL(IC,KC)*OVA2YL(IC,KC)*(STRVYL(IC,KC)+ACOUYL(IC,KC))

              WRHS(IC,JSTAL,KC) = WRHS(IC,JSTAL,KC)
     +                   - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRWYL(IC,KC)

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +               - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                              + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                              + OVGMYL(IC,KC))
              ENDIF
            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL
              IF((STRVYL(IC,KC).GT.ZERO).AND.(FLAG_PIO_YL==1))THEN

                FORNOW = BCLYYL(IC,KC,ISPEC)*STRDYL(IC,KC)

                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                            - FORNOW*STRHYL(IC,KC,ISPEC)

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     + - (BCL2YL(IC,KC)+BCL5YL(IC,KC)*OVA2YL(IC,KC))*STRYYL(IC,KC,ISPEC)
     +                                  - FORNOW


              ELSE
                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     +                 - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRYYL(IC,KC,ISPEC)
              ENDIF
              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCYL.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYL(IC,KC) = SORPYL(IC,KC)
     +        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L2Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                      *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
C             L1Y UNCHANGED
              BCL2YL(IC,KC) = -BCL2YL(IC,KC)
              BCL3YL(IC,KC) = -BCL3YL(IC,KC)
              BCL4YL(IC,KC) = -BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*SORPYL(IC,KC)
     +                     + COBCYL*ACOUYL(IC,KC)*(STRPYL(IC,KC)-PINFYL)
     +                      - BCL5YL(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYYL(IC,KC,ISPEC) = STRVYL(IC,KC)*BCLYYL(IC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
                BCLYYL(IC,KC,ISPEC) = RATEYL(IC,KC,ISPEC)/STRDYL(IC,KC)
     +                              - BCLYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

              URHS(IC,JSTAL,KC) = URHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRUYL(IC,KC)
     +                          - BCL3YL(IC,KC)*STRDYL(IC,KC)
     +                      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRUYL(IC,KC)

              VRHS(IC,JSTAL,KC) = VRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRVYL(IC,KC)
     +      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*(STRVYL(IC,KC)+ACOUYL(IC,KC))

              WRHS(IC,JSTAL,KC) = WRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRWYL(IC,KC)
     +                          - BCL4YL(IC,KC)*STRDYL(IC,KC)
     +                      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRWYL(IC,KC)

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STREYL(IC,KC)
     +                      - BCL3YL(IC,KC)*STRDYL(IC,KC)*STRUYL(IC,KC)
     +                      - BCL4YL(IC,KC)*STRDYL(IC,KC)*STRWYL(IC,KC)
     +                      - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                                     + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                                     + OVGMYL(IC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYYL(IC,KC,ISPEC)*STRDYL(IC,KC)

                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                            - FORNOW*STRHYL(IC,KC,ISPEC)

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     + - (BCL2YL(IC,KC)+BCL5YL(IC,KC)*OVA2YL(IC,KC))*STRYYL(IC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYL.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYL(IC,KC) = ZERO
              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SYDTYL(IC,KC) = SYDTYL(IC,KC)
     +                        + DYDTYL(IC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPYL(IC,KC) = SORPYL(IC,KC)
     +                        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYL(IC,KC) = SYDTYL(IC,KC)/STRRYL(IC,KC)
              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y,L2Y,L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                      *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL5YL(IC,KC) = BCL1YL(IC,KC)
     +                      - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
     +                      - BCL5YL(IC,KC)
              BCL2YL(IC,KC) = GAM1YL(IC,KC)*OVA2YL(IC,KC)
     +                       *(BCL1YL(IC,KC)+BCL5YL(IC,KC))
     +                      + STRDYL(IC,KC)*(DTDTYL(IC,KC)/STRTYL(IC,KC)
     +                      - SORPYL(IC,KC)/STRPYL(IC,KC)
     +                      + SYDTYL(IC,KC))
     +                      - BCL2YL(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCYL.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                    *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              FORNOW = BCL1YL(IC,KC)
     +               - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
              BCL2YL(IC,KC) = -DDDTYL(IC,KC)
     +                      - OVA2YL(IC,KC)*(BCL1YL(IC,KC)+FORNOW)
     +                      - BCL2YL(IC,KC)
              BCL3YL(IC,KC) = -DUDTYL(IC,KC) - BCL3YL(IC,KC)
              BCL4YL(IC,KC) = -DWDTYL(IC,KC) - BCL4YL(IC,KC)
              BCL5YL(IC,KC) = FORNOW - BCL5YL(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STREYL(IC,KC)
     +                  - BCL3YL(IC,KC)*STRDYL(IC,KC)*STRUYL(IC,KC)
     +                  - BCL4YL(IC,KC)*STRDYL(IC,KC)*STRWYL(IC,KC)
     +                  - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                                 + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                                 + OVGMYL(IC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                BCLYYL(IC,KC,ISPEC) = RATEYL(IC,KC,ISPEC)/STRDYL(IC,KC)
     +                              - DYDTYL(IC,KC,ISPEC)
     +                              - STRVYL(IC,KC)*BCLYYL(IC,KC,ISPEC)

                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +          - BCLYYL(IC,KC,ISPEC)*STRDYL(IC,KC)*STRHYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

C        IF(NSBCYL.EQ.NSBCW1)THEN
C
CC         WALL BOUNDARY CONDITION No 1
CC         NO-SLIP WALL - ADIABATIC
C
CC         ALL VELOCITY COMPONENTS IMPOSED
CC         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
CC         SET IN SUBROUTINE BOUNDT
C
CC         SPECIFY L's AS REQUIRED
CC         L1Y,L3Y-L5Y
C          DO KC = KSTAL,KSTOL
C            DO IC = ISTAL,ISTOL
C           
CC             OLD VALUE OF L's
C              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
C              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
C     +                      *(BCL5YL(IC,KC)-FORNOW)
C              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
C              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
C              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
C     +                      *(BCL5YL(IC,KC)+FORNOW)
C
CC             SUBTRACT FROM NEW VALUE OF L's
CC             L1Y,L2Y UNCHANGED
C              BCL3YL(IC,KC) = -DUDTYL(IC,KC) - BCL3YL(IC,KC)
C              BCL4YL(IC,KC) = -DWDTYL(IC,KC) - BCL4YL(IC,KC)
C              BCL5YL(IC,KC) = BCL1YL(IC,KC)
C     +                      - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
C     +                      - BCL5YL(IC,KC)
C
C            ENDDO
C          ENDDO
C
CC         ADD TO CONSERVATIVE SOURCE TERMS
C          DO KC = KSTAL,KSTOL
C            DO IC = ISTAL,ISTOL
C
C              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
C     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)
C
C              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
C     +                  - BCL3YL(IC,KC)*STRDYL(IC,KC)*STRUYL(IC,KC)
C     +                  - BCL4YL(IC,KC)*STRDYL(IC,KC)*STRWYL(IC,KC)
C     +                  - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
C     +                                 + STRVYL(IC,KC)/ACOUYL(IC,KC)
C     +                                 + OVGMYL(IC,KC))
C
C            ENDDO
C          ENDDO
C
C          DO ISPEC = 1,NSPEC
C
C            DO KC = KSTAL,KSTOL
C              DO IC = ISTAL,ISTOL
C
C                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
C     +                 - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRYYL(IC,KC,ISPEC)
C
C              ENDDO
C            ENDDO
C
C          ENDDO
C
C        ENDIF  

C       VM & NC CHANGING ADIABATIC WALL TREATMENT
        IF(NSBCYL.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ADIABATIC

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

          DO KC=KSTAL,KSTOL
            DO IC=ISTAL,ISTOL
              STRTYL(IC,KC)=(48.0*TRUN(IC,JSTAL+1,KC)
     +                      -36.0*TRUN(IC,JSTAL+2,KC)
     +                      +16.0*TRUN(IC,JSTAL+3,KC)
     +                      -3.0*TRUN(IC,JSTAL+4,KC))/25.0
              TRUN(IC,JSTAL,KC)=STRTYL(IC,KC)
            ENDDO
          ENDDO

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYL(IC,KC) = SORPYL(IC,KC)
     +                        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                    *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL3YL(IC,KC) = ZERO
              BCL4YL(IC,KC) = ZERO
C              BCL3YL(IC,KC) = -DUDTYL(IC,KC) - BCL3YL(IC,KC)
C              BCL4YL(IC,KC) = -DWDTYL(IC,KC) - BCL4YL(IC,KC)
              BCL5YL(IC,KC) = BCL1YL(IC,KC)
     +                      - BCL5YL(IC,KC)
C              BCL5YL(IC,KC) = BCL1YL(IC,KC)
C     +                      - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
C     +                      - BCL5YL(IC,KC)
              BCL2YL(IC,KC) = ZERO
C              BCL2YL(IC,KC) = GAM1YL(IC,KC)*OVA2YL(IC,KC)
C     +                       *(BCL1YL(IC,KC)+BCL5YL(IC,KC))
C     +                      + STRDYL(IC,KC)*(DTDTYL(IC,KC)/STRTYL(IC,KC)
C     +                      - SORPYL(IC,KC)/STRPYL(IC,KC))
C     +                      - BCL2YL(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYY
                BCLYYL(IC,KC,ISPEC) = ZERO
C                BCLYYL(IC,KC,ISPEC) = STRVYL(IC,KC)*BCLYYL(IC,KC,ISPEC)

C               UPDATE L2Y
C                BCL2YL(IC,KC) = BCL2YL(IC,KC)
C     +                        + (RATEYL(IC,KC,ISPEC)
C     +                         - STRDYL(IC,KC)*BCLYYL(IC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRYL(IC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)
              URHS(IC,JSTAL,KC) = ZERO
              VRHS(IC,JSTAL,KC) = ZERO
              WRHS(IC,JSTAL,KC) = ZERO
              ERHS(IC,JSTAL,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     + - (BCL2YL(IC,KC)+BCL5YL(IC,KC)*OVA2YL(IC,KC))*STRYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCYL.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYL(IC,KC) = SORPYL(IC,KC)
     +                        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                    *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL3YL(IC,KC) = ZERO
              BCL4YL(IC,KC) = ZERO
C              BCL3YL(IC,KC) = -DUDTYL(IC,KC) - BCL3YL(IC,KC)
C              BCL4YL(IC,KC) = -DWDTYL(IC,KC) - BCL4YL(IC,KC)
              BCL5YL(IC,KC) = BCL1YL(IC,KC)
     +                      - BCL5YL(IC,KC)
C              BCL5YL(IC,KC) = BCL1YL(IC,KC)
C     +                      - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
C     +                      - BCL5YL(IC,KC)
              BCL2YL(IC,KC) = ZERO
C              BCL2YL(IC,KC) = GAM1YL(IC,KC)*OVA2YL(IC,KC)
C     +                       *(BCL1YL(IC,KC)+BCL5YL(IC,KC))
C     +                      + STRDYL(IC,KC)*(DTDTYL(IC,KC)/STRTYL(IC,KC)
C     +                      - SORPYL(IC,KC)/STRPYL(IC,KC))
C     +                      - BCL2YL(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYY
                BCLYYL(IC,KC,ISPEC) = ZERO
C                BCLYYL(IC,KC,ISPEC) = STRVYL(IC,KC)*BCLYYL(IC,KC,ISPEC)

C               UPDATE L2Y
C                BCL2YL(IC,KC) = BCL2YL(IC,KC)
C     +                        + (RATEYL(IC,KC,ISPEC)
C     +                         - STRDYL(IC,KC)*BCLYYL(IC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRYL(IC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)
              URHS(IC,JSTAL,KC) = ZERO
              VRHS(IC,JSTAL,KC) = ZERO
              WRHS(IC,JSTAL,KC) = ZERO
              ERHS(IC,JSTAL,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     + - (BCL2YL(IC,KC)+BCL5YL(IC,KC)*OVA2YL(IC,KC))*STRYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Y-DIRECTION LEFT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Y-DIRECTION RIGHT-HAND END
C     --------------------------
      IF(FYRCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUYR = PRIMITIVE U-VELOCITY COMPONENT
C       STRVYR = PRIMITIVE V-VELOCITY COMPONENT
C       STRWYR = PRIMITIVE W-VELOCITY COMPONENT
C       STRPYR = PRESSURE
C       STRDYR = DENSITY
C       STRTYR = TEMPERATURE
C       STREYR = INTERNAL ENERGY
C       STRGYR = MIXTURE CP
C       STRRYR = MIXTURE SPECIFIC GAS CONSTANT
C       STRYYR(ISPEC) = SPECIES MASS FRACTION
C       RATEYR(ISPEC) = SPECIES REACTION RATE
C       STRHYR(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1YR = DVDY
C       BCL2YR = DRHODY
C       BCL3YR = DUDY
C       BCL4YR = DWDY
C       BCL5YR = DPDY
C       BCLYYR(ISPEC) = DYDY

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRHYR(IC,KC,ISPEC) = STRHYR(IC,KC,ISPEC)
     +         - STRGYR(IC,KC)*STRTYR(IC,KC)*RGSPEC(ISPEC)/STRRYR(IC,KC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            GAM1YR(IC,KC) = STRGYR(IC,KC) - STRRYR(IC,KC)
            STREYR(IC,KC) = STREYR(IC,KC) - GAM1YR(IC,KC)*STRTYR(IC,KC)

            GAM1YR(IC,KC) = STRRYR(IC,KC)/GAM1YR(IC,KC)
            OVGMYR(IC,KC) = ONE/GAM1YR(IC,KC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGYR(IC,KC)*GAM1YR(IC,KC)*STRTYR(IC,KC)
            ACOUYR(IC,KC) = SQRT(FORNOW)
            OVA2YR(IC,KC) = ONE/FORNOW

C     THESE ARE INTRODUCED FOR LODATO'S BC- NC

            TT1YR(IC,KC)=T51BYR(IC,KC)+
     +                   T52BYR(IC,KC)*(GAM1YR(IC,KC)+1.0)-
     +        STRDYR(IC,KC)* ACOUYR(IC,KC)*T2BYR(IC,KC)         

            TT2YR(IC,KC)=ACOUYR(IC,KC)*ACOUYR(IC,KC)*
     +      T1BYR(IC,KC)-T51BYR(IC,KC)-(GAM1YR(IC,KC)+1.0)*
     +      T52BYR(IC,KC)

            TT3YR(IC,KC)=T3BYR(IC,KC)
            TT4YR(IC,KC)=T4BYR(IC,KC)
            
            TT5YR(IC,KC)=T51BYR(IC,KC)+
     +                   T52BYR(IC,KC)*(GAM1YR(IC,KC)+1.0)+
     +        STRDYR(IC,KC)* ACOUYR(IC,KC)*T2BYR(IC,KC)   

            DO IS=1,NSPEC
            TT6YR(IC,KC,IS)=T6BYR(IC,KC,IS)
            END DO

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCYR.EQ.NSBCO1)THEN
          FLAG_PIO_YR=NYRPRM(1)!0
          FLAG_BET_YR=NYRPRM(2)!1

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYR(IC,KC) = SORPYR(IC,KC)
     +        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L1Y AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                      *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
C             OLD VALUE OF L1Y
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +        *(BCL5YR(IC,KC)-STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC))

C             SUBTRACT FROM NEW VALUE OF L1Y
C     TERMS INTRODUCED FOR LODATO'S BC- NC
              IF(FLAG_BET_YR==1) THEN
                BET=STRUYR(IC,KC)*STRUYR(IC,KC)+
     +              STRVYR(IC,KC)*STRVYR(IC,KC)+
     +              STRWYR(IC,KC)*STRWYR(IC,KC)
                BET=SQRT(BET)/ACOUYR(IC,KC)         
              END IF

              IF((STRVYR(IC,KC).LT.ZERO).AND.(FLAG_PIO_YR==1))THEN
              BCL2YR(IC,KC) =
     +                       -BCL2YR(IC,KC)
     +                       -OVA2YR(IC,KC)*SORPYR(IC,KC)

              BCL3YR(IC,KC) =
     +                        0.1*(STRUYR(IC,KC)-0.0D0)
     +                       -BCL3YR(IC,KC)

              BCL4YR(IC,KC) =
     +                        0.1*(STRWYR(IC,KC)-0.0D0)
     +                       -BCL4YR(IC,KC)
              BCL1YR(IC,KC)= HALF*SORPYR(IC,KC)
     +                     + COBCYR*ACOUYR(IC,KC)*(STRPYR(IC,KC)-PINFYR)
     ++0.5*(1.0-BET)*TT1YR(IC,KC)- BCL1YR(IC,KC)
C              17 JAN 2024 VM: change making negative
C     +                     + 100.0*STRDYR(IC,KC)*(ONE/YGDLEN)!ASK NC: CHANGED TO YGDLEN
     +                     - 100.0*STRDYR(IC,KC)*(ONE/YGDLEN)!ASK NC: CHANGED TO YGDLEN
     +               *(ACOUYR(IC,KC)**2.0-STRVYR(IC,KC)**2.0)
     +               *(STRVYR(IC,KC)-0.0)
              ELSE
              BCL1YR(IC,KC)= HALF*SORPYR(IC,KC)
     +                     + COBCYR*ACOUYR(IC,KC)*(STRPYR(IC,KC)-PINFYR)
     ++0.5*(1.0-BET)*TT1YR(IC,KC)- BCL1YR(IC,KC)
              ENDIF
            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
              IF((STRVYR(IC,KC).LT.ZERO).AND.(FLAG_PIO_YR==1))THEN

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)

              URHS(IC,JSTOL,KC) = URHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRUYR(IC,KC)
     +                          - BCL2YR(IC,KC)*STRUYR(IC,KC)
     +                          - BCL3YR(IC,KC)*STRDYR(IC,KC)

              VRHS(IC,JSTOL,KC) = VRHS(IC,JSTOL,KC)
     +      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STRVYR(IC,KC)

              WRHS(IC,JSTOL,KC) = WRHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRWYR(IC,KC)
     +                          - BCL2YR(IC,KC)*STRWYR(IC,KC)
     +                          - BCL4YR(IC,KC)*STRDYR(IC,KC)

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
C              16 JAN 2024: VM+NC sign change
C     +                                     + STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                     - STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                     + OVGMYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STREYR(IC,KC)
     +                      - BCL3YR(IC,KC)*STRDYR(IC,KC)*STRUYR(IC,KC)
     +                      - BCL4YR(IC,KC)*STRDYR(IC,KC)*STRWYR(IC,KC)

              ELSE
              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)

              URHS(IC,JSTOL,KC) = URHS(IC,JSTOL,KC)
     +                   - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRUYR(IC,KC)

              VRHS(IC,JSTOL,KC) = VRHS(IC,JSTOL,KC)
     +   - BCL1YR(IC,KC)*OVA2YR(IC,KC)*(STRVYR(IC,KC)-ACOUYR(IC,KC))

              WRHS(IC,JSTOL,KC) = WRHS(IC,JSTOL,KC)
     +                   - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRWYR(IC,KC)

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +               - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
     +                              - STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                              + OVGMYR(IC,KC))
              ENDIF
            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL
              IF((STRVYR(IC,KC).LT.ZERO).AND.(FLAG_PIO_YR==1))THEN

                FORNOW = BCLYYR(IC,KC,ISPEC)*STRDYR(IC,KC)

                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                            - FORNOW*STRHYR(IC,KC,ISPEC)

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     + - (BCL2YR(IC,KC)+BCL1YR(IC,KC)*OVA2YR(IC,KC))*STRYYR(IC,KC,ISPEC)
     +                                  - FORNOW

              ELSE
                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     +                 - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRYYR(IC,KC,ISPEC)
              ENDIF

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCYR.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYR(IC,KC) = SORPYR(IC,KC)
     +        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L4Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                      *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
C             L5Y UNCHANGED
              BCL1YR(IC,KC) = HALF*SORPYR(IC,KC)
     +                     + COBCYR*ACOUYR(IC,KC)*(STRPYR(IC,KC)-PINFYR)
     +                      - BCL1YR(IC,KC)
              BCL2YR(IC,KC) = -BCL2YR(IC,KC)
              BCL3YR(IC,KC) = -BCL3YR(IC,KC)
              BCL4YR(IC,KC) = -BCL4YR(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYYR(IC,KC,ISPEC) = STRVYR(IC,KC)*BCLYYR(IC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
                BCLYYR(IC,KC,ISPEC) = RATEYR(IC,KC,ISPEC)/STRDYR(IC,KC)
     +                              - BCLYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)

              URHS(IC,JSTOL,KC) = URHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRUYR(IC,KC)
     +                          - BCL2YR(IC,KC)*STRUYR(IC,KC)
     +                          - BCL3YR(IC,KC)*STRDYR(IC,KC)

              VRHS(IC,JSTOL,KC) = VRHS(IC,JSTOL,KC)
     +      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STRVYR(IC,KC)

              WRHS(IC,JSTOL,KC) = WRHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRWYR(IC,KC)
     +                          - BCL2YR(IC,KC)*STRWYR(IC,KC)
     +                          - BCL4YR(IC,KC)*STRDYR(IC,KC)

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
C                     16 JAN 2024: VM+NC sign change
C     +                                     + STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                     - STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                     + OVGMYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STREYR(IC,KC)
     +                      - BCL3YR(IC,KC)*STRDYR(IC,KC)*STRUYR(IC,KC)
     +                      - BCL4YR(IC,KC)*STRDYR(IC,KC)*STRWYR(IC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYYR(IC,KC,ISPEC)*STRDYR(IC,KC)

                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                            - FORNOW*STRHYR(IC,KC,ISPEC)

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     + - (BCL2YR(IC,KC)+BCL1YR(IC,KC)*OVA2YR(IC,KC))*STRYYR(IC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYR.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYR(IC,KC) = ZERO
              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SYDTYR(IC,KC) = SYDTYR(IC,KC)
     +                        + DYDTYR(IC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPYR(IC,KC) = SORPYR(IC,KC)
     +                        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYR(IC,KC) = SYDTYR(IC,KC)/STRRYR(IC,KC)
              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y,L2Y,L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                      *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Y UNCHANGED
              BCL1YR(IC,KC) = BCL5YR(IC,KC)
     +                      + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
     +                      - BCL1YR(IC,KC)
              BCL2YR(IC,KC) = GAM1YR(IC,KC)*OVA2YR(IC,KC)
     +                       *(BCL1YR(IC,KC)+BCL5YR(IC,KC))
     +                      + STRDYR(IC,KC)*(DTDTYR(IC,KC)/STRTYR(IC,KC)
     +                      - SORPYR(IC,KC)/STRPYR(IC,KC)
     +                      + SYDTYR(IC,KC))
     +                      - BCL2YR(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCYR.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                    *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Y UNCHANGED
              FORNOW = BCL5YR(IC,KC)
     +               + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
              BCL1YR(IC,KC) = FORNOW - BCL1YR(IC,KC)
              BCL2YR(IC,KC) = -DDDTYR(IC,KC)
     +                      - OVA2YR(IC,KC)*(BCL1YR(IC,KC)+FORNOW)
     +                      - BCL2YR(IC,KC)
              BCL3YR(IC,KC) = -DUDTYR(IC,KC) - BCL3YR(IC,KC)
              BCL4YR(IC,KC) = -DWDTYR(IC,KC) - BCL4YR(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                  - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
     +                                 + STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                 + OVGMYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STREYR(IC,KC)
     +                  - BCL3YR(IC,KC)*STRDYR(IC,KC)*STRUYR(IC,KC)
     +                  - BCL4YR(IC,KC)*STRDYR(IC,KC)*STRWYR(IC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                BCLYYR(IC,KC,ISPEC) = RATEYR(IC,KC,ISPEC)/STRDYR(IC,KC)
     +                              - DYDTYR(IC,KC,ISPEC)
     +                              - STRVYR(IC,KC)*BCLYYR(IC,KC,ISPEC)

                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +          - BCLYYR(IC,KC,ISPEC)*STRDYR(IC,KC)*STRHYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

C        IF(NSBCYR.EQ.NSBCW1)THEN
C
CC         WALL BOUNDARY CONDITION No 1
CC         NO-SLIP WALL - ADIABATIC
C
CC         ALL VELOCITY COMPONENTS IMPOSED
CC         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
CC         SET IN SUBROUTINE BOUNDT
C
CC         SPECIFY L's AS REQUIRED
CC         L1Y,L3Y-L5Y
C          DO KC = KSTAL,KSTOL
C            DO IC = ISTAL,ISTOL
C           
CC             OLD VALUE OF L's
C              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
C              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
C     +                      *(BCL5YR(IC,KC)-FORNOW)
C              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
C              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
C              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
C     +                      *(BCL5YR(IC,KC)+FORNOW)
C
CC             SUBTRACT FROM NEW VALUE OF L's
CC             L2Y,L5Y UNCHANGED
C              BCL1YR(IC,KC) = BCL5YR(IC,KC)
C     +                      + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
C     +                      - BCL1YR(IC,KC)
C              BCL3YR(IC,KC) = -DUDTYR(IC,KC) - BCL3YR(IC,KC)
C              BCL4YR(IC,KC) = -DWDTYR(IC,KC) - BCL4YR(IC,KC)
C
C            ENDDO
C          ENDDO
C
CC         ADD TO CONSERVATIVE SOURCE TERMS
C          DO KC = KSTAL,KSTOL
C            DO IC = ISTAL,ISTOL
C
C              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
C     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
C
C              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
C     +                  - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
C     +                                 + STRVYR(IC,KC)/ACOUYR(IC,KC)
C     +                                 + OVGMYR(IC,KC))
C     +                  - BCL3YR(IC,KC)*STRDYR(IC,KC)*STRUYR(IC,KC)
C     +                  - BCL4YR(IC,KC)*STRDYR(IC,KC)*STRWYR(IC,KC)
C
C            ENDDO
C          ENDDO
C
C          DO ISPEC = 1,NSPEC
C
C            DO KC = KSTAL,KSTOL
C              DO IC = ISTAL,ISTOL
C
C                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
C     +                 - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRYYR(IC,KC,ISPEC)
C
C              ENDDO
C            ENDDO
C
C          ENDDO
C
C        ENDIF  
C       VM & NC CHANGING ADIABATIC WALL TREATMENT
        IF(NSBCYR.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ADIABATIC

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

          DO KC=KSTAL,KSTOL
            DO IC=ISTAL,ISTOL
              STRTYR(IC,KC)=(48.0*TRUN(IC,JSTOL-1,KC)
     +                      -36.0*TRUN(IC,JSTOL-2,KC)
     +                      +16.0*TRUN(IC,JSTOL-3,KC)
     +                      -3.0*TRUN(IC,JSTOL-4,KC))/25.0
              TRUN(IC,JSTOL,KC)=STRTYR(IC,KC)
            ENDDO
          ENDDO

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYR(IC,KC) = SORPYR(IC,KC)
     +                        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                    *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Y UNCHANGED
              BCL1YR(IC,KC) = BCL5YR(IC,KC)
     +                      - BCL1YR(IC,KC)
C              BCL1YR(IC,KC) = BCL5YR(IC,KC)
C     +                      + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
C     +                      - BCL1YR(IC,KC)
              BCL3YR(IC,KC) = ZERO
              BCL4YR(IC,KC) = ZERO
              BCL2YR(IC,KC) = ZERO
C              BCL3YR(IC,KC) = -DUDTYR(IC,KC) - BCL3YR(IC,KC)
C              BCL4YR(IC,KC) = -DWDTYR(IC,KC) - BCL4YR(IC,KC)
C              BCL2YR(IC,KC) = GAM1YR(IC,KC)*OVA2YR(IC,KC)
C     +                       *(BCL1YR(IC,KC)+BCL5YR(IC,KC))
C     +                      + STRDYR(IC,KC)*(DTDTYR(IC,KC)/STRTYR(IC,KC)
C     +                      - SORPYR(IC,KC)/STRPYR(IC,KC))
C     +                      - BCL2YR(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYY
                BCLYYR(IC,KC,ISPEC) = ZERO
C                BCLYYR(IC,KC,ISPEC) = STRVYR(IC,KC)*BCLYYR(IC,KC,ISPEC)

C               UPDATE L2Y
C                BCL2YR(IC,KC) = BCL2YR(IC,KC)
C     +                        + (RATEYR(IC,KC,ISPEC)
C     +                         - STRDYR(IC,KC)*BCLYYR(IC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRYR(IC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)
              URHS(IC,JSTOL,KC) = ZERO
              VRHS(IC,JSTOL,KC) = ZERO
              WRHS(IC,JSTOL,KC) = ZERO
              ERHS(IC,JSTOL,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     + - (BCL2YR(IC,KC)+BCL1YR(IC,KC)*OVA2YR(IC,KC))*STRYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCYR.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYR(IC,KC) = SORPYR(IC,KC)
     +                        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                    *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Y UNCHANGED
              BCL1YR(IC,KC) = BCL5YR(IC,KC)
     +                      - BCL1YR(IC,KC)
C              BCL1YR(IC,KC) = BCL5YR(IC,KC)
C     +                      + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
C     +                      - BCL1YR(IC,KC)
              BCL3YR(IC,KC) = ZERO
              BCL4YR(IC,KC) = ZERO
              BCL2YR(IC,KC) = ZERO
C              BCL3YR(IC,KC) = -DUDTYR(IC,KC) - BCL3YR(IC,KC)
C              BCL4YR(IC,KC) = -DWDTYR(IC,KC) - BCL4YR(IC,KC)
C              BCL2YR(IC,KC) = GAM1YR(IC,KC)*OVA2YR(IC,KC)
C     +                       *(BCL1YR(IC,KC)+BCL5YR(IC,KC))
C     +                      + STRDYR(IC,KC)*(DTDTYR(IC,KC)/STRTYR(IC,KC)
C     +                      - SORPYR(IC,KC)/STRPYR(IC,KC))
C     +                      - BCL2YR(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYY
                BCLYYR(IC,KC,ISPEC) = ZERO
C                BCLYYR(IC,KC,ISPEC) = STRVYR(IC,KC)*BCLYYR(IC,KC,ISPEC)

C               UPDATE L2Y
C                BCL2YR(IC,KC) = BCL2YR(IC,KC)
C     +                        + (RATEYR(IC,KC,ISPEC)
C     +                         - STRDYR(IC,KC)*BCLYYR(IC,KC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRYR(IC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)
              URHS(IC,JSTOL,KC) = ZERO
              VRHS(IC,JSTOL,KC) = ZERO
              WRHS(IC,JSTOL,KC) = ZERO
              ERHS(IC,JSTOL,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     + - (BCL2YR(IC,KC)+BCL1YR(IC,KC)*OVA2YR(IC,KC))*STRYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Y-DIRECTION RIGHT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Z-DIRECTION LEFT-HAND END
C     -------------------------
      IF(FZLCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUZL = PRIMITIVE U-VELOCITY COMPONENT
C       STRVZL = PRIMITIVE V-VELOCITY COMPONENT
C       STRWZL = PRIMITIVE W-VELOCITY COMPONENT
C       STRPZL = PRESSURE
C       STRDZL = DENSITY
C       STRTZL = TEMPERATURE
C       STREZL = INTERNAL ENERGY
C       STRGZL = MIXTURE CP
C       STRRZL = MIXTURE SPECIFIC GAS CONSTANT
C       STRYZL(ISPEC) = SPECIES MASS FRACTION
C       RATEZL(ISPEC) = SPECIES REACTION RATE
C       STRHZL(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1ZL = DWDZ
C       BCL2ZL = DRHODZ
C       BCL3ZL = DUDZ
C       BCL4ZL = DVDZ
C       BCL5ZL = DPDZ
C       BCLYZL(ISPEC) = DYDZ

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRHZL(IC,JC,ISPEC) = STRHZL(IC,JC,ISPEC)
     +         - STRGZL(IC,JC)*STRTZL(IC,JC)*RGSPEC(ISPEC)/STRRZL(IC,JC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            GAM1ZL(IC,JC) = STRGZL(IC,JC) - STRRZL(IC,JC)
            STREZL(IC,JC) = STREZL(IC,JC) - GAM1ZL(IC,JC)*STRTZL(IC,JC)

            GAM1ZL(IC,JC) = STRRZL(IC,JC)/GAM1ZL(IC,JC)
            OVGMZL(IC,JC) = ONE/GAM1ZL(IC,JC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGZL(IC,JC)*GAM1ZL(IC,JC)*STRTZL(IC,JC)
            ACOUZL(IC,JC) = SQRT(FORNOW)
            OVA2ZL(IC,JC) = ONE/FORNOW

C     THESE ARE INTRODUCED FOR LODATO'S BC- NC

            TT1ZL(IC,JC)=T51BZL(IC,JC)+
     +                   T52BZL(IC,JC)*(GAM1ZL(IC,JC)+1.0)-
     +        STRDZL(IC,JC)* ACOUZL(IC,JC)*T2BZL(IC,JC)         

            TT2ZL(IC,JC)=ACOUZL(IC,JC)*ACOUZL(IC,JC)*
     +      T1BZL(IC,JC)-T51BZL(IC,JC)-(GAM1ZL(IC,JC)+1.0)*
     +      T52BZL(IC,JC)

            TT3ZL(IC,JC)=T3BZL(IC,JC)
            TT4ZL(IC,JC)=T4BZL(IC,JC)
            
            TT5ZL(IC,JC)=T51BZL(IC,JC)+
     +                   T52BZL(IC,JC)*(GAM1ZL(IC,JC)+1.0)+
     +        STRDZL(IC,JC)* ACOUZL(IC,JC)*T2BZL(IC,JC)   

            DO IS=1,NSPEC
            TT6ZL(IC,JC,IS)=T6BZL(IC,JC,IS)
            END DO

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCZL.EQ.NSBCO1)THEN
          FLAG_PIO_ZL=NZLPRM(1)!0
          FLAG_BET_ZL=NZLPRM(2)!1

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZL(IC,JC) = SORPZL(IC,JC)
     +        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L5Z AS REQUIRED
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                      *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)

C             OLD VALUE OF L5Z
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +        *(BCL5ZL(IC,JC)+STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC))

C             SUBTRACT FROM NEW VALUE OF L5Z
C             TERMS INTRODUCED FOR LODATO'S BC- NC
              IF(FLAG_BET_ZL==1) THEN
                BET=STRUZL(IC,JC)*STRUZL(IC,JC)+
     +              STRVZL(IC,JC)*STRVZL(IC,JC)+
     +              STRWZL(IC,JC)*STRWZL(IC,JC)
                BET=SQRT(BET)/ACOUZL(IC,JC)
              END IF
              IF((STRWZL(IC,JC).GT.ZERO).AND.(FLAG_PIO_ZL==1))THEN
              BCL2ZL(IC,JC) =
     +                       -BCL2ZL(IC,JC)
     +                       -OVA2ZL(IC,JC)*SORPZL(IC,JC)

              BCL3ZL(IC,JC) =
     +                        0.1*(STRUZL(IC,JC)-0.0D0)
     +                       -BCL3ZL(IC,JC)

              BCL4ZL(IC,JC) =
     +                        0.1*(STRVZL(IC,JC)-0.0D0)
     +                       -BCL4ZL(IC,JC)
              BCL5ZL(IC,JC)= HALF*SORPZL(IC,JC)
     +                     + COBCZL*ACOUZL(IC,JC)*(STRPZL(IC,JC)-PINFZL)
     ++0.5*(1.0-BET)*TT5ZL(IC,JC)- BCL5ZL(IC,JC)
     +                     + 100.0*STRDZL(IC,JC)*(ONE/ZGDLEN)
     +               *(ACOUZL(IC,JC)**2.0-STRWZL(IC,JC)**2.0)
     +               *(STRWZL(IC,JC)-0.0)
              ELSE
              BCL5ZL(IC,JC)= HALF*SORPZL(IC,JC)
     +                     + COBCZL*ACOUZL(IC,JC)*(STRPZL(IC,JC)-PINFZL)
     ++0.5*(1.0-BET)*TT5ZL(IC,JC)- BCL5ZL(IC,JC)
              ENDIF
            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
              IF((STRWZL(IC,JC).GT.ZERO).AND.(FLAG_PIO_ZL==1))THEN
              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

              URHS(IC,JC,KSTAL) = URHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRUZL(IC,JC)
     +                          - BCL3ZL(IC,JC)*STRDZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRUZL(IC,JC)

              VRHS(IC,JC,KSTAL) = VRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRVZL(IC,JC)
     +                          - BCL4ZL(IC,JC)*STRDZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRVZL(IC,JC)

              WRHS(IC,JC,KSTAL) = WRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRWZL(IC,JC)
     +      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*(STRWZL(IC,JC)+ACOUZL(IC,JC))

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STREZL(IC,JC)
     +                      - BCL3ZL(IC,JC)*STRDZL(IC,JC)*STRUZL(IC,JC)
     +                      - BCL4ZL(IC,JC)*STRDZL(IC,JC)*STRVZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                                     + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                                     + OVGMZL(IC,JC))
              ELSE
              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

              URHS(IC,JC,KSTAL) = URHS(IC,JC,KSTAL)
     +                   - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRUZL(IC,JC)

              VRHS(IC,JC,KSTAL) = VRHS(IC,JC,KSTAL)
     +                   - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRVZL(IC,JC)

              WRHS(IC,JC,KSTAL) = WRHS(IC,JC,KSTAL)
     +   - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*(STRWZL(IC,JC)+ACOUZL(IC,JC))

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +               - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                              + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                              + OVGMZL(IC,JC))
              ENDIF
            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL
              IF((STRWZL(IC,JC).GT.ZERO).AND.(FLAG_PIO_ZL==1))THEN
                FORNOW = BCLYZL(IC,JC,ISPEC)*STRDZL(IC,JC)

                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                            - FORNOW*STRHZL(IC,JC,ISPEC)

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     + - (BCL2ZL(IC,JC)+BCL5ZL(IC,JC)*OVA2ZL(IC,JC))*STRYZL(IC,JC,ISPEC)
     +                                  - FORNOW
              ELSE
                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     +                 - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRYZL(IC,JC,ISPEC)
              ENDIF
              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCZL.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZL(IC,JC) = SORPZL(IC,JC)
     +        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L2Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                      *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
C             L1Z UNCHANGED
              BCL2ZL(IC,JC) = -BCL2ZL(IC,JC)
              BCL3ZL(IC,JC) = -BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = -BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*SORPZL(IC,JC)
     +                     + COBCZL*ACOUZL(IC,JC)*(STRPZL(IC,JC)-PINFZL)
     +                      - BCL5ZL(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYZL(IC,JC,ISPEC) = STRWZL(IC,JC)*BCLYZL(IC,JC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
                BCLYZL(IC,JC,ISPEC) = RATEZL(IC,JC,ISPEC)/STRDZL(IC,JC)
     +                              - BCLYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

              URHS(IC,JC,KSTAL) = URHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRUZL(IC,JC)
     +                          - BCL3ZL(IC,JC)*STRDZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRUZL(IC,JC)

              VRHS(IC,JC,KSTAL) = VRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRVZL(IC,JC)
     +                          - BCL4ZL(IC,JC)*STRDZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRVZL(IC,JC)

              WRHS(IC,JC,KSTAL) = WRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRWZL(IC,JC)
     +      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*(STRWZL(IC,JC)+ACOUZL(IC,JC))

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STREZL(IC,JC)
     +                      - BCL3ZL(IC,JC)*STRDZL(IC,JC)*STRUZL(IC,JC)
     +                      - BCL4ZL(IC,JC)*STRDZL(IC,JC)*STRVZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                                     + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                                     + OVGMZL(IC,JC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYZL(IC,JC,ISPEC)*STRDZL(IC,JC)

                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                            - FORNOW*STRHZL(IC,JC,ISPEC)

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     + - (BCL2ZL(IC,JC)+BCL5ZL(IC,JC)*OVA2ZL(IC,JC))*STRYZL(IC,JC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZL.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZL(IC,JC) = ZERO
              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SYDTZL(IC,JC) = SYDTZL(IC,JC)
     +                        + DYDTZL(IC,JC,ISPEC)*RGSPEC(ISPEC)
                SORPZL(IC,JC) = SORPZL(IC,JC)
     +                        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZL(IC,JC) = SYDTZL(IC,JC)/STRRZL(IC,JC)
              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z,L2Z,L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                      *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Z UNCHANGED
              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
     +                      - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
     +                      - BCL5ZL(IC,JC)
              BCL2ZL(IC,JC) = GAM1ZL(IC,JC)*OVA2ZL(IC,JC)
     +                       *(BCL1ZL(IC,JC)+BCL5ZL(IC,JC))
     +                      + STRDZL(IC,JC)*(DTDTZL(IC,JC)/STRTZL(IC,JC)
     +                      - SORPZL(IC,JC)/STRPZL(IC,JC)
     +                      + SYDTZL(IC,JC))
     +                      - BCL2ZL(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZL.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                    *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Z UNCHANGED
              FORNOW = BCL1ZL(IC,JC)
     +               - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
              BCL2ZL(IC,JC) = -DDDTZL(IC,JC)
     +                      - OVA2ZL(IC,JC)*(BCL1ZL(IC,JC)+FORNOW)
     +                      - BCL2ZL(IC,JC)
              BCL3ZL(IC,JC) = -DUDTZL(IC,JC) - BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = -DVDTZL(IC,JC) - BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = FORNOW - BCL5ZL(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STREZL(IC,JC)
     +                  - BCL3ZL(IC,JC)*STRDZL(IC,JC)*STRUZL(IC,JC)
     +                  - BCL4ZL(IC,JC)*STRDZL(IC,JC)*STRVZL(IC,JC)
     +                  - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                                 + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                                 + OVGMZL(IC,JC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                BCLYZL(IC,JC,ISPEC) = RATEZL(IC,JC,ISPEC)/STRDZL(IC,JC)
     +                              - DYDTZL(IC,JC,ISPEC)
     +                              - STRWZL(IC,JC)*BCLYZL(IC,JC,ISPEC)

                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +          - BCLYZL(IC,JC,ISPEC)*STRDZL(IC,JC)*STRHZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

C        IF(NSBCZL.EQ.NSBCW1)THEN
C
CC         WALL BOUNDARY CONDITION No 1
CC         NO-SLIP WALL - ADIABATIC
C
CC         ALL VELOCITY COMPONENTS IMPOSED
CC         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
CC         SET IN SUBROUTINE BOUNDT
C
CC         SPECIFY L's AS REQUIRED
CC         L1Z,L3Z-L5Z
C          DO JC = JSTAL,JSTOL
C            DO IC = ISTAL,ISTOL
C           
CC             OLD VALUE OF L's
C              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
C              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
C     +                      *(BCL5ZL(IC,JC)-FORNOW)
C              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
C              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
C              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
C     +                      *(BCL5ZL(IC,JC)+FORNOW)
C
CC             SUBTRACT FROM NEW VALUE OF L's
CC             L1Z,L2Z UNCHANGED
C              BCL3ZL(IC,JC) = -DUDTZL(IC,JC) - BCL3ZL(IC,JC)
C              BCL4ZL(IC,JC) = -DVDTZL(IC,JC) - BCL4ZL(IC,JC)
C              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
C     +                      - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
C     +                      - BCL5ZL(IC,JC)
C
C            ENDDO
C          ENDDO
C
CC         ADD TO CONSERVATIVE SOURCE TERMS
C          DO JC = JSTAL,JSTOL
C            DO IC = ISTAL,ISTOL
C
C              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
C     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)
C
C              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
C     +                  - BCL3ZL(IC,JC)*STRDZL(IC,JC)*STRUZL(IC,JC)
C     +                  - BCL4ZL(IC,JC)*STRDZL(IC,JC)*STRVZL(IC,JC)
C     +                  - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
C     +                                 + STRWZL(IC,JC)/ACOUZL(IC,JC)
C     +                                 + OVGMZL(IC,JC))
C
C            ENDDO
C          ENDDO
C
C          DO ISPEC = 1,NSPEC
C
C            DO JC = JSTAL,JSTOL
C              DO IC = ISTAL,ISTOL
C
C                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
C     +                 - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRYZL(IC,JC,ISPEC)
C
C              ENDDO
C            ENDDO
C
C          ENDDO
C
C        ENDIF  
C       VM & NC CHANGING ADIABATIC WALL TREATMENT
        IF(NSBCZL.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ADIABATIC

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

          DO JC=JSTAL,JSTOL
            DO IC=ISTAL,ISTOL
              STRTZL(IC,JC)=(48.0*TRUN(IC,JC,KSTAL+1)
     +                      -36.0*TRUN(IC,JC,KSTAL+2)
     +                      +16.0*TRUN(IC,JC,KSTAL+3)
     +                      -3.0*TRUN(IC,JC,KSTAL+4))/25.0
              TRUN(IC,JC,KSTAL)=STRTZL(IC,JC)
            ENDDO
          ENDDO

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZL(IC,JC) = SORPZL(IC,JC)
     +                        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                    *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL3ZL(IC,JC) = ZERO
              BCL4ZL(IC,JC) = ZERO
C              BCL3ZL(IC,JC) = -DUDTZL(IC,JC) - BCL3ZL(IC,JC)
C              BCL4ZL(IC,JC) = -DVDTZL(IC,JC) - BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
     +                      - BCL5ZL(IC,JC)
C              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
C     +                      - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
C     +                      - BCL5ZL(IC,JC)
              BCL2ZL(IC,JC) = ZERO
C              BCL2ZL(IC,JC) = GAM1ZL(IC,JC)*OVA2ZL(IC,JC)
C     +                       *(BCL1ZL(IC,JC)+BCL5ZL(IC,JC))
C     +                      + STRDZL(IC,JC)*(DTDTZL(IC,JC)/STRTZL(IC,JC)
C     +                      - SORPZL(IC,JC)/STRPZL(IC,JC))
C     +                      - BCL2ZL(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYZ
                BCLYZL(IC,JC,ISPEC) = ZERO
C                BCLYZL(IC,JC,ISPEC) = STRWZL(IC,JC)*BCLYZL(IC,JC,ISPEC)

C               UPDATE L2Z
C                BCL2ZL(IC,JC) = BCL2ZL(IC,JC)
C     +                        + (RATEZL(IC,JC,ISPEC)
C     +                         - STRDZL(IC,JC)*BCLYZL(IC,JC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRZL(IC,JC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)
              URHS(IC,JC,KSTAL) = ZERO
              VRHS(IC,JC,KSTAL) = ZERO
              WRHS(IC,JC,KSTAL) = ZERO
              ERHS(IC,JC,KSTAL) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     + - (BCL2ZL(IC,JC)+BCL5ZL(IC,JC)*OVA2ZL(IC,JC))*STRYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZL.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZL(IC,JC) = SORPZL(IC,JC)
     +                        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                    *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL3ZL(IC,JC) = ZERO
              BCL4ZL(IC,JC) = ZERO
C              BCL3ZL(IC,JC) = -DUDTZL(IC,JC) - BCL3ZL(IC,JC)
C              BCL4ZL(IC,JC) = -DVDTZL(IC,JC) - BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
     +                      - BCL5ZL(IC,JC)
C              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
C     +                      - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
C     +                      - BCL5ZL(IC,JC)
              BCL2ZL(IC,JC) = ZERO
C              BCL2ZL(IC,JC) = GAM1ZL(IC,JC)*OVA2ZL(IC,JC)
C     +                       *(BCL1ZL(IC,JC)+BCL5ZL(IC,JC))
C     +                      + STRDZL(IC,JC)*(DTDTZL(IC,JC)/STRTZL(IC,JC)
C     +                      - SORPZL(IC,JC)/STRPZL(IC,JC))
C     +                      - BCL2ZL(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYZ
                BCLYZL(IC,JC,ISPEC) = ZERO
C                BCLYZL(IC,JC,ISPEC) = STRWZL(IC,JC)*BCLYZL(IC,JC,ISPEC)

C               UPDATE L2Z
C                BCL2ZL(IC,JC) = BCL2ZL(IC,JC)
C     +                        + (RATEZL(IC,JC,ISPEC)
C     +                         - STRDZL(IC,JC)*BCLYZL(IC,JC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRZL(IC,JC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)
              URHS(IC,JC,KSTAL) = ZERO
              VRHS(IC,JC,KSTAL) = ZERO
              WRHS(IC,JC,KSTAL) = ZERO
              ERHS(IC,JC,KSTAL) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     + - (BCL2ZL(IC,JC)+BCL5ZL(IC,JC)*OVA2ZL(IC,JC))*STRYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Z-DIRECTION LEFT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Z-DIRECTION RIGHT-HAND END
C     --------------------------
      IF(FZRCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUZR = PRIMITIVE U-VELOCITY COMPONENT
C       STRVZR = PRIMITIVE V-VELOCITY COMPONENT
C       STRWZR = PRIMITIVE W-VELOCITY COMPONENT
C       STRPZR = PRESSURE
C       STRDZR = DENSITY
C       STRTZR = TEMPERATURE
C       STREZR = INTERNAL ENERGY
C       STRGZR = MIXTURE CP
C       STRRZR = MIXTURE SPECIFIC GAS CONSTANT
C       STRYZR(ISPEC) = SPECIES MASS FRACTION
C       RATEZR(ISPEC) = SPECIES REACTION RATE
C       STRHZR(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1ZR = DWDR
C       BCL2ZR = DRHODR
C       BCL3ZR = DUDR
C       BCL4ZR = DVDR
C       BCL5ZR = DPDR
C       BCLYZR(ISPEC) = DYDR

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRHZR(IC,JC,ISPEC) = STRHZR(IC,JC,ISPEC)
     +         - STRGZR(IC,JC)*STRTZR(IC,JC)*RGSPEC(ISPEC)/STRRZR(IC,JC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            GAM1ZR(IC,JC) = STRGZR(IC,JC) - STRRZR(IC,JC)
            STREZR(IC,JC) = STREZR(IC,JC) - GAM1ZR(IC,JC)*STRTZR(IC,JC)

            GAM1ZR(IC,JC) = STRRZR(IC,JC)/GAM1ZR(IC,JC)
            OVGMZR(IC,JC) = ONE/GAM1ZR(IC,JC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGZR(IC,JC)*GAM1ZR(IC,JC)*STRTZR(IC,JC)
            ACOUZR(IC,JC) = SQRT(FORNOW)
            OVA2ZR(IC,JC) = ONE/FORNOW

C     THESE ARE INTRODUCED FOR LODATO'S BC- NC

            TT1ZR(IC,JC)=T51BZR(IC,JC)+
     +                   T52BZR(IC,JC)*(GAM1ZR(IC,JC)+1.0)-
     +        STRDZR(IC,JC)* ACOUZR(IC,JC)*T2BZR(IC,JC)         

            TT2ZR(IC,JC)=ACOUZR(IC,JC)*ACOUZR(IC,JC)*
     +      T1BZR(IC,JC)-T51BZR(IC,JC)-(GAM1ZR(IC,JC)+1.0)*
     +      T52BZR(IC,JC)

            TT3ZR(IC,JC)=T3BZR(IC,JC)
            TT4ZR(IC,JC)=T4BZR(IC,JC)
            
            TT5ZR(IC,JC)=T51BZR(IC,JC)+
     +                   T52BZR(IC,JC)*(GAM1ZR(IC,JC)+1.0)+
     +        STRDZR(IC,JC)* ACOUZR(IC,JC)*T2BZR(IC,JC)   

            DO IS=1,NSPEC
            TT6ZR(IC,JC,IS)=T6BZR(IC,JC,IS)
            END DO

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCZR.EQ.NSBCO1)THEN
          FLAG_PIO_ZR=NZRPRM(1)!0
          FLAG_BET_ZR=NZRPRM(2)!1

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZR(IC,JC) = SORPZR(IC,JC)
     +        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L1Z AS REQUIRED
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                      *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)
C             OLD VALUE OF L1Z
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +        *(BCL5ZR(IC,JC)-STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC))

C             SUBTRACT FROM NEW VALUE OF L1Z
C             TERMS INTRODUCED FOR LODATO'S BC- NC
              IF(FLAG_BET_ZR==1) THEN
                BET=STRUZR(IC,JC)*STRUZR(IC,JC)+
     +              STRVZR(IC,JC)*STRVZR(IC,JC)+
     +              STRWZR(IC,JC)*STRWZR(IC,JC)
                BET=SQRT(BET)/ACOUZR(IC,JC)
              END IF
              IF((STRWZR(IC,JC).LT.ZERO).AND.(FLAG_PIO_ZR==1))THEN
              BCL2ZR(IC,JC) =
     +                       -BCL2ZR(IC,JC)
     +                       -OVA2ZR(IC,JC)*SORPZR(IC,JC)

              BCL3ZR(IC,JC) =
     +                        0.1*(STRUZR(IC,JC)-0.0D0)
     +                       -BCL3ZR(IC,JC)

              BCL4ZR(IC,JC) =
     +                        0.1*(STRVZR(IC,JC)-0.0D0)
     +                       -BCL4ZR(IC,JC)
              BCL1ZR(IC,JC)= HALF*SORPZR(IC,JC)
     +                     + COBCZR*ACOUZR(IC,JC)*(STRPZR(IC,JC)-PINFZR)
     +                     + 0.5*(1.0-BET)*TT1ZR(IC,JC)- BCL1ZR(IC,JC)
C              19 JAN 2024 VM: change making negative
C     +                     + 100.0*STRDZR(IC,JC)*(ONE/ZGDLEN)
     +                     - 100.0*STRDZR(IC,JC)*(ONE/ZGDLEN)
     +               *(ACOUZR(IC,JC)**2.0-STRWZR(IC,JC)**2.0)
     +               *(STRWZR(IC,JC)-0.0)
              ELSE
              BCL1ZR(IC,JC)= HALF*SORPZR(IC,JC)
     +                     + COBCZR*ACOUZR(IC,JC)*(STRPZR(IC,JC)-PINFZR)
     +                     + 0.5*(1.0-BET)*TT1ZR(IC,JC)- BCL1ZR(IC,JC)
              ENDIF
            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
              IF((STRWZR(IC,JC).LT.ZERO).AND.(FLAG_PIO_ZR==1))THEN
              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)

              URHS(IC,JC,KSTOL) = URHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRUZR(IC,JC)
     +                          - BCL2ZR(IC,JC)*STRUZR(IC,JC)
     +                          - BCL3ZR(IC,JC)*STRDZR(IC,JC)

              VRHS(IC,JC,KSTOL) = VRHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRVZR(IC,JC)
     +                          - BCL2ZR(IC,JC)*STRVZR(IC,JC)
     +                          - BCL4ZR(IC,JC)*STRDZR(IC,JC)

              WRHS(IC,JC,KSTOL) = WRHS(IC,JC,KSTOL)
     +      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STRWZR(IC,JC)

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
C              16 JAN 2024: VM+NC sign change
C     +                                     + STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                     - STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                     + OVGMZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STREZR(IC,JC)
     +                      - BCL3ZR(IC,JC)*STRDZR(IC,JC)*STRUZR(IC,JC)
     +                      - BCL4ZR(IC,JC)*STRDZR(IC,JC)*STRVZR(IC,JC)
              ELSE
              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)

              URHS(IC,JC,KSTOL) = URHS(IC,JC,KSTOL)
     +                   - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRUZR(IC,JC)

              VRHS(IC,JC,KSTOL) = VRHS(IC,JC,KSTOL)
     +                   - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRVZR(IC,JC)

              WRHS(IC,JC,KSTOL) = WRHS(IC,JC,KSTOL)
     +   - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*(STRWZR(IC,JC)-ACOUZR(IC,JC))

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +               - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
     +                              - STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                              + OVGMZR(IC,JC))
              ENDIF
            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL
              IF((STRWZR(IC,JC).LT.ZERO).AND.(FLAG_PIO_ZR==1))THEN
                FORNOW = BCLYZR(IC,JC,ISPEC)*STRDZR(IC,JC)

                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                            - FORNOW*STRHZR(IC,JC,ISPEC)

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     + - (BCL2ZR(IC,JC)+BCL5ZR(IC,JC)*OVA2ZR(IC,JC))*STRYZR(IC,JC,ISPEC)
     +                                  - FORNOW

              ELSE
                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     +                 - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRYZR(IC,JC,ISPEC)
              ENDIF
              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCZR.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZR(IC,JC) = SORPZR(IC,JC)
     +        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L4Y
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                      *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
C             L5Z UNCHANGED
              BCL1ZR(IC,JC) = HALF*SORPZR(IC,JC)
     +                     + COBCZR*ACOUZR(IC,JC)*(STRPZR(IC,JC)-PINFZR)
     +                      - BCL1ZR(IC,JC)
              BCL2ZR(IC,JC) = -BCL2ZR(IC,JC)
              BCL3ZR(IC,JC) = -BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = -BCL4ZR(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYZR(IC,JC,ISPEC) = STRWZR(IC,JC)*BCLYZR(IC,JC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
                BCLYZR(IC,JC,ISPEC) = RATEZR(IC,JC,ISPEC)/STRDZR(IC,JC)
     +                              - BCLYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)

              URHS(IC,JC,KSTOL) = URHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRUZR(IC,JC)
     +                          - BCL2ZR(IC,JC)*STRUZR(IC,JC)
     +                          - BCL3ZR(IC,JC)*STRDZR(IC,JC)

              VRHS(IC,JC,KSTOL) = VRHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRVZR(IC,JC)
     +                          - BCL2ZR(IC,JC)*STRVZR(IC,JC)
     +                          - BCL4ZR(IC,JC)*STRDZR(IC,JC)

              WRHS(IC,JC,KSTOL) = WRHS(IC,JC,KSTOL)
     +      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STRWZR(IC,JC)

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
C              16 JAN 2024: VM+NC sign change
C     +                                     + STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                     - STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                     + OVGMZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STREZR(IC,JC)
     +                      - BCL3ZR(IC,JC)*STRDZR(IC,JC)*STRUZR(IC,JC)
     +                      - BCL4ZR(IC,JC)*STRDZR(IC,JC)*STRVZR(IC,JC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYZR(IC,JC,ISPEC)*STRDZR(IC,JC)

                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                            - FORNOW*STRHZR(IC,JC,ISPEC)

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     + - (BCL2ZR(IC,JC)+BCL5ZR(IC,JC)*OVA2ZR(IC,JC))*STRYZR(IC,JC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZR.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZR(IC,JC) = ZERO
              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SYDTZR(IC,JC) = SYDTZR(IC,JC)
     +                        + DYDTZR(IC,JC,ISPEC)*RGSPEC(ISPEC)
                SORPZR(IC,JC) = SORPZR(IC,JC)
     +                        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZR(IC,JC) = SYDTZR(IC,JC)/STRRZR(IC,JC)
              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z,L2Z,L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                      *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL5ZR(IC,JC) = HALF*(STRWZR(IC,JC)+ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Z UNCHANGED
              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
     +                      + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
     +                      - BCL1ZR(IC,JC)
              BCL2ZR(IC,JC) = GAM1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                       *(BCL1ZR(IC,JC)+BCL5ZR(IC,JC))
     +                      + STRDZR(IC,JC)*(DTDTZR(IC,JC)/STRTZR(IC,JC)
     +                      - SORPZR(IC,JC)/STRPZR(IC,JC)
     +                      + SYDTZR(IC,JC))
     +                      - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZR.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                    *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Z UNCHANGED
              FORNOW = BCL5ZR(IC,JC)
     +               + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
              BCL1ZR(IC,JC) = FORNOW - BCL1ZR(IC,JC)
              BCL2ZR(IC,JC) = -DDDTZR(IC,JC)
     +                      - OVA2ZR(IC,JC)*(BCL1ZR(IC,JC)+FORNOW)
     +                      - BCL2ZR(IC,JC)
              BCL3ZR(IC,JC) = -DUDTZR(IC,JC) - BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = -DVDTZR(IC,JC) - BCL4ZR(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                  - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
     +                                 + STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                 + OVGMZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STREZR(IC,JC)
     +                  - BCL3ZR(IC,JC)*STRDZR(IC,JC)*STRUZR(IC,JC)
     +                  - BCL4ZR(IC,JC)*STRDZR(IC,JC)*STRVZR(IC,JC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                BCLYZR(IC,JC,ISPEC) = RATEZR(IC,JC,ISPEC)/STRDZR(IC,JC)
     +                              - DYDTZR(IC,JC,ISPEC)
     +                              - STRWZR(IC,JC)*BCLYZR(IC,JC,ISPEC)

                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +          - BCLYZR(IC,JC,ISPEC)*STRDZR(IC,JC)*STRHZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

C        IF(NSBCZR.EQ.NSBCW1)THEN
C
CC         WALL BOUNDARY CONDITION No 1
CC         NO-SLIP WALL - ADIABATIC
C
CC         ALL VELOCITY COMPONENTS IMPOSED
CC         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
CC         SET IN SUBROUTINE BOUNDT
C
CC         SPECIFY L's AS REQUIRED
CC         L1Z,L3Z-L5Z
C          DO JC = JSTAL,JSTOL
C            DO IC = ISTAL,ISTOL
C           
CC             OLD VALUE OF L's
C              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
C              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
C     +                      *(BCL5ZR(IC,JC)-FORNOW)
C              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
C              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)
C              BCL5ZR(IC,JC) = HALF*(STRWZR(IC,JC)+ACOUZR(IC,JC))
C     +                      *(BCL5ZR(IC,JC)+FORNOW)
C
CC             SUBTRACT FROM NEW VALUE OF L's
CC             L2Z,L5Z UNCHANGED
C              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
C     +                      + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
C     +                      - BCL1ZR(IC,JC)
C              BCL3ZR(IC,JC) = -DUDTZR(IC,JC) - BCL3ZR(IC,JC)
C              BCL4ZR(IC,JC) = -DVDTZR(IC,JC) - BCL4ZR(IC,JC)
C
C            ENDDO
C          ENDDO
C
CC         ADD TO CONSERVATIVE SOURCE TERMS
C          DO JC = JSTAL,JSTOL
C            DO IC = ISTAL,ISTOL
C
C              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
C     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
C
C              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
C     +                  - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
C     +                                 + STRWZR(IC,JC)/ACOUZR(IC,JC)
C     +                                 + OVGMZR(IC,JC))
C     +                  - BCL3ZR(IC,JC)*STRDZR(IC,JC)*STRUZR(IC,JC)
C     +                  - BCL4ZR(IC,JC)*STRDZR(IC,JC)*STRVZR(IC,JC)
C
C            ENDDO
C          ENDDO
C
C          DO ISPEC = 1,NSPEC
C
C            DO JC = JSTAL,JSTOL
C              DO IC = ISTAL,ISTOL
C
C                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
C     +                 - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRYZR(IC,JC,ISPEC)
C
C              ENDDO
C            ENDDO
C
C          ENDDO
C
C        ENDIF  
C       VM & NC CHANGING ADIABATIC WALL TREATMENT
        IF(NSBCZR.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 1
C         NO-SLIP WALL - ADIABATIC

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

          DO JC=JSTAL,JSTOL
            DO IC=ISTAL,ISTOL
              STRTZR(IC,JC)=(48.0*TRUN(IC,JC,KSTOL-1)
     +                      -36.0*TRUN(IC,JC,KSTOL-2)
     +                      +16.0*TRUN(IC,JC,KSTOL-3)
     +                      -3.0*TRUN(IC,JC,KSTOL-4))/25.0
              TRUN(IC,JC,KSTOL)=STRTZR(IC,JC)
            ENDDO
          ENDDO

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZR(IC,JC) = SORPZR(IC,JC)
     +                        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                    *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)
              BCL5ZR(IC,JC) = HALF*(STRWZR(IC,JC)+ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Z UNCHANGED
              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
     +                      - BCL1ZR(IC,JC)
C              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
C     +                      + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
C     +                      - BCL1ZR(IC,JC)
              BCL3ZR(IC,JC) = ZERO
              BCL4ZR(IC,JC) = ZERO
              BCL2ZR(IC,JC) = ZERO
C              BCL3ZR(IC,JC) = -DUDTZR(IC,JC) - BCL3ZR(IC,JC)
C              BCL4ZR(IC,JC) = -DVDTZR(IC,JC) - BCL4ZR(IC,JC)
C              BCL2ZR(IC,JC) = GAM1ZR(IC,JC)*OVA2ZR(IC,JC)
C     +                       *(BCL1ZR(IC,JC)+BCL5ZR(IC,JC))
C     +                      + STRDZR(IC,JC)*(DTDTZR(IC,JC)/STRTZR(IC,JC)
C     +                      - SORPZR(IC,JC)/STRPZR(IC,JC))
C     +                      - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYZ
                BCLYZR(IC,JC,ISPEC) = ZERO
C                BCLYZR(IC,JC,ISPEC) = STRWZR(IC,JC)*BCLYZR(IC,JC,ISPEC)

C               UPDATE L2Z
C                BCL2ZR(IC,JC) = BCL2ZR(IC,JC)
C     +                        + (RATEZR(IC,JC,ISPEC)
C     +                         - STRDZR(IC,JC)*BCLYZR(IC,JC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRZR(IC,JC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)
              URHS(IC,JC,KSTOL) = ZERO
              VRHS(IC,JC,KSTOL) = ZERO
              WRHS(IC,JC,KSTOL) = ZERO
              ERHS(IC,JC,KSTOL) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     + - (BCL2ZR(IC,JC)+BCL1ZR(IC,JC)*OVA2ZR(IC,JC))*STRYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZR.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZR(IC,JC) = SORPZR(IC,JC)
     +                        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                    *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)
              BCL5ZR(IC,JC) = HALF*(STRWZR(IC,JC)+ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Z UNCHANGED
              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
     +                      - BCL1ZR(IC,JC)
C              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
C     +                      + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
C     +                      - BCL1ZR(IC,JC)
              BCL3ZR(IC,JC) = ZERO
              BCL4ZR(IC,JC) = ZERO
              BCL2ZR(IC,JC) = ZERO
C              BCL3ZR(IC,JC) = -DUDTZR(IC,JC) - BCL3ZR(IC,JC)
C              BCL4ZR(IC,JC) = -DVDTZR(IC,JC) - BCL4ZR(IC,JC)
C              BCL2ZR(IC,JC) = GAM1ZR(IC,JC)*OVA2ZR(IC,JC)
C     +                       *(BCL1ZR(IC,JC)+BCL5ZR(IC,JC))
C     +                      + STRDZR(IC,JC)*(DTDTZR(IC,JC)/STRTZR(IC,JC)
C     +                      - SORPZR(IC,JC)/STRPZR(IC,JC))
C     +                      - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYZ
                BCLYZR(IC,JC,ISPEC) = ZERO
C                BCLYZR(IC,JC,ISPEC) = STRWZR(IC,JC)*BCLYZR(IC,JC,ISPEC)

C               UPDATE L2Z
C                BCL2ZR(IC,JC) = BCL2ZR(IC,JC)
C     +                        + (RATEZR(IC,JC,ISPEC)
C     +                         - STRDZR(IC,JC)*BCLYZR(IC,JC,ISPEC))
C     +                          *RGSPEC(ISPEC)/STRRZR(IC,JC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)
              URHS(IC,JC,KSTOL) = ZERO
              VRHS(IC,JC,KSTOL) = ZERO
              WRHS(IC,JC,KSTOL) = ZERO
              ERHS(IC,JC,KSTOL) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     + - (BCL2ZR(IC,JC)+BCL1ZR(IC,JC)*OVA2ZR(IC,JC))*STRYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Z-DIRECTION RIGHT-HAND END

C     =========================================================================


      RETURN
      END
