      SUBROUTINE RHSCAL 

C     *************************************************************************
C
C     RHSCAL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     12-NOV-2002:  CREATED
C     26-OCT-2008:  RSC/TDD BUG FIX FZLCON
C     08-AUG-2012:  RSC EVALUATE ALL SPECIES
C     17-APR-2013:  RSC MIXTURE AVERAGED TRANSPORT
C     14-JUL-2013:  RSC RADIATION HEAT LOSS
C     08-JUN-2015:  RSC REMOVE Nth SPECIES TREATMENT
C     08-JUN-2015:  RSC UPDATED WALL BCS
C     01-APR-2022:  RSC REFITTED TRANSPORT COEFFICIENTS
C     18-SEP-2022:  RSC BC TRANSVERSE TERMS
C     19-MAY-2023:  RSC BUG FIX GRADIENTS OF LN(TEMPERATURE)
C     26-MAY-2023:  RSC BUG FIX THERMAL DIFFUSION RATIO
C     13-JUN-2023:  RSC/PJB BUG FIX J-INDEX
C     13-JUN-2023:  RSC/PJB BUG FIX J-INDEX
C     22-MAR-2026:  RSC BUG FIX THERMAL DIFFUSION RATIO
C     
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES RIGHT-HAND-SIDES FOR TIME INTEGRATION OF SCALAR PDEs
C     INCLUDES MULTIPLE SCALARS AND MULTI-STEP CHEMISTRY
C     ENERGY EQUATION REQUIRES PRESSURE-WORK AND VISCOUS WORK TERMS
C     COMPUTED IN RHSVEL
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION CTRANS(NSPCMX)
      DOUBLE PRECISION FORNOW
      DOUBLE PRECISION TRATIO,COMBO1,COMBO2
      INTEGER IC,JC,KC,ISPEC,JSPEC
      INTEGER ITINT,ICP,IINDEX,IPOWER,ICOEF1,ICOEF2


C     BEGIN
C     =====

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     EVALUATE THE TEMPERATURE
C     ------------------------
C     ALSO PRESSURE, MIXTURE CP AND MIXTURE GAS CONSTANT
      CALL TEMPER
C                                                               TRANSP = MIX CP
C                                                           STORE7 = RHO*MIX RG
C     =========================================================================

C     COLLECT MIXTURE CP AND GAS CONSTANT FOR BCs
C     -------------------------------------------
 
C     X-DIRECTION
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRGXL(JC,KC) = TRANSP(ISTAL,JC,KC)
            STRRXL(JC,KC) = STORE7(ISTAL,JC,KC)/DRHS(ISTAL,JC,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRGXR(JC,KC) = TRANSP(ISTOL,JC,KC)
            STRRXR(JC,KC) = STORE7(ISTOL,JC,KC)/DRHS(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF

C     Y-DIRECTION
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRGYL(IC,KC) = TRANSP(IC,JSTAL,KC)
            STRRYL(IC,KC) = STORE7(IC,JSTAL,KC)/DRHS(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRGYR(IC,KC) = TRANSP(IC,JSTOL,KC)
            STRRYR(IC,KC) = STORE7(IC,JSTOL,KC)/DRHS(IC,JSTOL,KC)

          ENDDO
        ENDDO
      ENDIF

C     Z-DIRECTION
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRGZL(IC,JC) = TRANSP(IC,JC,KSTAL)
            STRRZL(IC,JC) = STORE7(IC,JC,KSTAL)/DRHS(IC,JC,KSTAL)

          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRGZR(IC,JC) = TRANSP(IC,JC,KSTOL)
            STRRZR(IC,JC) = STORE7(IC,JC,KSTOL)/DRHS(IC,JC,KSTOL)

          ENDDO
        ENDDO
      ENDIF
C                                                               TRANSP = MIX CP
C                                                              ALL STORES CLEAR
C     =========================================================================

C     MASS FLUX DIVERGENCE
C     --------------------
C     URHS,VRHS,WRHS CONTAIN RHO U, RHO V, RHO W

      CALL DFBYDX(URHS,STORE1)
      CALL DFBYDY(VRHS,STORE2)
      CALL DFBYDZ(WRHS,STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            DIVM(IC,JC,KC) = STORE1(IC,JC,KC)
     +                     + STORE2(IC,JC,KC)
     +                     + STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               TRANSP = MIX CP
C                                                              ALL STORES CLEAR
C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     INTERNAL ENERGY EQUATION
C     ========================

C     CONVERT INTERNAL ENERGY
C     -----------------------

C     ERHS CONTAINS RHO E: CONVERT TO E
C     E IS PARALLEL
      DO KC = KSTALT,KSTOLT
        DO JC = JSTALT,JSTOLT
          DO IC = ISTALT,ISTOLT

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)/DRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               TRANSP = MIX CP
C                                                              ALL STORES CLEAR
C     =========================================================================

C     COLLECT INTERNAL ENERGY FOR BCs
C     -------------------------------
 
C     X-DIRECTION
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STREXL(JC,KC) = ERHS(ISTAL,JC,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STREXR(JC,KC) = ERHS(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF

C     Y-DIRECTION
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STREYL(IC,KC) = ERHS(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STREYR(IC,KC) = ERHS(IC,JSTOL,KC)

          ENDDO
        ENDDO
      ENDIF

C     Z-DIRECTION
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STREZL(IC,JC) = ERHS(IC,JC,KSTAL)

          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STREZR(IC,JC) = ERHS(IC,JC,KSTOL)

          ENDDO
        ENDDO
      ENDIF
C                                                               TRANSP = MIX CP
C                                                              ALL STORES CLEAR
C     =========================================================================

C     E EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF E DIV RHO U

C     COLLECT E DIV RHO U IN STORE4 FOR NOW
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE4(IC,JC,KC) = ERHS(IC,JC,KC)*DIVM(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               TRANSP = MIX CP
C                                                          STORE4 = E DIV RHO U
C     =========================================================================

C     E EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF DIV RHO U E

C     D/DX RHO U E
C     RHO U E IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = ERHS(IC,JC,KC)*URHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

      CALL DFBYDX(STORE7,STORE1)

C     D/DY RHO V E
C     RHO V E IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = ERHS(IC,JC,KC)*VRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
      CALL DFBYDY(STORE7,STORE2)

C     D/DZ RHO W E
C     RHO W E IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = ERHS(IC,JC,KC)*WRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
      CALL DFBYDZ(STORE7,STORE3)

C     COLLECT DIV RHO U E IN STORE4 FOR NOW
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE4(IC,JC,KC) = STORE4(IC,JC,KC)
     +                       + STORE1(IC,JC,KC)
     +                       + STORE2(IC,JC,KC)
     +                       + STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               TRANSP = MIX CP
C                                            STORE4 = E DIV RHO U + DIV RHO U E
C     =========================================================================

C     E EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF RHO U.DEL E
      
      CALL DFBYDX(ERHS,STORE1)
      CALL DFBYDY(ERHS,STORE2)
      CALL DFBYDZ(ERHS,STORE3)

C     COLLECT ALL CONVECTIVE TERMS IN ERHS
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = -HALF*(STORE4(IC,JC,KC)
     +                            + STORE1(IC,JC,KC)*URHS(IC,JC,KC)
     +                            + STORE2(IC,JC,KC)*VRHS(IC,JC,KC)
     +                            + STORE3(IC,JC,KC)*WRHS(IC,JC,KC))

          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------
C     E EQUATION: CONVECTIVE TERMS COMPLETE
C     -------------------------------------
C                                                               TRANSP = MIX CP
C                                                              ALL STORES CLEAR
C     =========================================================================

C     E-EQUATION: HEAT FLUX TERMS
C     ---------------------------

C     TEMPERATURE GRADIENTS
      CALL DFBYDX(TRUN,STORE1)
      CALL DFBYDY(TRUN,STORE2)
      CALL DFBYDZ(TRUN,STORE3)
C                                                               TRANSP = MIX CP
C                                                         STORE1,2,3 = DTDX,Y,Z
C     =========================================================================

C     COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs
C     ---------------------------------------------

C     X-DIRECTION
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRTXL(JC,KC) = TRUN(ISTAL,JC,KC)
            BCTXXL(JC,KC) = STORE1(ISTAL,JC,KC)
            BCTYXL(JC,KC) = STORE2(ISTAL,JC,KC)
            BCTZXL(JC,KC) = STORE3(ISTAL,JC,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRTXR(JC,KC) = TRUN(ISTOL,JC,KC)
            BCTXXR(JC,KC) = STORE1(ISTOL,JC,KC)
            BCTYXR(JC,KC) = STORE2(ISTOL,JC,KC)
            BCTZXR(JC,KC) = STORE3(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF

C     Y-DIRECTION
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRTYL(IC,KC) = TRUN(IC,JSTAL,KC)
            BCTXYL(IC,KC) = STORE1(IC,JSTAL,KC)
            BCTYYL(IC,KC) = STORE2(IC,JSTAL,KC)
            BCTZYL(IC,KC) = STORE3(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRTYR(IC,KC) = TRUN(IC,JSTOL,KC)
            BCTXYR(IC,KC) = STORE1(IC,JSTOL,KC)
            BCTYYR(IC,KC) = STORE2(IC,JSTOL,KC)
            BCTZYR(IC,KC) = STORE3(IC,JSTOL,KC)

          ENDDO
        ENDDO
      ENDIF

C     Z-DIRECTION
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRTZL(IC,JC) = TRUN(IC,JC,KSTAL)
            BCTXZL(IC,JC) = STORE1(IC,JC,KSTAL)
            BCTYZL(IC,JC) = STORE2(IC,JC,KSTAL)
            BCTZZL(IC,JC) = STORE3(IC,JC,KSTAL)

          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRTZR(IC,JC) = TRUN(IC,JC,KSTOL)
            BCTXZR(IC,JC) = STORE1(IC,JC,KSTOL)
            BCTYZR(IC,JC) = STORE2(IC,JC,KSTOL)
            BCTZZR(IC,JC) = STORE3(IC,JC,KSTOL)

          ENDDO
        ENDDO
      ENDIF
C                                                               TRANSP = MIX CP
C                                                         STORE1,2,3 = DTDX,Y,Z
C     =========================================================================

C     E-EQUATION: HEAT FLUX TERMS
C     ---------------------------

C     THERMAL CONDUCTIVITY
C     --------------------

C     MIXTURE AVERAGED TRANSPORT
C     RSC 17-APR-2013
C     RSC 01-APR-2022
      IF(FLMAVT)THEN

C       -----------------------------------------------------------------------

C       MIXTURE AVERAGED TRANSPORT

C       MIXTURE MOLAR MASS AND SPECIES MOLE FRACTIONS ARE PARALLEL
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              COMBO1 = ZERO
              DO ISPEC = 1, NSPEC
                CTRANS(ISPEC) = YRHS(IC,JC,KC,ISPEC)*OVWMOL(ISPEC)
                COMBO1 = COMBO1 + CTRANS(ISPEC)
              ENDDO
              COMBO1 = ONE/COMBO1
              WMOMIX(IC,JC,KC) = DRHS(IC,JC,KC)*COMBO1

              DO ISPEC = 1, NSPEC
                XMOLFR(ISPEC,IC,JC,KC) = CTRANS(ISPEC)*COMBO1
              ENDDO

            ENDDO
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       THERMAL CONDUCTIVITY IS PARALLEL
C       STORE CONDUCTIVITY IN STORE7 FOR NOW
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

C             CONDUCTIVITY FOR EACH SPECIES
              TRATIO = TRUN(IC,JC,KC)*TDIFGB
              DO ISPEC = 1, NSPEC
                FORNOW = CONDCO(NCOCON,ISPEC)
                DO ICP = NCOCM1,1,-1
                  FORNOW = FORNOW*TRATIO + CONDCO(ICP,ISPEC)
                ENDDO
                CTRANS(ISPEC) = FORNOW
              ENDDO

C             COMBINATION RULE FOR CONDUCTIVITY
C             RSC/GVN 08-MAR-2014 BUG FIX
              COMBO1 = ZERO
              COMBO2 = ZERO
              DO ISPEC = 1, NSPEC
                FORNOW = XMOLFR(ISPEC,IC,JC,KC)
                COMBO1 = COMBO1 + FORNOW*CTRANS(ISPEC)
                COMBO2 = COMBO2 + FORNOW/CTRANS(ISPEC)
              ENDDO
              STORE7(IC,JC,KC) = HALF*(COMBO1 + ONE/COMBO2) 

            ENDDO
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

      ELSE

C       -----------------------------------------------------------------------

C       CONSTANT LEWIS NUMBER TRANSPORT

C       THERMAL CONDUCTIVITY IS PARALLEL
C       STORE CONDUCTIVITY IN STORE7 FOR NOW
C       STORE CONDUCTIVITY/CP IN TRANSP FOR USE IN DIFFUSIVITY AND VISCOSITY
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              FORNOW = ALAMDA*EXP(RLAMDA*LOG(TRUN(IC,JC,KC)))
              STORE7(IC,JC,KC) = FORNOW*TRANSP(IC,JC,KC)
              TRANSP(IC,JC,KC) = FORNOW

            ENDDO
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

      ENDIF

C     CONDUCTIVITY GRADIENTS
      CALL DFBYDX(STORE7,STORE4)
      CALL DFBYDY(STORE7,STORE5)
      CALL DFBYDZ(STORE7,STORE6)

C     BOUNDARY CONDITIONS
C     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FXLCON)CALL ZEROXL(STORE4)
      IF(FXRCON)CALL ZEROXR(STORE4)
C     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FYLCON)CALL ZEROYL(STORE5)
      IF(FYRCON)CALL ZEROYR(STORE5)
C     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
C     RSC/TDD BUG FIX FZLCON
      IF(FZLCON)CALL ZEROZL(STORE6)
      IF(FZRCON)CALL ZEROZR(STORE6)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     +(STORE4(IC,JC,KC)*STORE1(IC,JC,KC)
     +                     + STORE5(IC,JC,KC)*STORE2(IC,JC,KC)
     +                     + STORE6(IC,JC,KC)*STORE3(IC,JC,KC))
          ENDDO
        ENDDO
      ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         STORE1,2,3 = DTDX,Y,Z
C                                                         STORE7 = CONDUCTIVITY
C     =========================================================================

C     E-EQUATION: HEAT FLUX TERMS
C     ---------------------------
C     WALL BC: THERMAL CONDUCTION TERMS
      IF(FXLCNW)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            FORNOW = ZERO
            DO IC = ISTAP1,ISTOW

              FORNOW = FORNOW
     +            + ACBCXL(IC-1)*STORE7(IC,JC,KC)*STORE1(IC,JC,KC)

            ENDDO
            ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC) + FORNOW
  
          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNW)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            FORNOW = ZERO
            DO IC = ISTAW,ISTOM1

              FORNOW = FORNOW
     +            + ACBCXR(ISTOL-IC)*STORE7(IC,JC,KC)*STORE1(IC,JC,KC)

            ENDDO
            ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC) + FORNOW
  
          ENDDO
        ENDDO
      ENDIF
      IF(FYLCNW)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = ZERO
            DO JC = JSTAP1,JSTOW

              FORNOW = FORNOW
     +            + ACBCYL(JC-1)*STORE7(IC,JC,KC)*STORE2(IC,JC,KC)

            ENDDO
            ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC) + FORNOW
  
          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNW)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = ZERO
            DO JC = JSTAW,JSTOM1

              FORNOW = FORNOW
     +            + ACBCYR(JSTOL-JC)*STORE7(IC,JC,KC)*STORE2(IC,JC,KC)

            ENDDO
            ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC) + FORNOW
  
          ENDDO
        ENDDO
      ENDIF
      IF(FZLCNW)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = ZERO
            DO KC = KSTAP1,KSTOW

              FORNOW = FORNOW
     +            + ACBCZL(KC-1)*STORE7(IC,JC,KC)*STORE3(IC,JC,KC)

            ENDDO
            ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL) + FORNOW
  
          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNW)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = ZERO
            DO KC = KSTAW,KSTOM1

              FORNOW = FORNOW
     +            + ACBCZR(KSTOL-KC)*STORE7(IC,JC,KC)*STORE3(IC,JC,KC)

            ENDDO
            ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL) + FORNOW
  
          ENDDO
        ENDDO
      ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         STORE1,2,3 = DTDX,Y,Z
C     =========================================================================

C     E-EQUATION: HEAT FLUX TERMS
C     ---------------------------
C     SECOND DERIVATIVE TERMS

C     TEMPERATURE SECOND DERIVATIVES
C     RSC 01-APR-2022 USE STORE4,5,6
      CALL D2FDX2(TRUN,STORE4)
      CALL D2FDY2(TRUN,STORE5)
      CALL D2FDZ2(TRUN,STORE6)

C     BOUNDARY CONDITIONS
C     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FXLCON)CALL ZEROXL(STORE4)
      IF(FXRCON)CALL ZEROXR(STORE4)
C     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FYLCON)CALL ZEROYL(STORE5)
      IF(FYRCON)CALL ZEROYR(STORE5)
C     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
C     RSC 28-JUN-2015 BUG FIX FZLCON
      IF(FZLCON)CALL ZEROZL(STORE6)
      IF(FZRCON)CALL ZEROZR(STORE6)

C     COLLECT CONDUCTIVITY TERMS
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     +(STORE4(IC,JC,KC)
     +                     + STORE5(IC,JC,KC)
     +                     + STORE6(IC,JC,KC))
     +                      *STORE7(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     ---------------------------------------------------
C     E-EQUATION: FURTHER HEAT FLUX TERMS EVALUATED BELOW
C     ---------------------------------------------------
C     E-EQUATION: PRESSURE-WORK AND VISCOUS WORK TERMS
C                 EVALUATED IN SUBROUTINE RHSVEL
C     ---------------------------------------------------
C                                                          TRANSP = MIX COND/CP
C     =========================================================================

C     E-EQUATION: RADIATION HEAT LOSS
C     -------------------------------
      IF(FLRADN)CALL RADCAL

C     =========================================================================

C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     SPECIES MASS FRACTION EQUATIONS
C     ===============================

C     REACTION RATE FOR ALL SPECIES
C     -----------------------------
      CALL CHRATE
C                                                          TRANSP = MIX COND/CP
C                                                          RATE = REACTION RATE
C     =========================================================================

C     COLLECT REACTION RATE FOR BCs
C     -----------------------------

C     X-DIRECTION
      IF(FXLCNV)THEN
        DO ISPEC = 1,NSPEC
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              RATEXL(JC,KC,ISPEC) = RATE(ISTAL,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO ISPEC = 1,NSPEC
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              RATEXR(JC,KC,ISPEC) = RATE(ISTOL,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDIF

C     Y-DIRECTION
      IF(FYLCNV)THEN
        DO ISPEC = 1,NSPEC
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              RATEYL(IC,KC,ISPEC) = RATE(IC,JSTAL,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO ISPEC = 1,NSPEC
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              RATEYR(IC,KC,ISPEC) = RATE(IC,JSTOL,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDIF

C     Z-DIRECTION
      IF(FZLCNV)THEN
        DO ISPEC = 1,NSPEC
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              RATEZL(IC,JC,ISPEC) = RATE(IC,JC,KSTAL,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO ISPEC = 1,NSPEC
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
  
              RATEZR(IC,JC,ISPEC) = RATE(IC,JC,KSTOL,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                          RATE = REACTION RATE
C     =========================================================================

C     ZERO THE ACCUMULATORS FOR THE DIFFUSION CORRECTION VELOCITY
C     AND ITS DIVERGENCE
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            UCOR(IC,JC,KC) = ZERO
            VCOR(IC,JC,KC) = ZERO
            WCOR(IC,JC,KC) = ZERO
            VTMP(IC,JC,KC) = ZERO

          ENDDO
        ENDDO
      ENDDO

C     ZERO THE ACCUMULATOR FOR THE MIXTURE ENTHALPY
C     MIXTURE H IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            WTMP(IC,JC,KC) = ZERO

          ENDDO
        ENDDO
      ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                          RATE = REACTION RATE
C                                                             U,V,WCOR CORR VEL
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C     =========================================================================

C     MIXTURE AVERAGED TRANSPORT
C     RSC 17-APR-2013
C     RSC 01-APR-2022
C     RSC 19-MAY-2023
C     EVALUATE FIRST AND SECOND DERIVATIVES
C     OF LN(MIXTURE MOLAR MASS), LN(PRESSURE) AND LN(TEMPERATURE)

C     MIXTURE MOLAR MASS
      IF(FLMIXW)THEN

        CALL DFBYDX(WMOMIX,WD1X)
        CALL DFBYDY(WMOMIX,WD1Y)
        CALL DFBYDZ(WMOMIX,WD1Z)
        CALL D2FDX2(WMOMIX,WD2X)
        CALL D2FDY2(WMOMIX,WD2Y)
        CALL D2FDZ2(WMOMIX,WD2Z)

C       RSC/PJB 13-JUN-2023 BUG FIX J-INDEX
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = ONE/WMOMIX(IC,JC,KC)
              WD1X(IC,JC,KC) = WD1X(IC,JC,KC)*FORNOW
              WD1Y(IC,JC,KC) = WD1Y(IC,JC,KC)*FORNOW
              WD1Z(IC,JC,KC) = WD1Z(IC,JC,KC)*FORNOW
              WD2X(IC,JC,KC) = WD2X(IC,JC,KC)*FORNOW
     +                       - WD1X(IC,JC,KC)*WD1X(IC,JC,KC)
              WD2Y(IC,JC,KC) = WD2Y(IC,JC,KC)*FORNOW
     +                       - WD1Y(IC,JC,KC)*WD1Y(IC,JC,KC)
              WD2Z(IC,JC,KC) = WD2Z(IC,JC,KC)*FORNOW
     +                       - WD1Z(IC,JC,KC)*WD1Z(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

      ENDIF

C     PRESSURE
      IF(FLMIXP)THEN

        CALL DFBYDX(PRUN,PD1X)
        CALL DFBYDY(PRUN,PD1Y)
        CALL DFBYDZ(PRUN,PD1Z)
        CALL D2FDX2(PRUN,PD2X)
        CALL D2FDY2(PRUN,PD2Y)
        CALL D2FDZ2(PRUN,PD2Z)

C       RSC/PJB 13-JUN-2023 BUG FIX J-INDEX
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = ONE/PRUN(IC,JC,KC)
              PD1X(IC,JC,KC) = PD1X(IC,JC,KC)*FORNOW
              PD1Y(IC,JC,KC) = PD1Y(IC,JC,KC)*FORNOW
              PD1Z(IC,JC,KC) = PD1Z(IC,JC,KC)*FORNOW
              PD2X(IC,JC,KC) = PD2X(IC,JC,KC)*FORNOW
     +                       - PD1X(IC,JC,KC)*PD1X(IC,JC,KC)
              PD2Y(IC,JC,KC) = PD2Y(IC,JC,KC)*FORNOW
     +                       - PD1Y(IC,JC,KC)*PD1Y(IC,JC,KC)
              PD2Z(IC,JC,KC) = PD2Z(IC,JC,KC)*FORNOW
     +                       - PD1Z(IC,JC,KC)*PD1Z(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

      ENDIF

C     TEMPERATURE
      IF(FLMIXT)THEN

C       RSC 19-MAY-2023 BUG FIX GRADIENTS OF LN(TEMPERATURE)
        CALL DFBYDX(TRUN,TD1X)
        CALL DFBYDY(TRUN,TD1Y)
        CALL DFBYDZ(TRUN,TD1Z)
        CALL D2FDX2(TRUN,TD2X)
        CALL D2FDY2(TRUN,TD2Y)
        CALL D2FDZ2(TRUN,TD2Z)

C       RSC/PJB 13-JUN-2023 BUG FIX J-INDEX
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = ONE/TRUN(IC,JC,KC)
              TD1X(IC,JC,KC) = TD1X(IC,JC,KC)*FORNOW
              TD1Y(IC,JC,KC) = TD1Y(IC,JC,KC)*FORNOW
              TD1Z(IC,JC,KC) = TD1Z(IC,JC,KC)*FORNOW
              TD2X(IC,JC,KC) = TD2X(IC,JC,KC)*FORNOW
     +                       - TD1X(IC,JC,KC)*TD1X(IC,JC,KC)
              TD2Y(IC,JC,KC) = TD2Y(IC,JC,KC)*FORNOW
     +                       - TD1Y(IC,JC,KC)*TD1Y(IC,JC,KC)
              TD2Z(IC,JC,KC) = TD2Z(IC,JC,KC)*FORNOW
     +                       - TD1Z(IC,JC,KC)*TD1Z(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

      ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                          RATE = REACTION RATE
C                                                             U,V,WCOR CORR VEL
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                              ALL STORES CLEAR
C     =========================================================================

C     RUN THROUGH ALL SPECIES
C     -----------------------
C     RSC 08-AUG-2012 EVALUATE ALL SPECIES
C     RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
      DO ISPEC = 1,NSPEC

C       =======================================================================

C       YRHS CONTAINS RHO Y: CONVERT TO Y
C       Y IS PARALLEL
        DO KC = KSTALT,KSTOLT
          DO JC = JSTALT,JSTOLT
            DO IC = ISTALT,ISTOLT

              YRHS(IC,JC,KC,ISPEC) = YRHS(IC,JC,KC,ISPEC)/DRHS(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

C       =======================================================================

C       Y EQUATION: CONVECTIVE TERMS
C       ----------------------------
C       HALF Y DIV RHO U

C       COLLECT Y SOURCE TERMS IN RATE FOR NOW
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC)
     +                        - HALF*YRHS(IC,JC,KC,ISPEC)*DIVM(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
C                                                         RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                              ALL STORES CLEAR
C       =======================================================================

C       Y EQUATION: CONVECTIVE TERMS
C       ----------------------------
C       HALF DIV RHO U Y

C       D/DX RHO U Y
C       RHO U Y IS PARALLEL
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              STORE7(IC,JC,KC) = YRHS(IC,JC,KC,ISPEC)*URHS(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
        CALL DFBYDX(STORE7,STORE1)

C       D/DY RHO V Y
C       RHO V Y IS PARALLEL
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              STORE7(IC,JC,KC) = YRHS(IC,JC,KC,ISPEC)*VRHS(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
        CALL DFBYDY(STORE7,STORE2)

C       D/DZ RHO W Y
C       RHO W Y IS PARALLEL
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              STORE7(IC,JC,KC) = YRHS(IC,JC,KC,ISPEC)*WRHS(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
        CALL DFBYDZ(STORE7,STORE3)

C       COLLECT DIV RHO U Y IN RATE FOR NOW
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC)
     +                       - HALF*(STORE1(IC,JC,KC)
     +                             + STORE2(IC,JC,KC)
     +                             + STORE3(IC,JC,KC))

            ENDDO
          ENDDO
        ENDDO
C                                       	          RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                              ALL STORES CLEAR
C       =======================================================================

C       SPECIES MASS FRACTION GRADIENT TERMS
C       ------------------------------------

C       SPECIES MASS FRACTION GRADIENTS
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              STORE7(IC,JC,KC) = YRHS(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
        CALL DFBYDX(STORE7,STORE1)
        CALL DFBYDY(STORE7,STORE2)
        CALL DFBYDZ(STORE7,STORE3)
C                                                         RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C       =======================================================================

C       COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs
C       -------------------------------------------------------
C       RSC 18-SEP-2022 BC TRANSVERSE TERMS

C       X-DIRECTION: DYDX
        IF(FXLCNV)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRYXL(JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
              BCYXXL(JC,KC,ISPEC) = STORE1(ISTAL,JC,KC)
              BCYYXL(JC,KC,ISPEC) = STORE2(ISTAL,JC,KC)
              BCYZXL(JC,KC,ISPEC) = STORE3(ISTAL,JC,KC)

            ENDDO
          ENDDO
        ENDIF
        IF(FXRCNV)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRYXR(JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
              BCYXXR(JC,KC,ISPEC) = STORE1(ISTOL,JC,KC)
              BCYYXR(JC,KC,ISPEC) = STORE2(ISTOL,JC,KC)
              BCYZXR(JC,KC,ISPEC) = STORE3(ISTOL,JC,KC)

            ENDDO
          ENDDO
        ENDIF

C       Y-DIRECTION: DYDY
        IF(FYLCNV)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRYYL(IC,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
              BCYXYL(IC,KC,ISPEC) = STORE1(IC,JSTAL,KC)
              BCYYYL(IC,KC,ISPEC) = STORE2(IC,JSTAL,KC)
              BCYZYL(IC,KC,ISPEC) = STORE3(IC,JSTAL,KC)

            ENDDO
          ENDDO
        ENDIF
        IF(FYRCNV)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRYYR(IC,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
              BCYXYR(IC,KC,ISPEC) = STORE1(IC,JSTOL,KC)
              BCYYYR(IC,KC,ISPEC) = STORE2(IC,JSTOL,KC)
              BCYZYR(IC,KC,ISPEC) = STORE3(IC,JSTOL,KC)

            ENDDO
          ENDDO
        ENDIF

C       Z-DIRECTION: DYDZ
        IF(FZLCNV)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRYZL(IC,JC,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
              BCYXZL(IC,JC,ISPEC) = STORE1(IC,JC,KSTAL)
              BCYYZL(IC,JC,ISPEC) = STORE2(IC,JC,KSTAL)
              BCYZZL(IC,JC,ISPEC) = STORE3(IC,JC,KSTAL)

            ENDDO
          ENDDO
        ENDIF
        IF(FZRCNV)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRYZR(IC,JC,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
              BCYXZR(IC,JC,ISPEC) = STORE1(IC,JC,KSTOL)
              BCYYZR(IC,JC,ISPEC) = STORE2(IC,JC,KSTOL)
              BCYZZR(IC,JC,ISPEC) = STORE3(IC,JC,KSTOL)

            ENDDO
          ENDDO
        ENDIF
C                                                         RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C       =======================================================================

C       Y EQUATION: CONVECTIVE TERMS
C       ----------------------------
C       HALF RHO U.DEL Y

C       COLLECT HALF RHO U.DEL Y IN RATE FOR NOW
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC)
     +                       - HALF*(STORE1(IC,JC,KC)*URHS(IC,JC,KC)
     +                             + STORE2(IC,JC,KC)*VRHS(IC,JC,KC)
     +                             + STORE3(IC,JC,KC)*WRHS(IC,JC,KC))

            ENDDO
          ENDDO
        ENDDO

C       ------------------------------------- 
C       Y-EQUATION: CONVECTIVE TERMS COMPLETE
C       ------------------------------------- 
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C       =======================================================================

C       Y-EQUATION: DIFFUSIVE TERMS
C       ---------------------------
C       E-EQUATION: FURTHER HEAT FLUX TERMS

C       MASS DIFFUSIVITY FOR Y
C       --------------------

C       MIXTURE AVERAGED TRANSPORT
C       RSC 17-APR-2013
C       RSC 01-APR-2022
        IF(FLMAVT)THEN

C         ---------------------------------------------------------------------

C         MASS DIFFUSIVITY FOR EACH SPECIES
C         RELATIVE TO CURRENT SPECIES
C         Y DIFFUSIVITY IS PARALLEL
C         STORE DIFFUSIVITY IN STORE7 FOR NOW
          DO KC = KSTAB,KSTOB
            DO JC = JSTAB,JSTOB
              DO IC = ISTAB,ISTOB

C               MASS DIFFUSIVITY FOR EACH SPECIES
                TRATIO = TRUN(IC,JC,KC)*TDIFGB
                DO JSPEC = 1, NSPEC
                  FORNOW = DIFFCO(NCODIF,JSPEC,ISPEC)
                  DO ICP = NCODM1,1,-1
                    FORNOW = FORNOW*TRATIO + DIFFCO(ICP,JSPEC,ISPEC)
                  ENDDO
                  CTRANS(JSPEC) = FORNOW
                ENDDO

C               COMBINATION RULE FOR MASS DIFFUSIVITY
                COMBO1 = ZERO
                DO JSPEC = 1, NSPEC
                  COMBO1 = COMBO1 + XMOLFR(JSPEC,IC,JC,KC)*CTRANS(JSPEC)
                ENDDO
                COMBO1 = COMBO1 - XMOLFR(ISPEC,IC,JC,KC)*CTRANS(ISPEC)
                COMBO1 = COMBO1*PRUN(IC,JC,KC)*PDIFGB
                COMBO1 = (ONE - YRHS(IC,JC,KC,ISPEC))/COMBO1
                DIFMIX(IC,JC,KC) = DRHS(IC,JC,KC)*COMBO1*TRATIO*TRATIO
                STORE7(IC,JC,KC) = DIFMIX(IC,JC,KC)
 
              ENDDO
            ENDDO
          ENDDO

C         ---------------------------------------------------------------------

        ELSE

C         ---------------------------------------------------------------------

C         CONSTANT LEWIS NUMBER TRANSPORT

C         Y DIFFUSIVITY IS PARALLEL
C         STORE DIFFUSIVITY IN STORE7 FOR NOW
          DO KC = KSTAB,KSTOB
            DO JC = JSTAB,JSTOB
              DO IC = ISTAB,ISTOB

                STORE7(IC,JC,KC) = TRANSP(IC,JC,KC)*OLEWIS(ISPEC)

              ENDDO
            ENDDO
          ENDDO

C         ---------------------------------------------------------------------

        ENDIF

C       -----------------------------------------------------------------------

C       MIXTURE AVERAGED TRANSPORT
C       RSC 17-APR-2013
C       RSC 01-APR-2022
C       RSC 26-MAY-2023 BUG FIX THERMAL DIFFUSION RATIO

C       THERMAL DIFFUSION RATIO FOR EACH SPECIES
C       RELATIVE TO CURRENT SPECIES
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              TDRMIX(IC,JC,KC) = ZERO

            ENDDO
          ENDDO
        ENDDO
  
        IF(FLMTDR(ISPEC))THEN

          DO JSPEC = 1, NSPEC

            DO KC = KSTAB,KSTOB
              DO JC = JSTAB,JSTOB
                DO IC = ISTAB,ISTOB

C                 THERMAL DIFFUSION RATIO FOR THIS SPECIES PAIR
                  TRATIO = TRUN(IC,JC,KC)*TDIFGB
                  FORNOW = TDRCCO(NCOTDR,JSPEC,ISPEC)
                  DO ICP = NCOTM1,1,-1
                    FORNOW = FORNOW*TRATIO + TDRCCO(ICP,JSPEC,ISPEC)
                  ENDDO

C                 COMBINATION RULE FOR THERMAL DIFFUSIION RATIO
                  TDRMIX(IC,JC,KC) = TDRMIX(IC,JC,KC)
     +                             + XMOLFR(JSPEC,IC,JC,KC)*FORNOW
  
                ENDDO
              ENDDO
            ENDDO
  
          ENDDO

        ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       DIFFUSION CORRECTION VELOCITY
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              UCOR(IC,JC,KC) = UCOR(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE1(IC,JC,KC)
              VCOR(IC,JC,KC) = VCOR(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE2(IC,JC,KC)
              WCOR(IC,JC,KC) = WCOR(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE3(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       SPECIES ENTHALPY
C       ----------------

C       TEMPERATURE INTERVAL INDEXING 
        IINDEX = 1 + (ISPEC-1)/NSPIMX
        IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
        ICOEF2 = NTBASE**IPOWER
        ICOEF1 = ICOEF2*NTBASE

C       SPECIES H IS PARALLEL
C       STORE SPECIES H IN UTMP FOR NOW
C       STORE MIXTURE H IN WTMP FOR NOW
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              ITINT = 1 + MOD(ITNDEX(IC,JC,KC,IINDEX),ICOEF1)/ICOEF2
              FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
              DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                FORNOW = FORNOW*TRUN(IC,JC,KC) + AMASCH(ICP,ITINT,ISPEC)
              ENDDO
              UTMP(IC,JC,KC) = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                       + FORNOW*TRUN(IC,JC,KC)

C             MIXTURE H
              WTMP(IC,JC,KC) = WTMP(IC,JC,KC)
     +                       + UTMP(IC,JC,KC)*YRHS(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
C                                                          TRANSP = MIX COND/CP
C                                       	          RATE = Y SOURCE TERMS
C                                       	               UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       COLLECT SPECIES H FOR BCs
C       -------------------------
 
C       X-DIRECTION
        IF(FXLCNV)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRHXL(JC,KC,ISPEC) = UTMP(ISTAL,JC,KC)

            ENDDO
          ENDDO
        ENDIF
        IF(FXRCNV)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRHXR(JC,KC,ISPEC) = UTMP(ISTOL,JC,KC)

            ENDDO
          ENDDO
        ENDIF

C       Y-DIRECTION
        IF(FYLCNV)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRHYL(IC,KC,ISPEC) = UTMP(IC,JSTAL,KC)

            ENDDO
          ENDDO
        ENDIF
        IF(FYRCNV)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRHYR(IC,KC,ISPEC) = UTMP(IC,JSTOL,KC)

            ENDDO
          ENDDO
        ENDIF

C       Z-DIRECTION
        IF(FZLCNV)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRHZL(IC,JC,ISPEC) = UTMP(IC,JC,KSTAL)

            ENDDO
          ENDDO
        ENDIF
        IF(FZRCNV)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRHZR(IC,JC,ISPEC) = UTMP(IC,JC,KSTOL)

            ENDDO
          ENDDO
        ENDIF
C                                                          TRANSP = MIX COND/CP
C                                       	          RATE = Y SOURCE TERMS
C                                       	               UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       MIXTURE AVERAGED TRANSPORT
C       RSC 23-APR-2013
C       ADD DUFOUR EFFECT TERMS TO SPECIES ENTHALPY
        IF(FLMDUF(ISPEC))THEN

          DO KC = KSTAB,KSTOB
            DO JC = JSTAB,JSTOB
              DO IC = ISTAB,ISTOB

                UTMP(IC,JC,KC) = UTMP(IC,JC,KC)
     +           + RGSPEC(ISPEC)*TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO

        ENDIF
C                                                          TRANSP = MIX COND/CP
C                                       	          RATE = Y SOURCE TERMS
C                                       	               UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                       	               WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       Y EQUATION: DIFFUSIVE TERMS
C       ---------------------------
C       E EQUATION: FURTHER HEAT FLUX TERMS
C       DIFFUSIVITY GRADIENT TERMS

C       DIFFUSIVITY GRADIENTS
        CALL DFBYDX(STORE7,STORE4)
        CALL DFBYDY(STORE7,STORE5)
        CALL DFBYDZ(STORE7,STORE6)

C       BOUNDARY CONDITIONS
C       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FXLDIF)CALL ZEROXL(STORE4)
        IF(FXRDIF)CALL ZEROXR(STORE4)
C       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FYLDIF)CALL ZEROYL(STORE5)
        IF(FYRDIF)CALL ZEROYR(STORE5)
C       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FZLDIF)CALL ZEROZL(STORE6)
        IF(FZRDIF)CALL ZEROZR(STORE6)

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = STORE4(IC,JC,KC)*STORE1(IC,JC,KC)
     +               + STORE5(IC,JC,KC)*STORE2(IC,JC,KC)
     +               + STORE6(IC,JC,KC)*STORE3(IC,JC,KC)
 
C             Y EQUATION
              RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C             DIFFUSION CORRECTION VELOCITY DIVERGENCE
              VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

            ENDDO
          ENDDO
        ENDDO

C       BOUNDARY CONDITIONS
C       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FXLADB)CALL ZEROXL(STORE4)
        IF(FXRADB)CALL ZEROXR(STORE4)
C       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FYLADB)CALL ZEROYL(STORE5)
        IF(FYRADB)CALL ZEROYR(STORE5)
C       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FZLADB)CALL ZEROZL(STORE6)
        IF(FZRADB)CALL ZEROZR(STORE6)

C       E EQUATION
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = STORE4(IC,JC,KC)*STORE1(IC,JC,KC)
     +               + STORE5(IC,JC,KC)*STORE2(IC,JC,KC)
     +               + STORE6(IC,JC,KC)*STORE3(IC,JC,KC)
 
              ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC) 
              
            ENDDO
          ENDDO
        ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       E-EQUATION: FURTHER HEAT FLUX TERMS
C       -----------------------------------
C       SPECIES ENTHALPY GRADIENT TERMS

C       SPECIES ENTHALPY GRADIENTS
        CALL DFBYDX(UTMP,STORE4)
        CALL DFBYDY(UTMP,STORE5)
        CALL DFBYDZ(UTMP,STORE6)

C       BOUNDARY CONDITIONS
C       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FXLDIF)CALL ZEROXL(STORE4)
        IF(FXRDIF)CALL ZEROXR(STORE4)
C       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FYLDIF)CALL ZEROYL(STORE5)
        IF(FYRDIF)CALL ZEROYR(STORE5)
C       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FZLDIF)CALL ZEROZL(STORE6)
        IF(FZRDIF)CALL ZEROZR(STORE6)

C       BOUNDARY CONDITIONS
C       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FXLADB)CALL ZEROXL(STORE4)
        IF(FXRADB)CALL ZEROXR(STORE4)
C       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FYLADB)CALL ZEROYL(STORE5)
        IF(FYRADB)CALL ZEROYR(STORE5)
C       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FZLADB)CALL ZEROZL(STORE6)
        IF(FZRADB)CALL ZEROZR(STORE6)

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = STORE4(IC,JC,KC)*STORE1(IC,JC,KC)
     +               + STORE5(IC,JC,KC)*STORE2(IC,JC,KC)
     +               + STORE6(IC,JC,KC)*STORE3(IC,JC,KC)
 
              ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*STORE7(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                         STORE1,2,3 = DYDX,Y,Z
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       Y-EQUATION: DIFFUSIVE TERMS
C       ---------------------------
C       WALL BC: MASS DIFFUSION TERMS
C       E-EQUATION: HEAT FLUX TERMS
C       WALL BC: ENTHALPY DIFFUSION TERMS
        IF(FXLDFW)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              FORNOW = ZERO
              DO IC = ISTAP1,ISTOW

                FORNOW = FORNOW
     +              + ACBCXL(IC-1)*STORE7(IC,JC,KC)*STORE1(IC,JC,KC)

              ENDDO
              RATE(ISTAL,JC,KC,ISPEC) = RATE(ISTAL,JC,KC,ISPEC)
     +                                + FORNOW
              VTMP(ISTAL,JC,KC) = VTMP(ISTAL,JC,KC) + FORNOW
              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                          + FORNOW*UTMP(ISTAL,JC,KC)
  
            ENDDO
          ENDDO
        ENDIF
        IF(FXRDFW)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              FORNOW = ZERO
              DO IC = ISTAW,ISTOM1

                FORNOW = FORNOW
     +              + ACBCXR(ISTOL-IC)*STORE7(IC,JC,KC)*STORE1(IC,JC,KC)

              ENDDO
              RATE(ISTOL,JC,KC,ISPEC) = RATE(ISTOL,JC,KC,ISPEC)
     +                                + FORNOW
              VTMP(ISTOL,JC,KC) = VTMP(ISTOL,JC,KC) + FORNOW
              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                          + FORNOW*UTMP(ISTOL,JC,KC)
  
            ENDDO
          ENDDO
        ENDIF
        IF(FYLDFW)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = ZERO
              DO JC = JSTAP1,JSTOW

                FORNOW = FORNOW
     +              + ACBCYL(JC-1)*STORE7(IC,JC,KC)*STORE2(IC,JC,KC)

              ENDDO
              RATE(IC,JSTAL,KC,ISPEC) = RATE(IC,JSTAL,KC,ISPEC)
     +                                + FORNOW
              VTMP(IC,JSTAL,KC) = VTMP(IC,JSTAL,KC) + FORNOW
              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                          + FORNOW*UTMP(IC,JSTAL,KC)
  
            ENDDO
          ENDDO
        ENDIF
        IF(FYRDFW)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = ZERO
              DO JC = JSTAW,JSTOM1

                FORNOW = FORNOW
     +              + ACBCYR(JSTOL-JC)*STORE7(IC,JC,KC)*STORE2(IC,JC,KC)

              ENDDO
              RATE(IC,JSTOL,KC,ISPEC) = RATE(IC,JSTOL,KC,ISPEC)
     +                                + FORNOW
              VTMP(IC,JSTOL,KC) = VTMP(IC,JSTOL,KC) + FORNOW
              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                          + FORNOW*UTMP(IC,JSTOL,KC)
  
            ENDDO
          ENDDO
        ENDIF
        IF(FZLDFW)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = ZERO
              DO KC = KSTAP1,KSTOW

                FORNOW = FORNOW
     +              + ACBCZL(KC-1)*STORE7(IC,JC,KC)*STORE3(IC,JC,KC)

              ENDDO
              RATE(IC,JC,KSTAL,ISPEC) = RATE(IC,JC,KSTAL,ISPEC)
     +                                + FORNOW
              VTMP(IC,JC,KSTAL) = VTMP(IC,JC,KSTAL) + FORNOW
              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                          + FORNOW*UTMP(IC,JC,KSTAL)
  
            ENDDO
          ENDDO
        ENDIF
        IF(FZRDFW)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = ZERO
              DO KC = KSTAW,KSTOM1

                FORNOW = FORNOW
     +              + ACBCZR(KSTOL-KC)*STORE7(IC,JC,KC)*STORE3(IC,JC,KC)

              ENDDO
              RATE(IC,JC,KSTOL,ISPEC) = RATE(IC,JC,KSTOL,ISPEC)
     +                                + FORNOW
              VTMP(IC,JC,KSTOL) = VTMP(IC,JC,KSTOL) + FORNOW
              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                          + FORNOW*UTMP(IC,JC,KSTOL)
  
            ENDDO
          ENDDO
        ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                          STORE7 = DIFFUSIVITY
C       =======================================================================

C       Y-EQUATION: DIFFUSIVE TERMS
C       ---------------------------
C       E-EQUATION: FURTHER HEAT FLUX TERMS
C       SECOND DERIVATIVE TERMS

C       SPECIES MASS FRACTION SECOND DERIVATIVES
C       MOVE DIFFUSIVITY TO STORE4
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STORE4(IC,JC,KC) = STORE7(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
C       MOVE MASS FRACTION TO STORE7
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              STORE7(IC,JC,KC) = YRHS(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
        CALL D2FDX2(STORE7,STORE1)
        CALL D2FDY2(STORE7,STORE2)
        CALL D2FDZ2(STORE7,STORE3)

C       BOUNDARY CONDITIONS
C       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FXLDIF)CALL ZEROXL(STORE1)
        IF(FXRDIF)CALL ZEROXR(STORE1)
C       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FYLDIF)CALL ZEROYL(STORE2)
        IF(FYRDIF)CALL ZEROYR(STORE2)
C       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FZLDIF)CALL ZEROZL(STORE3)
        IF(FZRDIF)CALL ZEROZR(STORE3)

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = (STORE1(IC,JC,KC)
     +               +  STORE2(IC,JC,KC)
     +               +  STORE3(IC,JC,KC))*STORE4(IC,JC,KC)

C             Y EQUATION
              RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C             DIFFUSION CORRECTION VELOCITY DIVERGENCE
              VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

            ENDDO
          ENDDO
        ENDDO

C       BOUNDARY CONDITIONS
C       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FXLADB)CALL ZEROXL(STORE1)
        IF(FXRADB)CALL ZEROXR(STORE1)
C       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FYLADB)CALL ZEROYL(STORE2)
        IF(FYRADB)CALL ZEROYR(STORE2)
C       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(FZLADB)CALL ZEROZL(STORE3)
        IF(FZRADB)CALL ZEROZR(STORE3)

C       E EQUATION
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = (STORE1(IC,JC,KC)
     +               +  STORE2(IC,JC,KC)
     +               +  STORE3(IC,JC,KC))*STORE4(IC,JC,KC)

              ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              ALL STORES CLEAR
C       =======================================================================

C       MIXTURE AVERAGED TRANSPORT
C       RSC 23-APR-2013
C       MOLAR MASS TERMS, PRESSURE TERMS, SORET EFFECT

C       MIXTURE MOLAR MASS TERMS
        IF(FLMIXW)THEN

C         FIRST AND SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED

          DO KC = KSTAB,KSTOB
            DO JC = JSTAB,JSTOB
              DO IC = ISTAB,ISTOB

                STORE7(IC,JC,KC) = DIFMIX(IC,JC,KC)*YRHS(IC,JC,KC,ISPEC)
  
              ENDDO
            ENDDO
          ENDDO
  
C         DIFFUSION CORRECTION VELOCITY
C         FIRST DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                UCOR(IC,JC,KC) = UCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*WD1X(IC,JC,KC)
                VCOR(IC,JC,KC) = VCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*WD1Y(IC,JC,KC)
                WCOR(IC,JC,KC) = WCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*WD1Z(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO


C         Y EQUATION: DIFFUSIVE TERMS
C         E EQUATION: FURTHER HEAT FLUX TERMS

C         DIFFUSIVITY GRADIENT TERMS

C         DIFFUSIVITY GRADIENTS
          CALL DFBYDX(STORE7,STORE1)
          CALL DFBYDY(STORE7,STORE2)
          CALL DFBYDZ(STORE7,STORE3)

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(STORE1)
          IF(FXRDIF)CALL ZEROXR(STORE1)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(STORE2)
          IF(FYRDIF)CALL ZEROYR(STORE2)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(STORE3)
          IF(FZRDIF)CALL ZEROZR(STORE3)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE1(IC,JC,KC)*WD1X(IC,JC,KC)
     +                 + STORE2(IC,JC,KC)*WD1Y(IC,JC,KC)
     +                 + STORE3(IC,JC,KC)*WD1Z(IC,JC,KC)
 
C               Y EQUATION
                RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C               DIFFUSION CORRECTION VELOCITY DIVERGENCE
                VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

              ENDDO
            ENDDO
          ENDDO

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(STORE1)
          IF(FXRADB)CALL ZEROXR(STORE1)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(STORE2)
          IF(FYRADB)CALL ZEROYR(STORE2)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(STORE3)
          IF(FZRADB)CALL ZEROZR(STORE3)

C         E EQUATION
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE1(IC,JC,KC)*WD1X(IC,JC,KC)
     +                 + STORE2(IC,JC,KC)*WD1Y(IC,JC,KC)
     +                 + STORE3(IC,JC,KC)*WD1Z(IC,JC,KC)
 
                ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC) 
              
              ENDDO
            ENDDO
          ENDDO


C         E-EQUATION: FURTHER HEAT FLUX TERMS
C         SPECIES ENTHALPY GRADIENT TERMS

C         SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(STORE4)
          IF(FXRDIF)CALL ZEROXR(STORE4)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(STORE5)
          IF(FYRDIF)CALL ZEROYR(STORE5)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(STORE6)
          IF(FZRDIF)CALL ZEROZR(STORE6)

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(STORE4)
          IF(FXRADB)CALL ZEROXR(STORE4)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(STORE5)
          IF(FYRADB)CALL ZEROYR(STORE5)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(STORE6)
          IF(FZRADB)CALL ZEROZR(STORE6)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE4(IC,JC,KC)*WD1X(IC,JC,KC)
     +                 + STORE5(IC,JC,KC)*WD1Y(IC,JC,KC)
     +                 + STORE6(IC,JC,KC)*WD1Z(IC,JC,KC)
 
              ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*STORE7(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C         =====================================================================

C         Y-EQUATION: DIFFUSIVE TERMS
C         ---------------------------
C         WALL BC: MOLAR MASS TERMS
C         E-EQUATION: HEAT FLUX TERMS
C         WALL BC: ENTHALPY DIFFUSION TERMS
          IF(FXLDFW)THEN
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = ZERO
                DO IC = ISTAP1,ISTOW

                  FORNOW = FORNOW
     +                + ACBCXL(IC-1)*STORE7(IC,JC,KC)*WD1X(IC,JC,KC)

                ENDDO
                RATE(ISTAL,JC,KC,ISPEC) = RATE(ISTAL,JC,KC,ISPEC)
     +                                  + FORNOW
                VTMP(ISTAL,JC,KC) = VTMP(ISTAL,JC,KC) + FORNOW
                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                            + FORNOW*UTMP(ISTAL,JC,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FXRDFW)THEN
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = ZERO
                DO IC = ISTAW,ISTOM1

                  FORNOW = FORNOW
     +                + ACBCXR(ISTOL-IC)*STORE7(IC,JC,KC)*WD1X(IC,JC,KC)

                ENDDO
                RATE(ISTOL,JC,KC,ISPEC) = RATE(ISTOL,JC,KC,ISPEC)
     +                                  + FORNOW
                VTMP(ISTOL,JC,KC) = VTMP(ISTOL,JC,KC) + FORNOW
                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                            + FORNOW*UTMP(ISTOL,JC,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FYLDFW)THEN
            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO JC = JSTAP1,JSTOW

                  FORNOW = FORNOW
     +                + ACBCYL(JC-1)*STORE7(IC,JC,KC)*WD1Y(IC,JC,KC)

                ENDDO
                RATE(IC,JSTAL,KC,ISPEC) = RATE(IC,JSTAL,KC,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JSTAL,KC) = VTMP(IC,JSTAL,KC) + FORNOW
                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                            + FORNOW*UTMP(IC,JSTAL,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FYRDFW)THEN
            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO JC = JSTAW,JSTOM1

                  FORNOW = FORNOW
     +                + ACBCYR(JSTOL-JC)*STORE7(IC,JC,KC)*WD1Y(IC,JC,KC)

                ENDDO
                RATE(IC,JSTOL,KC,ISPEC) = RATE(IC,JSTOL,KC,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JSTOL,KC) = VTMP(IC,JSTOL,KC) + FORNOW
                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                            + FORNOW*UTMP(IC,JSTOL,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FZLDFW)THEN
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO KC = KSTAP1,KSTOW

                  FORNOW = FORNOW
     +                + ACBCZL(KC-1)*STORE7(IC,JC,KC)*WD1Z(IC,JC,KC)

                ENDDO
                RATE(IC,JC,KSTAL,ISPEC) = RATE(IC,JC,KSTAL,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JC,KSTAL) = VTMP(IC,JC,KSTAL) + FORNOW
                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                            + FORNOW*UTMP(IC,JC,KSTAL)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FZRDFW)THEN
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO KC = KSTAW,KSTOM1

                  FORNOW = FORNOW
     +                + ACBCZR(KSTOL-KC)*STORE7(IC,JC,KC)*WD1Z(IC,JC,KC)

                ENDDO
                RATE(IC,JC,KSTOL,ISPEC) = RATE(IC,JC,KSTOL,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JC,KSTOL) = VTMP(IC,JC,KSTOL) + FORNOW
                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                            + FORNOW*UTMP(IC,JC,KSTOL)
  
              ENDDO
            ENDDO
          ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C         =====================================================================

C         Y-EQUATION: DIFFUSIVE TERMS
C         E-EQUATION: FURTHER HEAT FLUX TERMS
C         SECOND DERIVATIVE TERMS
C         SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(WD2X)
          IF(FXRDIF)CALL ZEROXR(WD2X)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(WD2Y)
          IF(FYRDIF)CALL ZEROYR(WD2Y)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(WD2Z)
          IF(FZRDIF)CALL ZEROZR(WD2Z)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = (WD2X(IC,JC,KC)
     +                 +  WD2Y(IC,JC,KC)
     +                 +  WD2Z(IC,JC,KC))*STORE7(IC,JC,KC)

C               Y EQUATION
                RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C               DIFFUSION CORRECTION VELOCITY DIVERGENCE
                VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

              ENDDO
            ENDDO
          ENDDO

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(WD2X)
          IF(FXRADB)CALL ZEROXR(WD2X)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(WD2Y)
          IF(FYRADB)CALL ZEROYR(WD2Y)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(WD2Z)
          IF(FZRADB)CALL ZEROZR(WD2Z)

C         E EQUATION
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = (WD2X(IC,JC,KC)
     +                 +  WD2Y(IC,JC,KC)
     +                 +  WD2Z(IC,JC,KC))*STORE7(IC,JC,KC)

                ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO

        ENDIF
C       MIXTURE MOLAR MASS TERMS
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C       =======================================================================

C       PRESSURE DIFFUSION TERMS
        IF(FLMIXP)THEN

C         FIRST AND SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED

          DO KC = KSTAB,KSTOB
            DO JC = JSTAB,JSTOB
              DO IC = ISTAB,ISTOB

                STORE7(IC,JC,KC) = DIFMIX(IC,JC,KC)*YRHS(IC,JC,KC,ISPEC)
     +                             *(ONE-WMOLAR(ISPEC)/WMOMIX(IC,JC,KC))
  
              ENDDO
            ENDDO
          ENDDO
  
C         DIFFUSION CORRECTION VELOCITY
C         FIRST DERIVATIVES OF LN(PRESSURE) ALREADY STORED
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                UCOR(IC,JC,KC) = UCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*PD1X(IC,JC,KC)
                VCOR(IC,JC,KC) = VCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*PD1Y(IC,JC,KC)
                WCOR(IC,JC,KC) = WCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*PD1Z(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO


C         Y EQUATION: DIFFUSIVE TERMS
C         E EQUATION: FURTHER HEAT FLUX TERMS

C         DIFFUSIVITY GRADIENT TERMS

C         DIFFUSIVITY GRADIENTS
          CALL DFBYDX(STORE7,STORE1)
          CALL DFBYDY(STORE7,STORE2)
          CALL DFBYDZ(STORE7,STORE3)

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(STORE1)
          IF(FXRDIF)CALL ZEROXR(STORE1)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(STORE2)
          IF(FYRDIF)CALL ZEROYR(STORE2)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(STORE3)
          IF(FZRDIF)CALL ZEROZR(STORE3)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE1(IC,JC,KC)*PD1X(IC,JC,KC)
     +                 + STORE2(IC,JC,KC)*PD1Y(IC,JC,KC)
     +                 + STORE3(IC,JC,KC)*PD1Z(IC,JC,KC)
 
C               Y EQUATION
                RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C               DIFFUSION CORRECTION VELOCITY DIVERGENCE
                VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

              ENDDO
            ENDDO
          ENDDO

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(STORE1)
          IF(FXRADB)CALL ZEROXR(STORE1)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(STORE2)
          IF(FYRADB)CALL ZEROYR(STORE2)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(STORE3)
          IF(FZRADB)CALL ZEROZR(STORE3)

C         E EQUATION
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE1(IC,JC,KC)*PD1X(IC,JC,KC)
     +                 + STORE2(IC,JC,KC)*PD1Y(IC,JC,KC)
     +                 + STORE3(IC,JC,KC)*PD1Z(IC,JC,KC)
 
                ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC) 
              
              ENDDO
            ENDDO
          ENDDO


C         E-EQUATION: FURTHER HEAT FLUX TERMS
C         SPECIES ENTHALPY GRADIENT TERMS

C         SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(STORE4)
          IF(FXRDIF)CALL ZEROXR(STORE4)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(STORE5)
          IF(FYRDIF)CALL ZEROYR(STORE5)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(STORE6)
          IF(FZRDIF)CALL ZEROZR(STORE6)

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(STORE4)
          IF(FXRADB)CALL ZEROXR(STORE4)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(STORE5)
          IF(FYRADB)CALL ZEROYR(STORE5)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(STORE6)
          IF(FZRADB)CALL ZEROZR(STORE6)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE4(IC,JC,KC)*PD1X(IC,JC,KC)
     +                 + STORE5(IC,JC,KC)*PD1Y(IC,JC,KC)
     +                 + STORE6(IC,JC,KC)*PD1Z(IC,JC,KC)
 
              ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*STORE7(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C         =====================================================================

C         Y-EQUATION: DIFFUSIVE TERMS
C         ---------------------------
C         WALL BC: PRESSURE TERMS
C         E-EQUATION: HEAT FLUX TERMS
C         WALL BC: ENTHALPY DIFFUSION TERMS
          IF(FXLDFW)THEN
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = ZERO
                DO IC = ISTAP1,ISTOW

                  FORNOW = FORNOW
     +                + ACBCXL(IC-1)*STORE7(IC,JC,KC)*PD1X(IC,JC,KC)

                ENDDO
                RATE(ISTAL,JC,KC,ISPEC) = RATE(ISTAL,JC,KC,ISPEC)
     +                                  + FORNOW
                VTMP(ISTAL,JC,KC) = VTMP(ISTAL,JC,KC) + FORNOW
                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                            + FORNOW*UTMP(ISTAL,JC,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FXRDFW)THEN
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = ZERO
                DO IC = ISTAW,ISTOM1

                  FORNOW = FORNOW
     +                + ACBCXR(ISTOL-IC)*STORE7(IC,JC,KC)*PD1X(IC,JC,KC)

                ENDDO
                RATE(ISTOL,JC,KC,ISPEC) = RATE(ISTOL,JC,KC,ISPEC)
     +                                  + FORNOW
                VTMP(ISTOL,JC,KC) = VTMP(ISTOL,JC,KC) + FORNOW
                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                            + FORNOW*UTMP(ISTOL,JC,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FYLDFW)THEN
            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO JC = JSTAP1,JSTOW

                  FORNOW = FORNOW
     +                + ACBCYL(JC-1)*STORE7(IC,JC,KC)*PD1Y(IC,JC,KC)

                ENDDO
                RATE(IC,JSTAL,KC,ISPEC) = RATE(IC,JSTAL,KC,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JSTAL,KC) = VTMP(IC,JSTAL,KC) + FORNOW
                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                            + FORNOW*UTMP(IC,JSTAL,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FYRDFW)THEN
            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO JC = JSTAW,JSTOM1

                  FORNOW = FORNOW
     +                + ACBCYR(JSTOL-JC)*STORE7(IC,JC,KC)*PD1Y(IC,JC,KC)

                ENDDO
                RATE(IC,JSTOL,KC,ISPEC) = RATE(IC,JSTOL,KC,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JSTOL,KC) = VTMP(IC,JSTOL,KC) + FORNOW
                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                            + FORNOW*UTMP(IC,JSTOL,KC)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FZLDFW)THEN
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO KC = KSTAP1,KSTOW

                  FORNOW = FORNOW
     +                + ACBCZL(KC-1)*STORE7(IC,JC,KC)*PD1Z(IC,JC,KC)

                ENDDO
                RATE(IC,JC,KSTAL,ISPEC) = RATE(IC,JC,KSTAL,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JC,KSTAL) = VTMP(IC,JC,KSTAL) + FORNOW
                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                            + FORNOW*UTMP(IC,JC,KSTAL)
  
              ENDDO
            ENDDO
          ENDIF
          IF(FZRDFW)THEN
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO KC = KSTAW,KSTOM1

                  FORNOW = FORNOW
     +                + ACBCZR(KSTOL-KC)*STORE7(IC,JC,KC)*PD1Z(IC,JC,KC)

                ENDDO
                RATE(IC,JC,KSTOL,ISPEC) = RATE(IC,JC,KSTOL,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JC,KSTOL) = VTMP(IC,JC,KSTOL) + FORNOW
                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                            + FORNOW*UTMP(IC,JC,KSTOL)
  
              ENDDO
            ENDDO
          ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C         =====================================================================

C         Y-EQUATION: DIFFUSIVE TERMS
C         E-EQUATION: FURTHER HEAT FLUX TERMS
C         SECOND DERIVATIVE TERMS
C         SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(PD2X)
          IF(FXRDIF)CALL ZEROXR(PD2X)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(PD2Y)
          IF(FYRDIF)CALL ZEROYR(PD2Y)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(PD2Z)
          IF(FZRDIF)CALL ZEROZR(PD2Z)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = (PD2X(IC,JC,KC)
     +                 +  PD2Y(IC,JC,KC)
     +                 +  PD2Z(IC,JC,KC))*STORE7(IC,JC,KC)

C               Y EQUATION
                RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C               DIFFUSION CORRECTION VELOCITY DIVERGENCE
                VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

              ENDDO
            ENDDO
          ENDDO

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(PD2X)
          IF(FXRADB)CALL ZEROXR(PD2X)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(PD2Y)
          IF(FYRADB)CALL ZEROYR(PD2Y)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(PD2Z)
          IF(FZRADB)CALL ZEROZR(PD2Z)

C         E EQUATION
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = (PD2X(IC,JC,KC)
     +                 +  PD2Y(IC,JC,KC)
     +                 +  PD2Z(IC,JC,KC))*STORE7(IC,JC,KC)

                ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO

        ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C       =======================================================================

C       SORET EFFECT (THERMAL DIFFUSION) TERMS
C       RSC 16-JUN-2023 BUG FIX THERMAL DIFFUSION RATIO
C       RSC 22-MAR-2026 BUG FIX THERMAL DIFFUSION RATIO
        IF(FLMSOR(ISPEC))THEN

C         FIRST AND SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED

          DO KC = KSTAB,KSTOB
            DO JC = JSTAB,JSTOB
              DO IC = ISTAB,ISTOB

                STORE7(IC,JC,KC) = DIFMIX(IC,JC,KC)*YRHS(IC,JC,KC,ISPEC)
     +                            *TDRMIX(IC,JC,KC)
C     +                          *TDRMIX(IC,JC,KC)*XMOLFR(ISPEC,IC,JC,KC)
  
              ENDDO
            ENDDO
          ENDDO
  
C         DIFFUSION CORRECTION VELOCITY
C         FIRST DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                UCOR(IC,JC,KC) = UCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*TD1X(IC,JC,KC)
                VCOR(IC,JC,KC) = VCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*TD1Y(IC,JC,KC)
                WCOR(IC,JC,KC) = WCOR(IC,JC,KC)
     +                         + STORE7(IC,JC,KC)*TD1Z(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO


C         Y EQUATION: DIFFUSIVE TERMS
C         E EQUATION: FURTHER HEAT FLUX TERMS

C         DIFFUSIVITY GRADIENT TERMS

C         DIFFUSIVITY GRADIENTS
          CALL DFBYDX(STORE7,STORE1)
          CALL DFBYDY(STORE7,STORE2)
          CALL DFBYDZ(STORE7,STORE3)

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(STORE1)
          IF(FXRDIF)CALL ZEROXR(STORE1)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(STORE2)
          IF(FYRDIF)CALL ZEROYR(STORE2)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(STORE3)
          IF(FZRDIF)CALL ZEROZR(STORE3)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE1(IC,JC,KC)*TD1X(IC,JC,KC)
     +                 + STORE2(IC,JC,KC)*TD1Y(IC,JC,KC)
     +                 + STORE3(IC,JC,KC)*TD1Z(IC,JC,KC)
 
C               Y EQUATION
                RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C               DIFFUSION CORRECTION VELOCITY DIVERGENCE
                VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

              ENDDO
            ENDDO
          ENDDO

C         SUBTRACT DUFOUR EFFECT TERMS TO RESTORE SPECIES ENTHALPY
C         RSC 08-JUN-2015 BUG FIX
          IF(FLMDUF(ISPEC))THEN
            DO KC = KSTAB,KSTOB
              DO JC = JSTAB,JSTOB
                DO IC = ISTAB,ISTOB

                  UTMP(IC,JC,KC) = UTMP(IC,JC,KC)
     +             - RGSPEC(ISPEC)*TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)

                ENDDO
              ENDDO
            ENDDO
          ENDIF

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(STORE1)
          IF(FXRADB)CALL ZEROXR(STORE1)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(STORE2)
          IF(FYRADB)CALL ZEROYR(STORE2)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(STORE3)
          IF(FZRADB)CALL ZEROZR(STORE3)

C         E EQUATION
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE1(IC,JC,KC)*TD1X(IC,JC,KC)
     +                 + STORE2(IC,JC,KC)*TD1Y(IC,JC,KC)
     +                 + STORE3(IC,JC,KC)*TD1Z(IC,JC,KC)
 
                ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC) 
              
              ENDDO
            ENDDO
          ENDDO


C         E-EQUATION: FURTHER HEAT FLUX TERMS
C         SPECIES ENTHALPY GRADIENT TERMS

C         EVALUATE SPECIES ENTHALPY GRADIENTS USING STORE4,5,6
          CALL DFBYDX(UTMP,STORE4)
          CALL DFBYDY(UTMP,STORE5)
          CALL DFBYDZ(UTMP,STORE6)

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(STORE4)
          IF(FXRDIF)CALL ZEROXR(STORE4)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(STORE5)
          IF(FYRDIF)CALL ZEROYR(STORE5)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(STORE6)
          IF(FZRDIF)CALL ZEROZR(STORE6)

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(STORE4)
          IF(FXRADB)CALL ZEROXR(STORE4)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(STORE5)
          IF(FYRADB)CALL ZEROYR(STORE5)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(STORE6)
          IF(FZRADB)CALL ZEROZR(STORE6)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = STORE4(IC,JC,KC)*TD1X(IC,JC,KC)
     +                 + STORE5(IC,JC,KC)*TD1Y(IC,JC,KC)
     +                 + STORE6(IC,JC,KC)*TD1Z(IC,JC,KC)
 
              ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*STORE7(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C         =====================================================================

C         Y-EQUATION: DIFFUSIVE TERMS
C         ---------------------------
C         WALL BC: SORET EFFECT TERMS
C         E-EQUATION: HEAT FLUX TERMS
C         WALL BC: ENTHALPY DIFFUSION TERMS
          IF(FXLDFW)THEN
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = ZERO
                DO IC = ISTAP1,ISTOW

                  FORNOW = FORNOW
     +                + ACBCXL(IC-1)*STORE7(IC,JC,KC)*TD1X(IC,JC,KC)

                ENDDO
                RATE(ISTAL,JC,KC,ISPEC) = RATE(ISTAL,JC,KC,ISPEC)
     +                                  + FORNOW
                VTMP(ISTAL,JC,KC) = VTMP(ISTAL,JC,KC) + FORNOW
                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                            + FORNOW*(UTMP(ISTAL,JC,KC)
     +            + RGSPEC(ISPEC)*TRUN(ISTAL,JC,KC)*TDRMIX(ISTAL,JC,KC))
  
              ENDDO
            ENDDO
          ENDIF
          IF(FXRDFW)THEN
            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = ZERO
                DO IC = ISTAW,ISTOM1

                  FORNOW = FORNOW
     +                + ACBCXR(ISTOL-IC)*STORE7(IC,JC,KC)*TD1X(IC,JC,KC)

                ENDDO
                RATE(ISTOL,JC,KC,ISPEC) = RATE(ISTOL,JC,KC,ISPEC)
     +                                  + FORNOW
                VTMP(ISTOL,JC,KC) = VTMP(ISTOL,JC,KC) + FORNOW
                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                            + FORNOW*(UTMP(ISTOL,JC,KC)
     +            + RGSPEC(ISPEC)*TRUN(ISTOL,JC,KC)*TDRMIX(ISTOL,JC,KC))
  
              ENDDO
            ENDDO
          ENDIF
          IF(FYLDFW)THEN
            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO JC = JSTAP1,JSTOW

                  FORNOW = FORNOW
     +                + ACBCYL(JC-1)*STORE7(IC,JC,KC)*TD1Y(IC,JC,KC)

                ENDDO
                RATE(IC,JSTAL,KC,ISPEC) = RATE(IC,JSTAL,KC,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JSTAL,KC) = VTMP(IC,JSTAL,KC) + FORNOW
                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                            + FORNOW*(UTMP(IC,JSTAL,KC)
     +            + RGSPEC(ISPEC)*TRUN(IC,JSTAL,KC)*TDRMIX(IC,JSTAL,KC))
  
              ENDDO
            ENDDO
          ENDIF
          IF(FYRDFW)THEN
            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO JC = JSTAW,JSTOM1

                  FORNOW = FORNOW
     +                + ACBCYR(JSTOL-JC)*STORE7(IC,JC,KC)*TD1Y(IC,JC,KC)

                ENDDO
                RATE(IC,JSTOL,KC,ISPEC) = RATE(IC,JSTOL,KC,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JSTOL,KC) = VTMP(IC,JSTOL,KC) + FORNOW
                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                            + FORNOW*(UTMP(IC,JSTOL,KC)
     +            + RGSPEC(ISPEC)*TRUN(IC,JSTOL,KC)*TDRMIX(IC,JSTOL,KC))
  
              ENDDO
            ENDDO
          ENDIF
          IF(FZLDFW)THEN
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO KC = KSTAP1,KSTOW

                  FORNOW = FORNOW
     +                + ACBCZL(KC-1)*STORE7(IC,JC,KC)*TD1Z(IC,JC,KC)

                ENDDO
                RATE(IC,JC,KSTAL,ISPEC) = RATE(IC,JC,KSTAL,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JC,KSTAL) = VTMP(IC,JC,KSTAL) + FORNOW
                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                            + FORNOW*(UTMP(IC,JC,KSTAL)
     +            + RGSPEC(ISPEC)*TRUN(IC,JC,KSTAL)*TDRMIX(IC,JC,KSTAL))
  
              ENDDO
            ENDDO
          ENDIF
          IF(FZRDFW)THEN
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = ZERO
                DO KC = KSTAW,KSTOM1

                  FORNOW = FORNOW
     +                + ACBCZR(KSTOL-KC)*STORE7(IC,JC,KC)*TD1Z(IC,JC,KC)

                ENDDO
                RATE(IC,JC,KSTOL,ISPEC) = RATE(IC,JC,KSTOL,ISPEC)
     +                                  + FORNOW
                VTMP(IC,JC,KSTOL) = VTMP(IC,JC,KSTOL) + FORNOW
                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                            + FORNOW*(UTMP(IC,JC,KSTOL)
     +            + RGSPEC(ISPEC)*TRUN(IC,JC,KSTOL)*TDRMIX(IC,JC,KSTOL))
  
              ENDDO
            ENDDO
          ENDIF

C         E-EQUATION: HEAT FLUX TERMS
C         WALL BC: SORET AND DUFOUR TERMS
          IF(FLMDUF(ISPEC))THEN
C           E-EQUATION: HEAT FLUX TERMS
C           WALL BC: ADIABATIC WALL
            IF(FXLCNW)THEN
              DO KC = KSTAL,KSTOL
                DO JC = JSTAL,JSTOL

                  FORNOW = ZERO
                  DO IC = ISTAP1,ISTOW

                  COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                  FORNOW = FORNOW
     +             + ACBCXL(IC-1)*COMBO1*STORE7(IC,JC,KC)*TD1X(IC,JC,KC)

                  ENDDO
                  ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FXRCNW)THEN
              DO KC = KSTAL,KSTOL
                DO JC = JSTAL,JSTOL

                  FORNOW = ZERO
                  DO IC = ISTAW,ISTOM1

                  COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                  FORNOW = FORNOW
     +        + ACBCXR(ISTOL-IC)*COMBO1*STORE7(IC,JC,KC)*TD1X(IC,JC,KC)

                  ENDDO
                  ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FYLCNW)THEN
              DO KC = KSTAL,KSTOL
                DO IC = ISTAL,ISTOL

                  FORNOW = ZERO
                  DO JC = JSTAP1,JSTOW

                  COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                  FORNOW = FORNOW
     +             + ACBCYL(JC-1)*COMBO1*STORE7(IC,JC,KC)*TD1Y(IC,JC,KC)

                  ENDDO
                  ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FYRCNW)THEN
              DO KC = KSTAL,KSTOL
                DO IC = ISTAL,ISTOL

                  FORNOW = ZERO
                  DO JC = JSTAW,JSTOM1

                  COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                  FORNOW = FORNOW
     +        + ACBCYR(JSTOL-JC)*COMBO1*STORE7(IC,JC,KC)*TD1Y(IC,JC,KC)

                  ENDDO
                  ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FZLCNW)THEN
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  FORNOW = ZERO
                  DO KC = KSTAP1,KSTOW

                  COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                  FORNOW = FORNOW
     +             + ACBCZL(KC-1)*COMBO1*STORE7(IC,JC,KC)*TD1Z(IC,JC,KC)

                  ENDDO
                  ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FZRCNW)THEN
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  FORNOW = ZERO
                  DO KC = KSTAW,KSTOM1

                  COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                  FORNOW = FORNOW
     +        + ACBCZR(KSTOL-KC)*COMBO1*STORE7(IC,JC,KC)*TD1Z(IC,JC,KC)

                  ENDDO
                  ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF

C           E-EQUATION: HEAT FLUX TERMS
C           WALL BC: ISOTHERMAL WALL
            IF(FXLADW)THEN
              DO KC = KSTAL,KSTOL
                DO JC = JSTAL,JSTOL

                  COMBO2 = TRUN(ISTAL,JC,KC)*TDRMIX(ISTAL,JC,KC)
                  COMBO2 = COMBO2*STORE7(ISTAL,JC,KC)*TD1X(ISTAL,JC,KC)
                  FORNOW = ZERO
                  DO IC = ISTAP1,ISTOW

                    COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                    COMBO1 = COMBO1*STORE7(IC,JC,KC)*TD1X(IC,JC,KC)
                    FORNOW = FORNOW
     +                     + ACBCXL(IC-1)*RGSPEC(ISPEC)*(COMBO1-COMBO2)

                  ENDDO
                  ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FXRADW)THEN
              DO KC = KSTAL,KSTOL
                DO JC = JSTAL,JSTOL

                  COMBO2 = TRUN(ISTOL,JC,KC)*TDRMIX(ISTOL,JC,KC)
                  COMBO2 = COMBO2*STORE7(ISTOL,JC,KC)*TD1X(ISTOL,JC,KC)
                  FORNOW = ZERO
                  DO IC = ISTAW,ISTOM1

                    COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                    COMBO1 = COMBO1*STORE7(IC,JC,KC)*TD1X(IC,JC,KC)
                    FORNOW = FORNOW
     +                 + ACBCXR(ISTOL-IC)*RGSPEC(ISPEC)*(COMBO2-COMBO1)

                  ENDDO
                  ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FYLADW)THEN
              DO KC = KSTAL,KSTOL
                DO IC = ISTAL,ISTOL

                  COMBO2 = TRUN(IC,JSTAL,KC)*TDRMIX(IC,JSTAL,KC)
                  COMBO2 = COMBO2*STORE7(IC,JSTAL,KC)*TD1Y(IC,JSTAL,KC)
                  FORNOW = ZERO
                  DO JC = JSTAP1,JSTOW

                    COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                    COMBO1 = COMBO1*STORE7(IC,JC,KC)*TD1Y(IC,JC,KC)
                    FORNOW = FORNOW
     +                     + ACBCYL(JC-1)*RGSPEC(ISPEC)*(COMBO1-COMBO2)

                  ENDDO
                  ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FYRADW)THEN
              DO KC = KSTAL,KSTOL
                DO IC = ISTAL,ISTOL

                  COMBO2 = TRUN(IC,JSTOL,KC)*TDRMIX(IC,JSTOL,KC)
                  COMBO2 = COMBO2*STORE7(IC,JSTOL,KC)*TD1Y(IC,JSTOL,KC)
                  FORNOW = ZERO
                  DO JC = JSTAW,JSTOM1

                    COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                    COMBO1 = COMBO1*STORE7(IC,JC,KC)*TD1Y(IC,JC,KC)
                    FORNOW = FORNOW
     +                 + ACBCYR(JSTOL-JC)*RGSPEC(ISPEC)*(COMBO2-COMBO1)

                  ENDDO
                  ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FZLADW)THEN
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  COMBO2 = TRUN(IC,JC,KSTAL)*TDRMIX(IC,JC,KSTAL)
                  COMBO2 = COMBO2*STORE7(IC,JC,KSTAL)*TD1Z(IC,JC,KSTAL)
                  FORNOW = ZERO
                  DO KC = KSTAP1,KSTOW

                    COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                    COMBO1 = COMBO1*STORE7(IC,JC,KC)*TD1Z(IC,JC,KC)
                    FORNOW = FORNOW
     +                     + ACBCZL(KC-1)*RGSPEC(ISPEC)*(COMBO1-COMBO2)

                  ENDDO
                  ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
            IF(FZRADW)THEN
              DO JC = JSTAL,JSTOL
                DO IC = ISTAL,ISTOL

                  COMBO2 = TRUN(IC,JC,KSTOL)*TDRMIX(IC,JC,KSTOL)
                  COMBO2 = COMBO2*STORE7(IC,JC,KSTOL)*TD1Z(IC,JC,KSTOL)
                  FORNOW = ZERO
                  DO KC = KSTAW,KSTOM1

                    COMBO1 = TRUN(IC,JC,KC)*TDRMIX(IC,JC,KC)
                    COMBO1 = COMBO1*STORE7(IC,JC,KC)*TD1Z(IC,JC,KC)
                    FORNOW = FORNOW
     +                 + ACBCZR(KSTOL-KC)*RGSPEC(ISPEC)*(COMBO2-COMBO1)

                  ENDDO
                  ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                              + RGSPEC(ISPEC)*FORNOW
  
                ENDDO
              ENDDO
            ENDIF
          ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              UTMP = SPECIES H
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              STORE7 = RHO D Y
C         =====================================================================

C         Y-EQUATION: DIFFUSIVE TERMS
C         E-EQUATION: FURTHER HEAT FLUX TERMS
C         SECOND DERIVATIVE TERMS
C         SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FXLDIF)CALL ZEROXL(TD2X)
          IF(FXRDIF)CALL ZEROXR(TD2X)
C         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FYLDIF)CALL ZEROYL(TD2Y)
          IF(FYRDIF)CALL ZEROYR(TD2Y)
C         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
          IF(FZLDIF)CALL ZEROZL(TD2Z)
          IF(FZRDIF)CALL ZEROZR(TD2Z)

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = (TD2X(IC,JC,KC)
     +                 +  TD2Y(IC,JC,KC)
     +                 +  TD2Z(IC,JC,KC))*STORE7(IC,JC,KC)

C               Y EQUATION
                RATE(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC) + FORNOW

C               DIFFUSION CORRECTION VELOCITY DIVERGENCE
                VTMP(IC,JC,KC) = VTMP(IC,JC,KC) + FORNOW

              ENDDO
            ENDDO
          ENDDO

C         BOUNDARY CONDITIONS
C         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FXLADB)CALL ZEROXL(TD2X)
          IF(FXRADB)CALL ZEROXR(TD2X)
C         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FYLADB)CALL ZEROYL(TD2Y)
          IF(FYRADB)CALL ZEROYR(TD2Y)
C         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
          IF(FZLADB)CALL ZEROZL(TD2Z)
          IF(FZRADB)CALL ZEROZR(TD2Z)

C         E EQUATION
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = (TD2X(IC,JC,KC)
     +                 +  TD2Y(IC,JC,KC)
     +                 +  TD2Z(IC,JC,KC))*STORE7(IC,JC,KC)

                ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC)

              ENDDO
            ENDDO
          ENDDO

        ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                           VTMP = DIV CORR VEL
C                                                              WTMP = MIXTURE H
C                                                              ALL STORES CLEAR
C       =======================================================================

C       ---------------------------------------------------------------- 
C       E-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
C       Y-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
C       ---------------------------------------------------------------- 

      ENDDO
C     RSC 08-AUG-2012 EVALUATE ALL SPECIES
C     END OF RUN THROUGH ALL SPECIES

C     =========================================================================

C     EVALUATE DIFFUSION CORRECTION VELOCITY TERMS
C     --------------------------------------------

C     E-EQUATION: FURTHER HEAT FLUX TERMS
C     -----------------------------------
C     MIXTURE ENTHALPY GRADIENTS
      CALL DFBYDX(WTMP,STORE1)
      CALL DFBYDY(WTMP,STORE2)
      CALL DFBYDZ(WTMP,STORE3)

C     BOUNDARY CONDITIONS
C     BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF(FXLDIF)CALL ZEROXL(STORE1)
      IF(FXRDIF)CALL ZEROXR(STORE1)
C     BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF(FYLDIF)CALL ZEROYL(STORE2)
      IF(FYRDIF)CALL ZEROYR(STORE2)
C     BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF(FZLDIF)CALL ZEROZL(STORE3)
      IF(FZRDIF)CALL ZEROZR(STORE3)

C     BOUNDARY CONDITIONS
C     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FXLADB)CALL ZEROXL(STORE1)
      IF(FXRADB)CALL ZEROXR(STORE1)
C     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FYLADB)CALL ZEROYL(STORE2)
      IF(FYRADB)CALL ZEROYR(STORE2)
C     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FZLADB)CALL ZEROZL(STORE3)
      IF(FZRADB)CALL ZEROZR(STORE3)

C     TRANSFER DIV CORR VEL TO TEMPORARY STORE
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE4(IC,JC,KC) = VTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FXLADB)CALL ZEROXL(STORE4)
      IF(FXRADB)CALL ZEROXR(STORE4)
C     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FYLADB)CALL ZEROYL(STORE4)
      IF(FYRADB)CALL ZEROYR(STORE4)
C     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF(FZLADB)CALL ZEROZL(STORE4)
      IF(FZRADB)CALL ZEROZR(STORE4)

C     DIV RHO VCORR HMIX
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - WTMP(IC,JC,KC)*STORE4(IC,JC,KC)
     +                     - STORE1(IC,JC,KC)*UCOR(IC,JC,KC)
     +                     - STORE2(IC,JC,KC)*VCOR(IC,JC,KC)
     +                     - STORE3(IC,JC,KC)*WCOR(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              ALL STORES CLEAR
C     =========================================================================

C     MIXTURE AVERAGED TRANSPORT
C     EVALUATE THE VISCOSITY
C     RSC 17-APR-2013
C     RSC 01-APR-2022
      IF(FLMAVT)THEN

C       MIXTURE AVERAGED TRANSPORT
C       VISCOSITY IS PARALLEL
C       STORE VISCOSITY IN DIFMIX FOR NOW
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

C             VISCOSITY FOR EACH SPECIES
              TRATIO = TRUN(IC,JC,KC)*TDIFGB
              DO JSPEC = 1, NSPEC
                FORNOW = VISCCO(NCOVIS,JSPEC)
                DO ICP = NCOVM1,1,-1
                  FORNOW = FORNOW*TRATIO + VISCCO(ICP,JSPEC)
                ENDDO
                CTRANS(JSPEC) = FORNOW
              ENDDO

C             COMBINATION RULE FOR VISCOSITY
              COMBO1 = ZERO
              DO ISPEC = 1, NSPEC
                COMBO2 = ZERO
                DO JSPEC = 1, NSPEC
                  FORNOW = SQRT(CTRANS(ISPEC)/CTRANS(JSPEC))
                  FORNOW = ONE + FORNOW*WILKO2(JSPEC,ISPEC)
                  FORNOW = WILKO1(JSPEC,ISPEC)*FORNOW*FORNOW
                  COMBO2 = COMBO2 + XMOLFR(JSPEC,IC,JC,KC)*FORNOW
                ENDDO
                FORNOW = CTRANS(ISPEC)/COMBO2
                COMBO1 = COMBO1 + XMOLFR(ISPEC,IC,JC,KC)*FORNOW

              ENDDO
              DIFMIX(IC,JC,KC) = COMBO1

            ENDDO
          ENDDO
        ENDDO

      ENDIF
C                                                          TRANSP = MIX COND/CP
C                                                         RATE = Y SOURCE TERMS
C                                                              ALL STORES CLEAR
C     =========================================================================

C     RUN THROUGH ALL SPECIES
C     -----------------------
C     RSC 08-AUG-2012 EVALUATE ALL SPECIES
C     RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
      DO ISPEC = 1,NSPEC

C       Y-EQUATION: DIFFUSIVE TERMS
C       ---------------------------
C       RECOMPUTE SPECIES MASS FRACTION GRADIENTS
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              STORE7(IC,JC,KC) = YRHS(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
        CALL DFBYDX(STORE7,STORE1)
        CALL DFBYDY(STORE7,STORE2)
        CALL DFBYDZ(STORE7,STORE3)

C       BOUNDARY CONDITIONS
C       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FXLDIF)CALL ZEROXL(STORE1)
        IF(FXRDIF)CALL ZEROXR(STORE1)
C       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FYLDIF)CALL ZEROYL(STORE2)
        IF(FYRDIF)CALL ZEROYR(STORE2)
C       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(FZLDIF)CALL ZEROZL(STORE3)
        IF(FZRDIF)CALL ZEROZR(STORE3)

C       DIV RHO VCORR Y
C       STORE Y SOURCE TERMS IN YRHS
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              YRHS(IC,JC,KC,ISPEC) = RATE(IC,JC,KC,ISPEC)
     +                             - YRHS(IC,JC,KC,ISPEC)*VTMP(IC,JC,KC)
     +                             - STORE1(IC,JC,KC)*UCOR(IC,JC,KC)
     +                             - STORE2(IC,JC,KC)*VCOR(IC,JC,KC)
     +                             - STORE3(IC,JC,KC)*WCOR(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     RSC 08-AUG-2012 EVALUATE ALL SPECIES
C     END OF RUN THROUGH ALL SPECIES
C                                                          TRANSP = MIX COND/CP
C                                                              ALL STORES CLEAR
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C     ------------------------------------------------ 
C     Y-EQUATION: SOURCE TERMS COMPLETE
C     ------------------------------------------------ 
C     E-EQUATION: PRESSURE-WORK AND VISCOUS WORK TERMS
C                 EVALUATED IN SUBROUTINE RHSVEL
C     ------------------------------------------------

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     COLLECT DENSITY AND ITS GRADIENTS FOR BCs
C     -----------------------------------------
C     RSC 18-SEP-2022 BC TRANSVERSE TERMS

CC     X-DIRECTION: DRHODX
C      IF(FXLCNV.OR.FXRCNV)THEN
      IF(FXLCNV.OR.FXRCNV.OR.FYLCNV.OR.FYRCNV.OR.FZLCNV.OR.FZRCNV)THEN

        CALL DFBYDX(DRHS,STORE1)
        CALL DFBYDY(DRHS,STORE2)
        CALL DFBYDZ(DRHS,STORE3)

        IF(FXLCNV)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRDXL(JC,KC) = DRHS(ISTAL,JC,KC)
              BCDXXL(JC,KC) = STORE1(ISTAL,JC,KC)
              BCDYXL(JC,KC) = STORE2(ISTAL,JC,KC)
              BCDZXL(JC,KC) = STORE3(ISTAL,JC,KC)

            ENDDO
          ENDDO
        ENDIF
        IF(FXRCNV)THEN
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRDXR(JC,KC) = DRHS(ISTOL,JC,KC)
              BCDXXR(JC,KC) = STORE1(ISTOL,JC,KC)
              BCDYXR(JC,KC) = STORE2(ISTOL,JC,KC)
              BCDZXR(JC,KC) = STORE3(ISTOL,JC,KC)

            ENDDO
          ENDDO
        ENDIF

C      ENDIF

CC     Y-DIRECTION: DRHODY
C      IF(FYLCNV.OR.FYRCNV)THEN
C
C        CALL DFBYDX(DRHS,STORE1)
C        CALL DFBYDY(DRHS,STORE2)
C        CALL DFBYDZ(DRHS,STORE3)

        IF(FYLCNV)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRDYL(IC,KC) = DRHS(IC,JSTAL,KC)
              BCDXYL(IC,KC) = STORE1(IC,JSTAL,KC)
              BCDYYL(IC,KC) = STORE2(IC,JSTAL,KC)
              BCDZYL(IC,KC) = STORE3(IC,JSTAL,KC)

            ENDDO
          ENDDO
        ENDIF
        IF(FYRCNV)THEN
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRDYR(IC,KC) = DRHS(IC,JSTOL,KC)
              BCDXYR(IC,KC) = STORE1(IC,JSTOL,KC)
              BCDYYR(IC,KC) = STORE2(IC,JSTOL,KC)
              BCDZYR(IC,KC) = STORE3(IC,JSTOL,KC)

            ENDDO
          ENDDO
        ENDIF

c      ENDIF

cC     Z-DIRECTION: DRHODZ
c      IF(FZLCNV.OR.FZRCNV)THEN
c
c        CALL DFBYDX(DRHS,STORE1)
c        CALL DFBYDY(DRHS,STORE2)
c        CALL DFBYDZ(DRHS,STORE3)

        IF(FZLCNV)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRDZL(IC,JC) = DRHS(IC,JC,KSTAL)
              BCDXZL(IC,JC) = STORE1(IC,JC,KSTAL)
              BCDYZL(IC,JC) = STORE2(IC,JC,KSTAL)
              BCDZZL(IC,JC) = STORE3(IC,JC,KSTAL)

            ENDDO
          ENDDO
        ENDIF
        IF(FZRCNV)THEN
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRDZR(IC,JC) = DRHS(IC,JC,KSTOL)
              BCDXZR(IC,JC) = STORE1(IC,JC,KSTOL)
              BCDYZR(IC,JC) = STORE2(IC,JC,KSTOL)
              BCDZZR(IC,JC) = STORE3(IC,JC,KSTOL)

            ENDDO
          ENDDO
        ENDIF

      ENDIF
C                                                          TRANSP = MIX COND/CP
C     =========================================================================


      RETURN
      END
