      SUBROUTINE RHSVEL 
 
C     *************************************************************************
C
C     RHSVEL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     08-DEC-2002:  CREATED
C     17-APR-2013:  RSC MIXTURE AVERAGED TRANSPORT
C     
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES RIGHT-HAND-SIDES FOR TIME INTEGRATION
C     OF CONTINUITY AND MOMENTUM EQUATIONS
C     EVALUATES PRESSURE WORK AND VISCOUS WORK TERMS IN ENERGY EQUATION
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FORNOW,PREFER
      INTEGER IC,JC,KC


C     BEGIN
C     =====

C     =========================================================================

C     CONVERT VELOCITIES
C     ------------------

C     U,V,WRHS CONTAIN RHO U,V,W: CONVERT TO U,V,W
C     U,V,W HELD IN U,V,WTMP THROUGHOUT THIS ROUTINE
C     U,V,W ARE PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            UTMP(IC,JC,KC) = URHS(IC,JC,KC)/DRHS(IC,JC,KC)
            VTMP(IC,JC,KC) = VRHS(IC,JC,KC)/DRHS(IC,JC,KC)
            WTMP(IC,JC,KC) = WRHS(IC,JC,KC)/DRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)


      CALL DFBYDY(UTMP,STORE1)
      CALL DFBYDZ(UTMP,STORE2)     

C     =========================================================================

C     COLLECT VELOCITY COMPONENTS FOR BCs
C     -----------------------------------

C     X-DIRECTION
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXL(JC,KC) = UTMP(ISTAL,JC,KC)
            STRVXL(JC,KC) = VTMP(ISTAL,JC,KC)
            STRWXL(JC,KC) = WTMP(ISTAL,JC,KC)
            
C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T2BXL(JC,KC)=-VTMP(ISTAL,JC,KC)*STORE1(ISTAL,JC,KC)
     +                       -WTMP(ISTAL,JC,KC)*STORE2(ISTAL,JC,KC)            
          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXR(JC,KC) = UTMP(ISTOL,JC,KC)
            STRVXR(JC,KC) = VTMP(ISTOL,JC,KC)
            STRWXR(JC,KC) = WTMP(ISTOL,JC,KC)

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T2BXR(JC,KC)=-VTMP(ISTOL,JC,KC)*STORE1(ISTOL,JC,KC)
     +                       -WTMP(ISTOL,JC,KC)*STORE2(ISTOL,JC,KC)
          ENDDO
        ENDDO
      ENDIF

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)


      CALL DFBYDX(VTMP,STORE1)
      CALL DFBYDZ(VTMP,STORE2)


C     Y-DIRECTION
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRUYL(IC,KC) = UTMP(IC,JSTAL,KC)
            STRVYL(IC,KC) = VTMP(IC,JSTAL,KC)
            STRWYL(IC,KC) = WTMP(IC,JSTAL,KC)

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T2BYL(IC,KC)=-UTMP(IC,JSTAL,KC)*STORE1(IC,JSTAL,KC)
     +                       -WTMP(IC,JSTAL,KC)*STORE2(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRUYR(IC,KC) = UTMP(IC,JSTOL,KC)
            STRVYR(IC,KC) = VTMP(IC,JSTOL,KC)
            STRWYR(IC,KC) = WTMP(IC,JSTOL,KC)

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T2BYR(IC,KC)=-UTMP(IC,JSTOL,KC)*STORE1(IC,JSTOL,KC)
     +                       -WTMP(IC,JSTOL,KC)*STORE2(IC,JSTOL,KC)
          ENDDO
        ENDDO
      ENDIF

      CALL DFBYDX(WTMP,STORE1)
      CALL DFBYDY(WTMP,STORE2)

C     Z-DIRECTION
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRUZL(IC,JC) = UTMP(IC,JC,KSTAL)
            STRVZL(IC,JC) = VTMP(IC,JC,KSTAL)
            STRWZL(IC,JC) = WTMP(IC,JC,KSTAL)

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T2BZL(IC,JC)=-UTMP(IC,JC,KSTAL)*STORE1(IC,JC,KSTAL)
     +                       -VTMP(IC,JC,KSTAL)*STORE2(IC,JC,KSTAL)
          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRUZR(IC,JC) = UTMP(IC,JC,KSTOL)
            STRVZR(IC,JC) = VTMP(IC,JC,KSTOL)
            STRWZR(IC,JC) = WTMP(IC,JC,KSTOL)

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T2BZR(IC,JC)=-UTMP(IC,JC,KSTOL)*STORE1(IC,JC,KSTOL)
     +                       -VTMP(IC,JC,KSTOL)*STORE2(IC,JC,KSTOL)
          ENDDO
        ENDDO
      ENDIF

      CALL DFBYDY(VTMP,STORE1)
      CALL DFBYDZ(VTMP,STORE2)

C     X-DIRECTION
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
           T3BXL(JC,KC)=-VTMP(ISTAL,JC,KC)*STORE1(ISTAL,JC,KC)
     +                  -WTMP(ISTAL,JC,KC)*STORE2(ISTAL,JC,KC)


          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T3BXR(JC,KC)=-VTMP(ISTOL,JC,KC)*STORE1(ISTOL,JC,KC)
     +                  -WTMP(ISTOL,JC,KC)*STORE2(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF

      CALL DFBYDX(WTMP,STORE1)
      CALL DFBYDZ(WTMP,STORE2)

C     Y-DIRECTION
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
          !NC and VM changed here: ycorrection
            T4BYL(IC,KC)=-UTMP(IC,JSTAL,KC)*STORE1(IC,JSTAL,KC)
     +                       -WTMP(IC,JSTAL,KC)*STORE2(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
          !NC and VM changed here: ycorrection
            T4BYR(IC,KC)=-UTMP(IC,JSTOL,KC)*STORE1(IC,JSTOL,KC)
     +                       -WTMP(IC,JSTOL,KC)*STORE2(IC,JSTOL,KC)
          ENDDO
        ENDDO
      ENDIF


      CALL DFBYDX(UTMP,STORE1)
      CALL DFBYDY(UTMP,STORE2)

C     Z-DIRECTION
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T3BZL(IC,JC)=-UTMP(IC,JC,KSTAL)*STORE1(IC,JC,KSTAL)
     +                       -VTMP(IC,JC,KSTAL)*STORE2(IC,JC,KSTAL)
          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T3BZR(IC,JC)=-UTMP(IC,JC,KSTOL)*STORE1(IC,JC,KSTOL)
     +                       -VTMP(IC,JC,KSTOL)*STORE2(IC,JC,KSTOL)
          ENDDO
        ENDDO
      ENDIF


      CALL DFBYDY(WTMP,STORE1)
      CALL DFBYDZ(WTMP,STORE2)

C     X-DIRECTION
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
           T4BXL(JC,KC)=-VTMP(ISTAL,JC,KC)*STORE1(ISTAL,JC,KC)
     +                  -WTMP(ISTAL,JC,KC)*STORE2(ISTAL,JC,KC)


          ENDDO
        ENDDO
      ENDIF

      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T4BXR(JC,KC)=-VTMP(ISTOL,JC,KC)*STORE1(ISTOL,JC,KC)
     +                  -WTMP(ISTOL,JC,KC)*STORE2(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF

      CALL DFBYDX(UTMP,STORE1)
      CALL DFBYDZ(UTMP,STORE2)

C     Y-DIRECTION
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
          !NC and VM changed here: ycorrection
            T3BYL(IC,KC)=-UTMP(IC,JSTAL,KC)*STORE1(IC,JSTAL,KC)
     +                       -WTMP(IC,JSTAL,KC)*STORE2(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
          !NC and VM changed here: ycorrection
            T3BYR(IC,KC)=-UTMP(IC,JSTOL,KC)*STORE1(IC,JSTOL,KC)
     +                       -WTMP(IC,JSTOL,KC)*STORE2(IC,JSTOL,KC)
          ENDDO
        ENDDO
      ENDIF

      CALL DFBYDX(VTMP,STORE1)
      CALL DFBYDY(VTMP,STORE2)

C     Z-DIRECTION
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T4BZL(IC,JC)=-UTMP(IC,JC,KSTAL)*STORE1(IC,JC,KSTAL)
     +                       -VTMP(IC,JC,KSTAL)*STORE2(IC,JC,KSTAL)
          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

C    EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
            T4BZR(IC,JC)=-UTMP(IC,JC,KSTOL)*STORE1(IC,JC,KSTOL)
     +                       -VTMP(IC,JC,KSTOL)*STORE2(IC,JC,KSTOL)
          ENDDO
        ENDDO
      ENDIF

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     MOMENTUM EQUATIONS: CONVECTIVE TERMS
C     ------------------------------------
C     RHO U U
C     RHO U U IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = URHS(IC,JC,KC)*UTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D/DX RHO U U
C     STRAIGHT INTO STORE4 FOR NOW
      CALL DFBYDX(STORE7,STORE4)


C     U-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     (HALF) D/DX RHO U U: ALREADY IN STORE4
C                                                   STORE4 = U CONVECTIVE TERMS
C                                                          U,V,WRHS = RHO U,V,W
C     =========================================================================

C     RHO U V
C     RHO U V IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = URHS(IC,JC,KC)*VTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D/DY RHO V U
C     D/DX RHO U V
      CALL DFBYDY(STORE7,STORE1)
      CALL DFBYDX(STORE7,STORE5)


C     U-EQUATION: CONVECTIVE TERMS
C     V-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     (HALF) D/DY RHO V U
C     (HALF) D/DX RHO U V: ALREADY IN STORE5
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE4(IC,JC,KC) = STORE4(IC,JC,KC) + STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                               STORE4,5 = U,V CONVECTIVE TERMS
C                                                          U,V,WRHS = RHO U,V,W
C     =========================================================================

C     RHO U W
C     RHO U W IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = URHS(IC,JC,KC)*WTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D/DZ RHO W U
C     D/DX RHO U W
      CALL DFBYDZ(STORE7,STORE1)
      CALL DFBYDX(STORE7,STORE6)


C     U-EQUATION: CONVECTIVE TERMS
C     W-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     (HALF) D/DZ RHO W U
C     (HALF) D/DX RHO U W: ALREADY IN STORE6
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE4(IC,JC,KC) = STORE4(IC,JC,KC) + STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
C                                                          U,V,WRHS = RHO U,V,W
C     =========================================================================

C     RHO V V
C     RHO V V IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = VRHS(IC,JC,KC)*VTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D/DY RHO V V
      CALL DFBYDY(STORE7,STORE1)


C     V-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     (HALF) D/DY RHO V V
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE5(IC,JC,KC) = STORE5(IC,JC,KC) + STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
C                                                          U,V,WRHS = RHO U,V,W
C     =========================================================================

C     RHO V W
C     RHO V W IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = VRHS(IC,JC,KC)*WTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D/DZ RHO W V
C     D/DY RHO V W
      CALL DFBYDZ(STORE7,STORE1)
      CALL DFBYDY(STORE7,STORE2)


C     V-EQUATION: CONVECTIVE TERMS
C     W-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     (HALF) D/DZ RHO W V
C     (HALF) D/DY RHO V W
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE5(IC,JC,KC) = STORE5(IC,JC,KC) + STORE1(IC,JC,KC)
            STORE6(IC,JC,KC) = STORE6(IC,JC,KC) + STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
C                                                          U,V,WRHS = RHO U,V,W
C     =========================================================================

C     RHO W W
C     RHO W W IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = WRHS(IC,JC,KC)*WTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D/DZ RHO W W
      CALL DFBYDZ(STORE7,STORE1)


C     W-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     (HALF) D/DZ RHO W W
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = STORE6(IC,JC,KC) + STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
C                                                          U,V,WRHS = RHO U,V,W
C     =========================================================================

C     VELOCITY NORMAL DERIVATIVES
C     ---------------------------
C     DUDX,DVDY,DWDZ
      CALL DFBYDX(UTMP,STORE1)
      CALL DFBYDY(VTMP,STORE2)
      CALL DFBYDZ(WTMP,STORE3)
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                           STORE4,5,6 = U,V,W CONVECTIVE TERMS
C                                                          U,V,WRHS = RHO U,V,W
C     =========================================================================

C     COLLECT VELOCITY DERIVATIVES FOR BCs
C     ------------------------------------

C     X-DIRECTION: DUDX
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            BCL1XL(JC,KC) = STORE1(ISTAL,JC,KC)            

          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            BCL1XR(JC,KC) = STORE1(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF

C     Y-DIRECTION: DVDY
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            BCL1YL(IC,KC) = STORE2(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            BCL1YR(IC,KC) = STORE2(IC,JSTOL,KC)

          ENDDO
        ENDDO
      ENDIF

C     Z-DIRECTION: DWDZ
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            BCL1ZL(IC,JC) = STORE3(IC,JC,KSTAL)

          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            BCL1ZR(IC,JC) = STORE3(IC,JC,KSTOL)

          ENDDO
        ENDDO
      ENDIF

C     =========================================================================

C     U-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     DIV RHO U U  + U DIV RHO U
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE7(IC,JC,KC) = STORE4(IC,JC,KC)
     +                       + UTMP(IC,JC,KC)*DIVM(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     RHO U DUDX
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE4(IC,JC,KC) = URHS(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     HALF DIV RHO U U + HALF RHO U DUDX + HALF U DIV RHO U
C     STORE IN URHS
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = -HALF*(STORE4(IC,JC,KC)
     +                            + STORE7(IC,JC,KC))

          ENDDO
        ENDDO
      ENDDO


C     V-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     DIV RHO U V + V DIV RHO U
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE7(IC,JC,KC) = STORE5(IC,JC,KC)
     +                       + VTMP(IC,JC,KC)*DIVM(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     RHO V DVDY
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE5(IC,JC,KC) = VRHS(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     HALF DIV RHO U V + HALF RHO V DVDY + HALF V DIV RHO U
C     STORE IN VRHS
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = -HALF*(STORE5(IC,JC,KC)
     +                            + STORE7(IC,JC,KC))

          ENDDO
        ENDDO
      ENDDO

  
C     W-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     DIV RHO U W + W DIV RHO U
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE7(IC,JC,KC) = STORE6(IC,JC,KC)
     +                       + WTMP(IC,JC,KC)*DIVM(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     RHO W DWDZ
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = WRHS(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     HALF DIV RHO U W + HALF RHO W DWDZ + HALF W DIV RHO U
C     STORE IN WRHS
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = -HALF*(STORE6(IC,JC,KC)
     +                            + STORE7(IC,JC,KC))

          ENDDO
        ENDDO
      ENDDO
C     U,V,WRHS CONTAIN U,V,W SOURCE TERMS THROUGHOUT REMAINDER OF THIS ROUTINE
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C     =========================================================================

C     PRESSURE GRADIENTS
C     ------------------

      PREFER = PRUN(IPREF,JPREF,KPREF)
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            STORE7(IC,JC,KC) = PRUN(IC,JC,KC) - PREFER

          ENDDO
        ENDDO
      ENDDO

C     DPDX,DPDY,DPDZ
      CALL DFBYDX(STORE7,STORE4)
      CALL DFBYDY(STORE7,STORE5)
      CALL DFBYDZ(STORE7,STORE6)
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                                   STORE4,5,6 = DPDX,DPDY,DPDZ
C     =========================================================================

C     COLLECT PRESSURE AND ITS GRADIENTS FOR BCs
C     ------------------------------------------

C     X-DIRECTION: P AND DPDX
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRPXL(JC,KC) = PRUN(ISTAL,JC,KC)
            BCL5XL(JC,KC) = STORE4(ISTAL,JC,KC)

C           BOUNDARY CONDITION DUE TO LODATO's FORMULATION
            T3BXL(JC,KC) = T3BXL(JC,KC)-
     +                     STORE5(ISTAL,JC,KC)/DRHS(ISTAL,JC,KC)
            T4BXL(JC,KC) = T4BXL(JC,KC)-
     +                     STORE6(ISTAL,JC,KC)/DRHS(ISTAL,JC,KC)
            T51BXL(JC,KC)=-VTMP(ISTAL,JC,KC)*STORE5(ISTAL,JC,KC)
     +                        -WTMP(ISTAL,JC,KC)*STORE6(ISTAL,JC,KC)
            T52BXL(JC,KC)=-PRUN(ISTAL,JC,KC)*(STORE2(ISTAL,JC,KC)+
     +                                            STORE3(ISTAL,JC,KC))
          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRPXR(JC,KC) = PRUN(ISTOL,JC,KC)
            BCL5XR(JC,KC) = STORE4(ISTOL,JC,KC)
C           BOUNDARY CONDITION DUE TO LODATO's FORMULATION
            T3BXR(JC,KC) = T3BXR(JC,KC)-
     +                     STORE5(ISTOL,JC,KC)/DRHS(ISTOL,JC,KC)
            T4BXR(JC,KC) = T4BXR(JC,KC)-
     +                     STORE6(ISTOL,JC,KC)/DRHS(ISTOL,JC,KC)
            T51BXR(JC,KC)=-VTMP(ISTOL,JC,KC)*STORE5(ISTOL,JC,KC)
     +                        -WTMP(ISTOL,JC,KC)*STORE6(ISTOL,JC,KC)
            T52BXR(JC,KC)=-PRUN(ISTOL,JC,KC)*(STORE2(ISTOL,JC,KC)+
     +                                            STORE3(ISTOL,JC,KC))
          ENDDO
        ENDDO
      ENDIF

C     Y-DIRECTION: P AND DPDY
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRPYL(IC,KC) = PRUN(IC,JSTAL,KC)
            BCL5YL(IC,KC) = STORE5(IC,JSTAL,KC)
C           BOUNDARY CONDITION DUE TO LODATO's FORMULATION
          !NC and VM changed here: ycorrection
            T4BYL(IC,KC) = T4BYL(IC,KC)-
     +                     STORE6(IC,JSTAL,KC)/DRHS(IC,JSTAL,KC)
            T3BYL(IC,KC) = T3BYL(IC,KC)-
     +                     STORE4(IC,JSTAL,KC)/DRHS(IC,JSTAL,KC)
            T51BYL(IC,KC)=-UTMP(IC,JSTAL,KC)*STORE4(IC,JSTAL,KC)
     +                        -WTMP(IC,JSTAL,KC)*STORE6(IC,JSTAL,KC)
            T52BYL(IC,KC)=-PRUN(IC,JSTAL,KC)*(STORE1(IC,JSTAL,KC)+
     +                                            STORE3(IC,JSTAL,KC))
          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRPYR(IC,KC) = PRUN(IC,JSTOL,KC)
            BCL5YR(IC,KC) = STORE5(IC,JSTOL,KC)
C           BOUNDARY CONDITION DUE TO LODATO's FORMULATION
          !NC and VM changed here: ycorrection
            T4BYR(IC,KC) = T4BYR(IC,KC)-
     +                     STORE6(IC,JSTOL,KC)/DRHS(IC,JSTOL,KC)
            T3BYR(IC,KC) = T3BYR(IC,KC)-
     +                     STORE4(IC,JSTOL,KC)/DRHS(IC,JSTOL,KC)
            T51BYR(IC,KC)=-UTMP(IC,JSTOL,KC)*STORE4(IC,JSTOL,KC)
     +                        -WTMP(IC,JSTOL,KC)*STORE6(IC,JSTOL,KC)
            T52BYR(IC,KC)=-PRUN(IC,JSTOL,KC)*(STORE1(IC,JSTOL,KC)+
     +                                            STORE3(IC,JSTOL,KC))
          ENDDO
        ENDDO
      ENDIF

C     Z-DIRECTION: P AND DPDZ
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRPZL(IC,JC) = PRUN(IC,JC,KSTAL)
            BCL5ZL(IC,JC) = STORE6(IC,JC,KSTAL)
C           BOUNDARY CONDITION DUE TO LODATO's FORMULATION
            T3BZL(IC,JC) = T3BZL(IC,JC)-
     +                     STORE4(IC,JC,KSTAL)/DRHS(IC,JC,KSTAL)
            T4BZL(IC,JC) = T4BZL(IC,JC)-
     +                     STORE5(IC,JC,KSTAL)/DRHS(IC,JC,KSTAL)
            T51BZL(IC,JC)=-UTMP(IC,JC,KSTAL)*STORE4(IC,JC,KSTAL)
     +                        -VTMP(IC,JC,KSTAL)*STORE5(IC,JC,KSTAL)
            T52BZL(IC,JC)=-PRUN(IC,JC,KSTAL)*(STORE1(IC,JC,KSTAL)+
     +                                            STORE2(IC,JC,KSTAL))
          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRPZR(IC,JC) = PRUN(IC,JC,KSTOL)
            BCL5ZR(IC,JC) = STORE6(IC,JC,KSTOL)
C           BOUNDARY CONDITION DUE TO LODATO's FORMULATION
            T3BZR(IC,JC) = T3BZR(IC,JC)-
     +                     STORE4(IC,JC,KSTOL)/DRHS(IC,JC,KSTOL)
            T4BZR(IC,JC) = T4BZR(IC,JC)-
     +                     STORE5(IC,JC,KSTOL)/DRHS(IC,JC,KSTOL)
            T51BZR(IC,JC)=-UTMP(IC,JC,KSTOL)*STORE4(IC,JC,KSTOL)
     +                        -VTMP(IC,JC,KSTOL)*STORE5(IC,JC,KSTOL)
            T52BZR(IC,JC)=-PRUN(IC,JC,KSTOL)*(STORE1(IC,JC,KSTOL)+
     +                                            STORE2(IC,JC,KSTOL))
          ENDDO
        ENDDO
      ENDIF

C     =========================================================================

C     U,V,W-EQUATIONS: PRESSURE GRADIENT TERMS
C     E-EQUATION: PRESSURE WORK TERMS
C     -------------------------------
C     DPDX,DPDY,DPDZ
C     U DPDX + V DPDY + W DPDZ

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC) - STORE4(IC,JC,KC)
            VRHS(IC,JC,KC) = VRHS(IC,JC,KC) - STORE5(IC,JC,KC)
            WRHS(IC,JC,KC) = WRHS(IC,JC,KC) - STORE6(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - UTMP(IC,JC,KC)*STORE4(IC,JC,KC)
     +                     - VTMP(IC,JC,KC)*STORE5(IC,JC,KC)
     +                     - WTMP(IC,JC,KC)*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C     =========================================================================

C     E-EQUATION: PRESSURE WORK TERMS
C     -------------------------------
C     P DIV U

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - PRUN(IC,JC,KC)*(STORE1(IC,JC,KC)
     +                                     + STORE2(IC,JC,KC)
     +                                     + STORE3(IC,JC,KC))

          ENDDO
        ENDDO
      ENDDO
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C     =========================================================================

C     VISCOSITY
C     ---------

C     VISCOSITY IS PARALLEL
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            TRANSP(IC,JC,KC) = TRANSP(IC,JC,KC)*PRANTL

          ENDDO
        ENDDO
      ENDDO
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C     -------------------------------------------------------------------------

C     MIXTURE-AVERAGED TRANSPORT
C     RSC 17-APR-2013
CC     DIAGNOSTICS
C      WRITE(6,*)'RHSVEL: visc: ',FLMAVT
      IF(FLMAVT)THEN


        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              TRANSP(IC,JC,KC) = DIFMIX(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

CC       DIAGNOSTICS
C        KC = 1
C        JC = 1
C        WRITE(6,'(4I5)')ITIME,IRKSTP,JC,KC
C        DO IC = ISTAB,ISTOB
C          WRITE(6,'(I5,2(1PE15.7))')IC,TRANSP(IC,JC,KC)
C        ENDDO

CC       DIAGNOSTICS
C        KC = 1
C        JC = 1
C        IC = 1
C        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
C        IC = 2
C        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
C        IC = 500
C        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
C        IC = 1000
C        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)
C        IC = 1001
C        WRITE(6,'(3I5,1PE15.7)')IC,JC,KC,TRANSP(IC,JC,KC)

      ENDIF

C     -------------------------------------------------------------------------

C     VISCOUS TERMS: TAUXXb,e,f
C     -------------

C     DVDY+DWDZ
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = STORE2(IC,JC,KC)+STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     TAUXXb,e,f
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = FTHD*STORE1(IC,JC,KC)
     +                       - TTHD*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXXb,e,f DUDX
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +    + TRANSP(IC,JC,KC)*STORE6(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                                           STORE6 = TAUXXb,e,f
C     =========================================================================

C     VISCOSITY GRADIENT: X COMPONENT
C     ------------------
      CALL DFBYDX(TRANSP,STORE4)
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                                                STORE4 = DMUDX
C                                                           STORE6 = TAUXXb,e,f
C     =========================================================================

C     U-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXX,Xb,e,f
C     U TAUXX,Xb,e,f

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = STORE6(IC,JC,KC)*STORE4(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN X: TAUXX,X TERM ZERO ON END POINTS
      IF(FXLVSN)CALL ZEROXL(STORE6)
      IF(FXRVSN)CALL ZEROXR(STORE6)
  
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC) + STORE6(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + UTMP(IC,JC,KC)*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                                                STORE4 = DMUDX
C     =========================================================================

C     VISCOUS TERMS: TAUYYb,e,f
C     -------------

C     DUDX+DWDZ
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = STORE1(IC,JC,KC)+STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     TAUYYb,e,f
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = FTHD*STORE2(IC,JC,KC)
     +                       - TTHD*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYYb,e,f DVDY
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +    + TRANSP(IC,JC,KC)*STORE6(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                                                STORE4 = DMUDX
C                                                           STORE6 = TAUYYb,e,f
C     =========================================================================

C     VISCOSITY GRADIENT: Y COMPONENT
C     ------------------
      CALL DFBYDY(TRANSP,STORE5)
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                                            STORE4,5 = DMUDX,Y
C                                                           STORE6 = TAUYYb,e,f
C     =========================================================================

C     V-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYY,Yb,e,f
C     V TAUYY,Yb,e,f

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE6(IC,JC,KC) = STORE6(IC,JC,KC)*STORE5(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Y: TAUYY,Y TERM ZERO ON END POINTS
      IF(FYLVSN)CALL ZEROYL(STORE6)
      IF(FYRVSN)CALL ZEROYR(STORE6)
  
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = VRHS(IC,JC,KC) + STORE6(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + VTMP(IC,JC,KC)*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                   STORE1,2,3 = DUDX,DVDY,DWDZ
C                                                            STORE4,5 = DMUDX,Y
C     =========================================================================

C     VISCOUS TERMS: TAUZZb,e,f
C     -------------

C     DUDX+DVDY
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE1(IC,JC,KC) = STORE1(IC,JC,KC)+STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     TAUZZb,e,f
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE1(IC,JC,KC) = FTHD*STORE3(IC,JC,KC)
     +                       - TTHD*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUZZb,e,f DWDZ
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +    + TRANSP(IC,JC,KC)*STORE1(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                                 STORE3 = DWDZ
C                                                            STORE4,5 = DMUDX,Y
C                                                           STORE1 = TAUZZb,e,f
C     =========================================================================

C     VISCOSITY GRADIENT: Z COMPONENT
C     ------------------
      CALL DFBYDZ(TRANSP,STORE6)

C                                                        STORE4,5,6 = DMUDX,Y,Z
C                                                           STORE1 = TAUZZb,e,f
C     =========================================================================

C     W-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUZZ,Zb,e,f
C     W TAUZZ,Zb,e,f

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Z: TAUZZ,Z TERM ZERO ON END POINTS
      IF(FZLVSN)CALL ZEROZL(STORE1)
      IF(FZRVSN)CALL ZEROZR(STORE1)
  
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = WRHS(IC,JC,KC) + STORE1(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + WTMP(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     VELOCITY CROSS-DERIVATIVES
C     --------------------------

C     DUDY
C     ----
      CALL DFBYDY(UTMP,STORE1)


C     COLLECT VELOCITY DERIVATIVE FOR BCs
C     -----------------------------------
C     Y-DIRECTION: DUDY
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            BCL3YL(IC,KC) = STORE1(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            BCL3YR(IC,KC) = STORE1(IC,JSTOL,KC)

          ENDDO
        ENDDO
      ENDIF
C                                                                 STORE1 = DUDY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     U-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF RHO V DUDY

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC)
     +             - HALF*DRHS(IC,JC,KC)*VTMP(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                                 STORE1 = DUDY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     DVDX
C     ----
      CALL DFBYDX(VTMP,STORE2)


C     COLLECT VELOCITY DERIVATIVES FOR BCs
C     ------------------------------------
C     X-DIRECTION: DVDX
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            BCL3XL(JC,KC) = STORE2(ISTAL,JC,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            BCL3XR(JC,KC) = STORE2(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF
C                                                                 STORE1 = DUDY
C                                                                 STORE2 = DVDX
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     V-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF RHO U DVDX

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = VRHS(IC,JC,KC)
     +             - HALF*DRHS(IC,JC,KC)*UTMP(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                                 STORE1 = DUDY
C                                                                 STORE2 = DVDX
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     VELOCITY SECOND CROSS DERIVATIVES
C     ---------------------------------

C     DUDY+DVDX
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE7(IC,JC,KC) = STORE1(IC,JC,KC) + STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D2UDXY
      CALL D2FDXY(UTMP,STORE1)

C     D2VDXY
      CALL D2FDXY(VTMP,STORE2)
C                                                               STORE1 = D2UDXY
C                                                               STORE2 = D2VDXY
C                                                          STORE7 = (DUDY+DVDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     U-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXX,Xc
C     U TAUXX,Xc

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TTHD*TRANSP(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN X: TAUXX,X TERMS ZERO ON END POINTS
      IF(FXLVSN)CALL ZEROXL(STORE3)
      IF(FXRVSN)CALL ZEROXR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC) - STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - UTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2UDXY
C                                                               STORE2 = D2VDXY
C                                                          STORE7 = (DUDY+DVDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     V-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYY,Yc
C     V TAUYY,Yc

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TTHD*TRANSP(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Y: TAUYY,Y TERMS ZERO ON END POINTS
      IF(FYLVSN)CALL ZEROYL(STORE3)
      IF(FYRVSN)CALL ZEROYR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = VRHS(IC,JC,KC) - STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - VTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2UDXY
C                                                               STORE2 = D2VDXY
C                                                          STORE7 = (DUDY+DVDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     D2UDY2
      CALL D2FDY2(UTMP,STORE3)

C     D2UDY2+D2VDXY
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = STORE3(IC,JC,KC) + STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     U-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXY,Ya,b,c
C     U TAUXY,Ya,b,c

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TRANSP(IC,JC,KC)*STORE3(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE5(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Y: TAUXY,Y TERM ZERO ON END POINTS
      IF(FYLVST)CALL ZEROYL(STORE3)
      IF(FYRVST)CALL ZEROYR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC) + STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + UTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2UDXY
C                                                          STORE7 = (DUDY+DVDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     D2VDX2
      CALL D2FDX2(VTMP,STORE3)

C     D2UDXY+D2VDX2
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = STORE3(IC,JC,KC) + STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     V-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYX,Xa,b,c
C     V TAUYX,Xa,b,c

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TRANSP(IC,JC,KC)*STORE3(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE4(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN X: TAUYX,X TERM ZERO ON END POINTS
      IF(FXLVST)CALL ZEROXL(STORE3)
      IF(FXRVST)CALL ZEROXR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = VRHS(IC,JC,KC) + STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + VTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                          STORE7 = (DUDY+DVDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXY(DUDY+DVDX)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +    + TRANSP(IC,JC,KC)*STORE7(IC,JC,KC)*STORE7(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     DUDZ
C     ----
      CALL DFBYDZ(UTMP,STORE1)


C     COLLECT VELOCITY DERIVATIVES FOR BCs
C     ------------------------------------
C     Z-DIRECTION: DUDZ
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            BCL3ZL(IC,JC) = STORE1(IC,JC,KSTAL)

          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            BCL3ZR(IC,JC) = STORE1(IC,JC,KSTOL)

          ENDDO
        ENDDO
      ENDIF
C                                                                 STORE1 = DUDZ
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     U-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF RHO W DUDZ

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC)
     +             - HALF*DRHS(IC,JC,KC)*WTMP(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                                 STORE1 = DUDZ
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     DWDX
C     ----
      CALL DFBYDX(WTMP,STORE2)


C     COLLECT VELOCITY DERIVATIVES FOR BCs
C     ------------------------------------
C     X-DIRECTION: DWDX
      IF(FXLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            BCL4XL(JC,KC) = STORE2(ISTAL,JC,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FXRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            BCL4XR(JC,KC) = STORE2(ISTOL,JC,KC)

          ENDDO
        ENDDO
      ENDIF
C                                                                 STORE1 = DUDZ
C                                                                 STORE2 = DWDX
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     W-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF RHO U DWDX

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = WRHS(IC,JC,KC)
     +             - HALF*DRHS(IC,JC,KC)*UTMP(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                                 STORE1 = DUDZ
C                                                                 STORE2 = DWDX
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================
C
C     VELOCITY SECOND CROSS DERIVATIVES
C     ---------------------------------

C     DUDZ+DWDX
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE7(IC,JC,KC) = STORE1(IC,JC,KC) + STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D2UDXZ
      CALL D2FDXZ(UTMP,STORE1)

C     D2WDXZ
      CALL D2FDXZ(WTMP,STORE2)
C                                                               STORE1 = D2UDXZ
C                                                               STORE2 = D2WDXZ
C                                                          STORE7 = (DUDZ+DWDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     U-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXX,Xd
C     U TAUXX,Xd

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TTHD*TRANSP(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN X: TAUXX,X TERMS ZERO ON END POINTS
      IF(FXLVSN)CALL ZEROXL(STORE3)
      IF(FXRVSN)CALL ZEROXR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC) - STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - UTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2UDXZ
C                                                               STORE2 = D2WDXZ
C                                                          STORE7 = (DUDZ+DWDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     W-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUZZ,Zc
C     W TAUZZ,Zc

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TTHD*TRANSP(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Z: TAUZZ,Z TERMS ZERO ON END POINTS
      IF(FZLVSN)CALL ZEROZL(STORE3)
      IF(FZRVSN)CALL ZEROZR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = WRHS(IC,JC,KC) - STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - WTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2UDXZ
C                                                               STORE2 = D2WDXZ
C                                                          STORE7 = (DUDZ+DWDX)
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     D2UDZ2
      CALL D2FDZ2(UTMP,STORE3)

C     D2UDZ2+D2WDXZ
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = STORE3(IC,JC,KC) + STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     U-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXZ,Za,b,c
C     U TAUXZ,Za,b,c

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TRANSP(IC,JC,KC)*STORE3(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Y: TAUXZ,Z TERM ZERO ON END POINTS
      IF(FZLVST)CALL ZEROZL(STORE3)
      IF(FZRVST)CALL ZEROZR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URHS(IC,JC,KC) = URHS(IC,JC,KC) + STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + UTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2UDXZ
C                                                          STORE7 = (DUDZ+DWDX)
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     D2WDX2
      CALL D2FDX2(WTMP,STORE3)

C     D2UDXZ+D2WDX2
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = STORE3(IC,JC,KC) + STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     W-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUZX,Xa,b,c
C     W TAUZX,Xa,b,c

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TRANSP(IC,JC,KC)*STORE3(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE4(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN X: TAUZX,X TERM ZERO ON END POINTS
      IF(FXLVST)CALL ZEROXL(STORE3)
      IF(FXRVST)CALL ZEROXR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = WRHS(IC,JC,KC) + STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + WTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                          STORE7 = (DUDZ+DWDX)
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXZ(DUDZ+DWDX)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +    + TRANSP(IC,JC,KC)*STORE7(IC,JC,KC)*STORE7(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     DVDZ
      CALL DFBYDZ(VTMP,STORE1)


C     COLLECT VELOCITY DERIVATIVES FOR BCs
C     ------------------------------------

C     Z-DIRECTION: DVDZ
      IF(FZLCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            BCL4ZL(IC,JC) = STORE1(IC,JC,KSTAL)

          ENDDO
        ENDDO
      ENDIF
      IF(FZRCNV)THEN
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            BCL4ZR(IC,JC) = STORE1(IC,JC,KSTOL)

          ENDDO
        ENDDO
      ENDIF
C                                                                 STORE1 = DVDZ
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     V-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF RHO W DVDZ

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = VRHS(IC,JC,KC)
     +             - HALF*DRHS(IC,JC,KC)*WTMP(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                                 STORE1 = DVDZ
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     DWDY
      CALL DFBYDY(WTMP,STORE2)


C     COLLECT VELOCITY DERIVATIVES FOR BCs
C     ------------------------------------

C     Y-DIRECTION: DWDY
      IF(FYLCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            BCL4YL(IC,KC) = STORE2(IC,JSTAL,KC)

          ENDDO
        ENDDO
      ENDIF
      IF(FYRCNV)THEN
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            BCL4YR(IC,KC) = STORE2(IC,JSTOL,KC)

          ENDDO
        ENDDO
      ENDIF
C                                                                 STORE1 = DVDZ
C                                                                 STORE2 = DWDY
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     W-EQUATION: CONVECTIVE TERMS
C     ----------------------------
C     HALF RHO V DWDY

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = WRHS(IC,JC,KC)
     +             - HALF*DRHS(IC,JC,KC)*VTMP(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                                 STORE1 = DVDZ
C                                                                 STORE2 = DWDY
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================
C
C     VELOCITY SECOND CROSS DERIVATIVES
C     ---------------------------------

C     DVDZ+DWDY
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE7(IC,JC,KC) = STORE1(IC,JC,KC) + STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     D2VDYZ
      CALL D2FDYZ(VTMP,STORE1)

C     D2WDYZ
      CALL D2FDYZ(WTMP,STORE2)
C
C                                                               STORE1 = D2VDYZ
C                                                               STORE2 = D2WDYZ
C                                                          STORE7 = (DVDZ+DWDY)
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     V-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYY,Yd
C     V TAUYY,Yd

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TTHD*TRANSP(IC,JC,KC)*STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Y: TAUYY,2 TERMS ZERO ON END POINTS
      IF(FYLVSN)CALL ZEROYL(STORE3)
      IF(FYRVSN)CALL ZEROYR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = VRHS(IC,JC,KC) - STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - VTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2VDYZ
C                                                               STORE2 = D2WDYZ
C                                                          STORE7 = (DVDZ+DWDY)
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     W-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUZZ,Zd
C     W TAUZZ,Zd

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TTHD*TRANSP(IC,JC,KC)*STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Z: TAUZZ,Z TERMS ZERO ON END POINTS
      IF(FZLVSN)CALL ZEROZL(STORE3)
      IF(FZRVSN)CALL ZEROZR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = WRHS(IC,JC,KC) - STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - WTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2VDYZ
C                                                               STORE2 = D2WDYZ
C                                                          STORE7 = (DVDZ+DWDY)
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     D2VDZ2
      CALL D2FDZ2(VTMP,STORE3)

C     D2VDZ2+D2WDYZ
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = STORE3(IC,JC,KC) + STORE2(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     V-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYZ,Za,b,c
C     V TAUYZ,Za,b,c

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TRANSP(IC,JC,KC)*STORE3(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE6(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Z: TAUYZ,Z TERM ZERO ON END POINTS
      IF(FZLVST)CALL ZEROZL(STORE3)
      IF(FZRVST)CALL ZEROZR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            VRHS(IC,JC,KC) = VRHS(IC,JC,KC) + STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + VTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                               STORE1 = D2VDYZ
C                                                          STORE7 = (DVDZ+DWDY)
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     D2WDY2
      CALL D2FDY2(WTMP,STORE3)

C     D2VDYZ+D2WDY2
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = STORE3(IC,JC,KC) + STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     W-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUZY,Ya,b,c
C     W TAUZY,Ya,b,c

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE3(IC,JC,KC) = TRANSP(IC,JC,KC)*STORE3(IC,JC,KC)
     +                       + STORE7(IC,JC,KC)*STORE5(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     BOUNDARY CONDITIONS
C     BC IN Y: TAUZY,Y TERM ZERO ON END POINTS
      IF(FYLVST)CALL ZEROYL(STORE3)
      IF(FYRVST)CALL ZEROYR(STORE3)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            WRHS(IC,JC,KC) = WRHS(IC,JC,KC) + STORE3(IC,JC,KC)
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     + WTMP(IC,JC,KC)*STORE3(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                          STORE7 = (DVDZ+DWDY)
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYZ(DVDZ+DWDY)

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +    + TRANSP(IC,JC,KC)*STORE7(IC,JC,KC)*STORE7(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
C                                                            TRANSP = VISCOSITY
C                                                        STORE4,5,6 = DMUDX,Y,Z
C     =========================================================================

C     VELOCITY SECOND NORMAL DERIVATIVE TERMS
C     ---------------------------------------
C     D2UDX2,D2VDY2,D2WDZ2
      CALL D2FDX2(UTMP,STORE1)
      CALL D2FDY2(VTMP,STORE2)
      CALL D2FDZ2(WTMP,STORE3)

C     BOUNDARY CONDITIONS
C     BC IN X: TAUXX,Xa TERM ZERO ON END POINTS
      IF(FXLVSN)CALL ZEROXL(STORE1)
      IF(FXRVSN)CALL ZEROXR(STORE1)
C     BC IN Y: TAUYY,Ya TERM ZERO ON END POINTS
      IF(FYLVSN)CALL ZEROYL(STORE2)
      IF(FYRVSN)CALL ZEROYR(STORE2)
C     BC IN Z: TAUZZ,Za TERM ZERO ON END POINTS
      IF(FZLVSN)CALL ZEROZL(STORE3)
      IF(FZRVSN)CALL ZEROZR(STORE3)


C     U-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUXX,Xa
C     U TAUXX,Xa

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = TRANSP(IC,JC,KC)*FTHD*STORE1(IC,JC,KC)
            URHS(IC,JC,KC) = URHS(IC,JC,KC) + FORNOW
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*UTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     V-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUYY,Ya
C     V TAUYY,Ya

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = TRANSP(IC,JC,KC)*FTHD*STORE2(IC,JC,KC)
            VRHS(IC,JC,KC) = VRHS(IC,JC,KC) + FORNOW
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*VTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO


C     W-EQUATION: VISCOUS STRESS TERMS
C     E-EQUATION: VISCOUS WORK TERMS
C     ------------------------------
C     TAUZZ,Za
C     W TAUZZ,Za

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = TRANSP(IC,JC,KC)*FTHD*STORE3(IC,JC,KC)
            WRHS(IC,JC,KC) = WRHS(IC,JC,KC) + FORNOW
            ERHS(IC,JC,KC) = ERHS(IC,JC,KC) + FORNOW*WTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     CONTINUITY EQUATION
C     -------------------
C     DIV RHO U 

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            DRHS(IC,JC,KC) = -DIVM(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
      
C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================


      RETURN
      END
