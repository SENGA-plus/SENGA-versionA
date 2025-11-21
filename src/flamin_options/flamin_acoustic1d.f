      SUBROUTINE FLAMIN 
 
C     *************************************************************************
C
C     FLAMIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     28-DEC-2003:  CREATED
C 
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     SETS INITIAL THERMOCHEMICAL FIELD
C     1D GAUSSIAN ACOUSTIC PULSE
C     
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
C     PRESSURE PULSE AMPLITUDE
      DOUBLE PRECISION PRPEAK
      PARAMETER(PRPEAK = 1.0D2)

C     PRESSURE PULSE LOCATION AND THICKNESS
      DOUBLE PRECISION PLOCAT,PTHICK
      PARAMETER(PLOCAT = 5.0D-3, PTHICK = 1.0D-4)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION DELTAG,XCOORD,ARGMNT,FORNOW
      DOUBLE PRECISION CPLOCL,CVLOCL,RGLOCL,GAMRAT
      INTEGER IC,JC,KC
      INTEGER ICPROC,IGOFST,IX
      INTEGER ISPEC,ITINT,ICP


C     BEGIN
C     =====

C     =========================================================================

C     SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
C     =========================================

C     GLOBAL COORDINATES
C     ------------------
      DELTAG = XGDLEN/(REAL(NXGLBL-1))

      IGOFST = 0
      DO ICPROC = 0, IXPROC-1
        IGOFST = IGOFST + NPMAPX(ICPROC)
      ENDDO


C     SET PRESSURE PROFILE
C     --------------------
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            IX = IGOFST + IC
            XCOORD = REAL(IX-1)*DELTAG
            ARGMNT = (XCOORD-PLOCAT)/PTHICK
            PRUN(IC,JC,KC) = PRIN + PRPEAK*EXP(-ARGMNT*ARGMNT)

          ENDDO
        ENDDO
      ENDDO


C     SET TEMPERATURE PROFILE ASSUMING CONSTANT ENTROPY
C     -------------------------------------------------
      CPLOCL = ZERO
      RGLOCL = ZERO
      DO ISPEC = 1,NSPEC

C       TEMPERATURE INTERVAL INDEX
        ITINT = 1
1010    CONTINUE
          IF(TRIN.GT.TINTHI(ITINT,ISPEC))THEN
            IF(ITINT.LT.NTINT(ISPEC))THEN
              ITINT = ITINT + 1
              GOTO 1010
            ENDIF
          ENDIF
C       END OF LOOP 1010

C       SPECIFIC HEAT CP
        FORNOW = AMASCP(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
        DO ICP = NCPOM1(ITINT,ISPEC),1,-1
          FORNOW = FORNOW*TRIN + AMASCP(ICP,ITINT,ISPEC)
        ENDDO
        CPLOCL = CPLOCL + FORNOW*YRIN(ISPEC)
 
C       MIXTURE GAS CONSTANT
        RGLOCL = RGLOCL + RGSPEC(ISPEC)*YRIN(ISPEC)

      ENDDO

C     SPECIFIC HEAT CP
      CVLOCL = CPLOCL - RGLOCL

C     GAMMA = CP/CV
      GAMRAT = CPLOCL/CVLOCL
      GAMRAT = (GAMRAT-ONE)/GAMRAT

C     SET TEMPERATURE
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            TRUN(IC,JC,KC) = TRIN*EXP(GAMRAT*LOG(PRUN(IC,JC,KC)/PRIN))

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
