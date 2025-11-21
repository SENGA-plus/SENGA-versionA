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
C     1D GAUSSIAN TEMPERATURE PULSE
C     
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
C     TEMPERATURE PULSE AMPLITUDE
      DOUBLE PRECISION TRPEAK
      PARAMETER(TRPEAK = 1.0D2)

C     TEMPERATURE PULSE LOCATION AND THICKNESS
      DOUBLE PRECISION TLOCAT,TTHICK
      PARAMETER(TLOCAT = 5.0D-3, TTHICK = 1.0D-4)


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


C     SET TEMPERATURE PROFILE
C     -----------------------
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

            IX = IGOFST + IC
            XCOORD = REAL(IX-1)*DELTAG
            ARGMNT = (XCOORD-TLOCAT)/TTHICK
            TRUN(IC,JC,KC) = TRIN + TRPEAK*EXP(-ARGMNT*ARGMNT)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
