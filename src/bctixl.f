      SUBROUTINE BCTIXL
 
C     *************************************************************************
C
C     BCTIXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-MAR-2006:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES TURBULENT INFLOW 
C
C     X-DIRECTION LEFT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
      CHARACTER*4 PNTIXL,PNTCXL
      CHARACTER*4 PNXDAT
      PARAMETER(PNTIXL='tixl', PNTCXL='tcxl',
     +          PNXDAT='.dat')


C     LOCAL DATA
C     ==========
      INTEGER ICPROC
      INTEGER IC,JC,KC
      INTEGER KXMODD,IXMODD
      INTEGER ITOTXL
      CHARACTER*6 PNPROC
      LOGICAL FXDUMP


C     BEGIN
C     =====

C     =========================================================================

C     BUILD THE FILENAMES FOR THE INLET TURBULENT VELOCITY FIELD
      WRITE(PNPROC,'(I6.6)')IPROC
C     RESTART FILE
      FNTIXL = PNTIXL//PNPROC//PNXDAT
C     COLD START FILE
      FNTCXL = PNTCXL//PNPROC//PNXDAT

      NCTIXL = 3

C     =========================================================================

C     INLET COLD START SWITCH
C     PARAMETER I2=0
      IF(NXLPRM(2).EQ.0)THEN

C       =======================================================================

C       INLET COLD START
C       ----------------

C       CHECK AND INITIALISE RESTART FILE
        INQUIRE(FILE=FNTIXL,EXIST=FXDUMP)
        IF(.NOT.FXDUMP)THEN
          OPEN(UNIT=NCTIXL,FILE=FNTIXL,STATUS='NEW',FORM='UNFORMATTED')
          CLOSE(NCTIXL)
        ENDIF

C       READ THE COLD START INLET TURBULENT VELOCITY FIELD
        OPEN(UNIT=NCTIXL,FILE=FNTCXL,STATUS='OLD',FORM='UNFORMATTED')
        READ(NCTIXL)STORE1,STORE4,STORE5,STORE6
        CLOSE(NCTIXL)

C       SET THE REAL PARTS
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URUN(IC,JC,KC) = STORE4(IC,JC,KC)
              VRUN(IC,JC,KC) = STORE5(IC,JC,KC)
              WRUN(IC,JC,KC) = STORE6(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

C       ZERO THE IMAGINARY PARTS
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              UTMP(IC,JC,KC) = ZERO
              VTMP(IC,JC,KC) = ZERO
              WTMP(IC,JC,KC) = ZERO

            ENDDO
          ENDDO
        ENDDO

C       PARTIAL (X-WISE) FOURIER TRANSFORM
        CALL BUFTXL

C       MEAN INLET VELOCITY
C       SCANNING PLANE LOCATION AND VELOCITY
C       PARAMETER R1=MEAN VEL, R2=EXTRA SCAN VEL, R3=INITIAL LOCATION
        BVELXL = RXLPRM(1)
        SVELXL = BVELXL + RXLPRM(2)
        ELOCXL = XGDLEN - RXLPRM(3)
 
C       =======================================================================

      ELSE

C       =======================================================================

C       INLET RESTART
C       -------------
C       READ THE RESTART INLET TURBULENT VELOCITY FIELD
        OPEN(UNIT=NCTIXL,FILE=FNTIXL,STATUS='OLD',FORM='UNFORMATTED')
        READ(NCTIXL)UFXL,VFXL,WFXL,ELOCXL,SVELXL,BVELXL
        CLOSE(NCTIXL)

C       =======================================================================

      ENDIF

C     =========================================================================

C     INITIALISE RUNNING LOCATION
      SLOCXL = ELOCXL

C     INITIALISE SCALE FACTORS
      TPOVXG = TWO*PI/XGDLEN
      SCAUXL = TWO/REAL(NXGLBL)
      SCDUXL = -TWO*PI*SCAUXL*REAL(NXGLBL-1)/REAL(NXGLBL)/XGDLEN

C     INITIALISE FLAGS AND INDICES FOR INLET PLANE DFT
      FLLIXL = .FALSE.
      FLTRXL = .FALSE.
      ISTAXL = 1
      ISTOXL = NPMAPX(IXPROC)

      KMINXL = 0
      DO ICPROC = 0, IXPROC-1
        KMINXL = KMINXL + NPMAPX(ICPROC)
      ENDDO
      KXMODD = MOD(KMINXL,2)
      ITOTXL = KMINXL + NPMAPX(IXPROC)

      IXMODD = MOD(ISTOXL,2)

      IF(KXMODD.EQ.0)THEN

C       EVEN NUMBER OF POINTS TO LH SIDE OF THIS PROCESSOR
C       LH PROCESSOR HAS TRAILING REAL VALUE
C       LOCAL PROCESSOR MUST HAVE LEADING IMAGINARY VALUE
C       EXCEPT FOR FIRST POINT
        IF(KMINXL.GT.0)FLLIXL = .TRUE.
        KMINXL = KMINXL/2
        ISTAXL = 2
        IF(IXMODD.EQ.0)THEN
C         EVEN NUMBER OF POINTS ON LOCAL PROCESSOR
C         LOCAL PROCESSOR MUST HAVE TRAILING REAL VALUE
          ISTOXL = ISTOXL - 1
          IF(ITOTXL.EQ.NXGLBL)THEN
C           END OF PENCIL: IGNORE TRAILING REAL VALUE
            FLTRXL = .FALSE.
          ELSE
            FLTRXL = .TRUE.
          ENDIF
        ELSE
C         ODD NUMBER OF POINTS ON LOCAL PROCESSOR
C         LOCAL PROCESSOR HAS NO TRAILING REAL VALUE
          FLTRXL = .FALSE.
        ENDIF

      ELSE

C       ODD NUMBER OF POINTS TO LH SIDE OF THIS PROCESSOR
C       LH PROCESSOR HAS NO TRAILING REAL VALUE
C       LOCAL PROCESSOR MUST HAVE NO LEADING IMAGINARY VALUE
        FLLIXL = .FALSE.
        KMINXL = KMINXL/2 + 1
        ISTAXL = 1
        IF(IXMODD.EQ.0)THEN
C         EVEN NUMBER OF POINTS ON LOCAL PROCESSOR
C         LOCAL PROCESSOR HAS NO TRAILING REAL VALUE
          FLTRXL = .FALSE.
        ELSE
C         ODD NUMBER OF POINTS ON LOCAL PROCESSOR
C         LOCAL PROCESSOR MUST HAVE TRAILING REAL VALUE
          ISTOXL = ISTOXL - 1
          IF(ITOTXL.EQ.NXGLBL)THEN
C           END OF PENCIL: IGNORE TRAILING REAL VALUE
            FLTRXL = .FALSE.
          ELSE
            FLTRXL = .TRUE.
          ENDIF
        ENDIF

      ENDIF

C     =========================================================================


      RETURN
      END
