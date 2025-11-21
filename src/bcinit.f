      SUBROUTINE BCINIT
 
C     *************************************************************************
C
C     BCINIT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     26-MAR-2003:  CREATED
C     26-OCT-2008:  RSC/TDD BUG FIX BCTIXL
C     23-JAN-2015:  RSC BUG FIX WALL MASS FLUX BC
C     19-MAY-2015:  RSC UPDATE WALL DIFFUSIVE BC
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     CARRIES OUT INITIALISATION OF BOUNDARY CONDITIONS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     BEGIN
C     =====

C     =========================================================================

C     APPLY GLOBAL BOUNDARY CONDITIONS TO LOCAL PROCESSOR
C     ===================================================

C     X-LEFT
      IF(IXPROC.EQ.0)THEN
        NSBCXL = NGBCXL
        NENDXL = NBOUND
        IF(NSBCXL.EQ.NSPERI)THEN
          NENDXL = NOBC
          IF(NXPROC.EQ.1)NENDXL = NPERI
        ENDIF
      ELSE
        NSBCXL = NSNULL
        NENDXL = NOBC
      ENDIF

C     X-RIGHT
      IF(IXPROC.EQ.NXPRM1)THEN
        NSBCXR = NGBCXR
        NENDXR = NBOUND
        IF(NSBCXR.EQ.NSPERI)THEN
          NENDXR = NOBC
          IF(NXPROC.EQ.1)NENDXR = NPERI
        ENDIF
      ELSE
        NSBCXR = NSNULL
        NENDXR = NOBC
      ENDIF

C     Y-LEFT
      IF(IYPROC.EQ.0)THEN
        NSBCYL = NGBCYL
        NENDYL = NBOUND
        IF(NSBCYL.EQ.NSPERI)THEN
          NENDYL = NOBC
          IF(NYPROC.EQ.1)NENDYL = NPERI
        ENDIF
      ELSE
        NSBCYL = NSNULL
        NENDYL = NOBC
      ENDIF

C     Y-RIGHT
      IF(IYPROC.EQ.NYPRM1)THEN
        NSBCYR = NGBCYR
        NENDYR = NBOUND
        IF(NSBCYR.EQ.NSPERI)THEN
          NENDYR = NOBC
          IF(NYPROC.EQ.1)NENDYR = NPERI
        ENDIF
      ELSE
        NSBCYR = NSNULL
        NENDYR = NOBC
      ENDIF

C     Z-LEFT
      IF(IZPROC.EQ.0)THEN
        NSBCZL = NGBCZL
        NENDZL = NBOUND
        IF(NSBCZL.EQ.NSPERI)THEN
          NENDZL = NOBC
          IF(NZPROC.EQ.1)NENDZL = NPERI
        ENDIF
      ELSE
        NSBCZL = NSNULL
        NENDZL = NOBC
      ENDIF

C     Z-RIGHT
      IF(IZPROC.EQ.NZPRM1)THEN
        NSBCZR = NGBCZR
        NENDZR = NBOUND
        IF(NSBCZR.EQ.NSPERI)THEN
          NENDZR = NOBC
          IF(NZPROC.EQ.1)NENDZR = NPERI
        ENDIF
      ELSE
        NSBCZR = NSNULL
        NENDZR = NOBC
      ENDIF

C     =========================================================================

C     INITIALISE PARALLEL TRANSFER BC FLAGS
C     =====================================
      PRGOXL = .FALSE.
      IF(NENDXL.EQ.NOBC)PRGOXL = .TRUE.
      PRGOXR = .FALSE.
      IF(NENDXR.EQ.NOBC)PRGOXR = .TRUE.
      PRGOYL = .FALSE.
      IF(NENDYL.EQ.NOBC)PRGOYL = .TRUE.
      PRGOYR = .FALSE.
      IF(NENDYR.EQ.NOBC)PRGOYR = .TRUE.
      PRGOZL = .FALSE.
      IF(NENDZL.EQ.NOBC)PRGOZL = .TRUE.
      PRGOZR = .FALSE.
      IF(NENDZR.EQ.NOBC)PRGOZR = .TRUE.

C     =========================================================================

C     INITIALISE BC FLAGS
C     ===================
C     FLAG = .TRUE. INDICATES THAT TERMS ARE TO BE SET TO ZERO ON THE BOUNDARY 
C     F--CNV CONVECTIVE TERMS
C     F--VSN NORMAL VISCOUS TERMS
C     F--VST TANGENTIAL VISCOUS TERMS
C     F--CON THERMAL CONDUCTION TERMS
C     F--ADB OTHER HEAT FLUX TERMS
C     F--DIF MASS DIFFUSION TERMS

      FXLCNV = .FALSE.
      FXLVSN = .FALSE.
      FXLVST = .FALSE.
      FXLCON = .FALSE.
      FXLADB = .FALSE.
      FXLDIF = .FALSE.
C     RSC 19-MAY-2015
      FXLDFW = .FALSE.
C     RSC 28-JUN-2015
      FXLCNW = .FALSE.
      FXLADW = .FALSE.

      FXRCNV = .FALSE.
      FXRVSN = .FALSE.
      FXRVST = .FALSE.
      FXRCON = .FALSE.
      FXRADB = .FALSE.
      FXRDIF = .FALSE.
C     RSC 19-MAY-2015
      FXRDFW = .FALSE.
C     RSC 28-JUN-2015
      FXRCNW = .FALSE.
      FXRADW = .FALSE.

      FYLCNV = .FALSE.
      FYLVSN = .FALSE.
      FYLVST = .FALSE.
      FYLCON = .FALSE.
      FYLADB = .FALSE.
      FYLDIF = .FALSE.
C     RSC 19-MAY-2015
      FYLDFW = .FALSE.
C     RSC 28-JUN-2015
      FYLCNW = .FALSE.
      FYLADW = .FALSE.

      FYRCNV = .FALSE.
      FYRVSN = .FALSE.
      FYRVST = .FALSE.
      FYRCON = .FALSE.
      FYRADB = .FALSE.
      FYRDIF = .FALSE.
C     RSC 19-MAY-2015
      FYRDFW = .FALSE.
C     RSC 28-JUN-2015
      FYRCNW = .FALSE.
      FYRADW = .FALSE.

      FZLCNV = .FALSE.
      FZLVSN = .FALSE.
      FZLVST = .FALSE.
      FZLCON = .FALSE.
      FZLADB = .FALSE.
      FZLDIF = .FALSE.
C     RSC 19-MAY-2015
      FZLDFW = .FALSE.
C     RSC 28-JUN-2015
      FZLCNW = .FALSE.
      FZLADW = .FALSE.

      FZRCNV = .FALSE.
      FZRVSN = .FALSE.
      FZRVST = .FALSE.
      FZRCON = .FALSE.
      FZRADB = .FALSE.
      FZRDIF = .FALSE.
C     RSC 19-MAY-2015
      FZRDFW = .FALSE.
C     RSC 28-JUN-2015
      FZRCNW = .FALSE.
      FZRADW = .FALSE.

C     =========================================================================

C     FLAGS FOR TURBULENT INFLOW VELOCITY FIELD
C     =========================================
      FXLTRB = .FALSE.
      FXRTRB = .FALSE.
      FYLTRB = .FALSE.
      FYRTRB = .FALSE.
      FZLTRB = .FALSE.
      FZRTRB = .FALSE.

C     =========================================================================

C     INITIALISE RK COUNTERS
C     ======================
      ISTALD = ISTAL
      ISTOLD = ISTOL
      JSTALD = JSTAL
      JSTOLD = JSTOL
      KSTALD = KSTAL
      KSTOLD = KSTOL

      ISTALU = ISTAL
      ISTOLU = ISTOL
      JSTALU = JSTAL
      JSTOLU = JSTOL
      KSTALU = KSTAL
      KSTOLU = KSTOL

      ISTALV = ISTAL
      ISTOLV = ISTOL
      JSTALV = JSTAL
      JSTOLV = JSTOL
      KSTALV = KSTAL
      KSTOLV = KSTOL

      ISTALW = ISTAL
      ISTOLW = ISTOL
      JSTALW = JSTAL
      JSTOLW = JSTOL
      KSTALW = KSTAL
      KSTOLW = KSTOL

      ISTALE = ISTAL
      ISTOLE = ISTOL
      JSTALE = JSTAL
      JSTOLE = JSTOL
      KSTALE = KSTAL
      KSTOLE = KSTOL

      ISTALY = ISTAL
      ISTOLY = ISTOL
      JSTALY = JSTAL
      JSTOLY = JSTOL
      KSTALY = KSTAL
      KSTOLY = KSTOL

      ISTALT = ISTAL-NHALOX
      ISTOLT = ISTOL+NHALOX
      JSTALT = JSTAL-NHALOY
      JSTOLT = JSTOL+NHALOY
      KSTALT = KSTAL-NHALOZ
      KSTOLT = KSTOL+NHALOZ

C     =========================================================================

C     SET BC FLAGS, RK COUNTERS AND PHYSICAL PARAMETERS
C     =================================================

C     =========================================================================

C     OUTFLOW BC No 1
C     ---------------
C     SUBSONIC NON-REFLECTING OUTFLOW
C     WITH OPTION TO SET PRESSURE AT INFINITY
C     EQUIVALENT TO NSCBC OUTFLOW 1
C     PARAMETERS: R1=P AT INFINITY, R2=COEFF, R3=MACH NO 
C     X-LEFT
      IF(NSBCXL.EQ.NSBCO1)THEN
        FXLCNV = .TRUE.
        FXLVST = .TRUE.
        FXLCON = .TRUE.
        FXLADB = .TRUE.
        FXLDIF = .TRUE.
C       RSC 19-MAY-2015
C       RSC 28-JUN-2015
C        FXLDFW = .TRUE.
        ISTALT = ISTAL
        COBCXL = RXLPRM(2)*(ONE-RXLPRM(3))/(TWO*XGDLEN)
        PINFXL = RXLPRM(1)
      ENDIF

C     X-RIGHT
      IF(NSBCXR.EQ.NSBCO1)THEN
        FXRCNV = .TRUE.
        FXRVST = .TRUE.
        FXRCON = .TRUE.
        FXRADB = .TRUE.
        FXRDIF = .TRUE.
C       RSC 19-MAY-2015
C       RSC 28-JUN-2015
C        FXRDFW = .TRUE.
        ISTOLT = ISTOL
        COBCXR = RXRPRM(2)*(ONE-RXRPRM(3))/(TWO*XGDLEN)
        PINFXR = RXRPRM(1)
      ENDIF

C     Y-LEFT
      IF(NSBCYL.EQ.NSBCO1)THEN
        FYLCNV = .TRUE.
        FYLVST = .TRUE.
        FYLCON = .TRUE.
        FYLADB = .TRUE.
        FYLDIF = .TRUE.
C       RSC 19-MAY-2015
C       RSC 28-JUN-2015
C        FYLDFW = .TRUE.
        JSTALT = JSTAL
        COBCYL = RYLPRM(2)*(ONE-RYLPRM(3))/(TWO*YGDLEN)
        PINFYL = RYLPRM(1)
      ENDIF

C     Y-RIGHT
      IF(NSBCYR.EQ.NSBCO1)THEN
        FYRCNV = .TRUE.
        FYRVST = .TRUE.
        FYRCON = .TRUE.
        FYRADB = .TRUE.
        FYRDIF = .TRUE.
C       RSC 19-MAY-2015
C       RSC 28-JUN-2015
C        FYRDFW = .TRUE.
        JSTOLT = JSTOL
        COBCYR = RYRPRM(2)*(ONE-RYRPRM(3))/(TWO*YGDLEN)
        PINFYR = RYRPRM(1)
      ENDIF

C     Z-LEFT
      IF(NSBCZL.EQ.NSBCO1)THEN
        FZLCNV = .TRUE.
        FZLVST = .TRUE.
        FZLCON = .TRUE.
        FZLADB = .TRUE.
        FZLDIF = .TRUE.
C       RSC 19-MAY-2015
C       RSC 28-JUN-2015
C        FZLDFW = .TRUE.
        KSTALT = KSTAL
        COBCZL = RZLPRM(2)*(ONE-RZLPRM(3))/(TWO*ZGDLEN)
        PINFZL = RZLPRM(1)
      ENDIF

C     Z-RIGHT
      IF(NSBCZR.EQ.NSBCO1)THEN
        FZRCNV = .TRUE.
        FZRVST = .TRUE.
        FZRCON = .TRUE.
        FZRADB = .TRUE.
        FZRDIF = .TRUE.
C       RSC 19-MAY-2015
C       RSC 28-JUN-2015
C        FZRDFW = .TRUE.
        KSTOLT = KSTOL
        COBCZR = RZRPRM(2)*(ONE-RZRPRM(3))/(TWO*ZGDLEN)
        PINFZR = RZRPRM(1)
      ENDIF

C     =========================================================================

C     INFLOW BC no 1
C     --------------
C     SUBSONIC NON-REFLECTING LAMINAR INFLOW
C     EQUIVALENT TO NSCBC INFLOW 4

C     X-LEFT
      IF(NSBCXL.EQ.NSBCI1)THEN
        FXLCNV = .TRUE.
        FXLVSN = .TRUE.
        ISTALT = ISTAL
      ENDIF

C     X-RIGHT
      IF(NSBCXR.EQ.NSBCI1)THEN
        FXRCNV = .TRUE.
        FXRVSN = .TRUE.
        ISTOLT = ISTOL
      ENDIF

C     Y-LEFT
      IF(NSBCYL.EQ.NSBCI1)THEN
        FYLCNV = .TRUE.
        FYLVSN = .TRUE.
        JSTALT = JSTAL
      ENDIF

C     Y-RIGHT
      IF(NSBCYR.EQ.NSBCI1)THEN
        FYRCNV = .TRUE.
        FYRVSN = .TRUE.
        JSTOLT = JSTOL
      ENDIF

C     Z-LEFT
      IF(NSBCZL.EQ.NSBCI1)THEN
        FZLCNV = .TRUE.
        FZLVSN = .TRUE.
        KSTALT = KSTAL
      ENDIF

C     Z-RIGHT
      IF(NSBCZR.EQ.NSBCI1)THEN
        FZRCNV = .TRUE.
        FZRVSN = .TRUE.
        KSTOLT = KSTOL
      ENDIF

C     =========================================================================

C     INFLOW BC No 2
C     --------------
C     SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
C     WITH OPTION FOR INFLOW TURBULENCE
C     EQUIVALENT TO NSCBC INFLOW 1
C     AS IMPLEMENTED BY SUTHERLAND+KENNEDY

C     X-LEFT

C     GLOBAL BC SUPPORT
C     TURBULENT INFLOW VELOCITY FIELD
      IF(NGBCXL.EQ.NSBCI2)THEN
        IF(NXLPRM(1).EQ.3)THEN
          FXLTRB = .TRUE.
          CALL BCTIXL
        ENDIF
      ENDIF

C     LOCAL BC SUPPORT
      IF(NSBCXL.EQ.NSBCI2)THEN

        FXLCNV = .TRUE.
        ISTALU = ISTAP1
        ISTALV = ISTAP1
        ISTALW = ISTAP1
        ISTALE = ISTAP1
        ISTALY = ISTAP1
        ISTALT = ISTAL
        FXLTRB = .FALSE.

      ENDIF

C     X-RIGHT
      IF(NSBCXR.EQ.NSBCI2)THEN

        FXRCNV = .TRUE.
        ISTOLU = ISTOM1
        ISTOLV = ISTOM1
        ISTOLW = ISTOM1
        ISTOLE = ISTOM1
        ISTOLY = ISTOM1
        ISTOLT = ISTOL

      ENDIF

C     Y-LEFT
      IF(NSBCYL.EQ.NSBCI2)THEN

        FYLCNV = .TRUE.
        JSTALU = JSTAP1
        JSTALV = JSTAP1
        JSTALW = JSTAP1
        JSTALE = JSTAP1
        JSTALY = JSTAP1
        JSTALT = JSTAL

      ENDIF

C     Y-RIGHT
      IF(NSBCYR.EQ.NSBCI2)THEN

        FYRCNV = .TRUE.
        JSTOLU = JSTOM1
        JSTOLV = JSTOM1
        JSTOLW = JSTOM1
        JSTOLE = JSTOM1
        JSTOLY = JSTOM1
        JSTOLT = JSTOL

      ENDIF

C     Z-LEFT
      IF(NSBCZL.EQ.NSBCI2)THEN

        FZLCNV = .TRUE.
        KSTALU = KSTAP1
        KSTALV = KSTAP1
        KSTALW = KSTAP1
        KSTALE = KSTAP1
        KSTALY = KSTAP1
        KSTALT = KSTAL

      ENDIF

C     Z-RIGHT
      IF(NSBCZR.EQ.NSBCI2)THEN

        FZRCNV = .TRUE.
        KSTOLU = KSTOM1
        KSTOLV = KSTOM1
        KSTOLW = KSTOM1
        KSTOLE = KSTOM1
        KSTOLY = KSTOM1
        KSTOLT = KSTOL

      ENDIF

C     RSC/TDD BUG FIX BCTIXL
CC     INITIALISE INLET TURBULENCE GENERATOR
C      IF(FXLTRB)CALL BCTIXL

C     =========================================================================

C     INFLOW BC No 3
C     --------------
C     SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
C     WITH OPTION FOR INFLOW TURBULENCE
C     EQUIVALENT TO NSCBC INFLOW 2

C     X-LEFT

C     GLOBAL BC SUPPORT
C     X-LEFT TURBULENT INFLOW VELOCITY FIELD
      IF(NGBCXL.EQ.NSBCI3)THEN
        IF(NXLPRM(1).EQ.3)THEN
          FXLTRB = .TRUE.
          CALL BCTIXL
        ENDIF
      ENDIF

C     LOCAL BC SUPPORT
      IF(NSBCXL.EQ.NSBCI3)THEN

        FXLCNV = .TRUE.
        FXLVSN = .TRUE.
        ISTALD = ISTAP1
        ISTALU = ISTAP1
        ISTALV = ISTAP1
        ISTALW = ISTAP1
        ISTALY = ISTAP1
        ISTALT = ISTAL
        FXLTRB = .FALSE.

      ENDIF

C     X-RIGHT
      IF(NSBCXR.EQ.NSBCI3)THEN

        FXRCNV = .TRUE.
        FXRVSN = .TRUE.
        ISTOLD = ISTOM1
        ISTOLU = ISTOM1
        ISTOLV = ISTOM1
        ISTOLW = ISTOM1
        ISTOLY = ISTOM1
        ISTOLT = ISTOL

      ENDIF

C     Y-LEFT
      IF(NSBCYL.EQ.NSBCI3)THEN

        FYLCNV = .TRUE.
        FYLVSN = .TRUE.
        JSTALD = JSTAP1
        JSTALU = JSTAP1
        JSTALV = JSTAP1
        JSTALW = JSTAP1
        JSTALY = JSTAP1
        JSTALT = JSTAL

      ENDIF

C     Y-RIGHT
      IF(NSBCYR.EQ.NSBCI3)THEN

        FYRCNV = .TRUE.
        FYRVSN = .TRUE.
        JSTOLD = JSTOM1
        JSTOLU = JSTOM1
        JSTOLV = JSTOM1
        JSTOLW = JSTOM1
        JSTOLY = JSTOM1
        JSTOLT = JSTOL

      ENDIF

C     Z-LEFT
      IF(NSBCZL.EQ.NSBCI3)THEN

        FZLCNV = .TRUE.
        FZLVSN = .TRUE.
        KSTALD = KSTAP1
        KSTALU = KSTAP1
        KSTALV = KSTAP1
        KSTALW = KSTAP1
        KSTALY = KSTAP1
        KSTALT = KSTAL

      ENDIF

C     Z-RIGHT
      IF(NSBCZR.EQ.NSBCI3)THEN

        FZRCNV = .TRUE.
        FZRVSN = .TRUE.
        KSTOLD = KSTOM1
        KSTOLU = KSTOM1
        KSTOLV = KSTOM1
        KSTOLW = KSTOM1
        KSTOLY = KSTOM1
        KSTOLT = KSTOL

      ENDIF

C     =========================================================================

C     INFLOW BC no 4
C     --------------
C     SUBSONIC NON-REFLECTING INFLOW 
C     FY - IMPLEMENTATION AS PER LODATO
C     SPECIFY TEMPERATURE, VELOCITY AND MASS FRACTIONS

C     X-LEFT
      IF(NSBCXL.EQ.NSBCI4)THEN
        FXLCNV = .TRUE.
        FXLVSN = .TRUE. 
        ISTALT = ISTAL
        FXLTRB = .FALSE.
      ENDIF

C     X-RIGHT
      IF(NSBCXR.EQ.NSBCI4)THEN
        FXRCNV = .TRUE.
        FXRVSN = .TRUE.
        ISTOLT = ISTOL
      ENDIF

C     Y-LEFT
      IF(NSBCYL.EQ.NSBCI4)THEN
        FYLCNV = .TRUE.
        FYLVSN = .TRUE.
        JSTALT = JSTAL
      ENDIF

C     Y-RIGHT
      IF(NSBCYR.EQ.NSBCI4)THEN
        FYRCNV = .TRUE.
        FYRVSN = .TRUE.
        JSTOLT = JSTOL
      ENDIF

C     Z-LEFT
      IF(NSBCZL.EQ.NSBCI4)THEN
        FZLCNV = .TRUE.
        FZLVSN = .TRUE.
        KSTALT = KSTAL
      ENDIF

C     Z-RIGHT
      IF(NSBCZR.EQ.NSBCI4)THEN
        FZRCNV = .TRUE.
        FZRVSN = .TRUE.
        KSTOLT = KSTOL
      ENDIF

C     =========================================================================

C     WALL BC No 1
C     ------------
C     NO-SLIP WALL - ADIABATIC

C     X-LEFT
      IF(NSBCXL.EQ.NSBCW1)THEN

        FXLCNV = .TRUE.
        FXLCON = .TRUE.
        FXLADB = .TRUE.
C       RSC 23-JAN-2015
        FXLDIF = .TRUE.
C       RSC 19-MAY-2015
        FXLCNW = .TRUE.
        FXLDFW = .TRUE.
        ISTOW  = ISTAP4
        ISTALU = ISTAP1
        ISTALV = ISTAP1
        ISTALW = ISTAP1
        ISTALT = ISTAL

      ENDIF

C     X-RIGHT
      IF(NSBCXR.EQ.NSBCW1)THEN

        FXRCNV = .TRUE.
        FXRCON = .TRUE.
        FXRADB = .TRUE.
C       RSC 23-JAN-2015
        FXRDIF = .TRUE.
C       RSC 19-MAY-2015
        FXRCNW = .TRUE.
        FXRDFW = .TRUE.
        ISTAW  = ISTOM4
        ISTOLU = ISTOM1
        ISTOLV = ISTOM1
        ISTOLW = ISTOM1
        ISTOLT = ISTOL

      ENDIF

C     Y-LEFT
      IF(NSBCYL.EQ.NSBCW1)THEN

        FYLCNV = .TRUE.
        FYLCON = .TRUE.
        FYLADB = .TRUE.
C       RSC 23-JAN-2015
        FYLDIF = .TRUE.
C       RSC 19-MAY-2015
        FYLCNW = .TRUE.
        FYLDFW = .TRUE.
        JSTOW  = JSTAP4
        JSTALU = JSTAP1
        JSTALV = JSTAP1
        JSTALW = JSTAP1
        JSTALT = JSTAL

      ENDIF

C     Y-RIGHT
      IF(NSBCYR.EQ.NSBCW1)THEN

        FYRCNV = .TRUE.
        FYRCON = .TRUE.
        FYRADB = .TRUE.
C       RSC 23-JAN-2015
        FYRDIF = .TRUE.
C       RSC 19-MAY-2015
        FYRCNW = .TRUE.
        FYRDFW = .TRUE.
        JSTAW  = JSTOM4
        JSTOLU = JSTOM1
        JSTOLV = JSTOM1
        JSTOLW = JSTOM1
        JSTOLT = JSTOL

      ENDIF

C     Z-LEFT
      IF(NSBCZL.EQ.NSBCW1)THEN

        FZLCNV = .TRUE.
        FZLCON = .TRUE.
        FZLADB = .TRUE.
C       RSC 23-JAN-2015
        FZLDIF = .TRUE.
C       RSC 19-MAY-2015
        FZLCNW = .TRUE.
        FZLDFW = .TRUE.
        KSTOW  = KSTAP4
        KSTALU = KSTAP1
        KSTALV = KSTAP1
        KSTALW = KSTAP1
        KSTALT = KSTAL

      ENDIF

C     Z-RIGHT
      IF(NSBCZR.EQ.NSBCW1)THEN

        FZRCNV = .TRUE.
        FZRCON = .TRUE.
        FZRADB = .TRUE.
C       RSC 23-JAN-2015
        FZRDIF = .TRUE.
C       RSC 19-MAY-2015
        FZRCNW = .TRUE.
        FZRDFW = .TRUE.
        KSTAW  = KSTOM4
        KSTOLU = KSTOM1
        KSTOLV = KSTOM1
        KSTOLW = KSTOM1
        KSTOLT = KSTOL

      ENDIF

C     =========================================================================

C     WALL BC No 2
C     ------------
C     NO-SLIP WALL - ISOTHERMAL

C     X-LEFT
      IF(NSBCXL.EQ.NSBCW2)THEN

        FXLCNV = .TRUE.
C       RSC 23-JAN-2015
        FXLDIF = .TRUE.
C       RSC 19-MAY-2015
        FXLADW = .TRUE.
        FXLDFW = .TRUE.
        ISTOW  = ISTAP4
        ISTALU = ISTAP1
        ISTALV = ISTAP1
        ISTALW = ISTAP1
        ISTALE = ISTAP1
        ISTALT = ISTAL

      ENDIF

C     X-RIGHT
      IF(NSBCXR.EQ.NSBCW2)THEN

        FXRCNV = .TRUE.
C       RSC 23-JAN-2015
        FXRDIF = .TRUE.
C       RSC 19-MAY-2015
        FXRADW = .TRUE.
        FXRDFW = .TRUE.
        ISTAW  = ISTOM4
        ISTOLU = ISTOM1
        ISTOLV = ISTOM1
        ISTOLW = ISTOM1
        ISTOLE = ISTOM1
        ISTOLT = ISTOL

      ENDIF

C     Y-LEFT
      IF(NSBCYL.EQ.NSBCW2)THEN

        FYLCNV = .TRUE.
C       RSC 23-JAN-2015
        FYLDIF = .TRUE.
C       RSC 19-MAY-2015
        FYLADW = .TRUE.
        FYLDFW = .TRUE.
        JSTOW  = JSTAP4
        JSTALU = JSTAP1
        JSTALV = JSTAP1
        JSTALW = JSTAP1
        JSTALE = JSTAP1
        JSTALT = JSTAL

      ENDIF

C     Y-RIGHT
      IF(NSBCYR.EQ.NSBCW2)THEN

        FYRCNV = .TRUE.
C       RSC 23-JAN-2015
        FYRDIF = .TRUE.
C       RSC 19-MAY-2015
        FYRADW = .TRUE.
        FYRDFW = .TRUE.
        JSTAW  = JSTOM4
        JSTOLU = JSTOM1
        JSTOLV = JSTOM1
        JSTOLW = JSTOM1
        JSTOLE = JSTOM1
        JSTOLT = JSTOL

      ENDIF

C     Z-LEFT
      IF(NSBCZL.EQ.NSBCW2)THEN

        FZLCNV = .TRUE.
C       RSC 23-JAN-2015
        FZLDIF = .TRUE.
C       RSC 19-MAY-2015
        FZLADW = .TRUE.
        FZLDFW = .TRUE.
        KSTOW  = KSTAP4
        KSTALU = KSTAP1
        KSTALV = KSTAP1
        KSTALW = KSTAP1
        KSTALE = KSTAP1
        KSTALT = KSTAL

      ENDIF

C     Z-RIGHT
      IF(NSBCZR.EQ.NSBCW2)THEN

        FZRCNV = .TRUE.
C       RSC 23-JAN-2015
        FZRDIF = .TRUE.
C       RSC 19-MAY-2015
        FZRADW = .TRUE.
        FZRDFW = .TRUE.
        KSTAW  = KSTOM4
        KSTOLU = KSTOM1
        KSTOLV = KSTOM1
        KSTOLW = KSTOM1
        KSTOLE = KSTOM1
        KSTOLT = KSTOL

      ENDIF

C     =========================================================================


      RETURN
      END
