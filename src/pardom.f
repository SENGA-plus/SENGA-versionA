      SUBROUTINE PARDOM
 
C     *************************************************************************
C
C     PARDOM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     01-OCT-1996:  CREATED
C     09-MAR-2003:  RSC UPDATED FOR SENGA2
C     27-SEP-2009:  RSC UPDATE INDEXING FOR FILTERS
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     CARRIES OUT PARALLEL DOMAIN DECOMPOSITION
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ===========
      INTEGER NEXTRA
      INTEGER ICPROC


C     BEGIN
C     =====

C     =========================================================================

C     INITIALISE PARALLEL MESSAGE PASSING
C     ===================================
C     INITIALISE
      CALL P_INIT

C     GET THE NUMBER OF PROCESSORS NPROC
C     AND THE LOCAL PROCESSOR ID IPROC: IPROC=0,NPROC-1
      CALL P_SIZE(NPROC,IPROC)

C     =========================================================================

C     GET XYZ COORDINATES OF LOCAL PROCESSOR
C     ======================================
C     IPROC = IXPROC+IYPROC*NXPROC+IZPROC*NXPROC*NYPROC
      IXPROC = MOD(IPROC,NXPROC)
      IYPROC = IPROC/NXPROC
      IZPROC = IYPROC/NYPROC
      IYPROC = MOD(IYPROC,NYPROC)

C     =========================================================================

C     ASSIGN GRID POINTS TO PROCESSORS
C     ================================
C     STATIC LOAD-BALANCING:
C     LOAD IS AUTOMATICALLY BALANCED IF NXGLBL IS EXACTLY DIVISIBLE BY NXPROC
C     OTHERWISE ANY REMAINDER IS SPLIT EQUALLY BETWEEN THE HIGHER-RANKED
C     PROCESSORS IN EACH DIRECTION
      NXNODE = NXGLBL/NXPROC
      NEXTRA = MOD(NXGLBL,NXPROC)
      NEXTRA = NXPROC-NEXTRA
      NXPRM1 = NXPROC-1
      DO ICPROC = 0,NXPRM1
        NPMAPX(ICPROC) = NXNODE
        IF(ICPROC.GE.NEXTRA)NPMAPX(ICPROC) = NXNODE+1
      ENDDO
      NXNODE = NPMAPX(IXPROC)
      NXNBIG = NXNODE+2*NHALOX

      NYNODE = NYGLBL/NYPROC
      NEXTRA = MOD(NYGLBL,NYPROC)
      NEXTRA = NYPROC-NEXTRA
      NYPRM1 = NYPROC-1
      DO ICPROC = 0,NYPRM1
        NPMAPY(ICPROC) = NYNODE
        IF(ICPROC.GE.NEXTRA)NPMAPY(ICPROC) = NYNODE+1
      ENDDO
      NYNODE = NPMAPY(IYPROC)
      NYNBIG = NYNODE+2*NHALOY

      NZNODE = NZGLBL/NZPROC
      NEXTRA = MOD(NZGLBL,NZPROC)
      NEXTRA = NZPROC-NEXTRA
      NZPRM1 = NZPROC-1
      DO ICPROC = 0,NZPRM1
        NPMAPZ(ICPROC) = NZNODE
        IF(ICPROC.GE.NEXTRA)NPMAPZ(ICPROC) = NZNODE+1
      ENDDO
      NZNODE = NPMAPZ(IZPROC)
      NZNBIG = NZNODE+2*NHALOZ

C     =========================================================================

C     MAP PROCESSORS ALONG CURRENT ROW
C     ================================
      DO ICPROC = 0,NXPRM1
        NPROCX(ICPROC) = ICPROC+NXPROC*(IYPROC+NYPROC*IZPROC)
      ENDDO
      DO ICPROC = 0,NYPRM1
        NPROCY(ICPROC) = IXPROC+NXPROC*(ICPROC+NYPROC*IZPROC)
      ENDDO
      DO ICPROC = 0,NZPRM1
        NPROCZ(ICPROC) = IXPROC+NXPROC*(IYPROC+NYPROC*ICPROC)
      ENDDO

C     =========================================================================

C     IDENTIFY NEAREST NEIGHBOURS FOR MESSAGE-PASSING
C     ===============================================
      IXPROM = IXPROC-1
      IF(IXPROM.LT.0)IXPROM = NXPROC-1
      IXPROP = IXPROC+1
      IF(IXPROP.GE.NXPROC)IXPROP = 0
      IXPROM = IXPROM+NXPROC*(IYPROC+IZPROC*NYPROC)
      IXPROP = IXPROP+NXPROC*(IYPROC+IZPROC*NYPROC)

      IYPROM = IYPROC-1
      IF(IYPROM.LT.0)IYPROM = NYPROC-1
      IYPROP = IYPROC+1
      IF(IYPROP.GE.NYPROC)IYPROP = 0
      IYPROM = IXPROC+NXPROC*(IYPROM+IZPROC*NYPROC)
      IYPROP = IXPROC+NXPROC*(IYPROP+IZPROC*NYPROC)

      IZPROM = IZPROC-1
      IF(IZPROM.LT.0)IZPROM = NZPROC-1
      IZPROP = IZPROC+1
      IF(IZPROP.GE.NZPROC)IZPROP = 0
      IZPROM = IXPROC+NXPROC*(IYPROC+IZPROM*NYPROC)
      IZPROP = IXPROC+NXPROC*(IYPROC+IZPROP*NYPROC)

C     =========================================================================

C     SET FLAGS AND TAGS FOR MESSAGE-PASSING
C     ======================================
C     LOCAL PROCESSOR EVEN/ODD FLAGS
      PRODDX = .TRUE.
      IF(MOD(IXPROC,2).EQ.0)PRODDX = .FALSE.
      PRODDY = .TRUE.
      IF(MOD(IYPROC,2).EQ.0)PRODDY = .FALSE.
      PRODDZ = .TRUE.
      IF(MOD(IZPROC,2).EQ.0)PRODDZ = .FALSE.

C     SET MESSAGE TAGS
C     EACH TAG UNIQUELY IDENTIFIES THE SENDING AND RECEIVING PROCESS
C     FOR EACH PROCESSOR FOR EACH DIRECTION THERE ARE FOUR TAGS
C     CORRESPONDING TO SEND AND RECEIVE LEFT AND RIGHT
      ITGXSL = IPROC*NPROC + IXPROM
      ITGXRL = IXPROM*NPROC + IPROC
      ITGXSR = IPROC*NPROC + IXPROP
      ITGXRR = IXPROP*NPROC + IPROC
      ITGYSL = IPROC*NPROC + IYPROM
      ITGYRL = IYPROM*NPROC + IPROC
      ITGYSR = IPROC*NPROC + IYPROP
      ITGYRR = IYPROP*NPROC + IPROC
      ITGZSL = IPROC*NPROC + IZPROM
      ITGZRL = IZPROM*NPROC + IPROC
      ITGZSR = IPROC*NPROC + IZPROP
      ITGZRR = IZPROP*NPROC + IPROC

C     =========================================================================

C     SET LIMITS ON SPATIAL COUNTERS
C     ==============================
C     ON THE LOCAL PROCESSOR

C     STANDARD SIZE ARRAYS
      ISTAL = 1
      ISTOL = NXNODE
      JSTAL = 1
      JSTOL = NYNODE
      KSTAL = 1
      KSTOL = NZNODE

      ISTAP1 = ISTAL+1
      ISTOM1 = ISTOL-1
      JSTAP1 = JSTAL+1
      JSTOM1 = JSTOL-1
      KSTAP1 = KSTAL+1
      KSTOM1 = KSTOL-1

      ISTAP2 = ISTAL+2
      ISTOM2 = ISTOL-2
      JSTAP2 = JSTAL+2
      JSTOM2 = JSTOL-2
      KSTAP2 = KSTAL+2
      KSTOM2 = KSTOL-2

      ISTAP3 = ISTAL+3
      ISTOM3 = ISTOL-3
      JSTAP3 = JSTAL+3
      JSTOM3 = JSTOL-3
      KSTAP3 = KSTAL+3
      KSTOM3 = KSTOL-3

      ISTAP4 = ISTAL+4
      ISTOM4 = ISTOL-4
      JSTAP4 = JSTAL+4
      JSTOM4 = JSTOL-4
      KSTAP4 = KSTAL+4
      KSTOM4 = KSTOL-4

      ISTAP5 = ISTAL+5
      ISTOM5 = ISTOL-5
      JSTAP5 = JSTAL+5
      JSTOM5 = JSTOL-5
      KSTAP5 = KSTAL+5
      KSTOM5 = KSTOL-5
      
      ISTAP6 = ISTAL+6
      ISTOM6 = ISTOL-6
      JSTAP6 = JSTAL+6
      JSTOM6 = JSTOL-6
      KSTAP6 = KSTAL+6
      KSTOM6 = KSTOL-6

      ISTAP7 = ISTAL+7
      ISTOM7 = ISTOL-7
      JSTAP7 = JSTAL+7
      JSTOM7 = JSTOL-7
      KSTAP7 = KSTAL+7
      KSTOM7 = KSTOL-7

      ISTAP8 = ISTAL+8
      ISTOM8 = ISTOL-8
      JSTAP8 = JSTAL+8
      JSTOM8 = JSTOL-8
      KSTAP8 = KSTAL+8
      KSTOM8 = KSTOL-8

      ISTAP9 = ISTAL+9
      ISTOM9 = ISTOL-9
      JSTAP9 = JSTAL+9
      JSTOM9 = JSTOL-9
      KSTAP9 = KSTAL+9
      KSTOM9 = KSTOL-9

      ISTAPA = ISTAL+10
      ISTOMA = ISTOL-10
      JSTAPA = JSTAL+10
      JSTOMA = JSTOL-10
      KSTAPA = KSTAL+10
      KSTOMA = KSTOL-10

      ISTAPB = ISTAL+11
      ISTOMB = ISTOL-11
      JSTAPB = JSTAL+11
      JSTOMB = JSTOL-11
      KSTAPB = KSTAL+11
      KSTOMB = KSTOL-11

C     BIGGER SIZE ARRAYS
      ISTAB = 1-NHALOX
      ISTOB = NXNODE+NHALOX
      JSTAB = 1-NHALOY
      JSTOB = NYNODE+NHALOY
      KSTAB = 1-NHALOZ
      KSTOB = NZNODE+NHALOZ

C     HALO DATA
      ISTALI = 1
      ISTOLI = NHALOX
      JSTALI = 1
      JSTOLI = NHALOY
      KSTALI = 1
      KSTOLI = NHALOZ
      ISTARI = NXNODE+1-NHALOX
      ISTORI = NXNODE
      JSTARI = NYNODE+1-NHALOY
      JSTORI = NYNODE
      KSTARI = NZNODE+1-NHALOZ
      KSTORI = NZNODE
      ISTALO = 1-NHALOX
      ISTOLO = 0
      JSTALO = 1-NHALOY
      JSTOLO = 0
      KSTALO = 1-NHALOZ
      KSTOLO = 0
      ISTARO = NXNODE+1
      ISTORO = NXNODE+NHALOX
      JSTARO = NYNODE+1
      JSTORO = NYNODE+NHALOY
      KSTARO = NZNODE+1
      KSTORO = NZNODE+NHALOZ

C     =========================================================================


      RETURN
      END
