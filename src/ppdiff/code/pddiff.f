      SUBROUTINE PDDIFF
 
C     *************************************************************************
C
C     PDDIFF
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-SEP-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     PREPROCESSES MOLECULAR DIFFUSION DATA FOR PPDIFF
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
      INTEGER NUMSTR
      PARAMETER(NUMSTR=12)


C     EXTERNAL FUNCTION
C     ==================
      LOGICAL CHKSTR
      EXTERNAL CHKSTR


C     LOCAL DATA
C     ==========
      INTEGER LENSTR(NUMSTR)
      INTEGER ISPEC
      INTEGER NPARSE
      INTEGER ISTR,ISTOTL
      INTEGER LENCHK
 
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR
      CHARACTER*10 SYMCHK

      LOGICAL SPCFLG(NSPCMX)
      LOGICAL VALIDS


C     BEGIN
C     =====

C     =========================================================================

C     READ THE DIFFUSION DATA
C     -----------------------
      OPEN(UNIT=NCDIFF,FILE=FNDIFF,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1010)

C     FIVE LINES OF COMMENT ARE IGNORED
C     ---------------------------------
      READ(NCDIFF,'(A)',END=1020)ALINE
      READ(NCDIFF,'(A)',END=1020)ALINE
      READ(NCDIFF,'(A)',END=1020)ALINE
      READ(NCDIFF,'(A)',END=1020)ALINE
      READ(NCDIFF,'(A)',END=1020)ALINE

C     INITIALISE THE SPECIES-DETECTED FLAGS
      DO ISPEC = 1, NSPEC
        SPCFLG(ISPEC) = .FALSE.
      ENDDO

C     IGNORE HEADER LINE
      READ(NCDIFF,'(A)',END=1020)ALINE

1000  CONTINUE

C       READ LINE AND INITIALISE PARSER
1100    CONTINUE
          READ(NCDIFF,'(A)',END=1020)ALINE
          NPARSE = LEN(ALINE)
C          WRITE(NCREPT,'(A)')ALINE
          CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
          IF(ISTOTL.EQ.0)THEN
            ERRSTR = 'blank line in diffusion data'
            CALL ERHAND(ERRSTR,IEWARN)
            GOTO 1100
          ENDIF
C       END OF LOOP 1100

C       CHECK FOR END OF LIST
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
        IF(VALIDS)THEN

C         END OF LIST DETECTED
C         CHECK THAT DATA HAS BEEN COLLECTED FOR ALL SPECIES
          ISPEC = 1
1200      CONTINUE
            IF(SPCFLG(ISPEC))THEN
              IF(ISPEC.LT.NSPEC)THEN
                ISPEC = ISPEC + 1
                GOTO 1200
              ENDIF
            ELSE
              ERRSTR = 'end of diffusion data list: data incomplete'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
C         END OF LOOP 1200

        ELSE

C         CHECK NUMBER OF ITEMS
          IF(NUMSTR.LT.7)THEN
            ERRSTR = 'invalid diffusion data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         CHECK SPECIES SYMBOL
          LENCHK = LENSTR(1)
          SYMCHK(1:LENCHK) = PARSTR(1)(1:LENCHK)

          ISPEC = 1
1300      CONTINUE
            IF(SYMCHK(1:LENCHK).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN

C             SPECIES SYMBOL RECOGNISED

C             SET SPECIES FLAG
              SPCFLG(ISPEC) = .TRUE.

C             PROCEED TO READ THE SPECIES DATA
              ERRSTR = 'invalid diffusion data format'

C             INTERPRET DATA LINE

C             ITEM 2: INTEGER = 0, 1 OR 2
              VALIDS = CHKSTR(PARSTR(2),LENSTR(2),ICHINT)
              IF(VALIDS)THEN
C               MOLECULAR GEOMETRY INDICATOR
                READ(PARSTR(2),*)IMGEOM(ISPEC)
                IF((IMGEOM(ISPEC).LT.0).OR.(IMGEOM(ISPEC).GT.2))
     +            CALL ERHAND(ERRSTR,IEFATL)
              ELSE
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

C             ITEM 3: REAL NUMBER IN F FORMAT
              VALIDS = CHKSTR(PARSTR(3),LENSTR(3),ICHFNO)
              IF(VALIDS)THEN
C               LENNARD-JONES POTENTIAL WELL DEPTH IN DEGREES KELVIN
                READ(PARSTR(3),*)EPSOKB(ISPEC)
              ELSE
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

C             ITEM 4: REAL NUMBER IN F FORMAT
              VALIDS = CHKSTR(PARSTR(4),LENSTR(4),ICHFNO)
              IF(VALIDS)THEN
C               LENNARD-JONES COLLISION DIAMETER IN ANGSTROMS
                READ(PARSTR(4),*)SIGMAD(ISPEC)
              ELSE
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

C             ITEM 5: REAL NUMBER IN F FORMAT
              VALIDS = CHKSTR(PARSTR(5),LENSTR(5),ICHFNO)
              IF(VALIDS)THEN
C               DIPOLE MOMENT IN DEBYE
                READ(PARSTR(5),*)DIMOMU(ISPEC)
              ELSE
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

C             ITEM 6: REAL NUMBER IN F FORMAT
              VALIDS = CHKSTR(PARSTR(6),LENSTR(6),ICHFNO)
              IF(VALIDS)THEN
C               POLARISABILITY IN CUBIC ANGSTROMS
                READ(PARSTR(6),*)POLALF(ISPEC)
              ELSE
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

C             ITEM 7: REAL NUMBER IN F FORMAT
              VALIDS = CHKSTR(PARSTR(7),LENSTR(7),ICHFNO)
              IF(VALIDS)THEN
C               ROTATIONAL RELAXATION COLLISION NUMBER AT 298K
                READ(PARSTR(7),*)ZEDROT(ISPEC)
              ELSE
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

            ELSE

C             SPECIES SYMBOL NOT RECOGNISED

              IF(ISPEC.LT.NSPEC)THEN

C               CHECK THE NEXT SPECIES SYMBOL
                ISPEC = ISPEC + 1
                GOTO 1300

              ELSE

C               ALL SPECIES CHECKED - NO MATCH
C               IGNORE THE SPECIES DATA
                CONTINUE

              ENDIF

C           END OF LOOP 1300

          ENDIF
C         CHECK SPECIES SYMBOL

C         READ ANOTHER LINE
          GOTO 1000

        ENDIF
C       CHECK FOR END OF FILE

C     END OF LOOP 1000

C     =========================================================================

C     READ THE CONVERSION FACTORS
C     ---------------------------
C     IGNORE TITLE LINE
      READ(NCDIFF,'(A)',END=1020)ALINE

C     READ CONVERSION FACTOR DATA
2000  CONTINUE
        READ(NCDIFF,'(A)',END=1020)ALINE

        NPARSE = LEN(ALINE)
C        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in conversion factor data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 2000
        ENDIF
C     END OF LOOP 2000

C     CHECK NUMBER OF ITEMS
      IF(NUMSTR.LT.4)THEN
        ERRSTR = 'invalid conversion factor data format'
        CALL ERHAND(ERRSTR,IEFATL)
      ENDIF

C     ITEM 1: REAL NUMBER IN E FORMAT
      VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHRNO)
      IF(VALIDS)THEN
C       CONVERSION FOR LENNARD-JONES POTENTIAL WELL DEPTH TO DEGREES KELVIN
        READ(PARSTR(1),*)CLJP2K
      ELSE
        CALL ERHAND(ERRSTR,IEFATL)
      ENDIF

C     ITEM 2: REAL NUMBER IN E FORMAT
      VALIDS = CHKSTR(PARSTR(2),LENSTR(2),ICHRNO)
      IF(VALIDS)THEN
C       CONVERSION FOR LENNARD-JONES COLLISION DIAMETER TO METRES
        READ(PARSTR(2),*)CLJD2M
      ELSE
        CALL ERHAND(ERRSTR,IEFATL)
      ENDIF

C     ITEM 3: REAL NUMBER IN E FORMAT
      VALIDS = CHKSTR(PARSTR(3),LENSTR(3),ICHRNO)
      IF(VALIDS)THEN
C       CONVERSION FOR DIPOLE MOMENT IN DEBYE
        READ(PARSTR(3),*)CDM2SI
      ELSE
        CALL ERHAND(ERRSTR,IEFATL)
      ENDIF

C     ITEM 4: REAL NUMBER IN E FORMAT
      VALIDS = CHKSTR(PARSTR(4),LENSTR(4),ICHRNO)
      IF(VALIDS)THEN
C       CONVERSION FOR POLARISABILITY TO CUBIC METRES
        READ(PARSTR(4),*)CPL2M3
      ELSE
        CALL ERHAND(ERRSTR,IEFATL)
      ENDIF

C     =========================================================================

C     CLOSE THE DIFFUSION DATA FILE
C     -----------------------------
      CLOSE(NCDIFF)

C     =========================================================================

C     DIAGNOSTICS
      WRITE(6,*)'PDDIFF: raw molecular data'
      WRITE(6,'(I5)')NSPEC
      WRITE(6,9000)
9000  FORMAT(17X,'  geom',' LJ well ','  col.diam ','dip.mom ',
     +           ' polariz.','  Zrot ')
      DO ISPEC = 1, NSPEC
        WRITE(6,'(I5,2X,A10,I5,5F9.3)')ISPEC,SPCSYM(ISPEC),
     +        IMGEOM(ISPEC),
     +        EPSOKB(ISPEC),SIGMAD(ISPEC),DIMOMU(ISPEC),
     +        POLALF(ISPEC),ZEDROT(ISPEC)
      ENDDO
      WRITE(6,*)

C     CONVERT UNITS TO SI
      DO ISPEC = 1, NSPEC

        EPSOKB(ISPEC) = EPSOKB(ISPEC)*CLJP2K
        SIGMAD(ISPEC) = SIGMAD(ISPEC)*CLJD2M
        DIMOMU(ISPEC) = DIMOMU(ISPEC)*CDM2SI
        POLALF(ISPEC) = POLALF(ISPEC)*CPL2M3

      ENDDO

C     =========================================================================

C     FILE ERROR HANDLER
      GOTO 1030
1010  ERRSTR = 'Error opening diffusion data file:'
      CALL ERHAND(ERRSTR,IEFATL)
1020  ERRSTR = 'unexpected end of diffusion data file'
      CALL ERHAND(ERRSTR,IEFATL)
1030  CONTINUE
C     END OF FILE ERROR HANDLER

C     =========================================================================

C     DIAGNOSTICS
      WRITE(6,*)'PDDIFF: molecular data in SI units'
      WRITE(6,'(I5)')NSPEC
      WRITE(6,9010)
9010  FORMAT(17X,'  LJ well   ','  col.diam  ','  dip.mom   ',
     +           '  polariz.')
      DO ISPEC = 1, NSPEC
        WRITE(6,'(I5,2X,A10,4(1PE12.4))')ISPEC,SPCSYM(ISPEC),
     +        EPSOKB(ISPEC),SIGMAD(ISPEC),DIMOMU(ISPEC),
     +        POLALF(ISPEC)
      ENDDO
      WRITE(6,*)

C     =========================================================================


      RETURN
      END
