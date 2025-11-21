      SUBROUTINE PDSPEC
 
C     *************************************************************************
C
C     PDSPEC
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
C     PREPROCESSES SPECIES DATA FOR PPDIFF
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
      PARAMETER(NUMSTR=4)


C     EXTERNAL FUNCTION
C     ==================
      LOGICAL CHKSTR
      EXTERNAL CHKSTR


C     LOCAL DATA
C     ==========
      INTEGER LENSTR(NUMSTR)
      INTEGER NPARSE
      INTEGER ISTR,JSTR,ISTOTL
      INTEGER LENCHK
      INTEGER ISPEC
 
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR
      CHARACTER*16 SYMCHK

      LOGICAL VALIDS


C     BEGIN
C     =====

C     DIAGNOSTICS
      WRITE(6,*)'PDSPEC: species data'
      WRITE(6,*)

C     =========================================================================

C     READ THE SPECIES DATA FILE
C     --------------------------
      OPEN(UNIT=NCSPEC,FILE=FNSPEC,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1000)

C     FIVE LINES OF COMMENT ARE IGNORED
C     ---------------------------------
      READ(NCSPEC,'(A)',END=1000)ALINE
      READ(NCSPEC,'(A)',END=1000)ALINE
      READ(NCSPEC,'(A)',END=1000)ALINE
      READ(NCSPEC,'(A)',END=1000)ALINE
      READ(NCSPEC,'(A)',END=1000)ALINE

C     =========================================================================

C     READ THE SPECIES LIST
C     ---------------------

C     READ AND IGNORE HEADER LINE
      READ(NCSPEC,'(A)',END=1000)ALINE
      WRITE(NCREPT,'(A)')ALINE

C     READ AND PARSE EACH LINE OF THE SPECIES LIST
      NSPEC = 0
      LENCHK = 0
2000  CONTINUE

C       INCREMENT SPECIES COUNTER AND CHECK
        NSPEC = NSPEC + 1
        IF(NSPEC.GT.NSPCMX)THEN
          ERRSTR = 'max no. of species exceeded'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF

C       READ LINE AND INITIALISE PARSER
        READ(NCSPEC,'(A)',END=1000)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE

        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)

C       CONDITION FOR END OF SPECIES LIST
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)

        IF(VALIDS)THEN

C         SPECIES LIST FINISHED
          NSPEC = NSPEC - 1

        ELSE

          IF(ISTOTL.EQ.0)THEN

            ERRSTR = 'blank line in species list'
            CALL ERHAND(ERRSTR,IEWARN)

          ELSE

C           CONVERT STRINGS TO SPECIES DATA

C           SPECIES NUMBER
C           CHECK FOR VALID NUMBERS
            DO JSTR = 1, LENSTR(1)
              IF((PARSTR(1)(JSTR:JSTR).LT.'0')
     +          .OR.(PARSTR(1)(JSTR:JSTR).GT.'9'))THEN
                ERRSTR = 'species number format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
            ENDDO
            READ(PARSTR(1)(1:LENSTR(1)),*)ISPEC
            IF(ISPEC.NE.NSPEC)THEN
              ERRSTR = 'species number out of order'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
            NUMSPC(NSPEC) = ISPEC

C           SPECIES SYMBOL
            LENCHK = LENSTR(2)
            SYMCHK(1:LENCHK) = PARSTR(2)(1:LENCHK)

C           CHECK SPECIES SYMBOL DOES NOT END WITH "D" OR "E"
C           REQUIRED TO AVOID CONFLICT WITH RATE DATA FORMAT
            IF((SYMCHK(LENCHK:LENCHK).EQ.'D')
     +      .OR.(SYMCHK(LENCHK:LENCHK).EQ.'E'))THEN
              ERRSTR = 'species name cannot end with D or E'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF

C           CHECK SPECIES SYMBOL IS NOT REPEATED
            ISPEC = 1
2030        CONTINUE
              IF(ISPEC.LT.NSPEC)THEN
              IF(SYMCHK(1:LENCHK).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN
                  ERRSTR = 'species symbol duplicated'
                  CALL ERHAND(ERRSTR,IEFATL)
                ELSE
                  ISPEC = ISPEC + 1
                  GOTO 2030
                ENDIF
              ELSE
                SPCSYM(NSPEC) = SYMCHK(1:LENCHK)
                LENSYM(NSPEC) = LENCHK
              ENDIF
C           END OF LOOP 2030

C           SPECIES MOLAR MASS
            VALIDS = CHKSTR(PARSTR(3),LENSTR(3),ICHRNO)
            IF(VALIDS)THEN
              READ(PARSTR(3)(1:LENSTR(3)),*)WMOLAR(NSPEC)
            ELSE
              ERRSTR = 'molar mass data format incorrect'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF

C           SPECIES LEWIS NUMBER
            VALIDS = CHKSTR(PARSTR(4),LENSTR(4),ICHRNO)
            IF(VALIDS)THEN
              READ(PARSTR(4)(1:LENSTR(4)),*)CLEWIS(NSPEC)
            ELSE
              ERRSTR = 'Lewis no. data format incorrect'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF

          ENDIF

C         READ ANOTHER LINE
          GOTO 2000

        ENDIF

C     END OF LOOP 2000

C     NSPEC NOW CONTAINS THE NUMBER OF SPECIES

C     =========================================================================

C     FILE ERROR HANDLER
      GOTO 1010
1000  WRITE(NCREPT,*)'Error in species data file:'
      WRITE(NCREPT,*)FNSPEC
      WRITE(NCREPT,*)'Premature end of file'
      WRITE(NCREPT,*)'- exiting from PPDIFF'
      WRITE(NCREPT,*)
      STOP
1010  CONTINUE
C     END OF FILE ERROR HANDLER

C     =========================================================================

C     READ THE LAST LINE
C     ------------------
      READ(NCSPEC,'(A)')ALINE
     
C     CLOSE THE SPECIES DATA FILE
C     ---------------------------
      CLOSE(NCSPEC)

C     =========================================================================

C     DIAGNOSTICS
      WRITE(6,*)
C      WRITE(6,*)'PDSPEC:'
C      WRITE(6,'(I5)')NSPEC
C      DO ISPEC = 1, NSPEC
C        WRITE(6,'(I5,2X,A10,2F5.2)')NUMSPC(ISPEC),
C     +    SPCSYM(ISPEC),WMOLAR(ISPEC),CLEWIS(ISPEC)
C      ENDDO

C     =========================================================================


      RETURN
      END
