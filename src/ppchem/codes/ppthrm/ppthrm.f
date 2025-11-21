      PROGRAM PPTHRM
 
C     *************************************************************************
C
C     PPTHRM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     10-SEP-2002:  CREATED
C     13-SEP-2009:  RSC UPDATES TO INPUT FORMATS
C      
C     DESCRIPTION
C     -----------
C     PREPROCESSES THERMODYNAMIC DATA FOR SENGA
C     USING A SPECIES FILE
C     AND A THERMODYNAMIC DATA SOURCE FILE IN CHEMKIN FORMAT
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ZERO, ONE
      PARAMETER(ZERO=0.0D0, ONE=1.0D0)
      DOUBLE PRECISION TINTDF
      PARAMETER(TINTDF=1.0D3)
      INTEGER NCINPT,NCREPT,NCDATA,NCOUTP
      PARAMETER(NCINPT=5, NCREPT=6, NCDATA=1, NCOUTP=2)
      INTEGER NUMSTR,LENMAX,LENSMX
      PARAMETER(NUMSTR=12, LENMAX=50, LENSMX=16)
      INTEGER NSPCMX,NSTPMX,NBDYMX
      PARAMETER(NSPCMX=100, NSTPMX=1000, NBDYMX=10)
      INTEGER NGIBMX,NLNDMX
      PARAMETER(NGIBMX=100, NLNDMX=10)
      INTEGER NLSMAX,NSSMAX
      PARAMETER(NLSMAX=20, NSSMAX=20)
      INTEGER NLGMAX,NLLMAX
      PARAMETER(NLGMAX=20, NLLMAX=20)
      INTEGER NTINMX,NCOFMX
      PARAMETER(NTINMX=2,NCOFMX=7)
      INTEGER IENULL,IEINFO,IEWARN,IERROR,IESEVR,IEFATL
      PARAMETER(IENULL=0,
     +          IEINFO=1, IEWARN=2, IERROR=3, IESEVR=4, IEFATL=5)
      INTEGER ICHNUL,ICHSTR,ICHEND,ICHINT,ICHSNT,ICHRNO
      PARAMETER(ICHNUL=0,
     +          ICHSTR=1, ICHEND=2, ICHINT=3, ICHSNT=4, ICHRNO=5)
      CHARACTER*50 BLANKS
      PARAMETER(
     +   BLANKS='                                                  ')


C     EXTERNAL FUNCTION
C     =================
      LOGICAL CHKSTR
      EXTERNAL CHKSTR


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ACOFCP(NCOFMX,NTINMX,NSPCMX)
      DOUBLE PRECISION TINTLO(NTINMX),TINTHI(NTINMX)
      DOUBLE PRECISION WMOLAR(NSPCMX),CLEWIS(NSPCMX)
      DOUBLE PRECISION PREFGB
      DOUBLE PRECISION CPOLY1,CPOLY2
      INTEGER LENSYM(NSPCMX)
      INTEGER NUMSPC(NSPCMX)
      INTEGER LENSTR(NUMSTR)
      INTEGER NCOEFF(NTINMX)
      INTEGER ISPEC,NSPEC
      INTEGER NPARSE
      INTEGER ISTR,JSTR,ISTOTL
      INTEGER LENCHK
      INTEGER ITINT,NTINT
      INTEGER IC
 
      CHARACTER*16 SPCSYM(NSPCMX)
      CHARACTER*80 FNSPEC,FNTHRM,FNOUTP
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR
      CHARACTER*16 SYMCHK
      CHARACTER*18 STR1
      CHARACTER*6  STR2
      CHARACTER*20 STR3
      CHARACTER*1  STR4
      CHARACTER*10 STR5
      CHARACTER*10 STR6
      CHARACTER*8  STR7
      CHARACTER*5  STR8
      CHARACTER*1  STR9

      LOGICAL FSPCNF(NSPCMX)
      LOGICAL VALIDS


C     BEGIN
C     =====

C     =========================================================================

C     ANNOUNCEMENT
C     ------------
      WRITE(NCREPT,*)'==========================================='
      WRITE(NCREPT,*)
      WRITE(NCREPT,*)'PPTHRM: thermo data preprocessor for PPCHEM'
      WRITE(NCREPT,*)
      WRITE(NCREPT,*)'==========================================='
      WRITE(NCREPT,*)

C     =========================================================================

C     GET THE FILENAMES AND CHECK THE FILES
C     -------------------------------------
      WRITE(NCREPT,*)'Please enter the name of the species data file:'
      READ(NCINPT,'(A)')FNSPEC
      WRITE(NCREPT,'(A)')FNSPEC
      WRITE(NCREPT,*)
      OPEN(UNIT=NCDATA,FILE=FNSPEC,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1000)
      CLOSE(NCDATA)

      WRITE(NCREPT,*)
     +      'Please enter the name of the thermo data source file:'
      READ(NCINPT,'(A)')FNTHRM
      WRITE(NCREPT,'(A)')FNTHRM
      WRITE(NCREPT,*)
      OPEN(UNIT=NCDATA,FILE=FNTHRM,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1010)
      CLOSE(NCDATA)

C     FILENAME ERROR HANDLER
      GOTO 1100
1000  WRITE(NCREPT,*)'Error opening species data file:'
      WRITE(NCREPT,*)FNSPEC
      WRITE(NCREPT,*)'- exiting from PPTHRM'
      WRITE(NCREPT,*)
      STOP
1010  WRITE(NCREPT,*)'Error opening thermo data source file:'
      WRITE(NCREPT,*)FNTHRM
      WRITE(NCREPT,*)'- exiting from PPTHRM'
      WRITE(NCREPT,*)
      STOP
1100  CONTINUE
C     END OF FILENAME ERROR HANDLER

      WRITE(NCREPT,*)'Please enter the reference pressure:'
      READ(NCINPT,*)PREFGB
      WRITE(NCREPT,'(1PE12.4)')PREFGB
      WRITE(NCREPT,*)

      WRITE(NCREPT,*)'Please enter the name of the output file:'
      READ(NCINPT,'(A)')FNOUTP
      WRITE(NCREPT,'(A)')FNOUTP
      WRITE(NCREPT,*)
      OPEN(UNIT=NCOUTP,FILE=FNOUTP,STATUS='UNKNOWN',FORM='FORMATTED',
     +     ERR=1020)

C     FILENAME ERROR HANDLER
      GOTO 1300
1020  WRITE(NCREPT,*)'Error opening output file:'
      WRITE(NCREPT,*)FNOUTP
      WRITE(NCREPT,*)'- exiting from PPCHEM'
      WRITE(NCREPT,*)
      STOP
1300  CONTINUE
C     END OF FILENAME ERROR HANDLER

C     HEADER
      WRITE(NCOUTP,*)'*****************************'
      WRITE(NCOUTP,*)'*                           *'
      WRITE(NCOUTP,*)'*  Output file from ppthrm  *'
      WRITE(NCOUTP,*)'*                           *'
      WRITE(NCOUTP,*)'*****************************'

C     REFERENCE PRESSURE
      WRITE(NCOUTP,'(1PE12.4)')PREFGB

C     =========================================================================

C     READ THE SPECIES DATA FILE
C     --------------------------
      OPEN(UNIT=NCDATA,FILE=FNSPEC,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1000)

C     FIVE LINES OF COMMENT ARE IGNORED
C     ---------------------------------
      READ(NCDATA,'(A)',END=1200)ALINE
      READ(NCDATA,'(A)',END=1200)ALINE
      READ(NCDATA,'(A)',END=1200)ALINE
      READ(NCDATA,'(A)',END=1200)ALINE
      READ(NCDATA,'(A)',END=1200)ALINE

C     =========================================================================

C     READ THE SPECIES LIST
C     ---------------------
C     PREINITIALISE FIRST SYMBOL

C     READ AND IGNORE HEADER LINE
      READ(NCDATA,'(A)',END=1200)ALINE
      WRITE(NCREPT,'(A)')ALINE

C     READ AND PARSE EACH LINE OF THE SPECIES LIST
      NSPEC = 0
      LENCHK = 0
2010  CONTINUE

C       INCREMENT SPECIES COUNTER AND CHECK
        NSPEC = NSPEC + 1
        IF(NSPEC.GT.NSPCMX)THEN
          ERRSTR = 'max no. of species exceeded'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF

C       READ LINE AND INITIALISE PARSER
        READ(NCDATA,'(A)',END=1200)ALINE
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

C           SET SPECIES FLAG
            FSPCNF(NSPEC) = .TRUE.

          ENDIF

C         READ ANOTHER LINE
          GOTO 2010

        ENDIF

C     END OF LOOP 2010

C     NSPEC NOW CONTAINS THE NUMBER OF SPECIES

C     =========================================================================

C     END-OF-FILE ERROR HANDLER
      GOTO 1210
1200  WRITE(NCREPT,*)'Error in species data file:'
      WRITE(NCREPT,*)FNSPEC
      WRITE(NCREPT,*)'Premature end of file'
      WRITE(NCREPT,*)'- exiting from PPCHEM'
      WRITE(NCREPT,*)
      STOP
1210  CONTINUE
C     END OF END-OF-FILE ERROR HANDLER

C     =========================================================================

C     READ THE LAST LINE
C     ------------------
      READ(NCDATA,'(A)')ALINE
     
C     CLOSE THE SPECIES DATA FILE
C     ---------------------------
      CLOSE(NCDATA)

C     =========================================================================

C     READ THE THERMO SOURCE FILE
C     ---------------------------
      OPEN(UNIT=NCDATA,FILE=FNTHRM,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1010)

C     =========================================================================

3000  CONTINUE

C       READ A LINE
        READ(NCDATA,'(A)',END=3010)ALINE
        NPARSE = LEN(ALINE)

        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)

C       CHECK FOR END OF FILE
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)

        IF(VALIDS)THEN

C         LAST LINE
          WRITE(NCREPT,*)'End of thermo data source file'

        ELSE

C         ASSUME DATA LINE AND INTERPRET USING CHEMKIN FIXED FORMAT
          STR1 = ALINE(1:18)
          STR2 = ALINE(19:24)
          STR3 = ALINE(25:44)
          STR4 = ALINE(45:45)
          STR5 = ALINE(46:55)
          STR6 = ALINE(56:65)
          STR7 = ALINE(66:73)
          STR8 = ALINE(74:78)
          STR9 = ALINE(80:80)

C         CHECK LINE NUMBER
          IF(STR9.NE.'1')THEN
            ERRSTR = 'thermo data source format incorrect'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         GET SPECIES NAME FROM FIXED FORMAT DATA
          ISTR = 1
3110      CONTINUE
            IF(STR1(ISTR:ISTR).EQ.' ')THEN

C             STRING IS FINISHED
              LENCHK = ISTR - 1

            ELSE

C             STRING IS NOT FINISHED
              IF(ISTR.LT.LENSMX)THEN
                SYMCHK(ISTR:ISTR) = STR1(ISTR:ISTR)
                ISTR = ISTR + 1
                GOTO 3110
              ELSE
                ERRSTR = 'thermo data source format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

            ENDIF

C         END OF LOOP 3110

C         COMPARE WITH SPECIES NAME FROM PARSER
          IF(SYMCHK(1:LENCHK).NE.PARSTR(1)(1:LENSTR(1)))THEN
            ERRSTR = 'thermo data source format incorrect'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         COMPARE WITH SPECIES LIST
          ISPEC = 1
3120      CONTINUE

            IF(SYMCHK(1:LENCHK).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN

C             SPECIES FOUND - READ AND PROCESS THE DATA
C             DIAGNOSTICS
              WRITE(6,*)'Species found'
              WRITE(6,'(A)')SYMCHK(1:LENCHK)
              FSPCNF(ISPEC) = .FALSE.

C             FIRST LINE: TEMPERATURE INTERVALS FOR CP DATA
              NTINT = 2
              DO ITINT = 1,NTINT
                NCOEFF(ITINT) = 7
              ENDDO
              READ(STR5,*)TINTLO(1)
              READ(STR6,*)TINTHI(2)
              IF(STR7.NE.'        ')THEN
                READ(STR7,*)TINTHI(1)
              ELSE
                TINTHI(1) = TINTDF
              ENDIF
              TINTLO(2) = TINTHI(1)

C             SECOND LINE
              READ(NCDATA,'(A)',END=3030)ALINE                 
              IF(ALINE(80:80).NE.'2')THEN
                ERRSTR = 'thermo data source format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
C             RSC 13-SEP-2009 UPDATE TO INPUT FORMAT
C              READ(ALINE,*)(ACOFCP(IC,2,ISPEC),IC=1,5)
              READ(ALINE,'(5(1PE15.8))')(ACOFCP(IC,2,ISPEC),IC=1,5)

C             THIRD LINE
              READ(NCDATA,'(A)',END=3030)ALINE                 
              IF(ALINE(80:80).NE.'3')THEN
                ERRSTR = 'thermo data source format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
C             RSC 13-SEP-2009 UPDATE TO INPUT FORMAT
C              READ(ALINE,*)ACOFCP(6,2,ISPEC),ACOFCP(7,2,ISPEC),
              READ(ALINE,'(5(1PE15.8))')
     +                     ACOFCP(6,2,ISPEC),ACOFCP(7,2,ISPEC),
     +                    (ACOFCP(IC,1,ISPEC),IC=1,3)

C             FOURTH LINE
              READ(NCDATA,'(A)',END=3030)ALINE                 
              IF(ALINE(80:80).NE.'4')THEN
                ERRSTR = 'thermo data source format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
C             RSC 13-SEP-2009 UPDATE TO INPUT FORMAT
C              READ(ALINE,*)(ACOFCP(IC,1,ISPEC),IC=4,7)
              READ(ALINE,'(5(1PE15.8))')(ACOFCP(IC,1,ISPEC),IC=4,7)

C             CHECK CONTINUITY OF THE CP POLYNOMIAL
              CPOLY1 = ACOFCP(5,1,ISPEC)
              DO IC = 1,4
                CPOLY1 = ACOFCP(5-IC,1,ISPEC) + CPOLY1*TINTHI(1)
              ENDDO
              CPOLY2 = ACOFCP(5,2,ISPEC)
              DO IC = 1,4
                CPOLY2 = ACOFCP(5-IC,2,ISPEC) + CPOLY2*TINTLO(2)
              ENDDO
              WRITE(6,'(2(1PE15.7))')CPOLY1,CPOLY2

C             CHECK CONTINUITY OF ENTHALPY
              CPOLY1 = ACOFCP(5,1,ISPEC)/5.0D0
              DO IC = 1,4
                CPOLY1 = ACOFCP(5-IC,1,ISPEC)/REAL(5-IC)
     +                 + CPOLY1*TINTHI(1)
              ENDDO
              CPOLY1 = CPOLY1*TINTHI(1) + ACOFCP(6,1,ISPEC)
              CPOLY2 = ACOFCP(5,2,ISPEC)/5.0D0
              DO IC = 1,4
                CPOLY2 = ACOFCP(5-IC,2,ISPEC)/REAL(5-IC)
     +                 + CPOLY2*TINTLO(2)
              ENDDO
              CPOLY2 = CPOLY2*TINTLO(2) + ACOFCP(6,2,ISPEC)
              WRITE(6,'(2(1PE15.7))')CPOLY1,CPOLY2

C             CHECK CONTINUITY OF ENTROPY
              CPOLY1 = ACOFCP(5,1,ISPEC)/4.0D0
              DO IC = 1,3
                CPOLY1 = ACOFCP(5-IC,1,ISPEC)/REAL(4-IC)
     +                 + CPOLY1*TINTHI(1)
              ENDDO
              CPOLY1 = CPOLY1*TINTHI(1) + ACOFCP(7,1,ISPEC)
     +               + ACOFCP(1,1,ISPEC)*LOG(TINTHI(1))
              CPOLY2 = ACOFCP(5,2,ISPEC)/4.0D0
              DO IC = 1,3
                CPOLY2 = ACOFCP(5-IC,2,ISPEC)/REAL(4-IC)
     +                 + CPOLY2*TINTLO(2)
              ENDDO
              CPOLY2 = CPOLY2*TINTLO(2) + ACOFCP(7,2,ISPEC)
     +               + ACOFCP(1,2,ISPEC)*LOG(TINTLO(2))
              WRITE(6,'(2(1PE15.7))')CPOLY1,CPOLY2

C             WRITE OUT THE CP DATA
C             ---------------------
              WRITE(NCOUTP,'(4X,A)')SPCSYM(ISPEC)(1:LENSYM(ISPEC))
              WRITE(NCOUTP,'(4X,1PE12.4)')WMOLAR(ISPEC)
              WRITE(NCOUTP,'(4X,1PE12.4)')CLEWIS(ISPEC)
              WRITE(NCOUTP,'(I5)')NTINT
              DO ITINT = 1,NTINT
                WRITE(NCOUTP,'(I5,2F8.2,I3)')
     +            ITINT,TINTLO(ITINT),TINTHI(ITINT),NCOEFF(ITINT)
                DO IC = 1,NCOEFF(ITINT)
                  WRITE(NCOUTP,'(1PE18.7))')ACOFCP(IC,ITINT,ISPEC)
                ENDDO
              ENDDO
              WRITE(NCOUTP,*)'    END'

            ELSE

C             SPECIES NOT FOUND
              IF(ISPEC.LT.NSPEC)THEN

C               KEEP LOOKING
                ISPEC = ISPEC + 1
                GOTO 3120

              ELSE

C               SPECIES NOT FOUND - READ AND IGNORE THE DATA
                READ(NCDATA,'(A)')ALINE
                READ(NCDATA,'(A)')ALINE
                READ(NCDATA,'(A)')ALINE

              ENDIF

            ENDIF
C         END OF LOOP 3120

C         GO BACK AND READ ANOTHER FIRST LINE
          GOTO 3000

        ENDIF
C       LAST LINE

C     END OF LOOP 3000

C     END-OF-FILE CONDITION HANDLER
      GOTO 3020
3010  WRITE(NCREPT,*)'End of thermo data source file:'
      WRITE(NCREPT,*)FNTHRM
3020  CONTINUE
C     END OF END-OF-FILE CONDITION HANDLER

C     END-OF-FILE ERROR HANDLER
      GOTO 3040
3030  WRITE(NCREPT,*)'Error: unexpected end of thermo data source file:'
      WRITE(NCREPT,*)FNTHRM
      WRITE(NCREPT,*)'- exiting from PPTHRM'
      WRITE(NCREPT,*)
      STOP
3040  CONTINUE
C     END OF END-OF-FILE ERROR HANDLER

C     SPECIES NOT PROCESSED ERROR HANDLER
      DO ISPEC = 1, NSPEC
        IF(FSPCNF(ISPEC))THEN
        WRITE(NCREPT,*)'Warning: species not found in thermo database:'
          WRITE(NCREPT,'(I5,4X,A)')ISPEC,SPCSYM(ISPEC)(1:LENSYM(ISPEC))
        ENDIF
      ENDDO

C     =========================================================================

C     CLOSE THE THERMO DATA SOURCE FILE
C     ---------------------------------
      CLOSE(NCDATA)

C     CLOSE THE OUTPUT FILE
C     ---------------------
      WRITE(NCOUTP,*)'END'
      CLOSE(NCOUTP)

C     =========================================================================

      STOP
      END
