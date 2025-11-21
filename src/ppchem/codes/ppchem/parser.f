      SUBROUTINE PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
 
C     *************************************************************************
C
C     PARSER
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     18-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     PARSES A LINE OF DATA FOR PPCHEM
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      CHARACTER*50 BLANKS
      PARAMETER(
     +   BLANKS='                                                  ')


C     ARGUMENTS
C     =========
      INTEGER NPARSE,NUMSTR,ISTR,ISTOTL
      INTEGER LENSTR(NUMSTR)
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)


C     LOCAL DATA
C     ==========
      INTEGER IPARSE,ILEN
      CHARACTER*1 PARCHR
      LOGICAL ONWORD


C     BEGIN
C     =====

C     INITIALISE PARSER
      ONWORD = .FALSE.
      IPARSE = 1
      ISTOTL = 0
      ISTR = 0
      ILEN = 0
1000  CONTINUE

C       PICK THE NEXT CHARACTER
        PARCHR = ALINE(IPARSE:IPARSE)

        IF(ONWORD)THEN

C         CURRENTLY BUILDING A STRING
          IF(PARCHR.EQ.' ')THEN

C           STRING IS FINISHED
            ONWORD = .FALSE.
            LENSTR(ISTR) = ILEN
            ISTOTL = ISTOTL + LENSTR(ISTR)
 
          ELSE

C           STRING IS NOT FINISHED
            ILEN = ILEN + 1
            PARSTR(ISTR)(ILEN:ILEN) = PARCHR

          ENDIF

        ELSE

C         NOT CURRENTLY BUILDING A STRING
          IF(PARCHR.EQ.' ')THEN

C           STRIP THE BLANK
            CONTINUE

          ELSE

C           BEGIN A NEW STRING
            ONWORD = .TRUE.
            ISTR = ISTR + 1
            ILEN = 1
            PARSTR(ISTR) = BLANKS
            PARSTR(ISTR)(1:1) = PARCHR

          ENDIF

        ENDIF

C       PREPARE FOR NEXT CHARACTER
        IPARSE = IPARSE + 1
        IF(IPARSE.LE.NPARSE)THEN
          GOTO 1000
        ELSE
C         END OF LINE - DROP THROUGH
          CONTINUE
        ENDIF
C             
C     END OF LOOP 1000


      RETURN
      END
