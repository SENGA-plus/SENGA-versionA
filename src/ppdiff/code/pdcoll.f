      SUBROUTINE PDCOLL
 
C     *************************************************************************
C
C     PDCOLL
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
C     PREPROCESSES COLLISION INTEGRAL DATA FOR PPDIFF
C     FIXED FILE FORMAT: SET BY MONCHICK AND MASON DATA
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
      INTEGER NPARSE
      INTEGER ISTR,JSTR,ISTOTL
      INTEGER IDELTA,IRTEMP
 
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR

      LOGICAL VALIDS


C     BEGIN
C     =====

C     DIAGNOSTICS
      WRITE(6,*)'PDCOLL: raw collision integral data'
      WRITE(6,*)

C     =========================================================================

C     READ THE COLLISION INTEGRAL DATA
C     --------------------------------
      OPEN(UNIT=NCCOLL,FILE=FNCOLL,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1010)

C     FIVE LINES OF COMMENT ARE IGNORED
C     ---------------------------------
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE

C     SEPARATOR LINES
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE

C     =========================================================================

C     OMEGA(1,1)*
C     ----------

C     IGNORE THE BLOCK HEADER
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE

C     FIRST LINE OF BLOCK: READ LINE AND INITIALISE PARSER
1000  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 1000
        ENDIF
C     END OF LOOP 1000

C     CHECK FIRST LINE OF BLOCK
      JSTR = 0
      IDELTA = 0
1100  CONTINUE
        JSTR = JSTR + 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN
          IDELTA = IDELTA + 1
          READ(PARSTR(JSTR),*)DELR11(IDELTA)
          IF(JSTR.LT.ISTR)GOTO 1100
        ELSE
          ERRSTR = 'invalid collision data format'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF
        NDEL11 = IDELTA
C     END OF LOOP 1100

C     REST OF BLOCK: READ LINE AND INITIALISE PARSER
      IRTEMP = 0
1200  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 1200
        ENDIF

C       CHECK FIRST ITEM
        JSTR = 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN

C         CHECK NUMBER OF ITEMS
          IF(ISTR.NE.NDEL11+1)THEN
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         CHECK REMAINING ITEMS
          IRTEMP = IRTEMP + 1
          READ(PARSTR(JSTR),*)TRED11(IRTEMP)
          IDELTA = 0
1210      CONTINUE
            JSTR = JSTR + 1
            VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
            IF(VALIDS)THEN
              IDELTA = IDELTA + 1
              READ(PARSTR(JSTR),*)OMEG11(IRTEMP,IDELTA)
              IF(JSTR.LT.ISTR)GOTO 1210
            ELSE
              ERRSTR = 'invalid collision data format'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
C         END OF LOOP 1210

          IF(IRTEMP.LT.NTRDMX)THEN
            GOTO 1200
          ELSE
            ERRSTR = 'too much collision data'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ELSE

C         CHECK FOR END OF BLOCK
          IF(PARSTR(1).EQ.'#')THEN
            WRITE(6,*)'End of block detected'
C           CHECK THAT SOME DATA HAS BEEN COLLECTED
            IF(IRTEMP.EQ.0)THEN
              ERRSTR = 'end of collision data: data incomplete'
              CALL ERHAND(ERRSTR,IEFATL)
            ELSE
              NTRD11 = IRTEMP
            ENDIF
          ELSE
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ENDIF
C     END OF LOOP 1100

C     =========================================================================

C     OMEGA(2,2)*
C     ----------

C     IGNORE THE BLOCK HEADER
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE

C     FIRST LINE OF BLOCK: READ LINE AND INITIALISE PARSER
2000  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 2000
        ENDIF
C     END OF LOOP 2000

C     CHECK FIRST LINE OF BLOCK
      JSTR = 0
      IDELTA = 0
2100  CONTINUE
        JSTR = JSTR + 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN
          IDELTA = IDELTA + 1
          READ(PARSTR(JSTR),*)DELR22(IDELTA)
          IF(JSTR.LT.ISTR)GOTO 2100
        ELSE
          ERRSTR = 'invalid collision data format'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF
        NDEL22 = IDELTA
C     END OF LOOP 2100

C     REST OF BLOCK: READ LINE AND INITIALISE PARSER
      IRTEMP = 0
2200  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 2200
        ENDIF

C       CHECK FIRST ITEM
        JSTR = 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN

C         CHECK NUMBER OF ITEMS
          IF(ISTR.NE.NDEL22+1)THEN
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         CHECK REMAINING ITEMS
          IRTEMP = IRTEMP + 1
          READ(PARSTR(JSTR),*)TRED22(IRTEMP)
          IDELTA = 0
2210      CONTINUE
            JSTR = JSTR + 1
            VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
            IF(VALIDS)THEN
              IDELTA = IDELTA + 1
              READ(PARSTR(JSTR),*)OMEG22(IRTEMP,IDELTA)
              IF(JSTR.LT.ISTR)GOTO 2210
            ELSE
              ERRSTR = 'invalid collision data format'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
C         END OF LOOP 2210

          IF(IRTEMP.LT.NTRDMX)THEN
            GOTO 2200
          ELSE
            ERRSTR = 'too much collision data'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ELSE

C         CHECK FOR END OF BLOCK
          IF(PARSTR(1).EQ.'#')THEN
            WRITE(6,*)'End of block detected'
C           CHECK THAT SOME DATA HAS BEEN COLLECTED
            IF(IRTEMP.EQ.0)THEN
              ERRSTR = 'end of collision data: data incomplete'
              CALL ERHAND(ERRSTR,IEFATL)
            ELSE
              NTRD22 = IRTEMP
            ENDIF
          ELSE
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ENDIF
C     END OF LOOP 2100

C     =========================================================================

C     A*
C     -

C     IGNORE THE BLOCK HEADER
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE

C     FIRST LINE OF BLOCK: READ LINE AND INITIALISE PARSER
3000  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 3000
        ENDIF
C     END OF LOOP 3000

C     CHECK FIRST LINE OF BLOCK
      JSTR = 0
      IDELTA = 0
3100  CONTINUE
        JSTR = JSTR + 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN
          IDELTA = IDELTA + 1
          READ(PARSTR(JSTR),*)DELRAA(IDELTA)
          IF(JSTR.LT.ISTR)GOTO 3100
        ELSE
          ERRSTR = 'invalid collision data format'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF
        NDELAA = IDELTA
C     END OF LOOP 3100

C     REST OF BLOCK: READ LINE AND INITIALISE PARSER
      IRTEMP = 0
3200  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 3200
        ENDIF

C       CHECK FIRST ITEM
        JSTR = 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN

C         CHECK NUMBER OF ITEMS
          IF(ISTR.NE.NDELAA+1)THEN
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         CHECK REMAINING ITEMS
          IRTEMP = IRTEMP + 1
          READ(PARSTR(JSTR),*)TREDAA(IRTEMP)
          IDELTA = 0
3210      CONTINUE
            JSTR = JSTR + 1
            VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
            IF(VALIDS)THEN
              IDELTA = IDELTA + 1
              READ(PARSTR(JSTR),*)ASTAR(IRTEMP,IDELTA)
              IF(JSTR.LT.ISTR)GOTO 3210
            ELSE
              ERRSTR = 'invalid collision data format'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
C         END OF LOOP 3210

          IF(IRTEMP.LT.NTRDMX)THEN
            GOTO 3200
          ELSE
            ERRSTR = 'too much collision data'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ELSE

C         CHECK FOR END OF BLOCK
          IF(PARSTR(1).EQ.'#')THEN
            WRITE(6,*)'End of block detected'
C           CHECK THAT SOME DATA HAS BEEN COLLECTED
            IF(IRTEMP.EQ.0)THEN
              ERRSTR = 'end of collision data: data incomplete'
              CALL ERHAND(ERRSTR,IEFATL)
            ELSE
              NTRDAA = IRTEMP
            ENDIF
          ELSE
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ENDIF
C     END OF LOOP 3100

C     =========================================================================

C     B*
C     -

C     IGNORE THE BLOCK HEADER
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE

C     FIRST LINE OF BLOCK: READ LINE AND INITIALISE PARSER
4000  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 4000
        ENDIF
C     END OF LOOP 4000

C     CHECK FIRST LINE OF BLOCK
      JSTR = 0
      IDELTA = 0
4100  CONTINUE
        JSTR = JSTR + 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN
          IDELTA = IDELTA + 1
          READ(PARSTR(JSTR),*)DELRBB(IDELTA)
          IF(JSTR.LT.ISTR)GOTO 4100
        ELSE
          ERRSTR = 'invalid collision data format'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF
        NDELBB = IDELTA
C     END OF LOOP 4100

C     REST OF BLOCK: READ LINE AND INITIALISE PARSER
      IRTEMP = 0
4200  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 4200
        ENDIF

C       CHECK FIRST ITEM
        JSTR = 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN

C         CHECK NUMBER OF ITEMS
          IF(ISTR.NE.NDELBB+1)THEN
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         CHECK REMAINING ITEMS
          IRTEMP = IRTEMP + 1
          READ(PARSTR(JSTR),*)TREDBB(IRTEMP)
          IDELTA = 0
4210      CONTINUE
            JSTR = JSTR + 1
            VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
            IF(VALIDS)THEN
              IDELTA = IDELTA + 1
              READ(PARSTR(JSTR),*)BSTAR(IRTEMP,IDELTA)
              IF(JSTR.LT.ISTR)GOTO 4210
            ELSE
              ERRSTR = 'invalid collision data format'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
C         END OF LOOP 4210

          IF(IRTEMP.LT.NTRDMX)THEN
            GOTO 4200
          ELSE
            ERRSTR = 'too much collision data'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ELSE

C         CHECK FOR END OF BLOCK
          IF(PARSTR(1).EQ.'#')THEN
            WRITE(6,*)'End of block detected'
C           CHECK THAT SOME DATA HAS BEEN COLLECTED
            IF(IRTEMP.EQ.0)THEN
              ERRSTR = 'end of collision data: data incomplete'
              CALL ERHAND(ERRSTR,IEFATL)
            ELSE
              NTRDBB = IRTEMP
            ENDIF
          ELSE
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ENDIF
C     END OF LOOP 4100

C     =========================================================================

C     C*
C     -

C     IGNORE THE BLOCK HEADER
      READ(NCCOLL,'(A)',END=1020)ALINE
      READ(NCCOLL,'(A)',END=1020)ALINE

C     FIRST LINE OF BLOCK: READ LINE AND INITIALISE PARSER
5000  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 5000
        ENDIF
C     END OF LOOP 5000

C     CHECK FIRST LINE OF BLOCK
      JSTR = 0
      IDELTA = 0
5100  CONTINUE
        JSTR = JSTR + 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN
          IDELTA = IDELTA + 1
          READ(PARSTR(JSTR),*)DELRCC(IDELTA)
          IF(JSTR.LT.ISTR)GOTO 5100
        ELSE
          ERRSTR = 'invalid collision data format'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF
        NDELCC = IDELTA
C     END OF LOOP 5100

C     REST OF BLOCK: READ LINE AND INITIALISE PARSER
      IRTEMP = 0
5200  CONTINUE
        READ(NCCOLL,'(A)',END=1020)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE
        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
        IF(ISTOTL.EQ.0)THEN
          ERRSTR = 'blank line in collision data'
          CALL ERHAND(ERRSTR,IEWARN)
          GOTO 5200
        ENDIF

C       CHECK FIRST ITEM
        JSTR = 1
        VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
        IF(VALIDS)THEN

C         CHECK NUMBER OF ITEMS
          IF(ISTR.NE.NDELCC+1)THEN
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

C         CHECK REMAINING ITEMS
          IRTEMP = IRTEMP + 1
          READ(PARSTR(JSTR),*)TREDCC(IRTEMP)
          IDELTA = 0
5210      CONTINUE
            JSTR = JSTR + 1
            VALIDS = CHKSTR(PARSTR(JSTR),LENSTR(JSTR),ICHRNO)
            IF(VALIDS)THEN
              IDELTA = IDELTA + 1
              READ(PARSTR(JSTR),*)CSTAR(IRTEMP,IDELTA)
              IF(JSTR.LT.ISTR)GOTO 5210
            ELSE
              ERRSTR = 'invalid collision data format'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
C         END OF LOOP 5210

          IF(IRTEMP.LT.NTRDMX)THEN
            GOTO 5200
          ELSE
            ERRSTR = 'too much collision data'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ELSE

C         CHECK FOR END OF BLOCK
          IF(PARSTR(1).EQ.'#')THEN
            WRITE(6,*)'End of block detected'
C           CHECK THAT SOME DATA HAS BEEN COLLECTED
            IF(IRTEMP.EQ.0)THEN
              ERRSTR = 'end of collision data: data incomplete'
              CALL ERHAND(ERRSTR,IEFATL)
            ELSE
              NTRDCC = IRTEMP
            ENDIF
          ELSE
            ERRSTR = 'invalid collision data format'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF

        ENDIF
C     END OF LOOP 5100

C     =========================================================================

C     CLOSE THE COLLISION INTEGRAL DATA FILE
C     --------------------------------------
      CLOSE(NCCOLL)

C     =========================================================================

C     FILE ERROR HANDLER
      GOTO 1030
1010  ERRSTR = 'Error opening collison integral data file:'
      CALL ERHAND(ERRSTR,IEFATL)
1020  ERRSTR = 'unexpected end of collison integral data file'
      CALL ERHAND(ERRSTR,IEFATL)
1030  CONTINUE
C     END OF FILE ERROR HANDLER

C     =========================================================================

C     DIAGNOSTICS
      WRITE(6,*)
      WRITE(6,*)'PDCOLL: replay collision integral data'
      WRITE(6,*)

      WRITE(6,*)'Omega11:'
      WRITE(6,'(10X,8F8.3)')(DELR11(IDELTA),IDELTA =1,NDEL11)
      DO IRTEMP = 1, NTRD11
        WRITE(6,'(F10.4,8F8.4)')TRED11(IRTEMP),
     +                         (OMEG11(IRTEMP,IDELTA),IDELTA =1,NDEL11)
      ENDDO
      WRITE(6,*)

      WRITE(6,*)'Omega22:'
      WRITE(6,'(10X,8F8.3)')(DELR22(IDELTA),IDELTA =1,NDEL22)
      DO IRTEMP = 1, NTRD22
        WRITE(6,'(F10.4,8F8.4)')TRED22(IRTEMP),
     +                         (OMEG22(IRTEMP,IDELTA),IDELTA =1,NDEL22)
      ENDDO
      WRITE(6,*)

      WRITE(6,*)'Astar:'
      WRITE(6,'(10X,8F8.3)')(DELRAA(IDELTA),IDELTA =1,NDELAA)
      DO IRTEMP = 1, NTRDAA
        WRITE(6,'(F10.4,8F8.4)')TREDAA(IRTEMP),
     +                         (ASTAR(IRTEMP,IDELTA),IDELTA =1,NDELAA)
      ENDDO
      WRITE(6,*)

      WRITE(6,*)'Bstar:'
      WRITE(6,'(10X,8F8.3)')(DELRBB(IDELTA),IDELTA =1,NDELBB)
      DO IRTEMP = 1, NTRDBB
        WRITE(6,'(F10.4,8F8.4)')TREDBB(IRTEMP),
     +                         (BSTAR(IRTEMP,IDELTA),IDELTA =1,NDELBB)
      ENDDO
      WRITE(6,*)

      WRITE(6,*)'Cstar:'
      WRITE(6,'(10X,8F8.3)')(DELRCC(IDELTA),IDELTA =1,NDELCC)
      DO IRTEMP = 1, NTRDCC
        WRITE(6,'(F10.4,8F8.4)')TREDCC(IRTEMP),
     +                         (CSTAR(IRTEMP,IDELTA),IDELTA =1,NDELCC)
      ENDDO
      WRITE(6,*)

C     =========================================================================


      RETURN
      END
