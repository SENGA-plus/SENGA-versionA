      SUBROUTINE OUTPUT
 
C     *************************************************************************
C
C     OUTPUT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT - CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     28-DEC-2003:  RSC MODIFIED FOR SENGA2
C     11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH; REPORT THE DUMP
C     29-AUG-2009:  RSC UPDATE NUMBER OF PROCESSORS
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     PROCESSES THE RESULTS
C
C     FILES: REPORT FILE       - UNIT NCREPT - NAME FNREPT  - FORMATTED
C            DUMP OUTPUT FILES - UNIT NCDMPO - NAMES FNDMPO - UNFORMATTED
C            STATISTICS FILE   - UNIT NCSTAT - NAME FNSTAT  - FORMATTED
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
#ifdef HDF5
      USE hdf5io
#endif
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL PARAMETERS
C     ================
C     DIAGNOSTICS
      CHARACTER*4 PNXRES
      PARAMETER(PNXRES = '.res')
      CHARACTER*2 PNXRYF
      PARAMETER(PNXRYF = 'yf')
      INTEGER NCDIAG
      PARAMETER(NCDIAG = 11)
      INTEGER IDDUMP

C     LOCAL DATA
C     ==========
C     DIAGNOSTICS
      DOUBLE PRECISION DELTAG,FORNOW
      DOUBLE PRECISION TTEMP(NXSIZE,NYSIZE,NZSIZE)
      DOUBLE PRECISION PTEMP(NXSIZE,NYSIZE,NZSIZE)
      DOUBLE PRECISION YTEMP(NXSIZE,NYSIZE,NZSIZE,NSPEC)

      INTEGER ISPEC
      INTEGER IC,JC,KC
      INTEGER IGOFST,IX
      CHARACTER*21 FNDIAG
C     RSC UPDATE NUMBER OF PROCESSORS
      CHARACTER*6 PNPROC
      CHARACTER*11 STRQTY
      CHARACTER*2 STRSPC
      LOGICAL BCFLAG

      CHARACTER*60 FNAME
      CHARACTER*4 PROC
      CHARACTER*5 IPDUMP

C     BEGIN
C     =====

C     =========================================================================

C     REPORT OUTPUT 
C     =============
      IF(MOD(ITIME,NTREPT).EQ.0)THEN

C       REPORT ON PROCESSOR NO.1 ONLY
C       ------
        IF(IPROC.EQ.0)THEN

          OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='OLD',FORM='FORMATTED')

C         GO TO EOF
1000      CONTINUE
            READ(NCREPT,9000,END=1010)
            GOTO 1000
1010      BACKSPACE(NCREPT)

          WRITE(NCREPT,9100)ITIME
          WRITE(NCREPT,9110)ETIME,TSTEP
          CLOSE(NCREPT)

        ENDIF

C       =======================================================================

C       DIAGNOSTICS
        JC = MAX(NYGLBL/2,1)
        KC = MAX(NZGLBL/2,1)

        WRITE(PNPROC,'(I6.6)')IPROC

C       GLOBAL INDEXING
        DELTAG = XGDLEN/(REAL(NXGLBL-1))

        IGOFST = 0
        DO IC = 0, IXPROC-1
          IGOFST = IGOFST + NPMAPX(IC)
        ENDDO

C        STRQTY = 'output/pres'
C        FNDIAG = STRQTY//PNPROC//PNXRES
C        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
CC       GO TO EOF
C8000    CONTINUE
C          READ(NCDIAG,9000,END=8001)
C          GOTO 8000
C8001    BACKSPACE(NCDIAG)
C        DO IC = ISTAL,ISTOL
C          IX = IGOFST + IC
C          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,PRUN(IC,JC,KC)
C        ENDDO
C        WRITE(NCDIAG,*)
C        CLOSE(NCDIAG)
C
C        STRQTY = 'output/uvel'
C        FNDIAG = STRQTY//PNPROC//PNXRES
C        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
CC       GO TO EOF
C8010    CONTINUE
C          READ(NCDIAG,9000,END=8011)
C          GOTO 8010
C8011    BACKSPACE(NCDIAG)
C        DO IC = ISTAL,ISTOL
C          IX = IGOFST + IC
C          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
C     +                      URUN(IC,JC,KC)/DRUN(IC,JC,KC)
C        ENDDO
C        WRITE(NCDIAG,*)
C        CLOSE(NCDIAG)
C
C        STRQTY = 'output/temp'
C        FNDIAG = STRQTY//PNPROC//PNXRES
C        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
CC       GO TO EOF
C8020    CONTINUE
C          READ(NCDIAG,9000,END=8021)
C          GOTO 8020
C8021    BACKSPACE(NCDIAG)
C        DO IC = ISTAL,ISTOL
C          IX = IGOFST + IC
C          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,TRUN(IC,JC,KC)
C        ENDDO
C        WRITE(NCDIAG,*)
C        CLOSE(NCDIAG)
C
C        STRQTY = 'output/dens'
C        FNDIAG = STRQTY//PNPROC//PNXRES
C        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
CC       GO TO EOF
C8030    CONTINUE
C          READ(NCDIAG,9000,END=8031)
C          GOTO 8030
C8031    BACKSPACE(NCDIAG)
C        DO IC = ISTAL,ISTOL
C          IX = IGOFST + IC
C          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,DRUN(IC,JC,KC)
C        ENDDO
C        WRITE(NCDIAG,*)
C        CLOSE(NCDIAG)
C
C        STRQTY = 'output/ener'
C        FNDIAG = STRQTY//PNPROC//PNXRES
C        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
CC       GO TO EOF
C8040    CONTINUE
C          READ(NCDIAG,9000,END=8041)
C          GOTO 8040
C8041    BACKSPACE(NCDIAG)
C        DO IC = ISTAL,ISTOL
C          IX = IGOFST + IC
C          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
C     +                      ERUN(IC,JC,KC)/DRUN(IC,JC,KC)
C        ENDDO
C        WRITE(NCDIAG,*)
C        CLOSE(NCDIAG)
C
C        DO ISPEC = 1, NSPEC
C
C          WRITE(STRSPC,'(I2.2)')ISPEC
C          STRQTY = 'output/'//PNXRYF//STRSPC
C          FNDIAG = STRQTY//PNPROC//PNXRES
C          OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
CC         GO TO EOF
C8050      CONTINUE
C            READ(NCDIAG,9000,END=8051)
C            GOTO 8050
C8051      BACKSPACE(NCDIAG)
C          DO IC = ISTAL,ISTOL
C            IX = IGOFST + IC
C              FORNOW = YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
C              FORNOW = MAX(FORNOW,1.0D-30)
C            WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,FORNOW
CC     +                        YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
C          ENDDO
C          WRITE(NCDIAG,*)
C          CLOSE(NCDIAG)
C
C        ENDDO
C
CC        ISPEC = NSPEC
CC        WRITE(STRSPC,'(I2.2)')ISPEC
CC        STRQTY = PNXRYF//STRSPC
CC        FNDIAG = STRQTY//PNPROC//PNXRES
CC        OPEN(UNIT=NCDIAG,FILE=FNDIAG,FORM='FORMATTED')
CCC       GO TO EOF
CC8050    CONTINUE
CC          READ(NCDIAG,9000,END=8051)
CC          GOTO 8050
CC8051    BACKSPACE(NCDIAG)
CC        DO IC = ISTAL,ISTOL
CC          IX = IGOFST + IC
CC          WRITE(NCDIAG,9300)REAL(IX-1)*DELTAG,
CC     +                      YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
CC        ENDDO
CC        WRITE(NCDIAG,*)
CC        CLOSE(NCDIAG)

      ENDIF

C     ========================================================================= 
C     UMOD START
C     DATA OUTPUT FOR POST-PROCESSING
C     ======
      IF(MOD(ITIME,NTDUMP).EQ.0)THEN
          IDDUMP=ITIME/NTDUMP
       WRITE(IPDUMP,'(I5.5)') IDDUMP
       WRITE(PROC,'(I4.4)') IPROC
 
       FNAME = 'output/out'//IPDUMP//PROC//PNXRES


        DO KC = 1,NZSIZE
         DO JC = 1,NYSIZE
            DO IC = 1,NXSIZE
              DO ISPEC =1,NSPEC
               YTEMP(IC,JC,KC,ISPEC)=YRUN(IC,JC,KC,ISPEC)/DRUN(IC,JC,KC)
               ENDDO 
               TTEMP(IC,JC,KC)=TRUN(IC,JC,KC)
               PTEMP(IC,JC,KC)=PRUN(IC,JC,KC)
            ENDDO
          ENDDO
        ENDDO

      OPEN(UNIT=16,FILE=TRIM(FNAME),FORM='UNFORMATTED',STATUS='NEW')
       WRITE(16)DRUN,URUN/DRUN,VRUN/DRUN,WRUN/DRUN,ERUN/DRUN,TTEMP,
     +          PTEMP,YTEMP,RRTE,ETIME
       CLOSE(16)
       ENDIF



C     UMOD END
C     =====
C     FULL DUMP OUTPUT 
C     ================
      IF(MOD(ITIME,NTDUMP).EQ.0)THEN

C       CARRY OUT A FULL DUMP
C       ---------------------
C       USE THE DUMP FILE INDICATED BY IDFLAG
C       RSC 11-JUL-2009 ADD A DUMP FORMAT SWITCH
#ifndef HDF5
        IF(NDOFMT.EQ.0)THEN
          IF (NXLPRM(1).EQ.4)THEN
            IF(IPROC.eq.0)THEN
             OPEN(UNIT=NCDMPO,FILE='output/intran.dat',STATUS='unknown',
     +           FORM='FORMATTED')
             WRITE(NCDMPO,*)INTRAN
             CLOSE(NCDMPO)
            ENDIF


            IF (IXPROC.EQ.0)THEN
              WRITE(FNAME,'(A,I4.4,A)')'output/inflow'//PNPROC//'.dat'
              OPEN(UNIT=NCDMPO,FILE=trim(FNAME),STATUS='unknown',
     +             FORM='FORMATTED')
              REWIND(NCDMPO)
              DO KC = 1,NZSIZE
                DO JC = 1,NYSIZE
                  WRITE(NCDMPO,*)UINF2(JC,KC),VINF2(JC,KC),WINF2(JC,KC),
     +                    (YINF2(JC,KC,ISPEC),ISPEC=1,NSPEC)
                ENDDO
              ENDDO
              CLOSE(NCDMPO)
            END IF
          ENDIF


C         UNFORMATTED DUMP OUTPUT
          OPEN(UNIT=NCDMPO,FILE=FNDMPO(IDFLAG+1),STATUS='OLD',
     +         FORM='UNFORMATTED')
          REWIND(NCDMPO)
          WRITE(NCDMPO)NXNODE,NYNODE,NZNODE,NSPEC,
     +                 DRUN,URUN,VRUN,WRUN,ERUN,YRUN,
     +                 ETIME,TSTEP,ERROLD,ERRLDR
          CLOSE(NCDMPO)

        ELSE

C         FORMATTED DUMP OUTPUT
          OPEN(UNIT=NCDMPO,FILE=FNDMPO(IDFLAG+1),STATUS='OLD',
     +         FORM='FORMATTED')
          REWIND(NCDMPO)
          WRITE(NCDMPO,*)NXNODE,NYNODE,NZNODE,NSPEC
          DO KC = 1,NZNODE
            DO JC = 1,NYNODE
              DO IC = 1,NXNODE
                WRITE(NCDMPO,*)DRUN(IC,JC,KC),
     +                     URUN(IC,JC,KC),VRUN(IC,JC,KC),WRUN(IC,JC,KC),
     +                         ERUN(IC,JC,KC),
     +                        (YRUN(IC,JC,KC,ISPEC),ISPEC=1,NSPEC)
              ENDDO
            ENDDO
          ENDDO
          WRITE(NCDMPO,*)ETIME,TSTEP,ERROLD,ERRLDR
          CLOSE(NCDMPO)

        ENDIF

#else
        CALL WRITE_H5_DUMPFILE
#endif

C       REPORT THE DUMP
C       RSC 11-JUL-2009
        IF(IPROC.EQ.0)THEN

          OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='OLD',FORM='FORMATTED')
3000      CONTINUE
            READ(NCREPT,9000,END=3010)
            GOTO 3000
3010      BACKSPACE(NCREPT)
          WRITE(NCREPT,9120)FNDMPO(IDFLAG+1)
          CLOSE(NCREPT)

        ENDIF

C       RESET THE DUMP FLAG
        IDFLAG = MOD(IDFLAG+1,2)

      ENDIF

C     =========================================================================

C     DUMP BC INFORMATION AS REQUIRED
C     ===============================
      IF(MOD(ITIME,NTDUMP).EQ.0)THEN

C       VM: BUG WHILE READING tixl FILES
        BCFLAG = (NSBCXL.EQ.NSBCI2).OR.(NSBCXL.EQ.NSBCI3)
        BCFLAG = BCFLAG.AND.(NXLPRM(1).EQ.3)

        IF(BCFLAG)THEN

C         DUMP THE INLET TURBULENT VELOCITY FIELD
          OPEN(UNIT=NCTIXL,FILE=FNTIXL,STATUS='OLD',
     +         FORM='UNFORMATTED')
          REWIND(NCTIXL)
          WRITE(NCTIXL)UFXL,VFXL,WFXL,SLOCXL,SVELXL,BVELXL
          CLOSE(NCTIXL)

        ENDIF

      ENDIF

C     =========================================================================

C     TIME STEP HISTORY
      IF(IPROC.EQ.0)THEN
        WRITE(*,'(I7,1PE12.4,I5)')ITIME,TSTEP,INDERR
      ENDIF

C     ========================================================================= 

C     STATISTICS ON THE FLY
C     =====================

C     STATISTICS MASTER SWITCH
C     ------------------------
      IF(NTSTAT.GE.0)THEN

C     ========================================================================= 

C       OUTPUT STATISTICS
C       =================

C       STATISTICS ON ONE PROCESSOR ONLY
C       ----------
        IF(IPROC.EQ.0)THEN

          IF(MOD(ITIME,NTSTAT).EQ.0)THEN

            OPEN(UNIT=NCSTAT,FILE=FNSTAT,STATUS='OLD',FORM='FORMATTED')

C           GO TO EOF
2000        CONTINUE
              READ(NCSTAT,9200,END=2010)
              GOTO 2000
2010        BACKSPACE(NCSTAT)

            WRITE(NCSTAT,9100)ITIME

            CLOSE(NCSTAT)

          ENDIF       

        ENDIF       

C       RESET STORAGE INDEX
        ITSTAT = 0
 
C     =========================================================================

C     STATISTICS MASTER SWITCH
      ENDIF

C     =========================================================================


      RETURN

9000  FORMAT(A)
9100  FORMAT('Time step number: ',I7)
9110  FORMAT('Elapsed time: ',1PE12.4,';',2X,'next time step:',1PE12.4)
9120  FORMAT('Dump completed: ',A)
9200  FORMAT(I5,/
     +       5X,6(1PE12.4)/
     +       5X,6(1PE12.4)/
     +       5X,6(1PE12.4)/
     +       5X,4(1PE12.4)/
     +       5X,3(1PE12.4)/
     +       5X,3(1PE12.4)/
     +       5X,3(1PE12.4)/
     +       5X,3(1PE12.4)/
     +       5X,3(1PE12.4)/
     +       5X,2(1PE12.4))
9300  FORMAT(2(1PE15.7))

      END
