      SUBROUTINE RADCIN 
 
C     *************************************************************************
C
C     RADCIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     14-JUL-2013:  CREATED
C      
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     RADIATION TREATMENT
C     USING OPTICALLY-THIN ASSUMPTION: Ju et al: JFM 342, 315-334, 1997.
C     INITIALISES DATA
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION PCOUNT
      INTEGER ISPEC,JSPEC,ICOEFF
      INTEGER NRDSPC,NCOUNT
      INTEGER IROOT
      CHARACTER*10 SPCRAD(NSPCMX)
 

C     BEGIN
C     =====

C     =========================================================================

C     CONTROL FLAGS
C     DEFAULTS: RADIATION TREATMENT NOT ACTIVATED
      FLRADN = .FALSE.

C     =========================================================================

C     CHECK SWITCH FOR RADIATION TREATMENT
      IF(NFRADN.EQ.1)THEN

C       =======================================================================

C       DATA IS READ BY THE LOWEST-RANKED PROCESSOR AND BROADCAST TO THE OTHERS
C       -----------------------------------------------------------------------

C       SET ID OF ROOT (BROADCASTING) PROCESSOR
        IROOT = 0

C       INITIALISE THE BROADCAST COUNTER
        PCOUNT = ZERO

C       INITIALISE THE BROADCAST ARRAY
        DO NCOUNT = 1, NPARAY
          PARRAY(NCOUNT) = ZERO
        ENDDO

C       =======================================================================

C       LOWEST-RANKED PROCESSOR
        IF(IPROC.EQ.0)THEN

C       =======================================================================

C         READ THE RADIATION DATA FILE
C         ----------------------------
C         ---------------------------------------------------------------------

          OPEN(UNIT=NCRADN,FILE=FNRADN,STATUS='OLD',FORM='FORMATTED')

C         ---------------------------------------------------------------------

C         READ AND IGNORE THE HEADER
          READ(NCRADN,*)
          READ(NCRADN,*)
          READ(NCRADN,*)
          READ(NCRADN,*)
          READ(NCRADN,*)

C         ---------------------------------------------------------------------

C         READ SPECIES LIST AND CHECK AGAINST CHEMICAL DATA
          READ(NCRADN,*)
          READ(NCRADN,'(I5)')NRDSPC
          IF(NRDSPC.NE.NSPEC)THEN
            WRITE(NCRADN,*)
            WRITE(NCRADN,*)'Warning: RADCIN: mismatch in no. of species'
            WRITE(NCRADN,*)'CHEMIN: ',NSPEC,'; RADCIN: ',NRDSPC
            WRITE(NCRADN,*)
          ENDIF

          DO JSPEC = 1, NSPEC
            READ(NCRADN,'(I5,3X,A)')ISPEC,SPCRAD(ISPEC)
            IF(SPCSYM(ISPEC).NE.SPCRAD(ISPEC))THEN
              WRITE(NCREPT,*)
            WRITE(NCREPT,*)'Warning: RADCIN: mismatch in species symbol'
              WRITE(NCREPT,*)'CHEMIN: ',SPCSYM(ISPEC),
     +                     '; RADCIN: ',SPCRAD(ISPEC)
              WRITE(NCREPT,*)
            ENDIF
          ENDDO

C         RADIATION DATA
          READ(NCRADN,*)
          READ(NCRADN,*)NSPRAD
          DO ISPEC = 1,NSPRAD
            READ(NCRADN,*)JSPEC,NKPRAD(ISPEC)
            READ(NCRADN,*)
     +          (AKPRAD(ICOEFF,ISPEC),ICOEFF=1,NKPRAD(ISPEC))
            NSPRID(ISPEC) = JSPEC
          ENDDO

C         ---------------------------------------------------------------------

C         READ END-OF-FILE LINE
          READ(NCRADN,*)
          CLOSE(NCRADN)

C         =====================================================================

C         COMPRESS THE DATA
C         -----------------

C         RADIATION DATA
          NCOUNT = 1
          DO ISPEC = 1,NSPEC
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NSPRAD)
            DO ICOEFF = 1, NSPRAD
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = AKPRAD(ICOEFF,ISPEC)
            ENDDO
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NSPRID(ISPEC))
          ENDDO

C         ---------------------------------------------------------------------

C         BROADCAST COUNTER
          PCOUNT = REAL(NCOUNT)

C         CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
          IF(NCOUNT.GT.NPARAY)THEN
          WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Warning: RADCIN: broadcast array size too small'
          WRITE(NCREPT,*)'Actual: ',NCOUNT,'; Available: ',NPARAY
          WRITE(NCREPT,*)
          ENDIF

C         =====================================================================

        ENDIF

C       =======================================================================

C       BROADCAST (OR RECEIVE) THE COUNTER
C       ----------------------------------
        CALL P_BCST(PCOUNT,1,1,IROOT)
        NCOUNT = NINT(PCOUNT)

C       BROADCAST (OR RECEIVE) THE DATA
C       -------------------------------
        CALL P_BCST(PARRAY,NPARAY,NCOUNT,IROOT)

C       =======================================================================

C       NOT THE LOWEST-RANKED PROCESSOR
        IF(IPROC.NE.0)THEN

C         =====================================================================

C         UNCOMPRESS THE DATA
C         -------------------

C         RADIATION DATA
          NCOUNT = 1
          DO ISPEC = 1,NSPEC
            NCOUNT = NCOUNT + 1
            NSPRAD = NINT(PARRAY(NCOUNT))
            DO ICOEFF = 1, NSPRAD
              NCOUNT = NCOUNT + 1
              AKPRAD(ICOEFF,ISPEC) = PARRAY(NCOUNT)
            ENDDO
            NCOUNT = NCOUNT + 1
            NSPRID(ISPEC) = NINT(PARRAY(NCOUNT))
          ENDDO

C         =====================================================================

        ENDIF

C       =======================================================================

C       EVALUATE DERIVED QUANTITIES
C       ===========================

C       COEFFICIENT COUNTERS
        DO ISPEC = 1,NSPEC
          NKPRM1(ISPEC) = NKPRAD(ISPEC) - 1
        ENDDO

C       =======================================================================

C       CONSTANTS
        FOURSB = FOUR*STEFBO
        TRFRTH = TREFRN*TREFRN*TREFRN*TREFRN

C       =======================================================================

C       CONTROL FLAGS
        FLRADN = .TRUE.

C       =======================================================================

      ENDIF

C     =========================================================================


      RETURN
      END
