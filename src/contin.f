      SUBROUTINE CONTIN 
 
C     *************************************************************************
C
C     CONTIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-DEC-2003:  CREATED
C     11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH
C     01-MAY-2013:  RSC MIXTURE AVERAGED TRANSPORT
C     21-JUL-2013:  RSC RADIATION TREATMENT
C      
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES RUN CONTROL DATA
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
      INCLUDE 'com_espect.h'
C     -------------------------------------------------------------------------


C     PARAMETER
C     =========
C     RSC 01-MAY-2013
      CHARACTER*32 STRMOL
      PARAMETER(STRMOL = 'Molecular transport control data')


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION PCOUNT
      INTEGER NCOUNT
      INTEGER ISPEC,JSPEC
      INTEGER IC
      INTEGER IROOT
C     RSC 01-MAY-2013
      CHARACTER*50 EOFLIN


C     BEGIN
C     =====

C     =========================================================================

C     DATA IS READ BY THE LOWEST-RANKED PROCESSOR AND BROADCAST TO THE OTHERS
C     -----------------------------------------------------------------------

C     SET ID OF ROOT (BROADCASTING) PROCESSOR
      IROOT = 0

C     INITIALISE THE BROADCAST COUNTER
      PCOUNT = ZERO

C     INITIALISE THE BROADCAST ARRAY
      DO NCOUNT = 1,NPARAY
        PARRAY(NCOUNT) = ZERO
      ENDDO

C     =========================================================================

C     LOWEST-RANKED PROCESSOR
      IF(IPROC.EQ.0)THEN

C       =======================================================================

C       READ THE RUN CONTROL DATA FILE
C       ------------------------------
C       -----------------------------------------------------------------------

        OPEN(UNIT=NCCONT,FILE=FNCONT,STATUS='OLD',FORM='FORMATTED')

C       -----------------------------------------------------------------------

C       READ AND IGNORE THE HEADER
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)
C       AND THE SUBSEQUENT BLANK LINE
        READ(NCCONT,*)

C       -----------------------------------------------------------------------

C       GRID SIZE DATA
C       --------------
C       REQUIRED GLOBAL GRID SIZE X,Y,Z
        READ(NCCONT,*)
        READ(NCCONT,*)XGDLEN,YGDLEN,ZGDLEN
        READ(NCCONT,*)

C       REQUIRED GLOBAL GRID SIZE NX,NY,NZ
        READ(NCCONT,*)
        READ(NCCONT,*)NXGREQ,NYGREQ,NZGREQ
        READ(NCCONT,*)

C       REQUIRED NUMBER OF PROCESSORS X,Y,Z
        READ(NCCONT,*)
        READ(NCCONT,*)NXPREQ,NYPREQ,NZPREQ
        READ(NCCONT,*)

C       TIME STEPPING DATA
C       ------------------
C       TIME STEP, START STEP, NO OF STEPS, STEP SWITCH
        READ(NCCONT,*)
        READ(NCCONT,*)TSTEP,NTIME1,NTIME,NSTPSW
        READ(NCCONT,*)

C       INTERVALS BETWEEN DUMPS,REPORTS,STATISTICS
C       ------------------------------------------
C       DUMP OUTPUT FORMAT
C       RSC 11-JUL-2009
        READ(NCCONT,*)
        READ(NCCONT,*)NTDUMP,NTREPT,NTSTAT,NDOFMT
        READ(NCCONT,*)

C       COLD START SWITCH (0=COLD START, 1=RESTART)
C       -----------------
C       DUMP INPUT FORMAT
C       RSC 11-JUL-2009
        READ(NCCONT,*)
        READ(NCCONT,*)NCDMPI,NDIFMT
        READ(NCCONT,*)

C       INITIAL TURBULENCE (0=INITIAL TURBULENCE OFF, 1=ON)
C       ------------------
C       RANDOM SEED, SPECTRUM PARAMETERS (FOUR)
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)INTURB,INSEED,(SPARAM(IC),IC=1,NSPARM)
        READ(NCCONT,*)

C       FLAME START OPTION (0=NO FLAME, 1=FLAME)
C       ------------------
        READ(NCCONT,*)
        READ(NCCONT,*)INFLAM
        READ(NCCONT,*)

C       DEFAULT INITIAL CONDITIONS
C       --------------------------
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)PRIN,TRIN,URIN,VRIN,WRIN
        READ(NCCONT,*)
        READ(NCCONT,*)NSPREQ
        DO ISPEC = 1, NSPREQ
          READ(NCCONT,*)JSPEC,YRIN(ISPEC)
        ENDDO
        READ(NCCONT,*)

C       GLOBAL BOUNDARY CONDITION TYPES
C       -------------------------------
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)
        READ(NCCONT,*)NGBCXL,(NXLPRM(IC),IC=1,NBCPRI),
     +                       (RXLPRM(IC),IC=1,NBCPRR)
C       FY - ADDED FOR READING OF EXTRA VARIABLES FOR NON-REFLECTING INFLOW
        IF(NGBCXL.EQ.NSBCI4)THEN
          READ(NCCONT,*)(RXLPRC(IC),IC=1,NBCPRC)
        ENDIF
        READ(NCCONT,*)NGBCXR,(NXRPRM(IC),IC=1,NBCPRI),
     +                       (RXRPRM(IC),IC=1,NBCPRR)
        READ(NCCONT,*)NGBCYL,(NYLPRM(IC),IC=1,NBCPRI),
     +                       (RYLPRM(IC),IC=1,NBCPRR)
        READ(NCCONT,*)NGBCYR,(NYRPRM(IC),IC=1,NBCPRI),
     +                       (RYRPRM(IC),IC=1,NBCPRR)
        READ(NCCONT,*)NGBCZL,(NZLPRM(IC),IC=1,NBCPRI),
     +                       (RZLPRM(IC),IC=1,NBCPRR)
        READ(NCCONT,*)NGBCZR,(NZRPRM(IC),IC=1,NBCPRI),
     +                       (RZRPRM(IC),IC=1,NBCPRR)
        READ(NCCONT,*)

C       -----------------------------------------------------------------------

C       RSC 01-MAY-2013
C       READ MOLECULAR TRANSPORT CONTROL DATA
C       RSC 21-JUL-2013
C       READ RADIATION CONTROL DATA
        NFMAVT = 0
        NFMIXW = 0
        NFMIXP = 0
        NFMSOR = 0
        NFMDUF = 0
        WMLTDR = 0.0D0
        NFRADN = 0
        TREFRN = 3.0D0
        READ(NCCONT,'(A)')EOFLIN
        IF(EOFLIN(1:32).EQ.STRMOL)THEN
          READ(NCCONT,*)
          READ(NCCONT,*)NFMAVT,NFMIXW,NFMIXP,NFMSOR,NFMDUF,WMLTDR
          READ(NCCONT,*)
          READ(NCCONT,*)
          READ(NCCONT,*)
          READ(NCCONT,*)NFRADN,TREFRN
          READ(NCCONT,*)
        ENDIF

C       READ END-OF-FILE LINE
        READ(NCCONT,'(A)')EOFLIN

        CLOSE(NCCONT)

C       =======================================================================

C       COMPRESS THE DATA
C       -----------------

C       REQUIRED GLOBAL DOMAIN SIZE (x,y,z)
        NCOUNT = 1
        PARRAY(NCOUNT) = XGDLEN
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = YGDLEN
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = ZGDLEN

C       TIME STEP, START STEP, NO OF STEPS, STEP SWITCH
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = TSTEP
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NTIME1)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NTIME)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NSTPSW)

C       INTERVALS BETWEEN DUMPS,REPORTS,STATISTICS; DUMP OUTPUT SWITCH
C       RSC 11-JUL-2009
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NTDUMP)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NTREPT)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NTSTAT)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NDOFMT)

C       COLD START SWITCH; DUMP INPUT SWITCH
C       RSC 11-JUL-2009
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NCDMPI)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NDIFMT)

C       INITIAL TURBULENCE
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(INTURB)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(INSEED)
        DO ISPEC = 1, NSPARM
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = SPARAM(ISPEC)
        ENDDO

C       FLAME START OPTION
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(INFLAM)

C       DEFAULT INITIAL CONDITIONS
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = PRIN
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = TRIN
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = URIN
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = VRIN
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = WRIN
        DO ISPEC = 1, NSPCMX
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = YRIN(ISPEC)
        ENDDO

C       GLOBAL BOUNDARY CONDITION TYPES
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NGBCXL)
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NXLPRM(IC))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = RXLPRM(IC)
        ENDDO
        IF(NGBCXL.EQ.NSBCI4)THEN ! FY - ADDED FOR NSBCI4
          DO IC = 1, NBCPRC
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = RXLPRC(IC)
          ENDDO
        ENDIF
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NGBCXR)
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NXRPRM(IC))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = RXRPRM(IC)
        ENDDO
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NGBCYL)
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NYLPRM(IC))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = RYLPRM(IC)
        ENDDO
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NGBCYR)
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NYRPRM(IC))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = RYRPRM(IC)
        ENDDO
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NGBCZL)
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NZLPRM(IC))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = RZLPRM(IC)
        ENDDO
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NGBCZR)
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NZRPRM(IC))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = RZRPRM(IC)
        ENDDO

C       MIXTURE AVERAGED TRANSPORT
C       RSC 05-MAY-2013
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NFMAVT)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NFMIXW)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NFMIXP)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NFMSOR)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NFMDUF)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = WMLTDR

C       RADIATION TREATMENT 
C       RSC 21-JUL-2013
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NFRADN)
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = TREFRN

C       BROADCAST COUNTER
        PCOUNT = REAL(NCOUNT)

C       CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
        IF(NCOUNT.GT.NPARAY)THEN
        WRITE(NCREPT,*)'Warning: INCONT: broadcast array size too small'
        ENDIF

C       =======================================================================

      ENDIF

C     =========================================================================

C     BROADCAST (OR RECEIVE) THE COUNTER
C     ----------------------------------
      CALL P_BCST(PCOUNT,1,1,IROOT)
      NCOUNT = NINT(PCOUNT)

C     BROADCAST (OR RECEIVE) THE DATA
C     -------------------------------
      CALL P_BCST(PARRAY,NPARAY,NCOUNT,IROOT)

C     =========================================================================

C     NOT THE LOWEST-RANKED PROCESSOR
      IF(IPROC.NE.0)THEN

C       =======================================================================

C       UNCOMPRESS THE DATA
C       -------------------

C       REQUIRED GLOBAL GRID SIZE X,Y,Z
        NCOUNT = 1
        XGDLEN = PARRAY(NCOUNT)
        NCOUNT = NCOUNT + 1
        YGDLEN = PARRAY(NCOUNT)
        NCOUNT = NCOUNT + 1
        ZGDLEN = PARRAY(NCOUNT)

C       TIME STEP, START STEP, NO OF STEPS, STEP SWITCH
        NCOUNT = NCOUNT + 1
        TSTEP = PARRAY(NCOUNT)
        NCOUNT = NCOUNT + 1
        NTIME1 = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NTIME = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NSTPSW = NINT(PARRAY(NCOUNT))

C       INTERVALS BETWEEN DUMPS,REPORTS,STATISTICS; DUMP OUTPUT SWITCH
C       RSC 11-JUL-2009
        NCOUNT = NCOUNT + 1
        NTDUMP = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NTREPT = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NTSTAT = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NDOFMT = NINT(PARRAY(NCOUNT))

C       COLD START SWITCH; DUMP INPUT SWITCH
C       RSC 11-JUL-2009
        NCOUNT = NCOUNT + 1
        NCDMPI = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NDIFMT = NINT(PARRAY(NCOUNT))

C       INITIAL TURBULENCE
        NCOUNT = NCOUNT + 1
        INTURB = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        INSEED = NINT(PARRAY(NCOUNT))
        DO ISPEC = 1, NSPARM
          NCOUNT = NCOUNT + 1
          SPARAM(ISPEC) = PARRAY(NCOUNT)
        ENDDO

C       FLAME START OPTION
        NCOUNT = NCOUNT + 1
        INFLAM = NINT(PARRAY(NCOUNT))

C       DEFAULT INITIAL CONDITIONS
        NCOUNT = NCOUNT + 1
        PRIN = PARRAY(NCOUNT)
        NCOUNT = NCOUNT + 1
        TRIN = PARRAY(NCOUNT)
        NCOUNT = NCOUNT + 1
        URIN = PARRAY(NCOUNT)
        NCOUNT = NCOUNT + 1
        VRIN = PARRAY(NCOUNT)
        NCOUNT = NCOUNT + 1
        WRIN = PARRAY(NCOUNT)
        DO ISPEC = 1, NSPCMX
          NCOUNT = NCOUNT + 1
          YRIN(ISPEC) = PARRAY(NCOUNT)
        ENDDO

C       GLOBAL BOUNDARY CONDITION TYPES
        NCOUNT = NCOUNT + 1
        NGBCXL = NINT(PARRAY(NCOUNT))
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          NXLPRM(IC) = NINT(PARRAY(NCOUNT))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          RXLPRM(IC) = PARRAY(NCOUNT)
        ENDDO
        IF(NGBCXL.EQ.NSBCI4)THEN ! FY - ADDED FOR NSBCI4
          DO IC = 1, NBCPRC
            NCOUNT = NCOUNT + 1
            RXLPRC(IC) = PARRAY(NCOUNT)
          ENDDO
        ENDIF
        NCOUNT = NCOUNT + 1
        NGBCXR = NINT(PARRAY(NCOUNT))
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          NXRPRM(IC) = NINT(PARRAY(NCOUNT))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          RXRPRM(IC) = PARRAY(NCOUNT)
        ENDDO
        NCOUNT = NCOUNT + 1
        NGBCYL = NINT(PARRAY(NCOUNT))
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          NYLPRM(IC) = NINT(PARRAY(NCOUNT))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          RYLPRM(IC) = PARRAY(NCOUNT)
        ENDDO
        NCOUNT = NCOUNT + 1
        NGBCYR = NINT(PARRAY(NCOUNT))
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          NYRPRM(IC) = NINT(PARRAY(NCOUNT))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          RYRPRM(IC) = PARRAY(NCOUNT)
        ENDDO
        NCOUNT = NCOUNT + 1
        NGBCZL = NINT(PARRAY(NCOUNT))
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          NZLPRM(IC) = NINT(PARRAY(NCOUNT))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          RZLPRM(IC) = PARRAY(NCOUNT)
        ENDDO
        NCOUNT = NCOUNT + 1
        NGBCZR = NINT(PARRAY(NCOUNT))
        DO IC = 1, NBCPRI
          NCOUNT = NCOUNT + 1
          NZRPRM(IC) = NINT(PARRAY(NCOUNT))
        ENDDO
        DO IC = 1, NBCPRR
          NCOUNT = NCOUNT + 1
          RZRPRM(IC) = PARRAY(NCOUNT)
        ENDDO

C       MIXTURE AVERAGED TRANSPORT
C       RSC 05-MAY-2013
        NCOUNT = NCOUNT + 1
        NFMAVT = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NFMIXW = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NFMIXP = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NFMSOR = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        NFMDUF = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        WMLTDR = PARRAY(NCOUNT)

C       RADIATION TREATMENT 
C       RSC 21-JUL-2013
        NCOUNT = NCOUNT + 1
        NFRADN = NINT(PARRAY(NCOUNT))
        NCOUNT = NCOUNT + 1
        TREFRN = PARRAY(NCOUNT)

C       =======================================================================

      ENDIF

C     =========================================================================


      RETURN
      END
