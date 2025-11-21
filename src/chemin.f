      SUBROUTINE CHEMIN 
 
C     *************************************************************************
C
C     CHEMIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     08-JAN-2003:  CREATED
C     29-DEC-2006:  RSC BUG FIX INDEXING (THANKS TO MDLLP)
C      
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES DATA FOR MULTI-SPECIES MULTI-STEP CHEMISTRY
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
      DOUBLE PRECISION GIBLET
      INTEGER ISPEC,ISTEP,IBODY,ILIND,ITROE,ISRIF
      INTEGER JSPEC,JSTEP,JBODY,JLIND,JTROE,JSRIF
      INTEGER ITINT,ICP
      INTEGER NCOUNT
      INTEGER IROOT
      

C     BEGIN
C     =====

C     =========================================================================

C     GLOBAL DATA
C     -----------
C     NSPEC total no of species
C     NSTEP total no of steps
C     NBODY total no of third bodies
C     NGIBB total no of Gibbs steps
C     NLIND total no of Lindemann steps
C     NTROE total no of Troe steps
C     NSRIF total no of SRI steps
C     NSPCMX max no of species
C     NSTPMX max no of steps
C     NSSMAX max length of step species-list
C     NRSMAX max length of step reactant- and product-lists
C     NBDYMX max no of third bodies
C     NLLMAX max no of Lindemann steps
C     RPARAM(1,NSTPMX) pre-exponential factor ln(A)
C     RPARAM(2,NSTPMX) temperature exponent n
C     RPARAM(3,NSTPMX) activation energy E/R0
C     NSSPEC(NSSMAX,NSTPMX) is the step species-list
C     NRSPEC(NRSMAX,NSTPMX) is the step reactant-list
C     NPSPEC(NRSMAX,NSTPMX) is the step product-list
C     NRCPEC(NRSMAX,NSTPMX) is the step coefficient reactant-list
C     NPCPEC(NRSMAX,NSTPMX) is the step coefficient product-list
C     NSSLEN(NSTPMX) contains the length of the species-list for each step
C     NRSLEN(NSTPMX) contains the length of the reactant-list for each step
C     NPSLEN(NSTPMX) contains the length of the product-list for each step
C     NRCLEN(NSTPMX) contains the length of the reactant coefficient-list for
C                    each step
C     NPCLEN(NSTPMX) contains the length of the reactant coefficient-list for
C                    each step
C     WMOLAR(NSPCMX) molar mass of species
C     OVWMOL(NSPCMX) reciprocal of molar mass of species
C     DIFFMU(NSSMAX,NSTPMX) is the species delta-mu for each step
C     DIFFMW(NSSMAX,NSTPMX) is the species delta-mu (times molar mass)
C                           for each step
C     CRSPEC(NSSMAX,NSTPMX) is the reactant stoichiometric coefficient-list
C                           for each step
C     CPSPEC(NSSMAX,NSTPMX) is the product stoichiometric coefficient-list
C                           for each step
C     MBLIST(NSTPMX) contains 0 if no third body in a step
C                    else contains the index number of the third body
C     MGLIST(NSTPMX) contains 0 if no Gibbs function evaluation in a step
C                    else contains the index number of the Gibbs step
C     MLLIST(NSTPMX) contains 0 if no Lindemann rate in a step
C                    else contains the index number of the Lindemann step
C     MTLIST(NSTPMX) contains 0 if no Troe form rate in a step
C                    else contains the index number of the Troe form step
C     MSLIST(NSTPMX) contains 0 if no SRI form rate in a step
C                    else contains the index number of the SRI form step
C     EFFY3B(NSPCMX,NBDYMX) contains the third-body efficiencies
C     RCLIND(1,MLLIST(NSTPMX)) Lindemann rate pre-exponential factor ln(A)
C     RCLIND(2,MLLIST(NSTPMX)) Lindemann rate temperature exponent n
C     RCLIND(3,MLLIST(NSTPMX)) Lindemann rate activation energy E/R0
C     RCLIND(4,MLLIST(NSTPMX)) Lindemann rate broadening factor
C     RCTROE(1,MTLIST(NSTPMX)) Troe form rate pre-exponential factor ln(A)
C     RCTROE(2,MTLIST(NSTPMX)) Troe form rate temperature exponent n
C     RCTROE(3,MTLIST(NSTPMX)) Troe form rate activation energy E/R0
C     RCTROE(4,MTLIST(NSTPMX)) Troe form rate temperature parameter alpha
C     RCTROE(5,MTLIST(NSTPMX)) Troe form rate temperature parameter no 1*
C     RCTROE(6,MTLIST(NSTPMX)) Troe form rate temperature parameter no 2*
C     RCTROE(7,MTLIST(NSTPMX)) Troe form rate temperature parameter no 3*
C     RCTROE(8,MTLIST(NSTPMX)) Troe form rate c-const no 1
C     RCTROE(9,MTLIST(NSTPMX)) Troe form rate c-const no 2
C     RCTROE(10,MTLIST(NSTPMX)) Troe form rate n-const no 1
C     RCTROE(11,MTLIST(NSTPMX)) Troe form rate n-const no 2
C     RCTROE(12,MTLIST(NSTPMX)) Troe form rate d-const
C     RCSRIF(1,MSLIST(NSTPMX)) SRI form rate pre-exponential factor ln(A)
C     RCSRIF(2,MSLIST(NSTPMX)) SRI form rate temperature exponent n
C     RCSRIF(3,MSLIST(NSTPMX)) SRI form rate activation energy E/R0
C     RCSRIF(4,MSLIST(NSTPMX)) SRI form rate parameter a
C     RCSRIF(5,MSLIST(NSTPMX)) SRI form rate parameter b
C     RCSRIF(6,MSLIST(NSTPMX)) SRI form rate parameter c
C     RCSRIF(7,MSLIST(NSTPMX)) SRI form rate parameter d
C     RCSRIF(8,MSLIST(NSTPMX)) SRI form rate parameter e

C     =========================================================================

C     DATA IS READ BY THE LOWEST-RANKED PROCESSOR AND BROADCAST TO THE OTHERS
C     -----------------------------------------------------------------------

C     SET ID OF ROOT (BROADCASTING) PROCESSOR
      IROOT = 0

C     INITIALISE THE BROADCAST COUNTER
      PCOUNT = ZERO

C     INITIALISE THE BROADCAST ARRAY
      DO NCOUNT = 1, NPARAY
        PARRAY(NCOUNT) = ZERO
      ENDDO

C     =========================================================================

C     LOWEST-RANKED PROCESSOR
      IF(IPROC.EQ.0)THEN

C       =======================================================================

C       READ THE CHEMICAL DATA FILE
C       ---------------------------
C       -----------------------------------------------------------------------

        OPEN(UNIT=NCCHEM,FILE=FNCHEM,STATUS='OLD',FORM='FORMATTED')

C       -----------------------------------------------------------------------

C       READ AND IGNORE THE HEADER
        READ(NCCHEM,*)
        READ(NCCHEM,*)
        READ(NCCHEM,*)
        READ(NCCHEM,*)
        READ(NCCHEM,*)

C       -----------------------------------------------------------------------

C       SPECIES LIST
        READ(NCCHEM,*)
        READ(NCCHEM,'(I5)')NSPEC
        DO JSPEC = 1, NSPEC
          READ(NCCHEM,'(I5,3X,A)')ISPEC,SPCSYM(ISPEC)
        ENDDO

C       SPECIES DATA
        READ(NCCHEM,*)
        READ(NCCHEM,*)PREFGB
        DO JSPEC = 1,NSPEC
          READ(NCCHEM,'(I5,1PE12.4)')ISPEC,WMOLAR(ISPEC)
          READ(NCCHEM,'(1PE12.4)')CLEWIS(ISPEC)
          READ(NCCHEM,'(I5)')NTINT(ISPEC)
          DO ITINT = 1,NTINT(ISPEC)
            READ(NCCHEM,'(2(1PE12.4),I5)')TINTLO(ITINT,ISPEC),
     +                    TINTHI(ITINT,ISPEC),NCOFCP(ITINT,ISPEC)
            DO ICP = 1,NCOFCP(ITINT,ISPEC)
              READ(NCCHEM,'(1PE15.7)')AMOLCP(ICP,ITINT,ISPEC)
            ENDDO
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       STEP RATE DATA
        READ(NCCHEM,*)
        READ(NCCHEM,'(I5)')NSTEP
        DO JSTEP = 1, NSTEP
        READ(NCCHEM,'(I5,3(1PE12.4))')ISTEP,(RPARAM(ICP,ISTEP),ICP=1,3)
        ENDDO

C       -----------------------------------------------------------------------

C       STEP SPECIES-LIST
        READ(NCCHEM,*)
        DO JSTEP = 1,NSTEP
          READ(NCCHEM,'(2I5)')ISTEP,NSSLEN(ISTEP)
          DO JSPEC = 1, NSSLEN(ISTEP)
            READ(NCCHEM,'(2I5,I8)')ISTEP,ISPEC,NSSPEC(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       STEP REACTANT-LIST
        READ(NCCHEM,*)
        DO JSTEP = 1,NSTEP
          READ(NCCHEM,'(2I5)')ISTEP,NRSLEN(ISTEP)
          DO JSPEC = 1, NRSLEN(ISTEP)
            READ(NCCHEM,'(2I5,I8)')ISTEP,ISPEC,NRSPEC(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       STEP PRODUCT-LIST
        READ(NCCHEM,*)
        DO JSTEP = 1,NSTEP
          READ(NCCHEM,'(2I5)')ISTEP,NPSLEN(ISTEP)
          DO JSPEC = 1, NPSLEN(ISTEP)
            READ(NCCHEM,'(2I5,I8)')ISTEP,ISPEC,NPSPEC(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       STEP REACTANT COEFFICIENT-LIST
        READ(NCCHEM,*)
        DO JSTEP = 1,NSTEP
          READ(NCCHEM,'(2I5)')ISTEP,NRCLEN(ISTEP)
          DO JSPEC = 1, NRCLEN(ISTEP)
            READ(NCCHEM,'(2I5,I8,1PE12.4)')ISTEP,ISPEC,
     +                      NRCPEC(ISPEC,ISTEP),CRSPEC(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       STEP PRODUCT COEFFICIENT-LIST
        READ(NCCHEM,*)
        DO JSTEP = 1,NSTEP
          READ(NCCHEM,'(2I5)')ISTEP,NPCLEN(ISTEP)
          DO JSPEC = 1, NPCLEN(ISTEP)
            READ(NCCHEM,'(2I5,I8,1PE12.4)')ISTEP,ISPEC,
     +                      NPCPEC(ISPEC,ISTEP),CPSPEC(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       SPECIES DELTA-LIST
        READ(NCCHEM,*)
        DO JSTEP = 1,NSTEP
          DO JSPEC = 1,NSSLEN(JSTEP)
            READ(NCCHEM,'(2I5,F5.1)')ISTEP,ISPEC,DIFFMU(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       THIRD-BODY LIST
        READ(NCCHEM,*)
        READ(NCCHEM,'(I5)')NBODY
        DO JBODY = 1, NBODY
          READ(NCCHEM,'(I5,3X,A)')IBODY,BDYSYM(IBODY)
        ENDDO

C       THIRD-BODY STEP-LIST
        READ(NCCHEM,*)
        IF(NBODY.GT.0)THEN
          DO JSTEP = 1, NSTEP
            READ(NCCHEM,'(I5,I8)')ISTEP,MBLIST(ISTEP)
          ENDDO
        ENDIF

C       THIRD-BODY EFFICIENCIES
        READ(NCCHEM,*)
        DO JBODY = 1,NBODY
          DO JSPEC = 1,NSPEC
            READ(NCCHEM,'(2I5,1PE12.4)')IBODY,ISPEC,EFFY3B(ISPEC,IBODY)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       GIBBS STEP-LIST
        READ(NCCHEM,*)
        READ(NCCHEM,'(I5)')NGIBB
        IF(NGIBB.GT.0)THEN
          DO JSTEP = 1, NSTEP
            READ(NCCHEM,'(I5,I8)')ISTEP,MGLIST(ISTEP)
          ENDDO
        ENDIF

C       -----------------------------------------------------------------------

C       LINDEMANN STEP-LIST
        READ(NCCHEM,*)
        READ(NCCHEM,'(I5)')NLIND
        IF(NLIND.GT.0)THEN
          DO JSTEP = 1, NSTEP
            READ(NCCHEM,'(I5,I8)')ISTEP,MLLIST(ISTEP)
          ENDDO
        ENDIF

C       LINDEMANN STEP RATE DATA
        READ(NCCHEM,*)
        DO JLIND = 1, NLIND
        READ(NCCHEM,'(I5,4(1PE12.4))')ILIND,(RCLIND(ICP,ILIND),ICP=1,4)
        ENDDO

C       -----------------------------------------------------------------------

C       TROE FORM STEP-LIST
        READ(NCCHEM,*)
        READ(NCCHEM,'(I5)')NTROE
        IF(NTROE.GT.0)THEN
          DO JSTEP = 1, NSTEP
            READ(NCCHEM,'(I5,I8)')ISTEP,MTLIST(ISTEP)
          ENDDO
        ENDIF

C       TROE FORM STEP RATE DATA
        READ(NCCHEM,*)
        DO JTROE = 1, NTROE
        READ(NCCHEM,'(I5,6(1PE12.4))')ITROE,(RCTROE(ICP,ITROE),ICP=1,6)
        READ(NCCHEM,'(I5,6(1PE12.4))')ITROE,(RCTROE(ICP,ITROE),ICP=7,12)
        ENDDO

C       -----------------------------------------------------------------------

C       SRI FORM STEP-LIST
        READ(NCCHEM,*)
        READ(NCCHEM,'(I5)')NSRIF
        IF(NSRIF.GT.0)THEN
          DO JSTEP = 1, NSTEP
            READ(NCCHEM,'(I5,I8)')ISTEP,MSLIST(ISTEP)
          ENDDO
        ENDIF

C       SRI FORM STEP RATE DATA
        READ(NCCHEM,*)
        DO JSRIF = 1, NSRIF
        READ(NCCHEM,'(I5,3(1PE12.4))')ISRIF,(RCSRIF(ICP,ISRIF),ICP=1,3)
        READ(NCCHEM,'(I5,5(1PE12.4))')ISRIF,(RCSRIF(ICP,ISRIF),ICP=4,8)
        ENDDO

C       -----------------------------------------------------------------------

C       READ END-OF-FILE LINE
        READ(NCCHEM,*)
        CLOSE(NCCHEM)

C       =======================================================================

C       COMPRESS THE DATA
C       -----------------

C       NO OF SPECIES
        NCOUNT = 1
        PARRAY(NCOUNT) = REAL(NSPEC)

C       SPECIES STRINGS
C       RSC 29-DEC-2006 BUG FIX INDEXING
        DO ISPEC = 1,NSPEC
          DO ICP = 1, NSPSTR
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(ICHAR(SPCSYM(ISPEC)(ICP:ICP)))
          ENDDO
        ENDDO

C       SPECIES DATA
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = PREFGB
        DO ISPEC = 1,NSPEC
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = WMOLAR(ISPEC)
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = CLEWIS(ISPEC)
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NTINT(ISPEC))
          DO ITINT = 1,NTINT(ISPEC)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = TINTLO(ITINT,ISPEC)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = TINTHI(ITINT,ISPEC)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NCOFCP(ITINT,ISPEC))
            DO ICP = 1,NCOFCP(ITINT,ISPEC)
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = AMOLCP(ICP,ITINT,ISPEC)
            ENDDO
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       NO OF STEPS
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NSTEP)

C       STEP RATE DATA
        DO ISTEP = 1, NSTEP
          DO ICP = 1,3
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = RPARAM(ICP,ISTEP)
          ENDDO
        ENDDO

C       STEP SPECIES-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NSSLEN(ISTEP))
          DO ISPEC = 1, NSSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NSSPEC(ISPEC,ISTEP))
          ENDDO
        ENDDO

C       STEP REACTANT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NRSLEN(ISTEP))
          DO ISPEC = 1, NRSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NRSPEC(ISPEC,ISTEP))
          ENDDO
        ENDDO

C       STEP PRODUCT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NPSLEN(ISTEP))
          DO ISPEC = 1, NPSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NPSPEC(ISPEC,ISTEP))
          ENDDO
        ENDDO

C       STEP REACTANT COEFFICIENT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NRCLEN(ISTEP))
          DO ISPEC = 1, NRCLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NRCPEC(ISPEC,ISTEP))
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = CRSPEC(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       STEP PRODUCT COEFFICIENT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = REAL(NPCLEN(ISTEP))
          DO ISPEC = 1, NPCLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NPCPEC(ISPEC,ISTEP))
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = CPSPEC(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       SPECIES DELTA-LIST
        DO ISTEP = 1,NSTEP
          DO ISPEC = 1, NSSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = DIFFMU(ISPEC,ISTEP)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       THIRD-BODY LIST
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NBODY)
        DO IBODY = 1, NBODY
          DO ICP = 1, NSPSTR
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(ICHAR(BDYSYM(IBODY)(ICP:ICP)))
          ENDDO
        ENDDO

C       THIRD-BODY STEP-LIST
        IF(NBODY.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(MBLIST(ISTEP))
          ENDDO
        ENDIF

C       THIRD-BODY EFFICIENCIES
        DO IBODY = 1, NBODY
          DO ISPEC = 1,NSPEC
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = EFFY3B(ISPEC,IBODY)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       GIBBS STEP-LIST
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NGIBB)
        IF(NGIBB.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(MGLIST(ISTEP))
          ENDDO
        ENDIF

C       -----------------------------------------------------------------------

C       LINDEMANN STEP-LIST
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NLIND)
        IF(NLIND.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(MLLIST(ISTEP))
          ENDDO
        ENDIF

C       LINDEMANN STEP RATE DATA
        DO ILIND = 1, NLIND
          DO ICP = 1,4
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = RCLIND(ICP,ILIND)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       TROE FORM STEP-LIST
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NTROE)
        IF(NTROE.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(MTLIST(ISTEP))
          ENDDO
        ENDIF

C       TROE FORM STEP RATE DATA
        DO ITROE = 1, NTROE
          DO ICP = 1,12
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = RCTROE(ICP,ITROE)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       SRI FORM STEP-LIST
        NCOUNT = NCOUNT + 1
        PARRAY(NCOUNT) = REAL(NSRIF)
        IF(NSRIF.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(MSLIST(ISTEP))
          ENDDO
        ENDIF

C       SRI FORM STEP RATE DATA
        DO ISRIF = 1, NSRIF
          DO ICP = 1,8
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = RCSRIF(ICP,ISRIF)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       BROADCAST COUNTER
        PCOUNT = REAL(NCOUNT)

C       CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
        IF(NCOUNT.GT.NPARAY)THEN
        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Warning: CHEMIN: broadcast array size too small'
        WRITE(NCREPT,*)'Actual: ',NCOUNT,'; Available: ',NPARAY
        WRITE(NCREPT,*)
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

C       NO OF SPECIES
        NCOUNT = 1
        NSPEC = NINT(PARRAY(NCOUNT))

C       SPECIES STRINGS
C       RSC 29-DEC-2006 BUG FIX INDEXING
        DO ISPEC = 1,NSPEC
          DO ICP = 1, NSPSTR
            NCOUNT = NCOUNT + 1
            SPCSYM(ISPEC)(ICP:ICP) = CHAR(NINT(PARRAY(NCOUNT)))
          ENDDO
        ENDDO

C       SPECIES DATA
        NCOUNT = NCOUNT + 1
        PREFGB = PARRAY(NCOUNT)
        DO ISPEC = 1,NSPEC
          NCOUNT = NCOUNT + 1
          WMOLAR(ISPEC) = PARRAY(NCOUNT)
          NCOUNT = NCOUNT + 1
          CLEWIS(ISPEC) = PARRAY(NCOUNT)
          NCOUNT = NCOUNT + 1
          NTINT(ISPEC) = NINT(PARRAY(NCOUNT))
          DO ITINT = 1,NTINT(ISPEC)
            NCOUNT = NCOUNT + 1
            TINTLO(ITINT,ISPEC) = PARRAY(NCOUNT)
            NCOUNT = NCOUNT + 1
            TINTHI(ITINT,ISPEC) = PARRAY(NCOUNT)
            NCOUNT = NCOUNT + 1
            NCOFCP(ITINT,ISPEC) = NINT(PARRAY(NCOUNT))
            DO ICP = 1,NCOFCP(ITINT,ISPEC)
              NCOUNT = NCOUNT + 1
              AMOLCP(ICP,ITINT,ISPEC) = PARRAY(NCOUNT)
            ENDDO
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       NO OF STEPS
        NCOUNT = NCOUNT + 1
        NSTEP = NINT(PARRAY(NCOUNT))

C       STEP RATE DATA
        DO ISTEP = 1, NSTEP
          DO ICP = 1,3
            NCOUNT = NCOUNT + 1
            RPARAM(ICP,ISTEP) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       STEP SPECIES-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          NSSLEN(ISTEP) = NINT(PARRAY(NCOUNT))
          DO ISPEC = 1, NSSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            NSSPEC(ISPEC,ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDDO

C       STEP REACTANT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          NRSLEN(ISTEP) = NINT(PARRAY(NCOUNT))
          DO ISPEC = 1, NRSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            NRSPEC(ISPEC,ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDDO

C       STEP PRODUCT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          NPSLEN(ISTEP) = NINT(PARRAY(NCOUNT))
          DO ISPEC = 1, NPSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            NPSPEC(ISPEC,ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDDO

C       STEP REACTANT COEFFICIENT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          NRCLEN(ISTEP) = NINT(PARRAY(NCOUNT))
          DO ISPEC = 1, NRCLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            NRCPEC(ISPEC,ISTEP) = NINT(PARRAY(NCOUNT))
            NCOUNT = NCOUNT + 1
            CRSPEC(ISPEC,ISTEP) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       STEP PRODUCT COEFFICIENT-LIST
        DO ISTEP = 1,NSTEP
          NCOUNT = NCOUNT + 1
          NPCLEN(ISTEP) = NINT(PARRAY(NCOUNT))
          DO ISPEC = 1, NPCLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            NPCPEC(ISPEC,ISTEP) = NINT(PARRAY(NCOUNT))
            NCOUNT = NCOUNT + 1
            CPSPEC(ISPEC,ISTEP) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       SPECIES DELTA-LIST
        DO ISTEP = 1,NSTEP
          DO ISPEC = 1, NSSLEN(ISTEP)
            NCOUNT = NCOUNT + 1
            DIFFMU(ISPEC,ISTEP) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       THIRD-BODY LIST
        NCOUNT = NCOUNT + 1
        NBODY = NINT(PARRAY(NCOUNT))
        DO IBODY = 1, NBODY
          DO ICP = 1, NSPSTR
            NCOUNT = NCOUNT + 1
            BDYSYM(IBODY)(ICP:ICP) = CHAR(NINT(PARRAY(NCOUNT)))
          ENDDO
        ENDDO

C       THIRD-BODY STEP-LIST
        IF(NBODY.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            MBLIST(ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDIF

C       THIRD-BODY EFFICIENCIES
        DO IBODY = 1, NBODY
          DO ISPEC = 1,NSPEC
            NCOUNT = NCOUNT + 1
            EFFY3B(ISPEC,IBODY) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       GIBBS STEP-LIST
        NCOUNT = NCOUNT + 1
        NGIBB = NINT(PARRAY(NCOUNT))
        IF(NGIBB.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            MGLIST(ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDIF

C       -----------------------------------------------------------------------

C       LINDEMANN STEP-LIST
        NCOUNT = NCOUNT + 1
        NLIND = NINT(PARRAY(NCOUNT))
        IF(NLIND.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            MLLIST(ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDIF

C       LINDEMANN STEP RATE DATA
        DO ILIND = 1, NLIND
          DO ICP = 1,4
            NCOUNT = NCOUNT + 1
            RCLIND(ICP,ILIND) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       TROE FORM STEP-LIST
        NCOUNT = NCOUNT + 1
        NTROE = NINT(PARRAY(NCOUNT))
        IF(NTROE.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            MTLIST(ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDIF

C       TROE FORM STEP RATE DATA
        DO ITROE = 1, NTROE
          DO ICP = 1,12
            NCOUNT = NCOUNT + 1
            RCTROE(ICP,ITROE) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       SRI FORM STEP-LIST
        NCOUNT = NCOUNT + 1
        NSRIF = NINT(PARRAY(NCOUNT))
        IF(NSRIF.GT.0)THEN
          DO ISTEP = 1, NSTEP
            NCOUNT = NCOUNT + 1
            MSLIST(ISTEP) = NINT(PARRAY(NCOUNT))
          ENDDO
        ENDIF

C       SRI FORM STEP RATE DATA
        DO ISRIF = 1, NSRIF
          DO ICP = 1,8
            NCOUNT = NCOUNT + 1
            RCSRIF(ICP,ISRIF) = PARRAY(NCOUNT)
          ENDDO
        ENDDO

C       =======================================================================

      ENDIF

C     =========================================================================

C     EVALUATE DERIVED QUANTITIES
C     ===========================

C     NO OF ACTIVE SPECIES
      NSPM1 = NSPEC - 1

C     -------------------------------------------------------------------------

C     CONVERT RATE PARAMETERS
      DO ISTEP = 1,NSTEP
        RPARAM(1,ISTEP) = LOG(RPARAM(1,ISTEP))
        RPARAM(3,ISTEP) = RPARAM(3,ISTEP)/RGUNIV
      ENDDO

C     LINDEMANN STEP RATE DATA
      DO ILIND = 1, NLIND
        RCLIND(1,ILIND) = LOG(RCLIND(1,ILIND))
        RCLIND(3,ILIND) = RCLIND(3,ILIND)/RGUNIV
      ENDDO

C     TROE FORM STEP RATE DATA
      DO ITROE = 1, NTROE
        RCTROE(1,ITROE) = LOG(RCTROE(1,ITROE))
        RCTROE(3,ITROE) = RCTROE(3,ITROE)/RGUNIV
        RCTROE(5,ITROE) = -ONE/RCTROE(5,ITROE)
C       VM: CHANGE based on Khalil
C        RCTROE(5,ITROE) = ONE/RCTROE(5,ITROE)
        RCTROE(7,ITROE) = -ONE/RCTROE(7,ITROE)
      ENDDO

C     SRI FORM STEP RATE DATA
      DO ISRIF = 1, NSRIF
        RCSRIF(1,ISRIF) = LOG(RCSRIF(1,ISRIF))
        RCSRIF(3,ISRIF) = RCSRIF(3,ISRIF)/RGUNIV
        RCSRIF(5,ISRIF) = -RCSRIF(5,ISRIF)
        RCSRIF(6,ISRIF) = -ONE/RCSRIF(6,ISRIF)
      ENDDO

C     -------------------------------------------------------------------------

C     STOICHIOMETRIC COEFFICIENTS TIMES MOLAR MASS
      DO ISTEP = 1,NSTEP
        DO ISPEC = 1,NSSLEN(ISTEP)
          JSPEC = NSSPEC(ISPEC,ISTEP)
          DIFFMW(ISPEC,ISTEP) = DIFFMU(ISPEC,ISTEP)*WMOLAR(JSPEC)
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------

C     RECIPROCAL OF MOLAR MASS
C     SPECIFIC GAS CONSTANT
      DO ISPEC = 1,NSPEC
        OVWMOL(ISPEC) = ONE/WMOLAR(ISPEC)
        RGSPEC(ISPEC) = RGUNIV*OVWMOL(ISPEC)
      ENDDO

C     -------------------------------------------------------------------------

C     SPECIFIC HEAT CAPACITY
C     ----------------------

C     NUMBER OF CP POLYNOMIAL COEFFICIENTS
C     NUMBER OF CP POLYNOMIAL COEFFICIENTS MINUS ONE
C     INDEX NUMBERS OF ENTHALPY AND ENTROPY COEFFICIENTS
      DO ISPEC = 1,NSPEC
        DO ITINT = 1,NTINT(ISPEC)
          NCPOLY(ITINT,ISPEC) = NCOFCP(ITINT,ISPEC)-2
          NCPOM1(ITINT,ISPEC) = NCPOLY(ITINT,ISPEC)-1
          NCENTH(ITINT,ISPEC) = NCOFCP(ITINT,ISPEC)-1
          NCENPY(ITINT,ISPEC) = NCOFCP(ITINT,ISPEC)
        ENDDO
      ENDDO

C     SPECIFIC HEAT CAPACITY CP PER UNIT MASS
      DO ISPEC = 1,NSPEC
        DO ITINT = 1,NTINT(ISPEC)
          DO ICP = 1,NCOFCP(ITINT,ISPEC)
            AMASCP(ICP,ITINT,ISPEC)
     +      = AMOLCP(ICP,ITINT,ISPEC)*RGUNIV*OVWMOL(ISPEC)
          ENDDO
        ENDDO
      ENDDO

C     COEFFICIENTS FOR TEMPERATURE
      DO ISPEC = 1,NSPEC
        DO ITINT = 1,NTINT(ISPEC)
          AMASCT(1,ITINT,ISPEC) = AMASCP(1,ITINT,ISPEC) - RGSPEC(ISPEC)
          DO ICP = 2,NCPOLY(ITINT,ISPEC)
            AMASCT(ICP,ITINT,ISPEC) = AMASCP(ICP,ITINT,ISPEC)/REAL(ICP)
          ENDDO
          AMASCT(NCENTH(ITINT,ISPEC),ITINT,ISPEC) = ZERO
          AMASCT(NCENPY(ITINT,ISPEC),ITINT,ISPEC) = ZERO
        ENDDO
      ENDDO

C     COEFFICIENTS FOR ENTHALPY PER UNIT MASS
      DO ISPEC = 1,NSPEC
        DO ITINT = 1,NTINT(ISPEC)
          AMASCH(1,ITINT,ISPEC) = AMASCP(1,ITINT,ISPEC)
          DO ICP = 2,NCPOLY(ITINT,ISPEC)
            AMASCH(ICP,ITINT,ISPEC) = AMASCP(ICP,ITINT,ISPEC)/REAL(ICP)
          ENDDO
          AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +    = AMASCP(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
          AMASCH(NCENPY(ITINT,ISPEC),ITINT,ISPEC) = ZERO
        ENDDO
      ENDDO

C     COEFFICIENTS FOR ENTROPY PER UNIT MASS
      DO ISPEC = 1,NSPEC
        DO ITINT = 1,NTINT(ISPEC)
          AMASCS(1,ITINT,ISPEC) = AMASCP(1,ITINT,ISPEC)
          DO ICP = 2,NCPOLY(ITINT,ISPEC)
            AMASCS(ICP,ITINT,ISPEC) = AMASCP(ICP,ITINT,ISPEC)
     +                               /REAL(ICP-1)
          ENDDO
          AMASCS(NCENTH(ITINT,ISPEC),ITINT,ISPEC) = ZERO
          AMASCS(NCENPY(ITINT,ISPEC),ITINT,ISPEC)
     +    = AMASCP(NCENPY(ITINT,ISPEC),ITINT,ISPEC)
        ENDDO
      ENDDO

C     COEFFICIENTS FOR GIBBS FUNCTION PER MOLE
C     ACTUALLY GIBBS/(R^0 T) WITH PRESSURE TERM
      GIBLET = LOG(PREFGB/RGUNIV)
      DO ISPEC = 1,NSPEC
        DO ITINT = 1,NTINT(ISPEC)
          AMOLGB(1,ITINT,ISPEC)= AMOLCP(NCENPY(ITINT,ISPEC),ITINT,ISPEC)
     +                         - AMOLCP(1,ITINT,ISPEC)
     +                         + GIBLET
          DO ICP = 2,NCPOLY(ITINT,ISPEC)
            AMOLGB(ICP,ITINT,ISPEC) = AMOLCP(ICP,ITINT,ISPEC)
     +                               /REAL(ICP*(ICP-1))
          ENDDO
          AMOLGB(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +         = AMOLCP(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
          AMOLGB(NCENPY(ITINT,ISPEC),ITINT,ISPEC)
     +         = AMOLCP(1,ITINT,ISPEC) - ONE
        ENDDO
      ENDDO

C     SPECIFIC HEAT CAPACITY CP PER MOLE
      DO ISPEC = 1,NSPEC
        DO ITINT = 1,NTINT(ISPEC)
          DO ICP = 1,NCOFCP(ITINT,ISPEC)
            AMOLCP(ICP,ITINT,ISPEC) = AMOLCP(ICP,ITINT,ISPEC)*RGUNIV
          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------

C     RECIPROCAL OF LEWIS NUMBER
      DO ISPEC = 1,NSPEC
        OLEWIS(ISPEC) = ONE/CLEWIS(ISPEC)
      ENDDO

C     -------------------------------------------------------------------------

C     CONDUCTIVITY COEFFICIENT
      ALAMDA = ALAMDC*EXP(-RLAMDA*LOG(TLAMDA))

C     =========================================================================


      RETURN
      END
