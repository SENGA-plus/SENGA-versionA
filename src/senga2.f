      PROGRAM SENGA2
 
C     *************************************************************************
C
C     SENGA2
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     12-NOV-2002:  CREATED
C     30-JUN-2006:  RSC SENGA2 VERSION 1.0 (beta)
C     29-DEC-2006:  RSC UPDATES TO INDEXING
C     24-SEP-2009:  RSC SENGA2 VERSION 1.1: UPDATES AND BUG FIXES
C     08-AUG-2012:  RSC SENGA2 VERSION 1.2: UPDATES AND BUG FIXES
C     08-MAY-2013:  RSC SENGA2 VERSION 2.0: UPDATES AND BUG FIXES
C     01-JUL-2015:  RSC SENGA2 VERSION 2.1: UPDATES AND BUG FIXES
C
C     BASED CLOSELY ON SENGA, CREATED 01-AUG-1996
C     01-AUG-1997:  RSC PARALLEL VERSION 1.0
C     1997-2002:    DR KARL JENKINS: PARALLEL DNS OF FLAME KERNELS
C     2001-2005:    DR NILAN CHAKRABORTY: INFLOW BOUNDARY CONDITIONS
C     2008-2009:    DR TOM DUNSTAN: BUG FIXES
C     2009-2012:    PRESSURE-DEPENDENT REACTION RATES
C     2012:         RSC/ZACK NICOLAOU BUG FIXES
C     2013:         RSC MIXTURE AVERAGED TRANSPORT
C     2015:         RSC/ALEX PHILPOTT UPDATED WALL BOUNDARY CONDITIONS
C
C     DESCRIPTION
C     -----------
C     THIS PROGRAM CARRIES OUT DIRECT NUMERICAL SIMULATION
C     OF COMPRESSSIBLE TURBULENT FLOW AND COMBUSTION
C     USING HIGH-ORDER METHODS IN TIME AND SPACE
C
C     PROVISION IS MADE FOR MULTIPLE SPECIES MASS FRACTIONS
C                           VARIABLE TRANSPORT PROPERTIES
C                           MULTIPLE STEP CHEMISTRY
C
C     REFERENCES
C     ==========
C     1) Kennedy, C.A., Carpenter, M.H., Lewis, R.M.: "Low-storage
C        Explicit Runge-Kutta Schemes for the Compressible Navier-Stokes
C        Equations", Appl. Numer. Math. 35 (3) 177-219, 2000.
C
C     2) Lele, S.K.:  "Compact Finite Difference Schemes with Spectral-like
C                      Resolution", J. Comput. Phys. 103, 16-29, 1992.
C
C     3) Poinsot, T.J., Lele, S.K.: "Boundary Conditions for Direct
C        Simulations of Compressible Viscous Flow", J. Comput. Phys. 101,
C        104-129, 1992.
C
C     4) Sutherland, J.C., Kennedy, C.A.: "Improved Boundary Conditions for
C        Viscous Reacting Compressible Flows", J. Comput. Phys. 2004.
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
C     RSC 29-DEC-2006 UPDATED INDEXING
      INTEGER JTIME,JRKSTP


C     BEGIN
C     =====

C     =========================================================================
 
C     INITIALISATION
C     ==============
C     PARALLEL DOMAIN DECOMPOSITION
      CALL PARDOM

C     INITIALISE THE DATA
      CALL INDATA

C     RECORD INITIAL CONDITIONS
      CALL OUTPUT

C     =========================================================================

C     TIME STEP LOOP
C     ==============

C     RSC 29-DEC-2006 UPDATED INDEXING
      DO JTIME = NTIME1,NTIME2

        ITIME = JTIME

        IF(NXLPRM(1).EQ.4) THEN
          CALL INFLOW
          CALL P_SYNC
        END IF

C       =======================================================================

C       RUNGE-KUTTA SUBSTEPS
C       ====================

C       STANDARD SUBSTEPS
C       -----------------
C       RSC 29-DEC-2006 UPDATED INDEXING
        DO JRKSTP = 1,NRKSM1

          IRKSTP = JRKSTP
        
C         APPLY BCS ON PRIMITIVE VARIABLES
          CALL BOUNDT

C         PARALLEL DATA TRANSFER
          CALL PARFER

C         EVALUATE RHS FOR SCALARS
          CALL RHSCAL

C         EVALUATE RHS FOR VELOCITIES
          CALL RHSVEL

C         APPLY BCS ON SOURCE TERMS
          CALL BOUNDS

C         RUNGE-KUTTA ADVANCEMENT
          CALL LINCOM

        ENDDO

C       =======================================================================
  
C       LAST SUBSTEP IS DIFFERENT
C       -------------------------
        IRKSTP = NRKSTP

C       APPLY BCS ON PRIMITIVE VARIABLES
        CALL BOUNDT

C       PARALLEL DATA TRANSFER
        CALL PARFER

C       EVALUATE RHS FOR SCALARS
        CALL RHSCAL

C       EVALUATE RHS FOR VELOCITIES
        CALL RHSVEL

C       APPLY BCS ON SOURCE TERMS
        CALL BOUNDS

C       RUNGE-KUTTA ADVANCEMENT
        CALL FINCOM

C       =======================================================================

C       UPDATE THE ELAPSED TIME
C       =======================
        ETIME = ETIME + TSTEP

C       =======================================================================

C       SYNCHRONISE THE TIME-DEPENDENT BCS
C       ==================================
        CALL BOUNTT

C       =======================================================================

C       FILTER THE SOLUTION
C       ===================
C       RSC 30-AUG-2009 HIGH ORDER FILTERING
C        CALL FLTREM

C       =======================================================================

C       ADJUST THE TIME STEP
C       ====================
        CALL ADAPTT

C       =======================================================================

C       PROCESS THE RESULTS
C       ===================
        CALL OUTPUT

C       =======================================================================

      ENDDO
C     END OF TIME STEP LOOP

C     =========================================================================
 
C     TERMINATION
C     ===========

C     TERMINATE THE PROGRAM
      CALL FINISH

C     =========================================================================


      STOP
      END
