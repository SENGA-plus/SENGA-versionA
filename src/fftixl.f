      SUBROUTINE FFTIXL(NPENCL,ILPROC,NLPROC,NLSIZE,NGSIZE,NPROCS)

C     *************************************************************************
C
C     FFTIXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT
C
C     CHANGE RECORD
C     -------------
C     07-APR-2006:  CREATED
C     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     FORWARD FFT FOR PARALLEL INVERSE DFT FOR INLET VELOCITY FIELD
C     X-LEFT
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETER
C     =========
      INTEGER IFORW
      PARAMETER(IFORW = 1)


C     ARGUMENTS
C     =========
      INTEGER NPENCL,ILPROC,NLPROC,NLSIZE,NGSIZE
      INTEGER NPROCS(0:NLPROC)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION PCOUNT
      INTEGER NSTOTL(0:NPRMAX)
      INTEGER NSIZE2,IOFSET
      INTEGER NLPRM1,NLSIZ2,NGSIZ2,NCOUNT
      INTEGER IRPROC,IRTAG
      INTEGER IPENCL,IDIREC,IX,ICOUNT,ICPROC


C     BEGIN
C     =====

C     =========================================================================

      NLPRM1 = NLPROC-1
      NLSIZ2 = NLSIZE*2
      NGSIZ2 = NGSIZE*2

C     =========================================================================

C     CHECK FOR LEFTMOST PROCESSOR
      IF(ILPROC.EQ.0) THEN

C       =======================================================================

C       LEFTMOST PROCESSOR OF THE ROW
C       -----------------------------

C       GATHER ALL PHYSICAL-SPACE DATA AND DO THE FFT
 
C       PHYSICAL-SPACE DATA ON LOCAL PROCESSOR
        DO IPENCL = 1, NPENCL
          DO IDIREC = 1,3
            DO IX = 1,NLSIZ2

              FFTROW(IX,IDIREC,IPENCL) = FTPART(IX,IDIREC,IPENCL)

            ENDDO
          ENDDO
        ENDDO
        NSTOTL(0) = NLSIZ2

C       RECEIVE PHYSICAL-SPACE DATA FROM REMOTE PROCESSORS
C       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
        IOFSET = 0
        DO ICPROC = 1,NLPRM1

          IRPROC = NPROCS(ICPROC)
          IRTAG = IRPROC*NPROC+IPROC

          CALL P_RECV(PCOUNT,1,IRPROC,IRTAG)
          CALL P_RECV(PARRAY,NPARAY,IRPROC,IRTAG)
 
          NCOUNT = INT(PCOUNT)
          NSIZE2 = NCOUNT/NPENCL/3
          IOFSET = IOFSET + NSTOTL(ICPROC-1)

          ICOUNT = 0
          DO IPENCL = 1,NPENCL
            DO IDIREC = 1,3
              DO IX = 1,NSIZE2

                ICOUNT = ICOUNT + 1
                FFTROW(IX+IOFSET,IDIREC,IPENCL) = PARRAY(ICOUNT)

              ENDDO
            ENDDO
          ENDDO
          NSTOTL(ICPROC) = NSIZE2

        ENDDO

C       -----------------------------------------------------------------------

C       DO THE FOURIER TRANSFORM
        CALL FFTGIN(NGSIZE,IFORW)

        DO IPENCL = 1,NPENCL
          DO IDIREC = 1,3

            DO IX = 1,NGSIZ2
              FFTINX(IX) = FFTROW(IX,IDIREC,IPENCL)
            ENDDO

            CALL FFTGEN(FFTINX,NGSIZE)

            DO IX = 1,NGSIZ2
              FFTROW(IX,IDIREC,IPENCL) = FFTINX(IX)
            ENDDO

          ENDDO
        ENDDO
         
C       -----------------------------------------------------------------------

C       COMPRESS FT BY REMOVING IMAGINARY PART OF ZEROTH FOURIER COEFFICIENT
        DO IPENCL = 1, NPENCL
          DO IDIREC = 1,3

            DO IX = 2,NGSIZE
              FFTROW(IX,IDIREC,IPENCL) = FFTROW(IX+1,IDIREC,IPENCL)
            ENDDO
            
          ENDDO
        ENDDO

C       -----------------------------------------------------------------------

C       FOURIER-SPACE DATA ON LOCAL PROCESSOR
        DO IPENCL = 1, NPENCL
          DO IDIREC = 1,3

            DO IX = 1,NLSIZE
              FTPART(IX,IDIREC,IPENCL) = FFTROW(IX,IDIREC,IPENCL)
            ENDDO
            
          ENDDO
        ENDDO

C       SEND FOURIER-SPACE DATA BACK TO REMOTE PROCESSORS
C       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
        IOFSET = 0
        DO ICPROC = 1,NLPRM1

          IOFSET = IOFSET + NSTOTL(ICPROC-1)/2
          NSIZE2 = NSTOTL(ICPROC)/2

          ICOUNT = 0
          DO IPENCL = 1,NPENCL
            DO IDIREC = 1,3
              DO IX = 1,NSIZE2

                ICOUNT = ICOUNT + 1
                PARRAY(ICOUNT) = FFTROW(IX+IOFSET,IDIREC,IPENCL)

              ENDDO
            ENDDO
          ENDDO

          NCOUNT = ICOUNT
          PCOUNT = REAL(NCOUNT)
          IRPROC = NPROCS(ICPROC)
          IRTAG = IPROC*NPROC+IRPROC
          CALL P_SEND(PCOUNT,1,1,IRPROC,IRTAG)
          CALL P_SEND(PARRAY,NPARAY,NCOUNT,IRPROC,IRTAG)

        ENDDO

C       =======================================================================

      ELSE

C       =======================================================================

C       NOT THE LEFTMOST PROCESSOR
C       --------------------------
C       IDENTIFY THE LEFTMOST PROCESSOR
        IRPROC = NPROCS(0)

C       SEND THE PHYSICAL-SPACE DATA TO THE LEFTMOST PROCRESSOR
C       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
        ICOUNT = 0
        DO IPENCL = 1,NPENCL
          DO IDIREC = 1,3
            DO IX = 1,NLSIZ2

              ICOUNT = ICOUNT + 1
              PARRAY(ICOUNT) = FTPART(IX,IDIREC,IPENCL)

            ENDDO
          ENDDO
        ENDDO

        NCOUNT = ICOUNT
        PCOUNT = REAL(ICOUNT)
        IRTAG = IPROC*NPROC+IRPROC
        CALL P_SEND(PCOUNT,1,1,IRPROC,IRTAG)
        CALL P_SEND(PARRAY,NPARAY,NCOUNT,IRPROC,IRTAG)

C       RECEIVE THE FOURIER-SPACE DATA FROM THE LEFTMOST PROCESSOR
C       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
        IRTAG = IRPROC*NPROC+IPROC
        CALL P_RECV(PCOUNT,1,IRPROC,IRTAG)
        CALL P_RECV(PARRAY,NPARAY,IRPROC,IRTAG)

        ICOUNT = 0
        DO IPENCL = 1,NPENCL
          DO IDIREC = 1,3
            DO IX = 1,NLSIZE

              ICOUNT = ICOUNT + 1
              FTPART(IX,IDIREC,IPENCL) = PARRAY(ICOUNT)

            ENDDO
          ENDDO
        ENDDO

C       =======================================================================

      ENDIF
C     LEFTMOST PROCESSOR

C     =========================================================================


      RETURN
      END
