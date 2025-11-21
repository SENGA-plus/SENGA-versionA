      SUBROUTINE BUFTXL 
 
C     *************************************************************************
C
C     BUFTXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     02-APR-2006:  CREATED
C 
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INLET TURBULENT VELOCITY FIELD
C     FORWARD FFT WITH FOLDING IN X
C     
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER JPEN(NPENMX),KPEN(NPENMX)
      INTEGER IC,JC,KC,IPENCL
      INTEGER IIM,IIC
      INTEGER NPENCL


C     BEGIN
C     =====

C     =========================================================================

C     CARRY OUT A FORWARD FOURIER TRANSFORM
C     =====================================

C     =========================================================================

C     X-DIRECTION
C     -----------

C     PENCIL COUNTER
      IPENCL = 0

C     LOOP OVER X-LINES
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL

C         PENCIL INDEXING
          IPENCL = IPENCL + 1
          JPEN(IPENCL) = JC
          KPEN(IPENCL) = KC

C         ASSEMBLE FFT DATA
          DO IC = ISTAL,ISTOL

            IIC = 2*IC
            IIM = IIC-1
            FTPART(IIM,1,IPENCL) = URUN(IC,JC,KC)        
            FTPART(IIC,1,IPENCL) = UTMP(IC,JC,KC)        
            FTPART(IIM,2,IPENCL) = VRUN(IC,JC,KC)        
            FTPART(IIC,2,IPENCL) = VTMP(IC,JC,KC)        
            FTPART(IIM,3,IPENCL) = WRUN(IC,JC,KC)        
            FTPART(IIC,3,IPENCL) = WTMP(IC,JC,KC)        

          ENDDO

C         CHECK FOR MAX NO OF PENCILS
          IF(IPENCL.EQ.NPENMX)THEN

C           DO THE FFT
            NPENCL = IPENCL
            CALL FFTIXL(NPENCL,IXPROC,NXPROC,NXNODE,NXGLBL,NPROCX)

C           RESTORE TRANSFORMED DATA
            DO IPENCL = 1,NPENCL

              DO IC = ISTAL,ISTOL

                UFXL(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                             = FTPART(IC,1,IPENCL)
                VFXL(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                             = FTPART(IC,2,IPENCL)
                WFXL(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                             = FTPART(IC,3,IPENCL)

              ENDDO

            ENDDO
            IPENCL = 0

          ENDIF

        ENDDO
      ENDDO
C     END OF LOOP OVER X-LINES

C     CHECK FOR ANY REMAINING STORED PENCILS
      IF(IPENCL.NE.0)THEN

C       DO THE FFT
        NPENCL = IPENCL
        CALL FFTIXL(NPENCL,IXPROC,NXPROC,NXNODE,NXGLBL,NPROCX)

C       RESTORE TRANSFORMED DATA
        DO IPENCL = 1,NPENCL

          DO IC = ISTAL,ISTOL

            UFXL(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IC,1,IPENCL)
            VFXL(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IC,2,IPENCL)
            WFXL(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IC,3,IPENCL)

          ENDDO

        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
