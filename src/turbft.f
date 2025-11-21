      SUBROUTINE TURBFT 
 
C     *************************************************************************
C
C     TURBFT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     28-MAR-1997:  CREATED
C     15-MAY-1999:  KWJ PARALLEL IMPLEMENTATION (INDEPENDENT)
C     20-OCT-1999:  KWJ PARALLEL FFT IMPLEMENTATION
C     24-MAY-2003:  RSC UPDATED FOR SENGA2
C 
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIAL TURBULENCE GENERATOR
C     INVERSE FFT WITH CONJUGATE ANTISYMMETRY
C     
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION SCALEX,SCALEY,SCALEZ
      INTEGER IPEN(NPENMX),JPEN(NPENMX),KPEN(NPENMX)
      INTEGER IC,JC,KC,IX,JX,KX,IPENCL,ICPROC
      INTEGER IIM,IIC,JJM,JJC,KKM,KKC
      INTEGER IGOFST,JGOFST,KGOFST,IGOFM1,JGOFM1,KGOFM1
      INTEGER NODBLX,NODBLY,NODBLZ
      INTEGER NDOSYM,NPENCL,NCOUNT


C     BEGIN
C     =====

C     =========================================================================

C     INDEXING
C     --------

C     SET ODDBALL WAVENUMBER INDICES
      NODBLX = NXGLBL/2
      NODBLY = NYGLBL/2
      NODBLZ = NZGLBL/2

C     PHYSICAL-SPACE GLOBAL INDEX OFFSETS
      IGOFST = 0
      DO ICPROC = 0, IXPROC-1
        IGOFST = IGOFST + NPMAPX(ICPROC)
      ENDDO
      IGOFM1 = IGOFST-1

      JGOFST = 0
      DO ICPROC = 0, IYPROC-1
        JGOFST = JGOFST + NPMAPY(ICPROC)
      ENDDO
      JGOFM1 = JGOFST-1

      KGOFST = 0
      DO ICPROC = 0, IZPROC-1
        KGOFST = KGOFST + NPMAPZ(ICPROC)
      ENDDO
      KGOFM1 = KGOFST-1

C     =========================================================================

C     SET SCALE FACTORS
      SCALEX = ONE
      SCALEY = ONE
      SCALEZ = ONE

C     =========================================================================

C     CARRY OUT AN INVERSE FOURIER TRANSFORM
C     ======================================

C     CHECK FFT PARALLEL DATA SIZE AGAINST PARALLEL ARRAY SIZE
      IF(IPROC.EQ.0)THEN
        NCOUNT = 2*NSZMAX*3*NPENMX
        IF(NCOUNT.GT.NPARAY)THEN
        WRITE(NCREPT,*)'Warning: TURBFT: parallel array size too small'
        ENDIF
      ENDIF

C     =========================================================================

C     X-DIRECTION
C     -----------

C     PENCIL COUNTER
      IPENCL = 0

C     LOOP OVER X-LINES
      DO KC = KSTAL,KSTOL

C       FOURIER-SPACE GLOBAL INDEXING
        KX = KGOFM1 + KC
        IF(KX.GT.NODBLZ)KX = KX-NZGLBL

        DO JC = JSTAL,JSTOL

C         FOURIER-SPACE GLOBAL INDEXING
          JX = JGOFM1 + JC
          IF(JX.GT.NODBLY)JX = JX-NYGLBL

C         TRANSFORM ONLY FOR UPPER-CENTRAL WAVENUMBERS
          IF(KX.GT.0)THEN

C           STANDARD UPPER-CENTRAL X-LINES

C           PENCIL INDEXING
            IPENCL = IPENCL + 1
            JPEN(IPENCL) = JC
            KPEN(IPENCL) = KC

C           ASSEMBLE FFT DATA
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

C           CHECK FOR MAX NO OF PENCILS
            IF(IPENCL.EQ.NPENMX)THEN

C             DO THE FFT (NO SYMMETRY IMPOSED)
              NDOSYM = 0
              NPENCL = IPENCL

          CALL FFTSYM(NDOSYM,NPENCL,IXPROC,NXPROC,NXNODE,NXGLBL,NPROCX)

C             RESTORE TRANSFORMED DATA
              DO IPENCL = 1,NPENCL

                DO IC = ISTAL,ISTOL

                  IIC = 2*IC
                  IIM = IIC-1
                  URUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                               = FTPART(IIM,1,IPENCL)*SCALEX
                  UTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                               = FTPART(IIC,1,IPENCL)*SCALEX
                  VRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                               = FTPART(IIM,2,IPENCL)*SCALEX
                  VTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                               = FTPART(IIC,2,IPENCL)*SCALEX
                  WRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                               = FTPART(IIM,3,IPENCL)*SCALEX
                  WTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                               = FTPART(IIC,3,IPENCL)*SCALEX

                ENDDO

              ENDDO
              IPENCL = 0

            ENDIF

          ELSE IF(KX.EQ.0)THEN

            IF(JX.GT.0)THEN

C             KX = 0 PLANE UPPER-CENTRAL X-LINES

C             PENCIL INDEXING
              IPENCL = IPENCL + 1
              JPEN(IPENCL) = JC
              KPEN(IPENCL) = KC

C             ASSEMBLE FFT DATA
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

C             CHECK FOR MAX NO OF PENCILS
              IF(IPENCL.EQ.NPENMX)THEN

C               DO THE FFT (NO SYMMETRY IMPOSED)
                NDOSYM = 0
                NPENCL = IPENCL

          CALL FFTSYM(NDOSYM,NPENCL,IXPROC,NXPROC,NXNODE,NXGLBL,NPROCX)

C               RESTORE TRANSFORMED DATA
                DO IPENCL = 1,NPENCL

                  DO IC = ISTAL,ISTOL

                    IIC = 2*IC
                    IIM = IIC-1
                    URUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIM,1,IPENCL)*SCALEX
                    UTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIC,1,IPENCL)*SCALEX
                    VRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIM,2,IPENCL)*SCALEX
                    VTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIC,2,IPENCL)*SCALEX
                    WRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIM,3,IPENCL)*SCALEX
                    WTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIC,3,IPENCL)*SCALEX

                  ENDDO

                ENDDO
                IPENCL = 0

              ENDIF

            ELSE IF(JX.EQ.0)THEN

C             SPECIAL CASE OF KY=KZ=0

C             CHECK FOR PRE-STORED PENCILS
              IF(IPENCL.NE.0)THEN

C               DO THE FFT (NO SYMMETRY IMPOSED)
                NDOSYM = 0
                NPENCL = IPENCL

           CALL FFTSYM(NDOSYM,NPENCL,IXPROC,NXPROC,NXNODE,NXGLBL,NPROCX)

C               RESTORE TRANSFORMED DATA
                DO IPENCL = 1,NPENCL

                  DO IC = ISTAL,ISTOL

                    IIC = 2*IC
                    IIM = IIC-1
                    URUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIM,1,IPENCL)*SCALEX
                    UTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIC,1,IPENCL)*SCALEX
                    VRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                = FTPART(IIM,2,IPENCL)*SCALEX
                    VTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIC,2,IPENCL)*SCALEX
                    WRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIM,3,IPENCL)*SCALEX
                    WTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                                 = FTPART(IIC,3,IPENCL)*SCALEX

                  ENDDO

                ENDDO

              ENDIF

C             ASSEMBLE FFT DATA FOR KY=KZ=0
              IPENCL = 1

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

C             DO THE FFT (SYMMETRY IMPOSED)
              NDOSYM = 1
              NPENCL = 1

          CALL FFTSYM(NDOSYM,NPENCL,IXPROC,NXPROC,NXNODE,NXGLBL,NPROCX)

C             RESTORE TRANSFORMED DATA
              DO IC = ISTAL,ISTOL

                IIC = 2*IC
                IIM = IIC-1

                URUN(IC,JC,KC) = FTPART(IIM,1,IPENCL)*SCALEX
                UTMP(IC,JC,KC) = FTPART(IIC,1,IPENCL)*SCALEX
                VRUN(IC,JC,KC) = FTPART(IIM,2,IPENCL)*SCALEX
                VTMP(IC,JC,KC) = FTPART(IIC,2,IPENCL)*SCALEX
                WRUN(IC,JC,KC) = FTPART(IIM,3,IPENCL)*SCALEX
                WTMP(IC,JC,KC) = FTPART(IIC,3,IPENCL)*SCALEX

              ENDDO
              IPENCL = 0

            ELSE

C             KX = 0 JX < 0 LOWER-CENTRAL WAVENUMBERS - NO ACTION
              CONTINUE

            ENDIF

          ELSE

C           KX < 0 LOWER-CENTRAL WAVENUMBERS - NO ACTION
            CONTINUE

          ENDIF

        ENDDO
      ENDDO
C     END OF LOOP OVER X-LINES

C     CHECK FOR ANY REMAINING STORED PENCILS
      IF(IPENCL.NE.0)THEN

C       DO THE FFT (NO SYMMETRY IMPOSED)
        NDOSYM = 0
        NPENCL = IPENCL

        CALL FFTSYM(NDOSYM,NPENCL,IXPROC,NXPROC,NXNODE,NXGLBL,NPROCX)

C       RESTORE TRANSFORMED DATA
        DO IPENCL = 1,NPENCL

          DO IC = ISTAL,ISTOL

            IIC = 2*IC
            IIM = IIC-1
            URUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IIM,1,IPENCL)*SCALEX
            UTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IIC,1,IPENCL)*SCALEX
            VRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IIM,2,IPENCL)*SCALEX
            VTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IIC,2,IPENCL)*SCALEX
            WRUN(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IIM,3,IPENCL)*SCALEX
            WTMP(IC,JPEN(IPENCL),KPEN(IPENCL))
     +                         = FTPART(IIC,3,IPENCL)*SCALEX

          ENDDO

        ENDDO

      ENDIF

C     =========================================================================

C     Y-DIRECTION
C     -----------

C     PENCIL COUNTER
      IPENCL = 0

      DO KC = KSTAL,KSTOL

C       FOURIER-SPACE GLOBAL INDEXING
        KX = KGOFM1 + KC
        IF(KX.GT.NODBLZ)KX = KX-NZGLBL

C       TRANSFORM ONLY FOR UPPER-CENTRAL WAVENUMBERS
        IF(KX.GT.0)THEN

          DO IC = ISTAL,ISTOL

C           FOURIER-SPACE GLOBAL INDEXING
            IX = IGOFM1 + IC
            IF(IX.GT.NODBLX)IX = IX-NXGLBL

C           PENCIL INDEXING
            IPENCL = IPENCL + 1
            IPEN(IPENCL) = IC
            KPEN(IPENCL) = KC

C           ASSEMBLE FFT DATA
            DO JC = JSTAL,JSTOL

              JJC = 2*JC
              JJM = JJC-1
              FTPART(JJM,1,IPENCL) = URUN(IC,JC,KC)        
              FTPART(JJC,1,IPENCL) = UTMP(IC,JC,KC)        
              FTPART(JJM,2,IPENCL) = VRUN(IC,JC,KC)        
              FTPART(JJC,2,IPENCL) = VTMP(IC,JC,KC)        
              FTPART(JJM,3,IPENCL) = WRUN(IC,JC,KC)        
              FTPART(JJC,3,IPENCL) = WTMP(IC,JC,KC)        

            ENDDO

C           CHECK FOR MAX NO OF PENCILS
            IF(IPENCL.EQ.NPENMX)THEN

C             DO THE FFT (NO SYMMETRY IMPOSED)
              NDOSYM = 0
              NPENCL = IPENCL

          CALL FFTSYM(NDOSYM,NPENCL,IYPROC,NYPROC,NYNODE,NYGLBL,NPROCY)

C             RESTORE TRANSFORMED DATA
              DO IPENCL = 1,NPENCL

                DO JC = JSTAL,JSTOL

                  JJC = 2*JC
                  JJM = JJC-1
                  URUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJM,1,IPENCL)*SCALEY
                  UTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJC,1,IPENCL)*SCALEY
                  VRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJM,2,IPENCL)*SCALEY
                  VTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJC,2,IPENCL)*SCALEY
                  WRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJM,3,IPENCL)*SCALEY
                  WTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJC,3,IPENCL)*SCALEY

                ENDDO

              ENDDO
              IPENCL = 0

            ENDIF

          ENDDO

        ELSE IF(KX.EQ.0)THEN         

C         SPECIAL CASE OF KZ=0

C         CHECK FOR PRE-STORED PENCILS
          IF(IPENCL.NE.0)THEN

C           DO THE FFT (NO SYMMETRY IMPOSED)
            NDOSYM = 0
            NPENCL = IPENCL

          CALL FFTSYM(NDOSYM,NPENCL,IYPROC,NYPROC,NYNODE,NYGLBL,NPROCY)

C           RESTORE TRANSFORMED DATA
            DO IPENCL = 1,NPENCL

              DO JC = JSTAL,JSTOL

                JJC = 2*JC
                JJM = JJC-1
                URUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                          = FTPART(JJM,1,IPENCL)*SCALEY
                UTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                          = FTPART(JJC,1,IPENCL)*SCALEY
                VRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                          = FTPART(JJM,2,IPENCL)*SCALEY
                VTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                          = FTPART(JJC,2,IPENCL)*SCALEY
                WRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                          = FTPART(JJM,3,IPENCL)*SCALEY
                WTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                          = FTPART(JJC,3,IPENCL)*SCALEY

              ENDDO

            ENDDO

          ENDIF

C         ASSEMBLE FFT DATA FOR KZ=0
          IPENCL = 0

          DO IC = ISTAL,ISTOL

C           FOURIER-SPACE GLOBAL INDEXING
            IX = IGOFM1 + IC
            IF(IX.GT.NODBLX)IX = IX-NXGLBL

C           PENCIL INDEXING
            IPENCL = IPENCL + 1
            IPEN(IPENCL) = IC
            KPEN(IPENCL) = KC

C           ASSEMBLE FFT DATA
            DO JC = JSTAL,JSTOL

              JJC = 2*JC
              JJM = JJC-1
              FTPART(JJM,1,IPENCL) = URUN(IC,JC,KC)        
              FTPART(JJC,1,IPENCL) = UTMP(IC,JC,KC)        
              FTPART(JJM,2,IPENCL) = VRUN(IC,JC,KC)        
              FTPART(JJC,2,IPENCL) = VTMP(IC,JC,KC)        
              FTPART(JJM,3,IPENCL) = WRUN(IC,JC,KC)        
              FTPART(JJC,3,IPENCL) = WTMP(IC,JC,KC)        

            ENDDO

C           CHECK FOR MAX NO OF PENCILS
            IF((IPENCL.EQ.NPENMX).OR.(IC.EQ.ISTOL))THEN

C             DO THE FFT (SYMMETRY IMPOSED)
              NDOSYM = 1
              NPENCL = IPENCL

          CALL FFTSYM(NDOSYM,NPENCL,IYPROC,NYPROC,NYNODE,NYGLBL,NPROCY)

C             RESTORE TRANSFORMED DATA
              DO IPENCL = 1,NPENCL

                DO JC = JSTAL,JSTOL

                  JJC = 2*JC
                  JJM = JJC-1
                  URUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJM,1,IPENCL)*SCALEY
                  UTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJC,1,IPENCL)*SCALEY
                  VRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJM,2,IPENCL)*SCALEY
                  VTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJC,2,IPENCL)*SCALEY
                  WRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJM,3,IPENCL)*SCALEY
                  WTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                            = FTPART(JJC,3,IPENCL)*SCALEY

                ENDDO

              ENDDO
              IPENCL = 0

            ENDIF

          ENDDO

        ELSE

C         LOWER-CENTRAL WAVENUMBERS - NO ACTION
          CONTINUE

        ENDIF

      ENDDO

C     CHECK FOR ANY REMAINING STORED PENCILS
      IF(IPENCL.NE.0)THEN

C       DO THE FFT (NO SYMMETRY IMPOSED)
        NDOSYM = 0
        NPENCL = IPENCL

        CALL FFTSYM(NDOSYM,NPENCL,IYPROC,NYPROC,NYNODE,NYGLBL,NPROCY)

C       RESTORE TRANSFORMED DATA
        DO IPENCL = 1,NPENCL

          DO JC = JSTAL,JSTOL

            JJC = 2*JC
            JJM = JJC-1
            URUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                      = FTPART(JJM,1,IPENCL)*SCALEY
            UTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                      = FTPART(JJC,1,IPENCL)*SCALEY
            VRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                      = FTPART(JJM,2,IPENCL)*SCALEY
            VTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                      = FTPART(JJC,2,IPENCL)*SCALEY
            WRUN(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                      = FTPART(JJM,3,IPENCL)*SCALEY
            WTMP(IPEN(IPENCL),JC,KPEN(IPENCL))
     +                      = FTPART(JJC,3,IPENCL)*SCALEY

          ENDDO

        ENDDO

      ENDIF

C     =========================================================================

C     Z-DIRECTION
C     -----------

C     PENCIL COUNTER
      IPENCL = 0

      DO JC = JSTAL,JSTOL

C       FOURIER-SPACE GLOBAL INDEXING
        JX = JGOFM1 + JC
        IF(JX.GT.NODBLY)JX = JX-NYGLBL

        DO IC = ISTAL,ISTOL

C         FOURIER-SPACE GLOBAL INDEXING
          IX = IGOFM1 + IC
          IF(IX.GT.NODBLX)IX = IX-NXGLBL

C         PENCIL INDEXING
          IPENCL = IPENCL + 1
          IPEN(IPENCL) = IC
          JPEN(IPENCL) = JC
     
C         ASSEMBLE FFT DATA
          DO KC = KSTAL,KSTOL

            KKC = 2*KC
            KKM = KKC-1
            FTPART(KKM,1,IPENCL) = URUN(IC,JC,KC)        
            FTPART(KKC,1,IPENCL) = UTMP(IC,JC,KC)        
            FTPART(KKM,2,IPENCL) = VRUN(IC,JC,KC)        
            FTPART(KKC,2,IPENCL) = VTMP(IC,JC,KC)        
            FTPART(KKM,3,IPENCL) = WRUN(IC,JC,KC)        
            FTPART(KKC,3,IPENCL) = WTMP(IC,JC,KC)        

          ENDDO

C         CHECK FOR MAX NO OF PENCILS
          IF(IPENCL.EQ.NPENMX)THEN

C           DO THE FFT (SYMMETRY IMPOSED)
            NDOSYM = 1
            NPENCL = IPENCL

          CALL FFTSYM(NDOSYM,NPENCL,IZPROC,NZPROC,NZNODE,NZGLBL,NPROCZ)

C           RESTORE TRANSFORMED DATA
            DO IPENCL = 1,NPENCL

              DO KC = KSTAL,KSTOL

                KKC = 2*KC
                KKM = KKC-1
                URUN(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                          = FTPART(KKM,1,IPENCL)*SCALEZ
                UTMP(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                          = FTPART(KKC,1,IPENCL)*SCALEZ
                VRUN(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                          = FTPART(KKM,2,IPENCL)*SCALEZ
                VTMP(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                          = FTPART(KKC,2,IPENCL)*SCALEZ
                WRUN(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                          = FTPART(KKM,3,IPENCL)*SCALEZ
                WTMP(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                          = FTPART(KKC,3,IPENCL)*SCALEZ

              ENDDO

            ENDDO
            IPENCL = 0

          ENDIF

        ENDDO

      ENDDO

C     CHECK FOR ANY REMAINING STORED PENCILS
      IF(IPENCL.NE.0)THEN

C       DO THE FFT (SYMMETRY IMPOSED)
        NDOSYM = 1
        NPENCL = IPENCL

        CALL FFTSYM(NDOSYM,NPENCL,IZPROC,NZPROC,NZNODE,NZGLBL,NPROCZ)

C       RESTORE TRANSFORMED DATA
        DO IPENCL = 1,NPENCL

          DO KC = KSTAL,KSTOL

            KKC = 2*KC
            KKM = KKC-1
            URUN(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                      = FTPART(KKM,1,IPENCL)*SCALEZ
            UTMP(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                      = FTPART(KKC,1,IPENCL)*SCALEZ
            VRUN(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                      = FTPART(KKM,2,IPENCL)*SCALEZ
            VTMP(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                      = FTPART(KKC,2,IPENCL)*SCALEZ
            WRUN(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                      = FTPART(KKM,3,IPENCL)*SCALEZ
            WTMP(IPEN(IPENCL),JPEN(IPENCL),KC)
     +                      = FTPART(KKC,3,IPENCL)*SCALEZ

          ENDDO

        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
