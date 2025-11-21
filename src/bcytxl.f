      SUBROUTINE BCYTXL
 
C     *************************************************************************
C
C     BCYTXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-DEC-2003:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
C     AND THEIR TIME DERIVATIVES
C
C     X-DIRECTION LEFT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER JC,KC
      INTEGER ISPEC
      double precision toty


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRYXL,DYDTXL
      DO ISPEC = 1,NSPEC

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C           SET MASS FRACTIONS TO CONSTANT (INITIAL) VALUES
            STRYXL(JC,KC,ISPEC) = YRIN(ISPEC)

C           SET MASS FRACTION TIME DERIVATIVES TO ZERO
            DYDTXL(JC,KC,ISPEC) = ZERO

          ENDDO
        ENDDO

      ENDDO

C     VM: SYNTHETIC SCALAR INFLOW
C     VM: NXLPRM(2)=1 IMPLIES THAT THE SCALAR SYTHETIC DIGITAL FILTERING
C     IS ON
      IF ((NXLPRM(2)==1).AND.(NXLPRM(1)==4).AND.(NGBCXL==12))THEN
        DO ISPEC=1,NSPEC
          DO KC=KSTAL,KSTOL
            DO JC=JSTAL,JSTOL
              STRYXL(JC,KC,ISPEC)=YRIN(ISPEC)+YINF2(JC,KC,ISPEC)
              if(stryxl(jc,kc,ispec).gt.1.0d0) then
                yinf2(jc,kc,ispec)=1.0d0-yrin(ispec)
                stryxl(jc,kc,ispec)=1.0d0
              endif
              if(stryxl(jc,kc,ispec).lt.0.0d0) then
                yinf2(jc,kc,ispec)=yrin(ispec)-0.0d0
                stryxl(jc,kc,ispec)=0.0d0
              endif
              DYDTXL(JC,KC,ISPEC)=(YINF2(JC,KC,ISPEC)-
     +                             YINF1(JC,KC,ISPEC))/TSTEP
            ENDDO
          ENDDO
        ENDDO
        do kc=kstal,kstol
          do jc=jstal,jstol
            toty=0.0d0
            do ispec=1,nspec-1
              toty = toty+stryxl(jc,kc,ispec)
            enddo
            stryxl(jc,kc,nspec)=1.0d0-toty
            stryxl(jc,kc,nspec)=max(0.0,stryxl(jc,kc,nspec))
            stryxl(jc,kc,nspec)=min(1.0,stryxl(jc,kc,nspec))
            dydtxl(jc,kc,nspec)=(yinf2(jc,kc,nspec)
     +                          -yinf1(jc,kc,nspec))/tstep
          enddo
        enddo
C        if(itime.eq.1) then
C          if (ixproc.eq.0) then
C            if (iyproc.eq.0)then
C              open(unit=1111,file='debug/rand0.dat',status='new',
C     +             form='formatted')
C              do kc=kstal,kstol
C                do jc=jstal,jstol
C                write(1111,*)(stryxl(jc,kc,ispc),ispc=1,nspec)
C                enddo
C              enddo
C              close(1111)
C            end if
C            if (iyproc.eq.1)then
C              open(unit=1112,file='debug/rand1.dat',status='new',
C     +             form='formatted')
C              do kc=kstal,kstol
C                do jc=jstal,jstol
C                write(1112,*)(stryxl(jc,kc,ispc),ispc=1,nspec)
C                enddo
C              enddo
C              close(1112)
C            end if
C            if (iyproc.eq.2)then
C              open(unit=1113,file='debug/rand2.dat',status='new',
C     +             form='formatted')
C              do kc=kstal,kstol
C                do jc=jstal,jstol
C                write(1113,*)(stryxl(jc,kc,ispc),ispc=1,nspec)
C                enddo
C              enddo
C              close(1113)
C            end if
C            if (iyproc.eq.3)then
C              open(unit=1114,file='debug/rand3.dat',status='new',
C     +             form='formatted')
C              do kc=kstal,kstol
C                do jc=jstal,jstol
C                write(1114,*)(stryxl(jc,kc,ispc),ispc=1,nspec)
C                enddo
C              enddo
C              close(1114)
C            end if
C          end if
C        end if
      ENDIF

C     =========================================================================


      RETURN
      END
