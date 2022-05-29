      SUBROUTINE SAMLM(X,N,XMOM,NMOM,ISORT,IRATIO)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. HOSKING <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 2.0    April 2013                                          *
C*                                                                     *
C*  Version 2.9    June 2020                                           *
C*  * Removed the "sanity check" that restricted higher-order          *
C*    L-moments to be no greater than l_2 in absoulte value.           *
C*                                                                     *
C***********************************************************************
C
C  Sample L-moments of a data array
C
C  Arguments of routine:
C  X      * input* Array of length N. Contains the data.
C  N      * input* Number of data values
C  XMOM   *output* Array of length NMOM. On exit, contains the sample
C                  L-moments, in the order l-1, l-2, l-3, l-4, ... .
C  NMOM   * input* Number of L-moments to be found.
C  ISORT  * input* If >0, data will be sorted.
C  IRATIO * input* If >0, L-moment ratios rather than L-moments will be
C                  returned.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        Note: Fortran 77 would require dimension of COEF to be known
C        at compile time.
      DOUBLE PRECISION X(N),XMOM(NMOM),COEF(2,NMOM)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/
C
C                          QSORT3 is R's internal sort routine
      IF (ISORT.GT.0) CALL QSORT3(X,1,N)
C
      DN=N
      DO 10 J=1,NMOM
   10 XMOM(J)=ZERO
      IF(NMOM.LE.2)GOTO 100
C
C         Unbiased estimates of L-moments -- the 'DO 30' loop
C         recursively calculates discrete Legendre polynomials, via
C         eq.(9) of Neuman and Schonbach (1974, Int.J.Num.Meth.Eng.)
C
      DO 20 J=3,NMOM
      TEMP=ONE/DFLOAT((J-1)*(N-J+1))
      COEF(1,J)=DFLOAT(J+J-3)*TEMP
      COEF(2,J)=DFLOAT((J-2)*(N+J-2))*TEMP
   20 CONTINUE
      TEMP=-DN-ONE
      CONST=ONE/(DN-ONE)
      NHALF=N/2
      DO 40 I=1,NHALF
      TEMP=TEMP+TWO
      XI=X(I)
      XII=X(N+1-I)
      TERMP=XI+XII
      TERMN=XI-XII
      XMOM(1)=XMOM(1)+TERMP
      S1=ONE
      S=TEMP*CONST
      XMOM(2)=XMOM(2)+S*TERMN
      DO 30 J=3,NMOM,2
      S2=S1
      S1=S
      S=COEF(1,J)*TEMP*S1-COEF(2,J)*S2
      XMOM(J)=XMOM(J)+S*TERMP
      IF(J.EQ.NMOM)GOTO 30
      JJ=J+1
      S2=S1
      S1=S
      S=COEF(1,JJ)*TEMP*S1-COEF(2,JJ)*S2
      XMOM(JJ)=XMOM(JJ)+S*TERMN
   30 CONTINUE
   40 CONTINUE
      IF(N.EQ.NHALF+NHALF)GOTO 60
      TERM=X(NHALF+1)
      S=ONE
      XMOM(1)=XMOM(1)+TERM
      DO 50 J=3,NMOM,2
      S=-COEF(2,J)*S
      XMOM(J)=XMOM(J)+S*TERM
   50 CONTINUE
   60 CONTINUE
C
C         Sanity check: higher-order L-moments cannot exceed l_2
C
C     DO 70 J=3,NMOM
C       IF (ABS(XMOM(J)).GT.XMOM(2)) XMOM(J)=SIGN(XMOM(2),XMOM(J))
C  70 CONTINUE
C
C         L-moments or L-moment ratios
C
      RN=ONE/DN
      IF (IRATIO.GT.0) THEN
        IF (XMOM(2).EQ.ZERO) THEN
          AMULT=ZERO
        ELSE
          AMULT=ONE/XMOM(2)
        END IF
      ELSE
        AMULT=RN
      END IF
      DO 80 J=3,NMOM
   80 XMOM(J)=XMOM(J)*AMULT
      XMOM(1)=XMOM(1)*RN
      XMOM(2)=XMOM(2)*RN
      RETURN
C
C         At most two L-moments
C
  100 CONTINUE
      SUM1=ZERO
      SUM2=ZERO
      TEMP=-DN+ONE
      DO 110 I=1,N
      SUM1=SUM1+X(I)
      SUM2=SUM2+X(I)*TEMP
      TEMP=TEMP+TWO
  110 CONTINUE
      XMOM(1)=SUM1/DN
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=SUM2/(DN*(DN-ONE))
C
C         Sanity check for positive data: l_2 cannot exceed l_1
C
      IF (ISORT.GT.0.AND.X(1).GE.ZERO.AND.XMOM(2).GT.XMOM(1))
     *  XMOM(2)=XMOM(1)
C
      RETURN
C
      END
