      DOUBLE PRECISION FUNCTION D1MACH(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (I.EQ.1) THEN
        D1MACH=TINY(0D0)
      ELSE IF (I.EQ.2) THEN
        D1MACH=HUGE(0D0)
      ELSE IF (I.EQ.4) THEN
        D1MACH=EPSILON(0D0)
      ELSE
C         should never happen in this program
        D1MACH=0D0
      ENDIF
      RETURN
      END
