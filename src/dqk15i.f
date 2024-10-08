      subroutine dqk15i(boun,inf,a,b,result,abserr,resabs,resasc,ienv)
c***begin prologue  dqk15i
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a2,h2a4a2
c***keywords  15-point transformed gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the original (infinite integration range is mapped
c            onto the interval (0,1) and (a,b) is a part of (0,1).
c            it is the purpose to compute
c            i = integral of transformed integrand over (a,b),
c            j = integral of abs(transformed integrand) over (a,b).
c***description
c
c           integration rule
c           standard fortran subroutine
c           double precision version
c
c           parameters
c            on entry
cC             f      - double precision
cC                      function subprogram defining the integrand
cC                      function f(x). the actual name for f needs to be
cC                      declared e x t e r n a l in the calling program.
c
C              F      - subroutine defining the integrand function.
C                       Should have 4 arguments:
C                         - fout - double precision
C                                  vector of dimension n
C                                  output values (function results).
C                         - fin  - double precision
C                                  vector of dimension n
C                                  input values (function arguments)
C                         - n    - integer
C                                  number of input values
C                         - ienv - integer
C                                  array, length 2.  Passed to/from C
C                                  code, not used in Fortran statements.
C                       F should be declared EXTERNAL in the calling
C                       program.
c
c              boun   - double precision
c                       finite bound of original integration
c                       range (set to zero if inf = +2)
c
c              inf    - integer
c                       if inf = -1, the original interval is
c                                   (-infinity,bound),
c                       if inf = +1, the original interval is
c                                   (bound,+infinity),
c                       if inf = +2, the original interval is
c                                   (-infinity,+infinity) and
c                       the integral is computed as the sum of two
c                       integrals, one over (-infinity,0) and one over
c                       (0,+infinity).
c
c              a      - double precision
c                       lower limit for integration over subrange
c                       of (0,1)
c
c              b      - double precision
c                       upper limit for integration over subrange
c                       of (0,1)
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule(resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule(resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral j
c
c              resasc - double precision
c                       approximation to the integral of
c                       abs((transformed integrand)-i/(b-a)) over (a,b)
C
C              ienv   - integer
C                       array, length 2. Passed to/from C code,
C                       not used in Fortran statements.
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15i
c
      double precision a,absc,absc1,absc2,abserr,b,boun,centr,dinf,
     *  d1mach,epmach,fc,fsum,fval1,fval2,fv1,fv2,hlgth,
     *  resabs,resasc,resg,resk,reskh,result,tabsc1,tabsc2,uflow,wg,wgk,
     *  xgk
      Double Precision fout1(15),fout2(15),fin(15)
      integer inf,j
      Integer ienv(2)
      external f
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(8)
c
c           the abscissae and weights are supplied for the interval
c           (-1,1).  because of symmetry only the positive abscissae and
c           their corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule, corresponding
c                    to the abscissae xgk(2), xgk(4), ...
c                    wg(1), wg(3), ... are set to zero.
c
      data wg(1) / 0.0d0 /
      data wg(2) / 0.1294849661 6886969327 0611432679 082d0 /
      data wg(3) / 0.0d0 /
      data wg(4) / 0.2797053914 8927666790 1467771423 780d0 /
      data wg(5) / 0.0d0 /
      data wg(6) / 0.3818300505 0511894495 0369775488 975d0 /
      data wg(7) / 0.0d0 /
      data wg(8) / 0.4179591836 7346938775 5102040816 327d0 /
c
      data xgk(1) / 0.9914553711 2081263920 6854697526 329d0 /
      data xgk(2) / 0.9491079123 4275852452 6189684047 851d0 /
      data xgk(3) / 0.8648644233 5976907278 9712788640 926d0 /
      data xgk(4) / 0.7415311855 9939443986 3864773280 788d0 /
      data xgk(5) / 0.5860872354 6769113029 4144838258 730d0 /
      data xgk(6) / 0.4058451513 7739716690 6606412076 961d0 /
      data xgk(7) / 0.2077849550 0789846760 0689403773 245d0 /
      data xgk(8) / 0.0000000000 0000000000 0000000000 000d0 /
c
      data wgk(1) / 0.0229353220 1052922496 3732008058 970d0 /
      data wgk(2) / 0.0630920926 2997855329 0700663189 204d0 /
      data wgk(3) / 0.1047900103 2225018383 9876322541 518d0 /
      data wgk(4) / 0.1406532597 1552591874 5189590510 238d0 /
      data wgk(5) / 0.1690047266 3926790282 6583426598 550d0 /
      data wgk(6) / 0.1903505780 6478540991 3256402421 014d0 /
      data wgk(7) / 0.2044329400 7529889241 4161999234 649d0 /
      data wgk(8) / 0.2094821410 8472782801 2999174891 714d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           tabsc* - transformed abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of the transformed
c                    integrand over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15i
      epmach = d1mach(4)
      uflow = d1mach(1)
      dinf = min0(1,inf)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)

      tabsc1 = boun+dinf*(0.1d+01-centr)/centr
      fin(8)=tabsc1
      Do 2 j=1,7
        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        tabsc1 = boun+dinf*(0.1d+01-absc1)/absc1
        tabsc2 = boun+dinf*(0.1d+01-absc2)/absc2
        fin(j)=tabsc1
        fin(16-j)=tabsc2
    2 Continue
      Call f(fout1,fin,15,ienv)
      If (inf.eq.2) Then
        Do 4 j=1,15
          fin(j)=-fin(j)
    4   Continue
        Call f(fout2,fin,15,ienv)
      End If
      fval1=fout1(8)
      If (inf.eq.2) fval1=fval1+fout2(8)
      fc = (fval1/centr)/centr
      resg = wg(8)*fc
      resk = wgk(8)*fc
      resabs = abs(resk)
      Do 10 j=1,7
        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1=fout1(j)
        fval2=fout1(16-j)
        If (inf.eq.2) Then
          fval1=fval1+fout2(j)
          fval2=fval2+fout2(16-j)
        End If
        fval1 = (fval1/absc1)/absc1
        fval2 = (fval2/absc2)/absc2
        fv1(j) = fval1
        fv2(j) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(j)*fsum
        resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
   10 Continue

C        tabsc1 = boun+dinf*(0.1d+01-centr)/centr
C        fval1 = f(tabsc1)
C        if(inf.eq.2) fval1 = fval1+f(-tabsc1)
C        fc = (fval1/centr)/centr
C  c
C  c           compute the 15-point kronrod approximation to
C  c           the integral, and estimate the error.
C  c
C        resg = wg(8)*fc
C        resk = wgk(8)*fc
C        resabs = abs(resk)
C        do 10 j=1,7
C          absc = hlgth*xgk(j)
C          absc1 = centr-absc
C          absc2 = centr+absc
C          tabsc1 = boun+dinf*(0.1d+01-absc1)/absc1
C          tabsc2 = boun+dinf*(0.1d+01-absc2)/absc2
C          fval1 = f(tabsc1)
C          fval2 = f(tabsc2)
C          if(inf.eq.2) fval1 = fval1+f(-tabsc1)
C          if(inf.eq.2) fval2 = fval2+f(-tabsc2)
C          fval1 = (fval1/absc1)/absc1
C          fval2 = (fval2/absc2)/absc2
C          fv1(j) = fval1
C          fv2(j) = fval2
C          fsum = fval1+fval2
C          resg = resg+wg(j)*fsum
C          resk = resk+wgk(j)*fsum
C          resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
C     10 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*abs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resasc = resasc*hlgth
      resabs = resabs*hlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d0) abserr = resasc*
     *  min(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach))
     *  abserr =  max((epmach*0.5d+02)*resabs,abserr)
      return
      end

