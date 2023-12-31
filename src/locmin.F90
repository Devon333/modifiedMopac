      subroutine locmin(m, xparam, nvar, p, ssq, alf, efs, ncount) 
      USE vast_kind_param, ONLY:  double 
      use molkst_C, only : numcal, keywrd, escf
      use chanel_C, only : iw
      use exchng_I
      use compfg_I
      implicit none
      integer , intent(in) :: m, nvar
      integer  :: ncount
      real(double) , intent(inout) :: ssq 
      real(double) , intent(inout) :: alf 
      real(double)  :: xparam(nvar) 
      real(double) , intent(inout) :: p(nvar) 
      real(double)  :: efs(nvar) 
      integer ::  left, right, center, mxcnt2, iprint, icalcn, i, mxct, ictr 
      real(double), dimension(nvar) ::  xstor, gstor 
      real(double), dimension(3) :: phi, vt 
      real(double) :: xmaxm, scale, eps, tee, ymaxst, xcrit, xminm, fin, t, &
        tscale, sqstor, energy = 0.d0, estor, alfs, tlast, flast, f, alpha, beta, &
        gamma, abg, s, amdis, sscale, sum, sum2 
      logical :: debug, lower      
      real(double), external :: ddot
      save  xmaxm, scale, eps, debug, tee, ymaxst, xcrit, mxcnt2, iprint, icalcn 
      data icalcn/ 0/  
!***********************************************************************
!
!    LOCMIN IS CALLED BY NLLSQ ONLY. IT IS A LINE-SEARCH PROCEDURE FOR
!    LOCATING A MINIMUM IN THE FUNCTION SPACE OF COMPFG.  SEE NLLSQ
!    FOR MORE DETAILS
!
!***********************************************************************
      if (icalcn /= numcal) then 
        icalcn = numcal 
        xmaxm = 1.D9 
        scale = 1.D0 
!
! THE ABOVE LINE IS TO TRY TO PREVENT OVERFLOW IN NLLSQ
!
        eps = 1.D-5 
        debug = index(keywrd,'LINMIN') /= 0 
        tee = 1.D-2 
        ymaxst = 0.005D0 
        xcrit = 0.0002D0 
        mxcnt2 = 30 
        iprint = 0 
        if (debug) iprint = -1 
      endif 
      xmaxm = 1.D-11 
      do i = 1, nvar 
        xmaxm = max(xmaxm,abs(p(i))) 
      end do 
      xminm = xmaxm*scale 
      xmaxm = ymaxst/xmaxm/scale 
      fin = ssq 
      lower = .FALSE. 
      t = alf 
      phi(1) = ssq 
      vt(1) = 0.0D0 
      vt(2) = t/4.0D0 
      vt(2) = min(xmaxm,vt(2)) 
      t = vt(2) 
      tscale = t*scale 
      xparam(:nvar) = xparam(:nvar) + p(:nvar)*tscale 
      call compfg (xparam, .TRUE., escf, .TRUE., efs, .TRUE.) 
      phi(2) = ddot(nvar,efs,1,efs,1) 
      call exchng (phi(2), sqstor, energy, estor, xparam, xstor, t, alfs, nvar) 
      gstor(:m) = efs(:m) 
      if (phi(1) <= phi(2)) then 
        vt(3) = -vt(2) 
        left = 3 
        center = 1 
        right = 2 
      else 
        vt(3) = 2.0D0*vt(2) 
        left = 1 
        center = 2 
        right = 3 
      endif 
      tlast = vt(3) 
      t = tlast - t 
      tscale = t*scale 
      xparam(:nvar) = xparam(:nvar) + p(:nvar)*tscale 
      flast = phi(2) 
      call compfg (xparam, .TRUE., escf, .TRUE., efs, .TRUE.) 
      f = ddot(nvar,efs,1,efs,1) 
      if (f < sqstor) call exchng (f, sqstor, energy, estor, xparam, xstor, t, alfs&
        , nvar) 
      gstor(:m) = efs(:m) 
      if (f < fin) lower = .TRUE. 
      ncount = ncount + 2 
      phi(3) = f 
      if (iprint < 0) then 
        write (iw, 310) vt(1), sqrt(phi(1)), vt(2), sqrt(phi(2)), vt(3), sqrt(&
          phi(3)) 
      endif 
      mxct = mxcnt2 
      do ictr = 3, mxct 
        xmaxm = xmaxm*3.D0 
        alpha = vt(2) - vt(3) 
        beta = vt(3) - vt(1) 
        gamma = vt(1) - vt(2) 
        if (alpha == 0.D0) alpha = 1.D-20 
        if (beta == 0.D0) beta = 1.D-20 
        if (gamma == 0.D0) gamma = 1.D-20 
        abg = -(phi(1)*alpha+phi(2)*beta+phi(3)*gamma)/alpha 
        abg = abg/beta 
        abg = abg/gamma 
        alpha = abg 
        beta = (phi(1)-phi(2))/gamma - alpha*(vt(1)+vt(2)) 
        if (alpha <= 0.D0) then 
          if (phi(right) <= phi(left)) then 
            t = 3.0D0*vt(right) - 2.0D0*vt(center) 
          else 
            t = 3.0D0*vt(left) - 2.0D0*vt(center) 
          endif 
          s = t - tlast 
          t = s + tlast 
        else 
          t = -beta/(2.0D0*alpha) 
          s = t - tlast 
          if (s <= 0.D0) then 
            if (s == 0.D0) exit  
            amdis = vt(left) - tlast - xmaxm 
          else 
            amdis = vt(right) - tlast + xmaxm 
          endif 
          if (abs(s) > abs(amdis)) s = amdis 
          t = s + tlast 
        endif 
        if (ictr>3 .and. abs(s*xminm)<xcrit) then 
          if (debug) write (iw, '('' EXIT DUE TO SMALL PROJECTED STEP'')') 
          exit  
        endif 
        t = s + tlast 
        sscale = s*scale 
        xparam(:nvar) = xparam(:nvar) + p(:nvar)*sscale 
        flast = f 
        call compfg (xparam, .TRUE., escf, .TRUE., efs, .TRUE.) 
        f = ddot(nvar,efs,1,efs,1) 
        if (f < sqstor) call exchng (f, sqstor, energy, estor, xparam, xstor, t, &
          alfs, nvar) 
        gstor(:m) = efs(:m) 
        if (f < fin) lower = .TRUE. 
        ncount = ncount + 1 
        if (iprint < 0) then 
          write (iw, 320) vt(left), sqrt(phi(left)), vt(center), sqrt(phi(&
            center)), vt(right), sqrt(phi(right)), t, sqrt(f) 
        endif 
!
!    TEST FOR EXCITED STATES AND POTHOLES
!
        if (abs(vt(center)) > 1.D-10) go to 200 
        if (abs(t)/(abs(vt(left))+1.D-15) > 0.3333D0) go to 200 
        if (2.5D0*f - phi(right) - phi(left) < 0.5D0*phi(center)) go to 200 
!
!   WE ARE STUCK ON A FALSE MINIMUM
!
        exit  
  200   continue 
!
! NOW FOR THE MAIN STOPPING TESTS.  LOCMIN WILL STOP IF:-
!     THE ERROR FUNCTION HAS BEEN REDUCED, AND
!     THE RATE OF DROP OF THE ERROR FUNCTION IS LESS THAN 0.5% PER STEP
!     AND
!     (A) THE RATIO OF THE PROPOSED STEP TO THE TOTAL STEP IS LESS THAN
!         EPS,   OR
!     (B) THE LAST DROP IN ERROR FUNCTION WAS LESS THAN 5%OFTHETOTALDROP
!         DURING THIS CALL TO LOCMIN.
!
        if (debug) write (iw, '('' F/FLAST'',F13.6)') f/flast 
        if (lower .and. f/flast>0.995D0) then 
          if (abs(t - tlast) <= eps*abs(t + tlast) + tee) then 
            if (debug) write (iw, '('' EXIT AS STEP IS ABSOLUTELY SMALL '')') 
            exit  
          endif 
          sum = min(min(abs(f - phi(1)),abs(f-phi(2))),abs(f-phi(3))) 
          sum2 = (fin - sqstor)*0.05D0 
          if (sum < sum2) then 
            if (debug) write (iw, '('' EXIT DUE TO HAVING REACHED BOTTOM'')') 
            exit  
          endif 
        endif 
        tlast = t 
        if (.not.(t>vt(right) .or. t>vt(center) .and. f<phi(center) .or. t>vt(&
          left) .and. t<vt(center) .and. f>phi(center))) then 
          vt(right) = t 
          phi(right) = f 
        else 
          vt(left) = t 
          phi(left) = f 
        endif 
        if (vt(center) >= vt(right)) then 
          i = center 
          center = right 
          right = i 
        endif 
        if (vt(left) >= vt(center)) then 
          i = left 
          left = center 
          center = i 
        endif 
        if (vt(center) < vt(right)) cycle  
        i = center 
        center = right 
        right = i 
      end do 
      call exchng (sqstor, f, estor, energy, xstor, xparam, alfs, t, nvar) 
      efs(:m) = gstor(:m) 
      ssq = f 
      alf = t 
      if (t < 0.D0) then 
        t = -t 
        p(:nvar) = -p(:nvar) 
      endif 
      alf = t 
      return  
  310 format(' ---LOCMIN'/,5x,'LEFT   ...',2f19.6,/,5x,'CENTER ...',2f19.6,/,5x&
        ,'RIGHT  ...',2f19.6,/,' ') 
  320 format(5x,'LEFT   ...',2f19.6,/,5x,'CENTER ...',2f19.6,/,5x,'RIGHT  ...',&
        2f19.6,/,5x,'NEW    ...',2f19.6,/,' ') 
      end subroutine locmin 
