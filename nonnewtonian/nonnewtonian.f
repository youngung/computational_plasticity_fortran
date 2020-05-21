c     non newtonian
      real function calc_f(sigma, eta0, alpha, edot)
      implicit none
      real sigma, et0, alpha, edot, eta0
      calc_f = sigma - eta0 * edot * exp(sigma / alpha)
      return
      end function

c     ---------------
      real function calc_df(sigma, eta0, alpha, edot)
      implicit none
      real sigma, eta0, alpha, edot
      calc_df = 1. - eta0*edot/alpha * exp(sigma/alpha)
      return
      end function

c     ---------------
      program main
      implicit none
      real s, tol, calc_f, calc_df, f, df, edot,alpha ,eta0
      integer kount
      parameter(tol=1e-5)

c     Input conditions
      eta0 = 13.
      alpha = 2.
      edot = 1e-3
c     -- File
      open(3,file='nr.txt',status='unknown')
c     --
      s = 1.                    ! initial guess
      f = tol * 2.              ! work-around
      kount = 0

c     -- Newton-Raphson loop
      do while(abs(f)>tol .and. kount < 10)
         f=calc_f(s,eta0,alpha,edot)
         df = calc_df(s,eta0,alpha,edot)
         write(3,'(i2.2,3e11.3)')kount, s, f, df
         s = s - f/df
         kount = kount + 1
      enddo

      close(3)
      end
