c     non newtonian
      real function calc_f(sigma, eta0, alpha, c, l)
      implicit none
      real sigma, et0, alpha, c,l, eta0
      calc_f = sigma - eta0  * exp(sigma / alpha) * c / l
      return
      end function

c     ---------------
      real function calc_df(sigma, eta0, alpha, c, l)
      implicit none
      real sigma, eta0, alpha, c,l
      calc_df = 1. - eta0/alpha * exp(sigma/alpha) *c/l
      return
      end function

c     ---
      program main
      implicit none
      real dt,alpha,eta0,vel,l,t,calc_f,calc_df,tol,f,df,dl,sigma
      integer kount, i
      character*12 cdt
      parameter(tol=1e-5)

c     input
      dt = 1.
      alpha=300.
      eta0=30.
      vel=0.0001

      do i=1,iargc()
         call getarg(i,cdt)
         read(cdt,'(e20.13)')dt
      enddo
c
      dl = vel* dt
c
      l = 10.                    ! initial length
      t = 0.
c
      sigma=0.                     ! the very initial guess on stress

c     file
      open(2,file='euler_nr.txt')

      do while(t<30.01)
c     solve the equation to obtain sigma
         f = tol *2.             ! work-around
         kount = 0
         do while(abs(f)>tol .and. kount < 10)
            f = calc_f(sigma,eta0,alpha,vel,l)
            df = calc_df(sigma,eta0,alpha,vel,l)
            sigma = sigma - f/df
            kount = kount + 1
         enddo
c         write(*,*) t, l, sigma, kount
         write(2,*) t, l, sigma, kount
         l=l+dl
         t=t+dt
      enddo

      close(2)
      end
