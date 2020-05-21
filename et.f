c     E as a function of time
      program et
      implicit none
      real l0, vel, dt, t, eps
c     inputs
      l0=1.                      ! [mm]
      vel=0.1                   ! [mm/s]
      dt=0.1
      t=0.
c     --
      open(1,file='et.txt')
      do while(t.lt.3)
         write(1,'(2f10.5)')t,eps(l0,vel,t)
         t=t+dt
      enddo
      close(1)
      end program
c     ---------------------------
c     analytical function eps as function of t
      real function eps(l0,vel,t)
      implicit none
      real vel,t,l0
      eps=log(l0+vel*t) - log (l0)
      return
      end function
