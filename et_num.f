c     Numerical solution to find E as function of time

      program et_num
      implicit none
      real dt, t, l, dl,vel, eps, deps
      dt=0.5

      t=0.
      l=1.
      vel=0.1
      eps=0.

      open(1,file='et_num.txt')
      do while(t.lt.3)
         write(1,'(2f10.5)')t,eps
         dl=vel*dt
         deps = dl/l
c        updates
         l=l+dl
         eps=eps+deps
         t=t+dt
      enddo

      close(1)

      end
