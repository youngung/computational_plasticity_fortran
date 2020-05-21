      program var_vel
      implicit none
      character*12 cdt
      integer i
      real eps, t, l, dt

c     initial condition
      l=5.
      t=0.
      eps=0.

      dt=0.01
      do i=1,iargc()
         call getarg(i,cdt)
         read(cdt,'(e20.13)')dt
      enddo

      write(*,*)'dt:',dt

      open(3,file='var_vel.txt')
      do while(t.lt.30)
         write(3,'(3f10.5)')  t,l,eps
         eps=eps+cos(t)/l * dt
         l = l + cos(t)*dt
         t = t + dt
      enddo
      close(3)

      end program
