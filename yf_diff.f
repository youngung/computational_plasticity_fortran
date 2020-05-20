c     program to numerically obtain derivatives of yield surface

      program ys_diff
      implicit none
      real s(6),sd(6),rsig(6),s_plus_delta(6),phi,phi1,phi2,delta,
     $     deltas(10)
      integer n,i


      s(:)=0.
      s(1)=1.
      s(2)=1.

      n = 1 ! degree of homogeneous function
      call vonmises(s,phi) ! phi tilde
      do i=1,6
         s(i)=s(i)*(1./phi)**(1./n)
      enddo
c     s is now on the VM yield surface with phi=1.
      call vonmises(s,phi)      ! phi should be one (if I have not screwed up so far)

c     Its derivative?
      delta=0.0001
      do i=1,6 ! loop over each component
         s_plus_delta(:)=s(:)   ! initialize
         s_plus_delta(i)=s_plus_delta(i)+delta ! only for specific component
         call vonmises(s_plus_delta,phi1)
         sd(i) = (phi1-phi)/delta       ! this is the finite difference
      enddo

      write(*,'(6f10.3)') (sd(i),i=1,6)

      end program
