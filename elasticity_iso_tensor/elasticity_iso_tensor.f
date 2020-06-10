      program iso_elast
      implicit none
      real e(3,3,3,3),eps(3,3), sig(3,3),mu,lambda
      real d(3,3),dum
      integer i, j,l,k

      lambda = 200.
      mu = 120.

c     kronecker delta

      do i=1,3
      do j=1,3
         if (i.eq.j) then
            d(i,j)=1.0
         else
            d(i,j)=0.0
         endif
      enddo
      enddo

      e(:,:,:,:)=0d0

      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         dum =lambda * d(i,j) * d(k,l)
     $        + mu*(d(i,j)*d(j,l)+d(i,l)*d(j,k))
         e(i,j,k,l)= e(k,j,k,l) + dum
      enddo
      enddo
      enddo
      enddo

c     examples.
      eps(:,:)=0d0
      eps(1,2)=0.001
      eps(2,1)=0.001

      call tensor_dot(e,eps,sig)

c      write(*,'(3f10.1)') ((sig(i,j),j=1,3),i=1,3)
      do i=1,3
         write(*,'(3e13.2)') (sig(i,j),j=1,3)
      enddo

      return
      end program

c     ----------------------------------------------------
      subroutine tensor_dot(a,b,c)
c     c_ij = a_ijkl b_kl
      implicit none
      real a(3,3,3,3),b(3,3),c(3,3)
      integer i,j,k,l
      c(:,:)=0d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         c(i,j)=c(i,j) + a(i,j,k,l)*b(k,l)
      enddo
      enddo
      enddo
      enddo

      return
      end subroutine
