      program iso_elast
      implicit none
      real e(3,3,3,3),eps(3,3), sig(3,3)
      integer d(3,3), i, j

      labmda = 200.
      mu = 120.

c     kronecker delta

      do i=1,3
      do j=1,3
         if (i.eq.j) then
            d(i,j)=1
         else
            d(i,j)=0
         endif
      enddo
      enddo

      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         e(i,j,k,l)= e(k,j,k,l) + lambda * d(i,j) * d(k,l)
     $        + mu*(d(i,j)*d(j,l)+d(i,l)*d(j,k))
      enddo
      enddo
      enddo
      enddo


c     examples.
      eps(:,:)=0d0
      eps(1,2)=1.


      return
      end program

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
