      function vector_sum(a,b)
      real a(3),b(3),vector_sum(3)
      integer i

      do i=1,3
         vector_sum(i) = a(i)+b(i)
      enddo

      return
      end function

      program main
      real a(3),b(3)
      real vector_sum(3)

      a(:)=0.
      b(:)=0.

      a(1)=3
      b(3)=2.

      write(*,*)vector_sum(a,b)
      end
