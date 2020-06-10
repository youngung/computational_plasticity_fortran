      subroutine gauss_nxn_swap(a,ndim,x)
      implicit none
      integer, intent(in)::ndim
      real*8, intent(inout):: a(ndim,ndim+1), x(ndim)
c     locals
      integer irow_orig(ndim),i,i0,j,j0
      logical iswap,iverb
      real*8 x_veri(ndim),dum,f,a_init(ndim,ndim+1)
      iverb=.false.
c      iverb=.true.
      iswap=.true.
c      iswap=.false.
      a_init(:ndim,:ndim+1)=a(:ndim,:ndim+1)
c     gauss elimination
      if (iverb) write(*,*)' ----------------------'

      do j0=1,ndim-1            ! This loop is for column (moves from left to right)
c        components along vertical line
         if (a(j0,j0).eq.0) then
            if (iverb) then
               write(*,*)'diagonal is zero'
               call prn_row(a(j0+1:ndim,j0),ndim-j0)
               write(*,*)'current row init:',j0
               write(*,*)'max row index:',i0+j0
               write(*,*)'before switching with ', j0 +1
               call prn(a,ndim)
            endif
            call get_max_ind(a(j0+1:ndim,j0),ndim-j0,i0)
            call swap_rows(a,ndim,j0,i0+j0)
            if (iverb) then
               write(*,*)'after switching with ', j0 +1
               call prn(a,ndim)
            endif
         endif

         if (j0.ne.ndim-1) then !switching not for the last
            if (iverb) call prn_row(a(j0+1:ndim,j0),ndim-j0)
            call get_max_ind(a(j0+1:ndim,j0),ndim-j0,i0)
            if (iverb) then
               write(*,*)'current row init:',j0+1
               write(*,*)'max row index:',i0
               write(*,*)'max row index:',i0+j0+1-1
               write(*,*)'before switching with ', j0 +1
               call prn(a,ndim)
            endif
            call swap_rows(a,ndim,j0+1,i0+j0+1-1)
            if (iverb) then
               write(*,*)'after'
               call prn(a,ndim)
            endif
         endif

         do j=j0+1,ndim         ! This loop is for row (moves from up to down)
            if (a(j,j0).ne.0d0) then
               f              = -a(j,j0)/a(j0,j0)
c               f              = -a(j0,j0)/a(j,j0) ! factor calculation
               if (iverb) then
                  write(*,'(a6,x,a1,i1,a1,i1,a1,a3,f6.2)')
     $                 '(j,j0)','(',j,',',j0,')','f:',f
                  write(*,'(a9,9f5.1)')'a1:',a(j,j0:ndim+1)
                  write(*,'(a9,9f5.1)')'a2:',a(j0,j0:ndim+1)
                  write(*,'(a9,9f5.1)')'f*a2:',f*a(j0,j0:ndim+1)
                  write(*,'(a9,9f5.1)')'a1+f*a2:',a(j,j0:ndim+1)
     $                 +f*a(j0,j0:ndim+1)
                  write(*,*)'before'
                  call prn(a,ndim)
               endif
               a(j,j0:ndim+1)=a(j,j0:ndim+1)+f*a(j0,j0:ndim+1)
               if (iverb) then
                  write(*,*)'after'
                  call prn(a,ndim)
               endif
            endif
         enddo
         if (iverb) then
            write(*,*)' ----------------------'
            call prn(a,ndim)
         endif
      enddo
c      a(:,:)=a0(:,:)
c     ----------------------
c     back substitute
      if (iverb) write(*,*)'back substitute'
      do i=ndim,1,-1
         dum = a(i,ndim+1)
         do j=ndim,i+1,-1
            dum = dum - a(i,j)*x(j)
         enddo
         x(i) = dum/a(i,i)
      enddo
      if (iverb)write(*,'(a9,50f12.2)')'x:',x(:ndim)
c     ----------------------------------
c     verification
      x_veri(:)=0d0
      do i=1,ndim
      do j=1,ndim
         x_veri(i)=x_veri(i)+a_init(i,j)*x(j)
      enddo
      enddo
c     ----------------------------------
c     x_veri is meant to be zero
      do i=1,ndim
         x_veri(i)=x_veri(i)-a_init(i,ndim+1)
      enddo
      dum=0d0
      do i=1,ndim
         dum=dum+x_veri(i)
      enddo
      if (iverb) write(*,'(a9,9e12.4)')'x_veri:',x_veri(:ndim)
      if (abs(dum).gt.1e-5) then
         write(*,*)'** something is wrong with the solution'
         stop -1
      endif

      return
      end subroutine
c     ----------------------------------------------------
      subroutine swap_rows(a,ndim,i0,j0)
c     a(i0,:) <-> a(j0,:)
      implicit none
      integer ndim,i0,j0,i,idum
      real*8 a(ndim,ndim+1),row_temp(ndim+1)
      logical iverb
      parameter (iverb=.false.)
c      parameter (iverb=.true.)
      if (iverb) write(*,'(10x,i2.2,x,a5,x,i2.2)')i0,'<--->',j0
      do i=1,ndim+1
         row_temp(i)=a(i0,i)
      enddo
      do i=1,ndim+1
         a(i0,i)=a(j0,i)
      enddo
      do i=1,ndim+1
         a(j0,i)=row_temp(i)
      enddo
      return
      end subroutine swap_rows
c     ----------------------------------------------------
      subroutine get_max_ind(a,ndim,ind)
      implicit none
      integer ind,ndim,i
      real*8 a(ndim),mx

      ind=1
      mx=abs(a(1))
      do i=2,ndim
         if (abs(a(i)).gt.mx) then
            mx=abs(a(i))
            ind=i
         endif
      enddo

      return
      end subroutine get_max_ind
c     ----------------------------------------------------
      subroutine prn(a,ndim)
      implicit none
      integer i,j,ndim
      real*8 a(ndim,ndim+1)
      do i=1,ndim
         write(*,'(40f7.2)',advance='no')(a(i,j),j=1,ndim+1)
         write(*,*)
      enddo
      return
      end subroutine prn
c     ----------------------------------------------------
      subroutine prnn(a,ndim)
      implicit none
      integer i,j,ndim
      real*8 a(ndim,ndim)
      do i=1,ndim
         write(*,'(40f7.2)',advance='no')(a(i,j),j=1,ndim)
         write(*,*)
      enddo
      return
      end subroutine prnn
c     ----------------------------------------------------
      subroutine prn_irow(a,ndim)
      implicit none
      integer i,ndim,a(ndim)
      write(*,'(40i2)')(a(i),i=1,ndim)
      return
      end subroutine prn_irow
c     ----------------------------------------------------
      subroutine prn_row(a,ndim)
      implicit none
      integer i,ndim
      real*8 a(ndim)
      write(*,'(40f7.2)')(a(i),i=1,ndim)
      return
      end subroutine prn_row
