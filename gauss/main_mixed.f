c     read input file and use gauss to solve the system of linear equation
      program main_mixed
      implicit none
      integer ndimx,ndim
      parameter(ndimx=50)
      real*8 a_init(ndimx,ndimx+1),a_new(ndimx,ndimx+1),x(ndimx),
     $     y(ndimx),b(ndimx),c(ndimx),c_new(ndimx),dum
      integer deltax(ndimx),deltay(ndimx),d0(ndimx),d1(ndimx),
     $     i,j,io,idum,kronecker(ndimx,ndimx)
      parameter(io=3)
      character*125 prosa
      real*8 dtime,time0
      logical iverb
c      parameter(iverb=.false.)
      parameter(iverb=.true.)
      call cpu_time(time0)
c     --- read inputs and preparation
      open(io,file='gauss_mix.inp',status='old')
      read(io,'(a)')prosa
      read(io,*)ndim
      kronecker(:,:)=0
      do i=1,ndim
         kronecker(i,i)=1
         d0(i)=0
         d1(i)=1
      enddo
      do i=1,ndim
         read(io,*)(a_init(i,j),j=1,ndim)
      enddo
      read(io,'(a)') prosa
      b(:)=0d0
      read(io,*) (deltax(i), b(i),i=1,ndim)
      read(io,'(a)')prosa
      read(io,*) (deltay(i), c(i),i=1,ndim)
      close(io)
c     --- verify bc
      idum=0
      do i=1,ndim
         idum=deltax(i)*deltay(i)
      enddo
      write(*,*)'idum:',idum
      if (idum.ne.0) then
         write(*,*)'**Error: ill-posed bc'
         stop
      endif
      do i=1,ndim
         idum=deltax(i)+deltay(i)
         if (idum.ne.1) then
            write(*,*)'**Error: ill-posed bc'
            stop
         endif
      enddo
c     ---
      if (iverb) then
         write(*,*)'Original system to solve'
         call prin_nxn_system(a_init(:ndim,:ndim),
     $        b(:ndim),c(:ndim),deltax(:ndim),
     $        deltay(:ndim),ndim)
      endif
c     --- checking if the problem is trivial
      idum=0
      do i=1,ndim
         idum=idum+deltax(i)
      enddo
      if (.true..and.idum.eq.ndim) then
         write(*,*)'** Warning: The problem is quite trivial'
         x(:ndim)=b(:ndim)
         do i=1,ndim
            y(i)=0
            do j=1,ndim
               y(i)=y(i)+a_init(i,j)*x(j)
            enddo
         enddo
         call prn_row(x(:ndim),ndim)
         call prn_row(y(:ndim),ndim)
         call cpu_time(dtime)
         dtime = dtime-time0
         write(*,'(a16,f6.1,x,a11)')'elapsed time:',
     $        dtime*1e6,'[micro sec]'
         write(*,*)'**fin.'
         stop
      endif

c     --- if the problem is not trivial, let's solve!
c     -- for mixed boundary condition, rearrange the system of
c     equations
      a_new(:,:)=0d0
      c_new(:)=0d0
      do i=1,ndim
         c_new(i)=(1-deltax(i))*c(i)
      do j=1,ndim
         a_new(i,j)=(1-deltax(j))*a_init(i,j) - kronecker(i,j)*deltax(j)
         c_new(i)=c_new(i) - a_init(i,j)*b(j)*deltax(j)
      enddo
      enddo
      if (iverb) then
         write(*,*)'New system to solve'
         call prin_nxn_system(a_new(:ndim,:ndim),
     $        b(:ndim),c_new(:ndim),d0(:ndim),d1(:ndim),ndim)
      endif
c     new augmented matrix.
      a_new(:ndim,ndim+1)=c_new(:ndim)
      if (iverb) then
         write(*,*)'A new matrix subjected to Gauss elimination'
         call prn(a_new(:ndim,:ndim+1),ndim)
      endif
c     -- solve the equation using gauss elimination method
      call gauss_nxn_swap(a_new(:ndim,:ndim+1),ndim,c_new(:ndim))
c     -- filling up initially unknown components
      do i=1,ndim
         if (deltax(i).eq.0) then
            x(i)=c_new(i)
            y(i)=c(i)
         else
            x(i)=b(i)
            y(i)=c_new(i)
         endif
      enddo
      write(*,*)'x and y'
      call prn_row(x(:ndim),ndim)
      call prn_row(y(:ndim),ndim)
c     -- verification
      c(:)=d0
      do i=1,ndim
         do j=1,ndim
            c(i)=c(i)+a_init(i,j)*x(j)
         enddo
         c(i)=c(i)-y(i)
      enddo
      dum=0d0
      do i=1,ndim
         dum=dum+abs(c(i))
      enddo
c
      write(*,'(a6,x,e9.2)')'tol.:',dum
      if (abs(dum).gt.1e-10) then
         write(*,*)'wrong'
         stop -1
      endif

      call cpu_time(dtime)
      dtime = dtime-time0
      write(*,'(a16,f6.1,x,a11)')'elapsed time:',dtime*1e6,'[micro sec]'
      write(*,*)'**fin.'

      stop
      end program
c     -----------------------------------------------------------
      subroutine prin_nxn_system(a,x,y,deltax,deltay,ndim)
      implicit none
      integer ndim
      real*8 a(ndim,ndim), x(ndim),y(ndim)
      integer deltax(ndim),deltay(ndim)
      integer i,j

      do i=1,ndim
         write(*,'(30f7.2)',advance='no') (a(i,j),j=1,ndim)
         if (deltax(i).eq.1) then
            write(*,'(f7.2)',advance='no') x(i)
         else
            write(*,'(a7)',advance='no') 'unkwn'
         endif
         if (deltay(i).eq.1) then
            write(*,'(f7.2)',advance='no') y(i)
         else
            write(*,'(a7)',advance='no') 'unkwn'
         endif
         write(*,*)
      enddo

      return
      end subroutine
