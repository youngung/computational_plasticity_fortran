c     read input file and use gauss to solve the system of linear equation
      program main
      implicit none
      integer ndimx,ndim
      parameter(ndimx=50)
      real*8 a_init(ndimx,ndimx+1),a(ndimx,ndimx+1),dum,
     $     a0(ndimx,ndimx+1),b(ndimx+1),f,x(ndimx),x_veri(ndimx)
      integer i,j,k,l,j0,i0
      open(3,file='gauss.inp')
      read(3,*)ndim
      do i=1,ndim
         read(3,*)a_init(i,:ndim+1)
      enddo
      close(3)
      call prn(a_init(:ndim,:ndim+1),ndim)
      call gauss_nxn_swap(a_init(:ndim,:ndim+1),ndim,x(:ndim))
      stop
      end program
c     ----------------------------------------------------
