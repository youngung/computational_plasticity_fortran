      program yf_plot
      implicit none
      real tsig(6),rsig(6),dangle,pi,angle,ang,rsig2(6)
      integer nangle
      parameter(nangle=360)
      integer i,j

      pi=3.141592
      dangle=2.*pi/(nangle-1)
      open(1,file='ys.txt',status='unknown')
      do i=1,nangle
         ang=-pi + dangle*(i-1.)
         tsig(1)=cos(ang)
         tsig(2)=sin(ang)
         tsig(3:)=0.

         call yf(tsig,rsig2,1)   ! von Mises
         call yf(tsig,rsig,2)   ! Hill 48
         write(1,'(9f10.4)')(tsig(j),j=1,3),(rsig(j),j=1,3),
     $        (rsig2(j),j=1,3)
      enddo
      close(1)
      end program
