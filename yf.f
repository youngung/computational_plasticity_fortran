      subroutine yf(tsig,rsig,iopt)
      real tsig(6),rsig(6),phitilde,phi
      integer n,iopt
      !! Von Mises
      if (iopt.eq.1) then
         n=1.
         call vonmises(tsig,phitilde)
         write(*,*)'phitilde:',phitilde
         do i=1,6
            rsig(i)=tsig(i)*(1./phitilde)**(1./n)
         enddo
         write(*,*)(rsig(i),i=1,3)
         call vonmises(rsig,phi)
         write(*,*)'phi (Von Mises):',phi
      elseif(iopt.eq.2) then
         n=1.
         call hill48(0.7,0.4,0.3,1.5,1.5,1.5,tsig,phitilde)
         do i=1,6
            rsig(i)=tsig(i)*(1./phitilde)**(1./n)
         enddo
         call hill48(0.5,0.5,0.5,0.5,0.5,0.5,rsig,phi)
      endif

      return
      end subroutine yf
