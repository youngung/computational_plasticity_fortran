      subroutine hill48(f,g,h,l,m,n,s,hill)
      implicit none
      real f,g,h,l,m,n
      real s(6),hill
      real s1,s2,s3,s4,s5,s6
      s1=s(1)
      s2=s(2)
      s3=s(3)
      s4=s(4)
      s5=s(5)
      s6=s(6)
      hill=f*(s2-s3)**2+g*(s3-s1)**2+h*(s1-s2)**2 + 2*l*s4**2+2*m*s5**2
     $     +2*n*s6**2
      hill = sqrt(hill)

      write(*,*)'hill:',hill
      write(*,'(6f10.2)')s1,s2,s3,s4,s5,s6
      return
      end subroutine hill48
