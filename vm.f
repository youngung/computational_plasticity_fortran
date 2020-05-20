      subroutine vonmises(s,vm)
      implicit none
      real s(6),vm
      real s1,s2,s3,s4,s5,s6

      s1=s(1)
      s2=s(2)
      s3=s(3)
      s4=s(4)
      s5=s(5)
      s6=s(6)

      vm = (s1-s2)**2+(s2-s3)**2+(s3-s1)**2 + 6*(s4**2+s5**2+s6**2)
      vm = sqrt(vm /2.)
      write(*,*)'vm:',vm
      return
      end subroutine
