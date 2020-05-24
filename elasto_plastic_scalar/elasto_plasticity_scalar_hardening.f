c     elasto-plasticity with Hollomon strain-hardening

c     scalar yield function f(s) = sqrt(s**2)
c     Homogeneous function of degree 1
      real function calc_yield_function(stress)
      implicit none
      real stress
      calc_yield_function = sqrt(stress**2.)
      return
      end function
c     Hardening based on Hollomon Equation
      real function hollomon(c,eps,n,k)
      real c,eps,n,k
      hollomon = c+k*(eps)**n
c      write(*,*)'hollomon:',hollomon
      return
      end function
c     --------------------------------------------
      program elasto_plasticity_scalar
      implicit none
      real calc_yield_function
      real dt, E, c, t, l, dl, eps, deps, deps_el, deps_pl, dsig,
     $     eps_el, eps_pl, stress, n, eeq_pl, K, H,hollomon
      real vel, f, tol
      integer kount,iplast,kount_outer
      parameter(tol=1e-8)

      open(3,file='elasto_plasticity_scalar_hardening.txt')
      dt = 0.02                 ! time increment
      E = 200000                ! elastic modulus
      c = 200.                  ! yield criterion
c     ---------
c     hollomon hardening model
      K=100.
      n=0.25
c     ---------
      stress = 0.               ! initial stress
      eps = 0.                  ! initial strain
      eps_el = 0.               ! initial el strain
      eps_pl = 0.               ! initial pl strain
      eeq_pl = 0.
      l = 1.                    ! initial length
      t = 0.                    ! initial time
      vel=0.01
      kount_outer=1

      write(*,'(3a9)')'sig+dsig','F-H','H'

      do while(t<10.0)
         if (eps>0.04) vel=-0.01
         dl = vel * dt
         deps = dl / l
c        initially assuming all strain is elastic
         deps_pl = 0.0
         deps_el = deps
c        guess on stress increment
         dsig = E * deps_el
         kount = 0
c        hardening, and yield surface size
         H = hollomon(c,eeq_pl+abs(deps_pl),n,K)
         f = calc_yield_function(stress + dsig) - H
         iplast=0
         do while (f.gt.tol.and.kount.lt.9) ! if exceeding plastic onset
            if (kount.eq.0) iplast=1
            if (eps.gt.0.04) vel=-0.01
            deps_pl = deps_pl -
     $           (f)/(-E*sign(1.,deps))
c     $            -n*K*(eeq_pl+abs(deps_pl)+0.0001)**(n-1))
c           New elastic strain increment
            deps_el = deps - deps_pl
c           New stress increment
            dsig = e * deps_el
c           new hardening, new F=f-H value.
            h = hollomon(c,eeq_pl+abs(deps_pl),n,K)
            f = calc_yield_function(stress+dsig)  - h
            kount = kount +1
         enddo

         write(3,'(2f9.4,f10.6,f10.2,2f10.6,2i2)')
     $        t,l,eps,stress,eps_el,eps_pl,iplast,kount ! write before update

         stress = stress + dsig
         eps = eps + deps
         eps_el = eps_el + deps_el
         eps_pl = eps_pl + deps_pl
         eeq_pl = eeq_pl + abs(deps_pl) ! In case of one-dimensional world.
         t = t + dt
         l = l + dl
         kount_outer=kount_outer+1
      enddo
      close(3)
      end program
