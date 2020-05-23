c     scalar yield function f(s) = sqrt(s**2)
      real function calc_yield_function(stress)
      implicit none
      real stress
      calc_yield_function = sqrt(stress**2.)
      return
      end function

c     scalar yield function df(s)/ds = 1.
      real function calc_df()
      implicit none
      calc_df = 1d0
      return
      end function

c     associated flow rule
      real function calc_deps_pl(stress,deps_pl)
      implicit none
      real stress,dlamb,deps_pl,dwork_pl,calc_yield_function,f

c     plastic work increment dw = s_ij de_ij
      dwork_pl = stress * deps_pl
c     We use yield function as g; (associated flow rule)
      f=calc_yield_function(stress) ! this also is sigma bar

c     Calculate dlambda
      dlamb = dwork_pl / f

c     calculate deps_pl_bar
      deps_pl = dlamb ! homogeneous yield function of degree 1

      return
      end function
c     ----------------------------
      program elasto_plasticity_scalar
      implicit none
      real calc_yield_function
      real calc_df
      real dt, E, c, t, l, dl, eps, deps, deps_el, deps_pl,dsig,
     $     stress
      real vel, f, tol
      integer kount
      parameter(tol=1e-6)

      open(3,file='elasto_plasticity_scalar.txt')

      dt = 0.01                 ! time increment
      E = 200000                ! elastic modulus
      c = 200.                  ! yield criterion


      stress = 0.               ! initial stress
      eps = 0.                  ! initial strain

      l = 1.                    ! initial length
      t = 0.
      vel = 0.01

      dl = vel * dt

      do while(t<0.5)
         deps = dl / l
c        initially assuming all strain is elastic
         deps_pl = 0.0
         deps_el = deps
c        guess on stress increment
         dsig = E * deps_el
         kount = 0
         f = calc_yield_function( stress+ dsig) - c
c         write(*,*)'f:',f
c         write(*,*)'stress+dsig:',stress+dsig
         do while (f.gt.tol.and.kount.lt.3) ! if exceeding plastic onset
c            write(*,*)'entering plasticity?'
            deps_el = deps - deps_pl
            dsig = e * deps_el
            f = calc_yield_function(stress+dsig) - c
c           estimate new plastic increment
            deps_pl = deps_pl - f/(-E)
c            write(*,*)'deps_pl',deps_pl
            kount = kount +1
         enddo
         write(3,*)t,l,kount,eps,stress
         stress = stress + dsig
         eps = eps + deps
         t = t + dt
         l = l + dl
      enddo


      close(3)

      end program
