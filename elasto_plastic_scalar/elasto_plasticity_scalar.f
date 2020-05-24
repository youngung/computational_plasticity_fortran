c     elasto-plasticity without strain-hardening

c     scalar yield function f(s) = sqrt(s**2)
c     Homogeneous function of degree 1
      real function calc_yield_function(stress)
      implicit none
      real stress
      calc_yield_function = sqrt(stress**2.)
      return
      end function
c     --------------------------------------------
      program elasto_plasticity_scalar
      implicit none
      real calc_yield_function
      real dt, E, c, t, l, dl, eps, deps, deps_el, deps_pl, dsig,
     $     eps_el, eps_pl, stress
      real vel, f, tol
      integer kount,iplast
      parameter(tol=1e-6)

      open(3,file='elasto_plasticity_scalar.txt')
      dt = 0.03                  ! time increment
      E = 200000                ! elastic modulus
      c = 200.                  ! yield criterion
      stress = 0.               ! initial stress
      eps = 0.                  ! initial strain
      eps_el = 0.               ! initial el strain
      eps_pl = 0.               ! initial pl strain
      l = 1.                    ! initial length
      t = 0.                    ! initial time
      do while(t<1.0)
c        Loading condition 1
         if (t.le.0.25) then
            vel = 0.01
         elseif (t.gt.0.25.and.t.le.0.55) then
            vel =-0.01
         else
            vel = 0.01
         endif
c        end of loading condition
         dl = vel * dt
         deps = dl / l
c        initially assuming all strain is elastic
         deps_pl = 0.0
         deps_el = deps
c        guess on stress increment
         dsig = E * deps_el
         kount = 0
         f = calc_yield_function( stress+ dsig) - c
         iplast=0
         do while (f.gt.tol.and.kount.lt.3) ! if exceeding plastic onset
            iplast=1
c           New plastic strain increment
            deps_pl = deps_pl - f/(-E)*sign(1.,deps)
c           New elastic strain increment
            deps_el = deps - deps_pl
c           New stress increment
            dsig = e * deps_el
            f = calc_yield_function(stress+dsig) - c
            kount = kount +1
         enddo
         write(3,'(2f9.4,f10.6,f10.2,2f10.6,2i2)')
     $        t,l,eps,stress,eps_el,eps_pl,iplast,kount ! write before update
         stress = stress + dsig
         eps = eps + deps
         eps_el = eps_el + deps_el
         eps_pl = eps_pl + deps_pl
         t = t + dt
         l = l + dl
      enddo
      close(3)
      end program




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
