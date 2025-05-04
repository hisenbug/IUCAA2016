! 4th Order symplectic Method (Forrest-Ruth Algorithm) for Kuzmin disc with 
! asymetric factor

program forest_ruth
    implicit none
    real(kind = 8), parameter :: year = 3.1556926d7, ms=1.989d30, me=5.973d24
    ! Seconds in years, Mass of sun, Mass of Earth
    real(kind = 8), parameter :: PI = 3.14159265358979323d0, G=6.67408d-11
    real(kind = 8), parameter :: lda = 1.35120719195965763d0, pc=3.085677581d16 
    ! SI4 constant=1.0/(2-2**(1.0/3.0)), Parsec in meter
	real(kind = 8), parameter :: rc = G*ms/256d6, vc = 16d3, tc = rc/vc
	! Scales to make quatities dimensionless
	! Taking $ G M_\odot t_c^2=r_c^3. \ and \ v_c = v/16d3 $
	real(kind = 8), parameter :: mr = ms/me, au = 1.496d11 !1.495978707d11
	! Mass Ratio, Astronomical Unit
	real(kind = 8), parameter :: leng = 2.5d-3*pc/rc, sma = 100*au, ec = 0.05d0
	! Length Scale, Semi-Major Axis, Eccentricity
	real(kind = 8), parameter :: ls2 = leng**2, mew = G*ms, cl1=rc/sma*1.0
	real(kind = 8), parameter :: md = 1d-1, om_p = 0, ep = 1d-2
	! Mass of Disc, Pattern Speed, Asymmetry parameter
	integer, parameter :: D1 = 2, D2 = 2 ! No. of bodies, Dimensions
	real(kind = 8), parameter :: tp = 2*PI*sqrt(sma**3/mew) !Time-period guess
	
	real(kind = 8) :: x(D1,D2), v(D1,D2), a(D1,D2)
	! Motion elements for Star(1), Planet(2). All =0 for Star
	real(kind = 8) :: dt=year*1d-2/tc, tt=1.01d4*tp/tc, t=0d0, tl(2)
	! Time step, Total Time, Intial time, Saving control
	real(kind = 8) :: prev(D2+2), st0, sqdci, r2, rin, ri3
	! Temproray variable to store or optimize
	integer :: i=0, j, samp=10000, flag=0 ! samp is sampling interval
	real(kind = 8) :: L02, L2, E02, e2, tor
	real(kind = 8) :: c(4), d(3), om, kp, pr, tol, rdi= 1d0, ratio
	! c, d: SI4 constants; 

	x(2,:) = [ sma*(1d0-ec)/rc, 0.0d0/rc ] ! Perigee
	x(1,:) = 0
	v(2,:) = [0d0, sqrt(mew*(2d0/(sma*(1d0-ec))-1d0/sma))/vc] !vis-viva eqn
	v(1,:) = 0 !-1d0*v(2,:)/mr

	c = [ 0.5d0*lda, 0.5d0*(1.0d0-lda), 0.5d0*(1.0d0-lda), 0.5d0*lda ]
	d = [ lda, 1.0d0-2.0d0*lda, lda ]
	tol = 2d-5/rc*sma ! Tolerance for SoS
	
	om = sqrt(mew*(sma**(-3) + md*sqrt((sma*sma+ls2*rc*rc)**(-3))))*year
	kp = sqrt(mew*(sma**(-3) + md*(sqrt(sma*sma+ls2*rc*rc)**(-5))*(4*ls2*rc*rc+sma*sma)))*year
	pr = om - kp

	write(*,*) "om=", om, 2*PI/om, "kp=", kp, 2*PI/kp, "pr=", pr, 2*PI/pr
	write(*,*) "v_c =", om/year*sma, st0

	open(11, file='fckr7-2e.05_1c.01.txt')
	! f = Forest-Ruth, k = KT Disc, r = rotating, e= eccentricity, _1=disc_mass, c=center_diff, .01=ep
	open(11, file='fckr6-3e.05_1c.txt') 
	open(12, file='fckr7-2e.05_1-2ssx.txt')
	tl(1) = 2.5d5*tp/tc
	tl(2) = tt-tl(1)
	r2 = sum((x(2,:))**2)
	rin = 1.0d0/sqrt(r2)
	!****************** cos(atan2(y,x)) = x/r **************************
	E02 = -rin + 0.5d0*sum(v(2,:)**2) - 0.5*om_p*om_p*r2 &
			- md*(1d0+ep*x(2,1)*rin)/sqrt(r2+ls2)
    L02 = x(2,1)*v(2,2) - x(2,2)*v(2,1)
    	
    write(*,*) "dt = ", dt, "rc = ", rc, "T_p = ", tp/year
	prev = [x(2,1), x(2,2)+1d-15, v(2,1), v(2,2)]
	st0 = x(2,2)

	do while(t<=tt)
		if(i==samp) i=0
		! ~ 60 points/orbit needed for good ellipse
		if( i==0 .and. (t<tl(1) .or. t>tl(2)) ) then
			! Only every {samp}^th value is recorded due to memory limitations
			r2 = sum(x(2,:)**2)
			rin = 1.0d0/sqrt(r2)
			e2 = -rin + 0.5d0*sum(v(2,:)**2) - 0.5*om_p*om_p*r2 &
				- md*(1d0+ep*x(2,1)*rin)/sqrt(r2+ls2)
            L2 = x(2,1)*v(2,2) - x(2,2)*v(2,1)
            ! th = atan2(x(2,2),x(2,1))
			r2 = sum(x(2,:)**2)
			rin = 1.0d0/sqrt(r2)
			ri3 = rin*rin*rin
			sqdci = 1d0/sqrt(r2+ls2)
			a(2,1) = -ri3*x(2,1) + ep*md*(x(2,2)**2)*ri3*sqdci &
					- md*x(2,1)*(1d0+ep*x(2,1)*rin)*(sqdci**3)
			a(2,2) = -ri3*x(2,2) - ep*md*x(2,2)*x(2,1)*ri3*sqdci &
					- md*x(2,2)*(1d0+ep*x(2,1)*rin)*(sqdci**3)
			tor = x(2,1)*a(2,2) - x(2,2)*a(2,1)
			
			write(11,'(10g13.5)') x(2,:)*cl1, t*tc/year, e2/E02-1.0, v(2,:)*vc, &
				L2/L02-1.0, rdi*x(2,1)*rin, rdi*x(2,2)*rin, tor
		end if
		
		!****This block finds the surface of section with y=0, v_y>0*****
		if( (x(2,2)*prev(2) .le. 0) .and. (v(2,2)>0) ) then
			if( (abs(x(2,2))+abs(prev(2))<tol) .and. (abs(prev(2)) > 0) ) then
				ratio = abs(x(2,2))/abs(prev(2))
				write(12,'(4g15.7)') (x(2,1)+ratio*prev(1))*cl1/(ratio+1), &
					(v(2,1)+ratio*prev(3))*vc/(ratio+1), x(2,2)*cl1, t*tc/year
			else if( abs(x(2,2))<tol/4 ) then
				write(12,'(4g15.7)') x(2,1)*cl1, v(2,1)*vc, x(2,2)*cl1, t*tc/year
			else if( abs(prev(2))<tol/4 ) then
				write(12,'(4g15.7)') prev(1)*cl1, prev(3)*vc, prev(2)*cl1, t*tc/year
			end if
		end if
		prev = [x(2,1), x(2,2), v(2,1), v(2,2)]

		!****************Main loop of SI4 algorithm**********************
		do j=1, 4, 1
			x(2,:) = x(2,:) + c(j)*dt*v(2,:)
			if(j==4) exit

			r2 = sum(x(2,:)**2)
			rin = 1.0d0/sqrt(r2)
			ri3 = rin*rin*rin
			sqdci = 1d0/sqrt(r2+ls2)
			a(2,1) = -ri3*x(2,1) + ep*md*(x(2,2)**2)*ri3*sqdci &
					- md*x(2,1)*(1d0+ep*x(2,1)*rin)*(sqdci**3)
			a(2,2) = -ri3*x(2,2) - ep*md*x(2,2)*x(2,1)*ri3*sqdci &
					- md*x(2,2)*(1d0+ep*x(2,1)*rin)*(sqdci**3)

			v(2,:) = v(2,:) + d(j)*dt*a(2,:)
		end do
		
		if( (flag<2) .and. (st0*x(2,2)<0) ) then ! Find orbital Time-Period
			flag = flag+1
			if(flag==2) write(*,*) "Time=", tc*t/year
		endif
		st0 = x(2,2)

		t = t+dt
		i = i+1
	end do
	
	close(11)
	close(12)
	write(*,*) "Final Time=", tc*t/year
end program forest_ruth

! * Compile/Execute Commands
! * gfortran fr_kt_fx_rot.f95 -O3 -ffast-math
! * time ./a.out
! * Gnuplot is used for plotting
