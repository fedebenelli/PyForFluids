Module read_params  !
Contains
    Subroutine read_parameters(no_file, vo_file, no, vo, &
                               Kr_file, nr_file, dr_file, tr_file, cr_file, &
                               Kpol, Kexp, nr, dr, tr, cr, &
                               Kij_file, betaij_file, dij_file, epsij_file, &
                               etaij_file, gammij_file, nij_file, tij_file, &
                               Fij_file, &
                               Kpolij, Kexpij, betaij, dij, epsij, etaij, &
                               gammij, nij, tij, Fij, &
                               reducing_file, Bv, Gv, Bt, Gt, &
                               critical_file, rho_c, T_c, M)

        Implicit None
        ! Ideal gas parameters
        character(len=100), intent(in):: no_file, vo_file
        real*8, dimension(21, 7):: no, vo
        ! Residual energy parameteres
        character(len=100), intent(in)::Kr_file, nr_file, dr_file, tr_file, cr_file
        real*8, intent(inout), dimension(21, 24):: nr, tr
        integer, intent(inout), dimension(21, 24):: cr, dr
        integer, intent(inout):: Kpol(21), Kexp(21)
        ! Departure function parameters
        character(len=100), intent(in):: Fij_file, Kij_file, betaij_file, dij_file, epsij_file, &
                                         etaij_file, gammij_file, nij_file, tij_file
        integer:: Kpolij(21, 21), Kexpij(21, 21)
        real*8, dimension(21, 21, 12), intent(inout):: betaij, epsij, etaij, gammij
        real*8, dimension(21, 21, 12), intent(inout):: nij, tij
        integer, intent(inout):: dij(21, 21, 12)
        real*8, dimension(21, 21):: Fij
        ! Reducing functions parameters
        character(len=100), intent(in):: reducing_file
        real*8, dimension(21, 21), intent(inout):: Bv, Gv, Bt, Gt
        ! Critical properties
        character(len=100):: critical_file
        real*8, dimension(21), intent(inout):: rho_c, T_c, M
        ! Common parameters
        integer:: i, j, k, io

        no = 0
        vo = 0
        Kpol = 0
        Kexp = 0
        nr = 0
        dr = 0
        tr = 0
        cr = 0
        Kpolij = 0
        Kexpij = 0
        betaij = 0
        dij = 0
        epsij = 0
        etaij = 0
        gammij = 0
        nij = 0
        tij = 0
        Fij = 0
        Bv = 0
        Gv = 0
        Bt = 0
        Gt = 0
        rho_c = 0
        T_c = 0
        M = 0

        ! Get ideal gas parameters
        ! ------------------------
        open (1, file=no_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, (no(i, j), j=1, 7)
        end do
        close (1)

        open (1, file=vo_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, (vo(i, j), j=4, 7)
        end do
        close (1)

        ! Get Residual parameters
        ! -----------------------
        open (1, file=Kr_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, Kpol(i), Kexp(i)
        end do
        close (1)
        open (1, file=nr_file)

        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, (nr(i, j), j=1, Kpol(i) + Kexp(i))
        end do
        close (1)

        open (1, file=dr_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, (dr(i, j), j=1, Kpol(i) + Kexp(i))
        end do
        close (1)

        open (1, file=cr_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, (cr(i, j), j=Kpol(i) + 1, Kpol(i) + Kexp(i))
        end do
        close (1)

        open (1, file=tr_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, (tr(i, j), j=1, Kpol(i) + Kexp(i))
        end do

        ! Get reducing function parameters
        ! --------------------------------
        open (1, file=reducing_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, Bv(i, j), Gv(i, j), Bt(i, j), Gt(i, j)
        end do
        close (1)

        ! Get departure function parameters
        ! ---------------------------------
        open (1, file=Kij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, Kpolij(i, j), Kexpij(i, j)
        end do
        close (1)
        open (1, file=Fij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, Fij(i, j)
        end do
        close (1)

        !! Polynomial parameters
        open (1, file=dij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, (dij(i, j, k), k=1, Kpolij(i, j) + Kexpij(i, j))
        end do
        close (1)
        open (1, file=nij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, (nij(i, j, k), k=1, Kpolij(i, j) + Kexpij(i, j))
        end do
        close (1)
        open (1, file=tij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, (tij(i, j, k), k=1, Kpolij(i, j) + Kexpij(i, j))
        end do
        close (1)

        !! Exponential parameters
        open (1, file=betaij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, (betaij(i, j, k), k=Kpolij(i, j) + 1, Kpolij(i, j) + Kexpij(i, j))
        end do
        close (1)
        open (1, file=epsij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, (epsij(i, j, k), k=Kpolij(i, j) + 1, Kpolij(i, j) + Kexpij(i, j))
        end do
        close (1)
        open (1, file=etaij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, (etaij(i, j, k), k=Kpolij(i, j) + 1, Kpolij(i, j) + Kexpij(i, j))
        end do
        close (1)
        open (1, file=gammij_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, j, (gammij(i, j, k), k=Kpolij(i, j) + 1, Kpolij(i, j) + Kexpij(i, j))
        end do
        close (1)

        ! Get critical properties
        ! -----------------------
        open (1, file=critical_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, rho_c(i), T_c(i), M(i)
        end do

        !call robbed_ideal(no, vo)
    End Subroutine read_parameters
End Module read_params

! ------------------------
! Thermodynamic Properties
! ------------------------
Module thermo_props  !
    Implicit None

Contains
    Subroutine zeta(delta, ar, z)
        real*8, intent(in):: delta, ar(3, 3)
        real*8, intent(inout):: z
        z = 1.d0 + delta*ar(2, 1)
    End Subroutine zeta

    Subroutine isochoric_heat(tau, R, Ao, Ar, cv) 
        real*8, intent(in):: tau, R, Ao(3, 3), Ar(3, 3)
        real*8, intent(inout)::cv
        
        cv = -tau**2.d0*(Ao(3, 2) + Ar(3, 2))
        cv = cv*R
    End Subroutine isochoric_heat

    Subroutine isobaric_heat(delta, tau, R, Ao, Ar, cp)
        real*8, intent(in):: delta, tau, R, Ao(3, 3), Ar(3, 3)
        real*8, intent(inout):: cp
        real*8:: up, down

        up = (1.d0 + delta*Ar(2, 1) - delta*tau*Ar(3, 3))**2
        down = 1.d0 + 2.d0*delta*Ar(2, 1) + delta**2*Ar(3, 1)

        cp = -tau**2.d0*(Ao(3, 2) + Ar(3, 2)) + up/down
        cp = cp*R
    End Subroutine isobaric_heat

    Subroutine sound_speed(delta, tau, R, T, M, Ao, Ar, w)
        real*8, intent(in):: delta, tau, R, T, M, Ao(3, 3), Ar(3, 3)
        real*8, intent(out):: w
        real*8:: up, down

        up = (1.d0 + delta*Ar(2, 1) - delta*tau*Ar(3, 3))**2
        down = (tau**2*(ao(3, 2) + ar(3, 2)))
        
        w = 1.d0 + 2*delta*ar(2, 1) + delta**2*ar(3, 1) - up/down
        w = sqrt(w*R*T/M)
    End Subroutine sound_speed

    Subroutine isothermal_thermal_coefficent(delta, tau, rho, Ar, delta_t)
        real*8, intent(in):: delta, Ar(3, 3), tau, rho
        real*8, intent(out):: delta_t
        real*8:: up, down

        up = 1 + delta*Ar(2, 1) - delta*tau*Ar(3, 3)
        down = 1 + 2*delta*Ar(2, 1) + delta**2*Ar(3, 1)
        
        delta_t = (1 - up/down)/rho
    End Subroutine isothermal_thermal_coefficent

    Subroutine pressure(delta, rho, R, T, Ar, p)
        real*8, intent(in):: delta, rho, R, T, Ar(3,3)
        real*8, intent(out):: p

        p = (1 + delta*Ar(2,1))*rho*R*T
    End Subroutine pressure
    
    Subroutine entropy(tau, R, Ao, Ar, s) 
        real*8, intent(in):: tau, R, Ao(3,3), Ar(3,3)
        real*8, intent(out):: s 

        s = (tau*(Ao(2,2) + Ar(2,2)) - Ao(1,1) - Ar(1,1))*R
    End Subroutine entropy

    Subroutine internal_energy(tau, R, T, Ao, Ar, u)
        real*8, intent(in):: tau, R, T, Ao(3,3), Ar(3,3)
        real*8, intent(out):: u

        u = tau*(Ao(2,2) + Ar(2,2))*R*T
    End Subroutine internal_energy

    Subroutine enthalpy(delta, tau, R, T, Ao, Ar, h)
        real*8, intent(in):: delta, tau, R, T, Ao(3,3), Ar(3,3)
        real*8, intent(out):: h

        h = ((1 + tau*(Ao(2,2) + Ar(2,2)) + delta*Ar(2,1)))*R*T
    End Subroutine enthalpy

    Subroutine gibbs_free_energy(delta, R, T, Ao, Ar, g) 
        real*8, intent(in):: delta, R, T, Ao(3,3), Ar(3,3)
        real*8, intent(out):: g

        g = (1 + Ao(1,1) + Ar(1,1) + delta*Ar(2,1))*R*T
    End Subroutine gibbs_free_energy

    Subroutine joule_thomson_coeff(delta, tau, rho, R, Ao, Ar, JT)
        real*8, intent(in):: delta, tau, rho, R, Ao(3,3), Ar(3,3)
        real*8, intent(out):: JT
        real*8:: up, down

        up = delta*Ar(2,1) + delta**2*Ar(3,1) + delta*tau*Ar(2,3)
        down = (1 + delta*Ar(2,1) - delta*tau*Ar(2,3))**2 - tau**2*(Ao(3,2) + Ar(3,2)*(1 + 2*delta*Ar(2,1) + delta**2*Ar(3,1)))

        JT = -up/down*R*rho
    End Subroutine joule_thomson_coeff

    Subroutine isentropic_exponent(delta, tau, Ao, Ar, k)
        real*8, intent(in):: delta, tau, Ao(3,3), Ar(3,3)
        real*8, intent(out):: k
        real*8:: up1, down1, up2, down2

        up1 = 1 + 2*delta*Ar(2,1) + delta**2*Ar(3,1)
        down1 = 1 + delta*Ar(2,1)
        up2 = (1 + delta*Ar(2,1) - delta*tau*Ar(2,3))**2
        down2 = tau**2*(Ao(3,2) + Ar(2,2))*(1 + 2*delta*Ar(2,1) + delta**2*Ar(3,2))

        k = (up1/down1)*(1 - up2/down2)
    End Subroutine isentropic_exponent


End Module thermo_props
