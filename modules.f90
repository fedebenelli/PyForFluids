Module get_params  !
Contains
    Subroutine get_parameters(no_file, vo_file, no, vo, &
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
        double precision, dimension(21, 7):: no, vo
        ! Residual energy parameteres
        character(len=100), intent(in)::Kr_file, nr_file, dr_file, tr_file, cr_file
        double precision, intent(inout), dimension(21, 24):: nr, dr, tr, cr
        integer, intent(inout):: Kpol(21), Kexp(21)
        ! Departure function parameters
        character(len=100), intent(in):: Fij_file, Kij_file, betaij_file, dij_file, epsij_file, &
                                         etaij_file, gammij_file, nij_file, tij_file
        integer:: Kpolij(21, 21), Kexpij(21, 21)
        double precision, dimension(21, 21, 12), intent(inout):: betaij, epsij, etaij, gammij
        double precision, dimension(21, 21, 12), intent(inout):: dij, nij, tij
        double precision, dimension(21, 21):: Fij
        ! Reducing functions parameters
        character(len=100), intent(in):: reducing_file
        double precision, dimension(21, 21), intent(inout):: Bv, Gv, Bt, Gt
        ! Critical properties
        character(len=100):: critical_file
        double precision, dimension(21), intent(inout):: rho_c, T_c, M
        ! Common parameters
        integer:: i, j, k, compound, io, N = 21

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

        ! Get ideal gas parameters
        ! ------------------------
        open (1, file=no_file)
        do i = 1, N
            read (1, *) compound, (no(i, j), j=1, 7)
        end do
        close (1)

        open (1, file=vo_file)
        do i = 1, N
            read (1, *) compound, (vo(i, j), j=4, 7)
        end do
        close (1)

        ! Get Residual parameters
        ! -----------------------
        open (1, file=Kr_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) compound, Kpol(compound), Kexp(compound)
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

        ! Get critical properties
        ! -----------------------

        open (1, file=critical_file)
        io = 0
        do while (io == 0)
            read (1, *, iostat=io) i, rho_c(i), T_c(i), M(i)
        end do
    End Subroutine get_parameters
End Module get_params

! ------------------------
! Thermodynamic Properties
! ------------------------
Module thermo_props  !
    Implicit None

Contains
    Subroutine zeta(delta, dArdD, z)
        double precision:: delta, dArdD
        double precision, intent(inout):: z
        z = 1 + delta*dArdD
    End Subroutine zeta

    Subroutine isobaric_heat(delta, tau, ao, ar, cp)
        double precision, intent(in):: delta, tau, ao(3, 3), ar(3, 3)
        double precision, intent(inout) :: cp
        cp = -tau**2*(ao(3, 2) + ar(3, 2)) &
             + (1 + delta*ar(2, 1) - delta*tau*ar(3, 3))**2 &
             /(1 + 2*delta*ar(2, 1) + delta**2*ar(3, 1))
    End Subroutine isobaric_heat

    Subroutine sound_speed(R, T, M, delta, tau, Ar, Ao, w)
        double precision, intent(in):: R, T, M, delta, tau, Ar(3, 3), Ao(3, 3)
        double precision, intent(inout):: w
        w = 1 + 2*delta*Ar(2, 1) + delta**2*Ar(3, 1)
        w = w - (1 + delta*Ar(2, 1) - delta*tau*Ar(3, 3))**2 &
            /(tau**2*(Ao(3, 2) + Ar(3, 2)))
        w = sqrt(w*R*T/M)
    End Subroutine sound_speed
End Module thermo_props

