module cubic
    implicit none
contains

        subroutine set_model(&
                nc, model, &
                critical_atractive, repulsive_parameter, delta_1_parameter, critical_k_parameter, &
                critical_temperature, critical_pressure, critical_density, accentric_factor, &
                mix_rule, temperature_dependence, kij0_matrix, kijinf_matrix, T_star_matrix, &
                lij_matrix, &
                volume_traslation, volume_shift &
            )
            ! ==================================================================
            !  Python input
            ! ------------------------------------------------------------------
            integer, intent(in) :: nc
            integer, intent(in) :: model

            ! Mixing rule
            integer, intent(in) :: mix_rule
            integer, intent(in) :: temperature_dependence

            ! Pure compound parameters
            real(8), intent(in) :: critical_atractive(nc)
            real(8), intent(in) :: repulsive_parameter(nc)
            real(8), intent(in), optional :: delta_1_parameter(nc)
            real(8), intent(in) :: critical_k_parameter(nc)

            real(8), intent(in) :: critical_temperature
            real(8), intent(in) :: critical_pressure
            real(8), intent(in) :: critical_density
            real(8), intent(in) :: accentric_factor
    
            ! - Kij values 
            real(8), intent(in) :: kij0_matrix(nc, nc)
            real(8), intent(in), optional :: kijinf_matrix(nc, nc)
            real(8), intent(in), optional :: T_star_matrix(nc, nc)

            ! - lij values
            real(8), intent(in) :: lij_matrix(nc, nc)

            ! Peneloux volume traslation
            integer, intent(in) :: volume_traslation
            real(8), intent(in), optional :: volume_shift(nc)
            ! ==================================================================
            
            integer :: i, j

            ! ==================================================================
            !  Original commons
            ! ------------------------------------------------------------------
            integer, parameter :: nco=64
            integer :: nmodel
            integer :: ncomb
            integer :: iVshift
            integer :: ntdep

            real(8) :: Vs
            real(8) :: ac(nco), b(nco), del1(nco), rk(nco), Kij(nco, nco), bij(nco, nco), lij(nco,nco)
            real(8) :: tc(nco), pc(nco), dceos(nco), om(nco)
            real(8) :: Kinf(nco, nco), Tstar(nco, nco)

            common /model/ nmodel
            common /rule/ ncomb
            common /Vshift/ iVshift, Vs(nco)
            common /components/ ac, b, del1, rk, Kij, NTDEP
            common /CRIT/ tc, pc, dceos, om
            common /bcross/ bij
            common /lforin/ lij
            common /Tdep/ Kinf, Tstar
            ! ==================================================================

            ! ==================================================================
            !  Assignation of variables
            ! ------------------------------------------------------------------
            
            ! Equation of state
            nmodel = model
            
            ! Pure compounds parameters
            ac = critical_atractive
            b = repulsive_parameter
            select case(nmodel)
            case(1)
                ! SRK
                del1 = 1.d0
            case(2)
                ! PR
                del1 = 1.d0 + dsqrt(2.d0)
            case(3)
                ! RKPR
                del1 = delta_1_parameter
            end select
            rk(:nc) = critical_k_parameter

            tc = critical_temperature
            pc = critical_pressure
            dceos = critical_density
            om = accentric_factor

            ! Mixing rule
            ncomb = mix_rule
            ntdep = temperature_dependence

            select case(ntdep)
            case(0)
                kij = kij0_matrix
            case (1)
                kij = kij0_matrix
                kinf = kijinf_matrix
                tstar = T_star_matrix
            end select

            select case(ncomb)
            case(0)
                do i = 1, nc
                    do j = i, nc
                       bij(i, j) = (1.d0 - lij(i, j))*(b(i) + b(j))/2.d0
                       bij(j, i) = bij(i, j)
                    end do
                end do
                lij = lij_matrix
            end select

            ! Volume translation
            iVshift = volume_traslation
            Vs = volume_shift
            ! ==================================================================
        end subroutine set_model

        !subroutine helmholtz(nc, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
        !    integer, intent(in) :: nc
        !    real(8), intent(in) :: rn(nc), V, T
        !    real(8), intent(out) :: Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2
        !    call ArVnder(nc, 1, 2, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
        !end subroutine helmholtz

        subroutine lnfug(nc, root_type, indicator, T, P, rn, V, PHILOG, DLPHIP, DLPHIT, FUGN)
            integer, intent(in) :: nc, root_type, indicator
            real(8), intent(in) :: T, P, rn(nc)

            real(8), intent(out) :: V
            real(8), intent(out) :: PHILOG(nc), DLPHIT(nc), DLPHIP(nc)
            real(8), intent(out) :: FUGN(nc)
            
            print *, ""

           call TERMO(nc, root_type, indicator, T, P, rn, V, PHILOG, DLPHIP, DLPHIT, FUGN)
        end subroutine lnfug

        subroutine pressure_calc(nc, concentrations, volume, temperature, pressure)
            integer, intent(in) :: nc
            real(8), intent(in) :: volume, temperature, concentrations(nc)
            real(8), intent(out) :: pressure

            real(8) :: Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2
            real(8) :: totn, RGAS
            integer :: nder, ntemp
            parameter(RGAS=0.08314472d0)
            
            print *, ""

            nder = 0
            ntemp = 1
            totn = sum(concentrations)

            call ArVnder(&
                nc, nder, ntemp, concentrations, volume, temperature, &
                Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2 &
            )
            
            pressure = totn*RGAS*temperature/volume - ArV
        end subroutine pressure_calc

end module cubic
