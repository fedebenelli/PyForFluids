!TODO:
!      - Derivatives wrt t of aij
!      - kij exponentialy dependant of t

module mixing_rules
    !! Module that contains the available mixing rules to be used.
    use constants, only: wp, R
    use properties
    use cubic_eos

    implicit none

    private
    public :: ClassicVdW
    public :: aij_classic
    public :: bij_classic

    type, abstract :: MixingRule
    contains
        procedure(abs_mix), deferred :: mix
    end type

    abstract interface
        !! Abstract interface for mixing rules
        subroutine abs_mix(self, compounds, d_mixture, b_mixture, d1_mixture)
            !! Abstract mixing routine
            use properties
            use cubic_eos
            import MixingRule
            class(MixingRule), intent(in out) :: self
            class(CubicEoS), intent(in out) :: compounds(:)

            type(scalar_property), intent(out) :: d_mixture
            type(scalar_property), intent(out) :: b_mixture
            type(scalar_property), intent(out) :: d1_mixture
        end subroutine abs_mix
    end interface

    type, extends(MixingRule) :: ClassicVdW
        real(wp), allocatable :: kij(:, :) !! Kij matrix
        real(wp), allocatable :: lij(:, :) !! lij matrix
    contains
        procedure :: mix => classic_mix
    end type

contains

    ! ==========================================================================
    !  Constant kij and lij
    ! --------------------------------------------------------------------------
    pure subroutine classic_mix(self, compounds, d_mixture, b_mixture, d1_mixture)
        !! Classic cuadratic mixing rules
        class(ClassicVdW), intent(in out) :: self
        class(CubicEoS), intent(in out) :: compounds(:)

        type(scalar_property), intent(out) :: d_mixture
        type(scalar_property), intent(out) :: b_mixture
        type(scalar_property), intent(out) :: d1_mixture

        type(scalar_property) :: a(size(compounds))
        type(scalar_property) :: b(size(compounds))
        type(scalar_property) :: d1(size(compounds))

        type(binary_property) :: aij, bij

        real(wp) :: moles(size(compounds))

        integer :: i, j, n

        n = size(compounds)

        ! =====================================================================
        !  Initialization of parameters
        ! ---------------------------------------------------------------------
        d_mixture = null_scalar_property(n)
        b_mixture = null_scalar_property(n)
        d1_mixture = null_scalar_property(n)

        aij = null_binary_property(n)
        bij = null_binary_property(n)

        moles = compounds%moles
        ! =====================================================================

        a = compounds%a_parameter()
        b = compounds%b_parameter()

        aij = aij_classic(n, a, self%kij)
        bij = bij_classic(n, b, self%lij)

        d_mixture = d_mix(n, aij, moles)
        b_mixture = b_mix(n, bij, moles)

        ! All del1 derivatives =0 for a two parameters eos
        d1_mixture%val = compounds(1)%del1
    end subroutine classic_mix

    pure function aij_classic(n, a, kij) result(aij)
        integer, intent(in) :: n
        type(scalar_property), intent(in) :: a(n)  !! Atractive parameters
        real(wp), intent(in) :: kij(n, n)          !! kij
        type(binary_property) :: aij               !! aij

        integer :: i, j

        aij = null_binary_property(n)

        do i = 1, n
            aij%val(i, i) = a(i)%val
            aij%dt(i, i) = a(i)%dt
            aij%dt2(i, i) = a(i)%dt2
            do j = 1, i-1
                aij%val(j, i) = sqrt(a(i)%val*a(j)%val) * (1 - kij(i, j))
                aij%dt(j, i) = (1 - kij(j, i)) * (sqrt(a(i)/a(j)) * a(j)%dt + sqrt(a(j)/a(i)) * a(i)%dT)/2
                aij%dt2(j, i) = (1 - kij(j, i))*(a(j)%dt * a(i)%dt / sqrt(a(i)*a(j)) &
                               + sqrt(a(i)/a(j))*(a(j)%dt2 - a(j)%dt**2/(2*a(j))) &
                               + sqrt(a(j)/a(i))*(a(i)%dt2 - a(i)%dt**2/(2*a(i))))/2

                aij%val(i, j) = aij%val(j, i)
                aij%dt(i, j) = aij%dt(j, i)
                aij%dt2(i, j) = aij%dt2(j, i)
            end do
        end do
    end function aij_classic

    pure function bij_classic(n, b, lij) result(bij)
        integer, intent(in) :: n
        type(scalar_property), intent(in) :: b(n)   !! Repulsive parameters
        real(wp), intent(in) :: lij(n, n)           !! lij matrix

        type(binary_property) :: bij

        integer :: i, j

        bij = null_binary_property(n)

        do i = 1, n
        do j = i, n
            bij%val(i, j) = (1 - lij(i, j)) * (b(i) + b(j))/2_wp
            bij%val(j, i) = bij%val(i, j)
        end do
        end do
    end function bij_classic

    pure function d_mix(n, aij, moles) result(d)
       !! Mixture attractive parameter times moles^2
        integer, intent(in) :: n
        type(binary_property), intent(in) :: aij
        real(wp), intent(in) :: moles(n)
        type(scalar_property) :: d

        real(wp) :: aux, aux2
        integer :: i, j

        d = null_scalar_property(n)
        do i = 1, n
            aux = 0.0_wp
            aux2 = 0.0_wp
            d%dn(i) = 0.0_wp
            d%dtn(i) = 0.0_wp
            do j = 1, n
                d%dn(i) = d%dn(i) + 2*moles(j)*aij%val(i, j)
                d%dtn(i) = d%dtn(i) + 2*moles(j)*aij%dt(i, j)
                d%dn2(i, j) = 2*aij%val(i, j)

                aux = aux + moles(j)*aij%val(i, j)
                aux2 = aux2 + moles(j)*aij%dt2(i, j)
            end do
            d%val = d%val + moles(i)*aux
            d%dt = d%dt + moles(i)*d%dtn(i)/2
            d%dt2 = d%dt2 + moles(i)*aux2
        end do
    end function d_mix

    pure function b_mix(n, bij, moles) result(b)
        !! Mixture repulsive parameter
        integer, intent(in) :: n                     !! Number of components
        type(binary_property), intent(in) :: bij  !! bij matrix
        real(wp), intent(in) :: moles(n)             !! Concentrations vector
        type(scalar_property) :: b                !! Mixture repulsive parameter

        real(wp) :: totn ! Total number of moles
        real(wp) :: aux(n)
        integer :: i, j

        b = null_scalar_property(n)

        totn = sum(moles)
        aux = 0.0_wp

        do i = 1, n
            do j = 1, n
                aux(i) = aux(i) + moles(j)*bij%val(i, j)
            end do
            b%val = b%val + moles(i)*aux(i)
        end do
        b%val = b%val/totn
        do i = 1, n
            b%dn(i) = (2*aux(i) - b%val)/totn
            do j = 1, i
                b%dn2(i, j) = (2*bij%val(i, j) - b%dn(i) - b%dn(j))/totn
                b%dn2(j, i) = b%dn2(i, j)
            end do
        end do
    end function b_mix

    pure function delta1_mix(n, delta1, moles) result(d1)
        integer, intent(in) :: n
        real(wp), intent(in) :: delta1(n)
        real(wp), intent(in) :: moles(n)

        type(scalar_property) :: d1

        real(wp) :: totn
        integer :: i, j

        totn = sum(moles)
        d1 = null_scalar_property(n)

        do i = 1, n
            d1%val = d1%val + moles(i)*delta1(i)
        end do
        d1%val = d1%val/totn

        do i = 1, n
            d1%dn(i) = (delta1(i) - d1%val)/totn
            do j = 1, n
                d1%dn2(i, j) = (2.0_wp*d1%val - delta1(i) - delta1(j))/totn**2
            end do
        end do
    end function delta1_mix
    ! ==========================================================================
end module mixing_rules
