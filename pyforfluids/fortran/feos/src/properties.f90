! TODO:
!   - Fully support operations oveload
!   - Add tests
module properties
    !! Datatype to represent different kind of properties with the objective
    !! to ease accesibility to derivatives and don't use an excess of arguments.
    use constants

    implicit none

    private
    public :: scalar_property
    public :: binary_property
    public :: null_scalar_property
    public :: null_binary_property

    public :: operator(+), operator(-), operator(*), &
              operator(**), operator(/), assignment(=)
    public :: print_scalar

    type :: scalar_property
        !! Scalar property
        real(wp) :: val                    !! Value
        real(wp) :: dt                     !! First deritvative with temperature
        real(wp) :: dt2                    !! Second derivative with temperature
        real(wp) :: dv                     !! First deritvative with volume
        real(wp) :: dv2                    !! Second derivative with volume
        real(wp) :: dtv                    !! Cross derivative with temperature and volume
        real(wp), allocatable :: dn(:)     !! First derivative with composition
        real(wp), allocatable :: dn2(:, :) !! Second derivative with composition
        real(wp), allocatable :: dvn(:)    !! Cross derivative with composition and volume
        real(wp), allocatable :: dtn(:)    !! Cross derivative with composition and temperature
    end type scalar_property

    type :: binary_property
        !! Binary property
        real(wp), allocatable :: val(:, :)  !! Value
        real(wp), allocatable :: dt(:, :)   !! First deritvative with temperature
        real(wp), allocatable :: dt2(:, :)  !! Second derivative with temperature
        real(wp), allocatable :: dv(:, :)   !! First deritvative with volume
        real(wp), allocatable :: dv2(:, :)  !! Second derivative with volume
        real(wp), allocatable :: dtv(:, :)  !! Crossed derivative with temperature and volume
    end type

    interface write
        module procedure :: print_scalar
    end interface

    interface operator(+)
        module procedure :: add_integer_scalar_left
        module procedure :: add_integer_scalar_right
        module procedure :: add_real_scalar_left
        module procedure :: add_real_scalar_right
        module procedure :: add_self_scalar
    end interface

    interface operator(-)
        module procedure :: sub_integer_scalar_left
        module procedure :: sub_integer_scalar_right
        module procedure :: sub_real_scalar_left
        module procedure :: sub_real_scalar_right
        module procedure :: sub_self_scalar
    end interface

    interface operator(*)
        module procedure :: product_real_scalar_left
        module procedure :: product_real_scalar_right
        module procedure :: product_integer_scalar_left
        module procedure :: product_integer_scalar_right
        module procedure :: product_self_scalar
    end interface

    interface operator(**)
        module procedure :: exponent_real_scalar_left
        module procedure :: exponent_real_scalar_right
        module procedure :: exponent_integer_scalar_left
        module procedure :: exponent_integer_scalar_right
    end interface

    interface operator(/)
        module procedure :: division_real_scalar_left
        module procedure :: division_real_scalar_right
        module procedure :: division_integer_scalar_left
        module procedure :: division_integer_scalar_right
        module procedure :: division_scalar_scalar
    end interface

    interface assignment(=)
        module procedure set_value_scalar
    end interface

contains
    ! =============================================================================
    !  Initializers
    ! -----------------------------------------------------------------------------
    pure function null_scalar_property(n) result(nsp)
        integer, intent(in) :: n
        type(scalar_property) :: nsp
        real(wp) :: null_1darr(n), null_2darr(n, n)

        allocate (nsp%dn(n), nsp%dvn(n), nsp%dtn(n), nsp%dn2(n, n))

        null_1darr = 0
        null_2darr = 0

        nsp%val = 0
        nsp%dt = 0
        nsp%dt2 = 0
        nsp%dv = 0
        nsp%dv2 = 0
        nsp%dtv = 0
        nsp%dn = null_1darr
        nsp%dn2 = null_2darr
        nsp%dvn = null_1darr
        nsp%dtn = null_1darr
    end function

    pure function null_binary_property(n) result(nbp)
        integer, intent(in) :: n
        type(binary_property) :: nbp

        allocate (nbp%val(n, n))
        allocate (nbp%dt(n, n))
        allocate (nbp%dt2(n, n))
        allocate (nbp%dv(n, n))
        allocate (nbp%dv2(n, n))
        allocate (nbp%dtv(n, n))

        nbp%val = 0
        nbp%dt = 0
        nbp%dt2 = 0
        nbp%dv = 0
        nbp%dv2 = 0
        nbp%dtv = 0
    end function
    ! =============================================================================

    ! =========================================================================
    !  Addition procedures
    !  addition overloads
    ! -------------------------------------------------------------------------

    elemental function add_integer_scalar_left(self, to_add) result(add)
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_add
        real(wp) :: add

        add = self%val + to_add
    end function add_integer_scalar_left

    elemental function add_integer_scalar_right(to_add, self) result(add)
        integer, intent(in) :: to_add
        class(scalar_property), intent(in) :: self
        real(wp) :: add
        add = self%val + to_add
    end function add_integer_scalar_right

    elemental function add_real_scalar_left(self, to_add) result(add)
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_add
        real(wp) :: add

        add = self%val + to_add
    end function add_real_scalar_left

    elemental function add_real_scalar_right(to_add, self) result(add)
        real(wp), intent(in) :: to_add
        class(scalar_property), intent(in) :: self
        real(wp) :: add

        add = self%val + to_add
    end function add_real_scalar_right

    elemental function add_self_scalar(self, to_add) result(add)
        class(scalar_property), intent(in) :: self
        class(scalar_property), intent(in) :: to_add
        real(wp) :: add

        add = self%val + to_add%val
    end function add_self_scalar
    ! =========================================================================

    ! =========================================================================
    !  Substraction procedures
    ! -------------------------------------------------------------------------
    elemental function sub_integer_scalar_left(self, to_sub) result(sub)
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_sub
        real(wp) :: sub

        sub = self%val - to_sub
    end function sub_integer_scalar_left

    elemental function sub_integer_scalar_right(to_sub, self) result(sub)
        integer, intent(in) :: to_sub
        class(scalar_property), intent(in) :: self
        real(wp) :: sub
        sub = self%val - to_sub
    end function sub_integer_scalar_right

    elemental function sub_real_scalar_left(self, to_sub) result(sub)
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_sub
        real(wp) :: sub

        sub = self%val - to_sub
    end function sub_real_scalar_left

    elemental function sub_real_scalar_right(to_sub, self) result(sub)
        real(wp), intent(in) :: to_sub
        class(scalar_property), intent(in) :: self
        real(wp) :: sub

        sub =  to_sub - self%val
    end function sub_real_scalar_right

    elemental function sub_self_scalar(self, to_sub) result(sub)
        class(scalar_property), intent(in) :: self
        class(scalar_property), intent(in) :: to_sub
        real(wp) :: sub

        sub = self%val - to_sub%val
    end function sub_self_scalar
    ! =========================================================================

    ! =========================================================================
    !  Product procedures
    ! -------------------------------------------------------------------------
    elemental function product_real_scalar_left(self, to_add) result(prod)
        ! Producto entre un real y una propiedad escalar
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_add
        real(wp) :: prod

        prod = self%val*to_add
    end function product_real_scalar_left

    elemental function product_real_scalar_right(to_add, self) result(prod)
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_add
        real(wp) :: prod

        prod = to_add*self%val
    end function product_real_scalar_right

    elemental function product_integer_scalar_left(self, to_add) result(prod)
        ! Producto entre un integer y una propiedad escalar
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_add
        real(wp) :: prod

        prod = self%val*to_add
    end function product_integer_scalar_left

    elemental function product_integer_scalar_right(to_add, self) result(prod)
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_add
        real(wp) :: prod

        prod = to_add*self%val
    end function product_integer_scalar_right
    
    elemental function product_self_scalar(self, other) result(prod)
        class(scalar_property), intent(in) :: self
        class(scalar_property), intent(in) :: other
        real(wp) :: prod

        prod = self%val * other%val
    end function product_self_scalar
    ! =========================================================================

    ! =========================================================================
    !  Exponent procedures
    ! -------------------------------------------------------------------------
    elemental function exponent_real_scalar_left(self, to_add) result(res)
        ! Producto entre un real y una propiedad escalar
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_add
        real(wp) :: res

        res = self%val**to_add
    end function exponent_real_scalar_left

    elemental function exponent_real_scalar_right(to_add, self) result(res)
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_add
        real(wp) :: res

        res = to_add**self%val
    end function exponent_real_scalar_right

    elemental function exponent_integer_scalar_left(self, to_add) result(res)
        ! Producto entre un integer y una propiedad escalar
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_add
        real(wp) :: res

        res = self%val**to_add
    end function exponent_integer_scalar_left

    elemental function exponent_integer_scalar_right(to_add, self) result(res)
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_add
        real(wp) :: res

        res = to_add**self%val
    end function exponent_integer_scalar_right
    ! =========================================================================

    ! =============================================================================
    !  Division procedures
    ! -----------------------------------------------------------------------------
    elemental function division_real_scalar_left(self, to_add) result(res)
        ! division entre un real y una propiedad escalar
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_add
        real(wp) :: res
        res = self%val/to_add
    end function division_real_scalar_left

    elemental function division_real_scalar_right(to_add, self) result(res)
        class(scalar_property), intent(in) :: self
        real(wp), intent(in) :: to_add
        real(wp) :: res

        res = to_add/self%val
    end function division_real_scalar_right

    elemental function division_integer_scalar_left(self, to_add) result(res)
        ! Producto entre un integer y una propiedad escalar
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_add
        real(wp) :: res

        res = self%val/to_add
    end function division_integer_scalar_left

    elemental function division_integer_scalar_right(to_add, self) result(res)
        class(scalar_property), intent(in) :: self
        integer, intent(in) :: to_add
        real(wp) :: res

        res = to_add/self%val
    end function division_integer_scalar_right
    
    elemental function division_scalar_scalar(self, other) result(res)
        class(scalar_property), intent(in) :: self
        class(scalar_property), intent(in) :: other
        real(wp) :: res

        res = self%val/other%val
    end function division_scalar_scalar

    ! =========================================================================

    ! =========================================================================
    !  Assignation of routines
    ! -------------------------------------------------------------------------
    elemental subroutine set_value_scalar(self, value)
        class(scalar_property), intent(in out) :: self
        real(wp), intent(in) :: value
        self%val = value
    end subroutine
    ! =========================================================================

    ! =========================================================================
    !  repr?
    ! -------------------------------------------------------------------------
    subroutine print_scalar(self)
        class(scalar_property), intent(in) :: self
        integer :: i
        print *, "shape", shape(self%dn2)
        print *, "val:", self%val
        print *, "dt:", self%dt
        print *, "dt2:", self%dt2
        print *, "dv:", self%dv
        print *, "dv2:", self%dv2
        print *, "dtv:", self%dtv
        print *, "dn:", self%dn
        do i = 1, 2
            write (*, *) "dn2: ", self%dn2(:, i)
        end do
        print *, "dvn:", self%dvn
        print *, "dtn:", self%dtn
    end subroutine print_scalar
    ! =========================================================================

end module properties
