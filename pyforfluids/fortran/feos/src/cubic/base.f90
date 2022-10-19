module base_ceos
   implicit none

   private
   public :: CEOS

   type, abstract :: CEOS
   contains
      procedure(a_parameter), deferred :: a_parameter
      procedure(b_parameter), deferred :: b_parameter
   end type

   abstract interface
      elemental function a_parameter(self) result(a)
         use properties, only: scalar_property
         import CEOS
         class(CEOS), intent(in) :: self
         type(scalar_property) :: a
      end function
      elemental function b_parameter(self) result(b)
         use properties, only: scalar_property
         import CEOS
         class(CEOS), intent(in) :: self
         type(scalar_property) :: b
      end function
   end interface
end module base_ceos
