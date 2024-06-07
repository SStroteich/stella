!> Define abstract type to use as a base type for derived types used to represent
!> the input configuration of a module (actual or logical). In other words it can
!> be used to represent a namelist state.
module abstract_config
   implicit none
   private

   !Types
   public :: abstract_config_type, CONFIG_MAX_NAME_LEN

   integer, parameter :: CONFIG_MAX_NAME_LEN = 64
   character(len=CONFIG_MAX_NAME_LEN), parameter :: default_name = "Undefined name."

   !> Controls the alignment of key/value pairs when writing config
   !> namelists. Key is left justified, align = a minimum of
   !> `key_align_width` characters to right of key start.
   integer, parameter :: key_align_width = 40

   !> Define abstract type to represent the config type,
   !> allows common code to be moved here
   type, abstract :: abstract_config_type
      logical, private :: initialised = .false. !< Have we initialised the instance?
      logical :: exist = .false. !< Does the related namelist exist in the target input file?
      logical, private :: requires_index = .false. !< Is this a numbered namelist with an index
      integer :: index = 0 !< Used to hold the specific index of numbered namelists
      character(len=CONFIG_MAX_NAME_LEN), private :: name = default_name !< The name of the namelist that we represent.
      logical :: skip_read = .false. !< Do we want to skip the read step in init?
      logical :: skip_broadcast = .false. !< Do we want to skip the broadcast step in init?
   contains
      procedure :: is_initialised => is_initialised_generic
      procedure :: init => init_generic
      procedure :: setup => setup_generic
      procedure :: write_namelist_header
      procedure :: get_name => get_name_generic
      procedure :: get_requires_index => get_requires_index_generic
      procedure, nopass :: write_namelist_footer
      procedure, private, nopass :: write_key_val_string
      procedure, private, nopass :: write_key_val_real
      procedure, private :: write_key_val_real_array
      procedure, private, nopass :: write_key_val_complex
      procedure, private :: write_key_val_complex_array
      procedure, private, nopass :: write_key_val_integer
      procedure, private :: write_key_val_integer_array
      procedure, private, nopass :: write_key_val_logical
      generic :: write_key_val => write_key_val_string, write_key_val_real, write_key_val_complex, &
         write_key_val_integer, write_key_val_logical, write_key_val_real_array, &
         write_key_val_complex_array, write_key_val_integer_array
      procedure(read_interface), deferred :: read
      procedure(write_interface), deferred :: write
      procedure(reset_interface), deferred :: reset
      procedure(broadcast_interface), deferred :: broadcast
      procedure(get_default_name_interface), deferred, nopass :: get_default_name
      procedure(get_default_requires_index_interface), deferred, nopass :: get_default_requires_index
   end type abstract_config_type

   interface
      subroutine read_interface(self)
         import
         class(abstract_config_type), intent(in out) :: self
      end subroutine read_interface

      subroutine write_interface(self, unit)
         import
         class(abstract_config_type), intent(in) :: self
         integer, intent(in), optional :: unit
      end subroutine write_interface

      subroutine reset_interface(self)
         import
         class(abstract_config_type), intent(in out) :: self
      end subroutine reset_interface

      subroutine broadcast_interface(self)
         import
         class(abstract_config_type), intent(in out) :: self
      end subroutine broadcast_interface

      function get_default_name_interface() result(default_name)
         import
         character(len=CONFIG_MAX_NAME_LEN) :: default_name
      end function get_default_name_interface

      function get_default_requires_index_interface() result(default_requires_index)
         logical :: default_requires_index
      end function get_default_requires_index_interface
   end interface
contains

   !> Is this instance initialised?
   function is_initialised_generic(self)
      class(abstract_config_type), intent(in) :: self
      logical :: is_initialised_generic
      is_initialised_generic = self%initialised
   end function is_initialised_generic

   !> Fully initialise the config object
   subroutine init_generic(self, name, requires_index, index)
      implicit none
      class(abstract_config_type), intent(inout) :: self
      character(len=*), intent(in), optional :: name
      logical, intent(in), optional :: requires_index
      integer, intent(in), optional :: index
      if (self%is_initialised()) return
      call self%setup(name, requires_index, index)
      if (.not. self%skip_read) call self%read()
      if (.not. self%skip_broadcast) call self%broadcast()
      self%initialised = .true.
   end subroutine init_generic

   !> Do some standard setup/checking
   subroutine setup_generic(self, name, requires_index, index)
      implicit none
      class(abstract_config_type), intent(inout) :: self
      character(len=*), intent(in), optional :: name
      logical, intent(in), optional :: requires_index
      integer, intent(in), optional :: index
      if (present(name)) then
         self%name = name
      else
         !Set the default name if not passed
         self%name = self%get_default_name()
      end if

      if (present(requires_index)) then
         self%requires_index = requires_index
      else
         !Set the default requires index if not passed
         self%requires_index = self%get_default_requires_index()
      end if

      if (present(index)) then
         self%index = index
      else
         !Set a default index number if required and not set
         if (self%get_requires_index()) self%index = 1
      end if

      self%initialised = .true.
   end subroutine setup_generic

   !> Returns the namelist name. Not very useful at the moment
   !> but may want to do more interesting things in the future
   function get_name_generic(self)
      implicit none
      class(abstract_config_type), intent(in) :: self
      character(len=CONFIG_MAX_NAME_LEN) :: get_name_generic
      if (self%is_initialised()) then
         get_name_generic = self%name
      else
         get_name_generic = self%get_default_name()
      end if
   end function get_name_generic

   !> Returns the requires_index value. Allows access whilst keeping
   !> the variable private
   function get_requires_index_generic(self)
      implicit none
      class(abstract_config_type), intent(in) :: self
      logical :: get_requires_index_generic
      if (self%is_initialised()) then
         get_requires_index_generic = self%requires_index
      else
         get_requires_index_generic = self%get_default_requires_index()
      end if
   end function get_requires_index_generic

   !> Write the namelist header for this instance
   subroutine write_namelist_header(self, unit)
      implicit none
      class(abstract_config_type), intent(in) :: self
      integer, intent(in) :: unit

      !Decide if we should include the index or not
      if (self%get_requires_index()) then
         write (unit, '("&",A,"_",I0)') trim(self%get_name()), self%index
      else
         write (unit, '("&",A)') trim(self%get_name())
      end if
   end subroutine write_namelist_header

   !> Write the namelist footer
   subroutine write_namelist_footer(unit)
      implicit none
      integer, intent(in) :: unit
      write (unit, '("/")')
      write (unit, '()')
   end subroutine write_namelist_footer

   !> Writes a {key,val} pair where the value is of type character
   subroutine write_key_val_string(key, val, unit)
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: val
      integer, intent(in) :: unit
      write (unit, '("  ",A," = "," ",A)') format_key(key), '"'//trim(val)//'"'
   end subroutine write_key_val_string

   !> Writes a {key,val} pair where the value is of type real
   subroutine write_key_val_real(key, val, unit)
      character(len=*), intent(in) :: key
      real, intent(in) :: val
      integer, intent(in) :: unit
      write (unit, '("  ",A," = ",e18.11)') format_key(key), val
   end subroutine write_key_val_real

   !> Writes a {key,val} pair where the value is of type real array
   subroutine write_key_val_real_array(self, key, val, unit)
      class(abstract_config_type), intent(in) :: self
      character(len=*), intent(in) :: key
      real, dimension(:), intent(in) :: val
      integer, intent(in) :: unit
      integer :: i
      character(len=12) :: subscript
      do i = 1, size(val)
         write (subscript, '("(",I0,")")') i
         call self%write_key_val(trim(key)//subscript, val(i), unit)
      end do
   end subroutine write_key_val_real_array

   !> Writes a {key,val} pair where the value is of type complex
   subroutine write_key_val_complex(key, val, unit)
      character(len=*), intent(in) :: key
      complex, intent(in) :: val
      integer, intent(in) :: unit
      write (unit, '("  ",A," = (",e18.11,", ",e18.11,")")') format_key(key), real(val), aimag(val)
   end subroutine write_key_val_complex

   !> Writes a {key,val} pair where the value is of type complex array
   subroutine write_key_val_complex_array(self, key, val, unit)
      class(abstract_config_type), intent(in) :: self
      character(len=*), intent(in) :: key
      complex, dimension(:), intent(in) :: val
      integer, intent(in) :: unit
      integer :: i
      character(len=12) :: subscript
      do i = 1, size(val)
         write (subscript, '("(",I0,")")') i
         call self%write_key_val(trim(key)//subscript, val(i), unit)
      end do
   end subroutine write_key_val_complex_array

   !> Writes a {key,val} pair where the value is of type integer
   subroutine write_key_val_integer(key, val, unit)
      character(len=*), intent(in) :: key
      integer, intent(in) :: val
      integer, intent(in) :: unit
      if (val < 0) then
         write (unit, '("  ",A," = ",I0)') format_key(key), val
      else
         write (unit, '("  ",A," = "," ",I0)') format_key(key), val
      end if
   end subroutine write_key_val_integer

   !> Writes a {key,val} pair where the value is of type integer array
   subroutine write_key_val_integer_array(self, key, val, unit)
      class(abstract_config_type), intent(in) :: self
      character(len=*), intent(in) :: key
      integer, dimension(:), intent(in) :: val
      integer, intent(in) :: unit
      integer :: i
      character(len=12) :: subscript
      do i = 1, size(val)
         write (subscript, '("(",I0,")")') i
         call self%write_key_val(trim(key)//subscript, val(i), unit)
      end do
   end subroutine write_key_val_integer_array

   !> Writes a {key,val} pair where the value is of type logical
   subroutine write_key_val_logical(key, val, unit)
      character(len=*), intent(in) :: key
      logical, intent(in) :: val
      integer, intent(in) :: unit
      write (unit, '("  ",A," = "," ",L1)') format_key(key), val
   end subroutine write_key_val_logical

   !> Takes a given key and formats it in a consistent style.
   !> Currently that style is left justified in a character variable
   !> of minimum length `key_align_width`.
   pure function format_key(key)
      implicit none
      character(len=:), allocatable :: format_key
      character(len=*), intent(in) :: key
      integer :: length
      length = max(key_align_width, len_trim(adjustl(key)))
      allocate (character(len=length)::format_key)
      write (format_key, '(A)') trim(adjustl(key))
   end function format_key

end module abstract_config
