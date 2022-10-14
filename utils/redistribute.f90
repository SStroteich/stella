! Modifications for optimised local copy in c_redist_22 and c_redist_32
! (and their inverse routines):
! (c) The Numerical Algorithms Group (NAG) Ltd, 2012
! on behalf of EPSRC for the HECToR project
module redistribute
!
! Redistribute distributed (integer, real, complex or logical)
! (1, 2, 3, or 4) dimensional arrays into two dimensional arrays with
! first index on local processor, and vice versa.
!
! The first operation is called 'gather' and the second is called 'scatter.'
!
! One can also do a 'fill' operation.  This consists of copying
! values from a (2, 3, or 4) dimensional array of
! (integer, real, complex, or logical ) values into
! another array with the same number of dimensions.
!
! One can also do a three index to four index redistribution for complex numbers.
!
   implicit none
   private

   public :: index_list_type, delete_list
   public :: redist_type, delete_redist
! TT>
   public :: report_map_property
   public :: gather_count, scatter_count, time_redist
! <TT

   public :: init_redist, gather, scatter
   public :: set_redist_character_type

   public :: parallel_scatter_complex
   public :: parallel_scatter_complex_commsplit


   interface gather
      module procedure c_redist_35, r_redist_35
   end interface

   interface scatter
      module procedure c_redist_35_inv, r_redist_35_inv
   end interface

! TT>

   integer :: gather_count = 0, scatter_count = 0
   real, save :: time_redist(2) = 0.
! <TT


   type :: index_map
      integer :: nn
      integer, dimension(:), pointer :: k => null()
      integer, dimension(:), pointer :: l => null()
      integer, dimension(:), pointer :: m => null()
      integer, dimension(:), pointer :: n => null()
      integer, dimension(:), pointer :: o => null()
   end type index_map

   ! TT: want to add map name, from_layout and to_layout
   type :: redist_type
      private
      integer, dimension(5) :: to_low, from_low, to_high, from_high
      type(index_map), dimension(:), pointer :: to => null()
      type(index_map), dimension(:), pointer :: from => null()
      complex, dimension(:), pointer :: complex_buff => null()
      real, dimension(:), pointer :: real_buff => null()
      integer, dimension(:), pointer :: integer_buff => null()
      logical, dimension(:), pointer :: logical_buff => null()
      character(len=3) :: redistname = ""
   end type redist_type

   type :: index_list_type
      integer, dimension(:), pointer :: first => null()
      integer, dimension(:), pointer :: second => null()
      integer, dimension(:), pointer :: third => null()
      integer, dimension(:), pointer :: fourth => null()
      integer, dimension(:), pointer :: fifth => null()
   end type index_list_type

contains

   subroutine set_redist_character_type(r, chartype)

      type(redist_type), intent(inout) :: r
      character(3), intent(in) :: chartype

      r%redistname = chartype

   end subroutine set_redist_character_type

   subroutine init_redist(r, char, to_low, to_high, to_list, &
                          from_low, from_high, from_list, ierr)

      use mp, only: iproc, nproc, proc0
      type(redist_type), intent(inout) :: r
      character(1), intent(in) :: char
      type(index_list_type), dimension(0:nproc - 1), intent(in) :: to_list, from_list
      integer, dimension(:), intent(in) :: from_low, to_high, from_high, to_low

      integer :: j, ip, n_to, n_from, buff_size
      integer, optional, intent(out) :: ierr

      allocate (r%to(0:nproc - 1), r%from(0:nproc - 1))

      if (present(ierr)) ierr = 0
      buff_size = 0

      do ip = 0, nproc - 1
         if (associated(to_list(ip)%first)) then
            n_to = size(to_list(ip)%first)
            r%to(ip)%nn = n_to

            allocate (r%to(ip)%k(n_to))
            allocate (r%to(ip)%l(n_to))

            r%to(ip)%k = to_list(ip)%first
            r%to(ip)%l = to_list(ip)%second

            if (associated(to_list(ip)%third)) then
               allocate (r%to(ip)%m(n_to))
               r%to(ip)%m = to_list(ip)%third
            end if

            if (associated(to_list(ip)%fourth)) then
               allocate (r%to(ip)%n(n_to))
               r%to(ip)%n = to_list(ip)%fourth
            end if

            if (associated(to_list(ip)%fifth)) then
               allocate (r%to(ip)%o(n_to))
               r%to(ip)%o = to_list(ip)%fifth
            end if

            if (ip /= iproc) buff_size = max(buff_size, n_to)
         else
            r%to(ip)%nn = 0
         end if
      end do

      do j = 1, size(from_low)
         r%from_low(j) = from_low(j)
      end do

      do j = 1, size(from_high)
         r%from_high(j) = from_high(j)
      end do

      do j = 1, size(to_high)
         r%to_high(j) = to_high(j)
      end do

      do j = 1, size(to_low)
         r%to_low(j) = to_low(j)
      end do

      do ip = 0, nproc - 1
         if (associated(from_list(ip)%first)) then

            n_from = size(from_list(ip)%first)
            r%from(ip)%nn = n_from

            allocate (r%from(ip)%k(n_from))
            allocate (r%from(ip)%l(n_from))

            r%from(ip)%k = from_list(ip)%first
            r%from(ip)%l = from_list(ip)%second

            if (associated(from_list(ip)%third)) then
               allocate (r%from(ip)%m(n_from))
               r%from(ip)%m = from_list(ip)%third
            end if

            if (associated(from_list(ip)%fourth)) then
               allocate (r%from(ip)%n(n_from))
               r%from(ip)%n = from_list(ip)%fourth
            end if

            if (associated(from_list(ip)%fifth)) then
               allocate (r%from(ip)%o(n_from))
               r%from(ip)%o = from_list(ip)%fifth
            end if

            if (ip /= iproc) buff_size = max(buff_size, n_from)
         else
            r%from(ip)%nn = 0
         end if
      end do

      select case (char)
      case ('c')
         if (buff_size > 0) allocate (r%complex_buff(buff_size))
      case ('r')
         if (buff_size > 0) allocate (r%real_buff(buff_size))
      case ('i')
         if (buff_size > 0) allocate (r%integer_buff(buff_size))
      case ('l')
         if (buff_size > 0) allocate (r%logical_buff(buff_size))
      case default
         if (proc0) then
            write (*, *) 'Type to be redistributed invalid.  Must stop.'
            write (*, *) char
         end if
         stop
      end select

   end subroutine init_redist


   subroutine delete_redist(r)

      use mp, only: nproc
      type(redist_type), intent(in out) :: r

      integer :: i

      if (associated(r%to)) then
         do i = 0, nproc - 1
            if (associated(r%to(i)%k)) deallocate (r%to(i)%k)
            if (associated(r%to(i)%l)) deallocate (r%to(i)%l)
            if (associated(r%to(i)%m)) deallocate (r%to(i)%m)
            if (associated(r%to(i)%n)) deallocate (r%to(i)%n)
         end do
         deallocate (r%to)
      end if

      if (associated(r%from)) then
         do i = 0, nproc - 1
            if (associated(r%from(i)%k)) deallocate (r%from(i)%k)
            if (associated(r%from(i)%l)) deallocate (r%from(i)%l)
            if (associated(r%from(i)%m)) deallocate (r%from(i)%m)
            if (associated(r%from(i)%n)) deallocate (r%from(i)%n)
         end do
         deallocate (r%from)
      end if

      if (associated(r%complex_buff)) deallocate (r%complex_buff)
      if (associated(r%real_buff)) deallocate (r%real_buff)
      if (associated(r%integer_buff)) deallocate (r%integer_buff)
      if (associated(r%logical_buff)) deallocate (r%logical_buff)

   end subroutine delete_redist

   subroutine delete_list(list)
      use mp, only: nproc
! TT> caused a problem on PGI compiler
!    type (index_list_type), dimension(0:) :: list
      type(index_list_type), dimension(0:nproc - 1), intent(inout) :: list
! <TT

      integer :: ip

      do ip = 0, nproc - 1
         if (associated(list(ip)%first)) deallocate (list(ip)%first)
         if (associated(list(ip)%second)) deallocate (list(ip)%second)
         if (associated(list(ip)%third)) deallocate (list(ip)%third)
         if (associated(list(ip)%fourth)) deallocate (list(ip)%fourth)
         if (associated(list(ip)%fifth)) deallocate (list(ip)%fifth)
      end do

   end subroutine delete_list

   subroutine r_redist_35(r, from_here, to_here)

      use mp, only: iproc, nproc, send, receive
      type(redist_type), intent(in out) :: r

      real, dimension(r%from_low(1):, &
                      r%from_low(2):, &
                      r%from_low(3):), intent(in) :: from_here

      real, dimension(r%to_low(1):, &
                      r%to_low(2):, &
                      r%to_low(3):, &
                      r%to_low(4):, &
                      r%to_low(5):), intent(in out) :: to_here

      integer :: i, idp, ipto, ipfrom, iadp

      ! redistribute from local processor to local processor
      do i = 1, r%from(iproc)%nn
         to_here(r%to(iproc)%k(i), &
                 r%to(iproc)%l(i), &
                 r%to(iproc)%m(i), &
                 r%to(iproc)%n(i), &
                 r%to(iproc)%o(i)) &
            = from_here(r%from(iproc)%k(i), &
                        r%from(iproc)%l(i), &
                        r%from(iproc)%m(i))
      end do

      ! redistribute to idpth next processor from idpth preceding processor
      ! or redistribute from idpth preceding processor to idpth next processor
      ! to avoid deadlocks
      do idp = 1, nproc - 1
         ipto = mod(iproc + idp, nproc)
         ipfrom = mod(iproc + nproc - idp, nproc)
         iadp = min(idp, nproc - idp)
         ! avoid deadlock AND ensure mostly parallel resolution
         if (mod(iproc / iadp, 2) == 0) then

            ! send to idpth next processor
            if (r%from(ipto)%nn > 0) then
               do i = 1, r%from(ipto)%nn
                  r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                             r%from(ipto)%l(i), &
                                             r%from(ipto)%m(i))
               end do
               call send(r%real_buff(1:r%from(ipto)%nn), ipto, idp)
            end if

            ! receive from idpth preceding processor
            if (r%to(ipfrom)%nn > 0) then
               call receive(r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%to(ipfrom)%nn
                  to_here(r%to(ipfrom)%k(i), &
                          r%to(ipfrom)%l(i), &
                          r%to(ipfrom)%m(i), &
                          r%to(ipfrom)%n(i), &
                          r%to(ipfrom)%o(i)) &
                     = r%real_buff(i)
               end do
            end if
         else
            ! receive from idpth preceding processor
            if (r%to(ipfrom)%nn > 0) then
               call receive(r%real_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%to(ipfrom)%nn
                  to_here(r%to(ipfrom)%k(i), &
                          r%to(ipfrom)%l(i), &
                          r%to(ipfrom)%m(i), &
                          r%to(ipfrom)%n(i), &
                          r%to(ipfrom)%o(i)) &
                     = r%real_buff(i)
               end do
            end if

            ! send to idpth next processor
            if (r%from(ipto)%nn > 0) then
               do i = 1, r%from(ipto)%nn
                  r%real_buff(i) = from_here(r%from(ipto)%k(i), &
                                             r%from(ipto)%l(i), &
                                             r%from(ipto)%m(i))
               end do
               call send(r%real_buff(1:r%from(ipto)%nn), ipto, idp)
            end if
         end if
      end do

   end subroutine r_redist_35

   subroutine r_redist_35_inv(r, from_here, to_here)

      use mp, only: iproc, nproc, send, receive
      type(redist_type), intent(in out) :: r

      real, dimension(r%to_low(1):, &
                      r%to_low(2):, &
                      r%to_low(3):, &
                      r%to_low(4):, &
                      r%to_low(5):), intent(in) :: from_here

      real, dimension(r%from_low(1):, &
                      r%from_low(2):, &
                      r%from_low(3):), intent(in out) :: to_here

      integer :: i, idp, ipto, ipfrom, iadp

      ! redistribute from local processor to local processor
      do i = 1, r%to(iproc)%nn
         to_here(r%from(iproc)%k(i), &
                 r%from(iproc)%l(i), &
                 r%from(iproc)%m(i)) &
            = from_here(r%to(iproc)%k(i), &
                        r%to(iproc)%l(i), &
                        r%to(iproc)%m(i), &
                        r%to(iproc)%n(i), &
                        r%to(iproc)%o(i))
      end do

      ! redistribute to idpth next processor from idpth preceding processor
      ! or redistribute from idpth preceding processor to idpth next processor
      ! to avoid deadlocks
      do idp = 1, nproc - 1
         ipto = mod(iproc + idp, nproc)
         ipfrom = mod(iproc + nproc - idp, nproc)
         iadp = min(idp, nproc - idp)
         ! avoid deadlock AND ensure mostly parallel resolution
         if (mod(iproc / iadp, 2) == 0) then

            ! send to idpth next processor
            if (r%to(ipto)%nn > 0) then
               do i = 1, r%to(ipto)%nn
                  r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                             r%to(ipto)%l(i), &
                                             r%to(ipto)%m(i), &
                                             r%to(ipto)%n(i), &
                                             r%to(ipto)%o(i))
               end do
               call send(r%real_buff(1:r%to(ipto)%nn), ipto, idp)
            end if

            ! receive from idpth preceding processor
            if (r%from(ipfrom)%nn > 0) then
               call receive(r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%from(ipfrom)%nn
                  to_here(r%from(ipfrom)%k(i), &
                          r%from(ipfrom)%l(i), &
                          r%from(ipfrom)%m(i)) &
                     = r%real_buff(i)
               end do
            end if
         else
            ! receive from idpth preceding processor
            if (r%from(ipfrom)%nn > 0) then
               call receive(r%real_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%from(ipfrom)%nn
                  to_here(r%from(ipfrom)%k(i), &
                          r%from(ipfrom)%l(i), &
                          r%from(ipfrom)%m(i)) &
                     = r%real_buff(i)
               end do
            end if

            ! send to idpth next processor
            if (r%to(ipto)%nn > 0) then
               do i = 1, r%to(ipto)%nn
                  r%real_buff(i) = from_here(r%to(ipto)%k(i), &
                                             r%to(ipto)%l(i), &
                                             r%to(ipto)%m(i), &
                                             r%to(ipto)%n(i), &
                                             r%to(ipto)%o(i))
               end do
               call send(r%real_buff(1:r%to(ipto)%nn), ipto, idp)
            end if

         end if
      end do

   end subroutine r_redist_35_inv

   subroutine c_redist_35(r, from_here, to_here)

      use mp, only: iproc, nproc, send, receive
      type(redist_type), intent(in out) :: r

      complex, dimension(r%from_low(1):, &
                         r%from_low(2):, &
                         r%from_low(3):), intent(in) :: from_here

      complex, dimension(r%to_low(1):, &
                         r%to_low(2):, &
                         r%to_low(3):, &
                         r%to_low(4):, &
                         r%to_low(5):), intent(in out) :: to_here

      integer :: i, idp, ipto, ipfrom, iadp

      ! redistribute from local processor to local processor
      do i = 1, r%from(iproc)%nn
         to_here(r%to(iproc)%k(i), &
                 r%to(iproc)%l(i), &
                 r%to(iproc)%m(i), &
                 r%to(iproc)%n(i), &
                 r%to(iproc)%o(i)) &
            = from_here(r%from(iproc)%k(i), &
                        r%from(iproc)%l(i), &
                        r%from(iproc)%m(i))
      end do

      ! redistribute to idpth next processor from idpth preceding processor
      ! or redistribute from idpth preceding processor to idpth next processor
      ! to avoid deadlocks
      do idp = 1, nproc - 1
         ipto = mod(iproc + idp, nproc)
         ipfrom = mod(iproc + nproc - idp, nproc)
         iadp = min(idp, nproc - idp)
         ! avoid deadlock AND ensure mostly parallel resolution
         if (mod(iproc / iadp, 2) == 0) then

            ! send to idpth next processor
            if (r%from(ipto)%nn > 0) then
               do i = 1, r%from(ipto)%nn
                  r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                                r%from(ipto)%l(i), &
                                                r%from(ipto)%m(i))
               end do
               call send(r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
            end if

            ! receive from idpth preceding processor
            if (r%to(ipfrom)%nn > 0) then
               call receive(r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%to(ipfrom)%nn
                  to_here(r%to(ipfrom)%k(i), &
                          r%to(ipfrom)%l(i), &
                          r%to(ipfrom)%m(i), &
                          r%to(ipfrom)%n(i), &
                          r%to(ipfrom)%o(i)) &
                     = r%complex_buff(i)
               end do
            end if
         else
            ! receive from idpth preceding processor
            if (r%to(ipfrom)%nn > 0) then
               call receive(r%complex_buff(1:r%to(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%to(ipfrom)%nn
                  to_here(r%to(ipfrom)%k(i), &
                          r%to(ipfrom)%l(i), &
                          r%to(ipfrom)%m(i), &
                          r%to(ipfrom)%n(i), &
                          r%to(ipfrom)%o(i)) &
                     = r%complex_buff(i)
               end do
            end if

            ! send to idpth next processor
            if (r%from(ipto)%nn > 0) then
               do i = 1, r%from(ipto)%nn
                  r%complex_buff(i) = from_here(r%from(ipto)%k(i), &
                                                r%from(ipto)%l(i), &
                                                r%from(ipto)%m(i))
               end do
               call send(r%complex_buff(1:r%from(ipto)%nn), ipto, idp)
            end if
         end if
      end do

   end subroutine c_redist_35

   subroutine c_redist_35_inv(r, from_here, to_here)

      use mp, only: iproc, nproc, send, receive
      type(redist_type), intent(in out) :: r

      complex, dimension(r%to_low(1):, &
                         r%to_low(2):, &
                         r%to_low(3):, &
                         r%to_low(4):, &
                         r%to_low(5):), intent(in) :: from_here

      complex, dimension(r%from_low(1):, &
                         r%from_low(2):, &
                         r%from_low(3):), intent(in out) :: to_here

      integer :: i, idp, ipto, ipfrom, iadp

      ! redistribute from local processor to local processor
      do i = 1, r%to(iproc)%nn
         to_here(r%from(iproc)%k(i), &
                 r%from(iproc)%l(i), &
                 r%from(iproc)%m(i)) &
            = from_here(r%to(iproc)%k(i), &
                        r%to(iproc)%l(i), &
                        r%to(iproc)%m(i), &
                        r%to(iproc)%n(i), &
                        r%to(iproc)%o(i))
      end do

      ! redistribute to idpth next processor from idpth preceding processor
      ! or redistribute from idpth preceding processor to idpth next processor
      ! to avoid deadlocks
      do idp = 1, nproc - 1
         ipto = mod(iproc + idp, nproc)
         ipfrom = mod(iproc + nproc - idp, nproc)
         iadp = min(idp, nproc - idp)
         ! avoid deadlock AND ensure mostly parallel resolution
         if (mod(iproc / iadp, 2) == 0) then

            ! send to idpth next processor
            if (r%to(ipto)%nn > 0) then
               do i = 1, r%to(ipto)%nn
                  r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                                r%to(ipto)%l(i), &
                                                r%to(ipto)%m(i), &
                                                r%to(ipto)%n(i), &
                                                r%to(ipto)%o(i))
               end do
               call send(r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
            end if

            ! receive from idpth preceding processor
            if (r%from(ipfrom)%nn > 0) then
               call receive(r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%from(ipfrom)%nn
                  to_here(r%from(ipfrom)%k(i), &
                          r%from(ipfrom)%l(i), &
                          r%from(ipfrom)%m(i)) &
                     = r%complex_buff(i)
               end do
            end if
         else
            ! receive from idpth preceding processor
            if (r%from(ipfrom)%nn > 0) then
               call receive(r%complex_buff(1:r%from(ipfrom)%nn), ipfrom, idp)
               do i = 1, r%from(ipfrom)%nn
                  to_here(r%from(ipfrom)%k(i), &
                          r%from(ipfrom)%l(i), &
                          r%from(ipfrom)%m(i)) &
                     = r%complex_buff(i)
               end do
            end if

            ! send to idpth next processor
            if (r%to(ipto)%nn > 0) then
               do i = 1, r%to(ipto)%nn
                  r%complex_buff(i) = from_here(r%to(ipto)%k(i), &
                                                r%to(ipto)%l(i), &
                                                r%to(ipto)%m(i), &
                                                r%to(ipto)%n(i), &
                                                r%to(ipto)%o(i))
               end do
               call send(r%complex_buff(1:r%to(ipto)%nn), ipto, idp)
            end if

         end if
      end do

   end subroutine c_redist_35_inv

  subroutine get_master_and_node_id( master_size, node_size, world_to_local, world_to_local_master, &
                                     world_addr, master_id, node_id ) 
    !integer, dimension(node_size) :: local_to_world
    integer :: master_size, node_size
    integer, dimension(master_size*node_size) :: world_to_local_master
    integer, dimension(master_size*node_size) :: world_to_local
    integer :: world_addr 
    integer, intent(out) :: master_id, node_id

    master_id = world_to_local_master(world_addr+1)
    node_id = world_to_local(world_addr+1)

  end subroutine get_master_and_node_id

  subroutine parallel_scatter_complex_commsplit_impl( r, from_here, to_here, world_rank, world_size, world_comm, &
                                                   master_rank, master_size, master_comm, node_rank, node_size, node_comm )
     use mpi
     use mp, only: iproc, send, receive, waitall, waitany, barrier, mp_comm, mpicmplx, &
                   comm_send_complex_array, comm_receive_complex_array
     type(redist_type), intent(in out) :: r

     complex, dimension(r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):), intent(in) :: from_here

     complex, dimension(r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent(in out) :: to_here

     integer :: i,j, idp, ipto, ipfrom, iadp

     !complex, dimension(master_size*node_size*10000) :: send_buff
     !complex, dimension(node_size*master_size*node_size*10000) :: gather_buff

     complex, dimension(master_size*node_size*10000) :: send_buff
     complex, dimension(master_size*node_size*10000) :: receive_buff
     !complex, dimension(node_size*(master_size*node_size*10000)) :: gather_buff ! TODO only on 0
     complex, dimension(:), allocatable :: gather_buff ! TODO only on 0
     complex, dimension(:), allocatable :: gather_gather_buff ! TODO only on 0
     integer, dimension(world_size) :: send_requests
     integer, dimension(world_size) :: receive_requests
     integer, dimension(MPI_STATUS_SIZE) :: statuses
     integer :: idx
     integer :: ierror

     integer :: base 
     integer :: offset
     integer :: info
     integer :: node_comm
     integer :: node_rank, node_size
     integer :: master_comm
     integer :: master_rank, master_size
     integer :: world_rank, world_size, world_comm

     integer, dimension(node_size) :: local_to_world
     integer, dimension(node_size) :: local_to_master
     integer, dimension(master_size*node_size) :: world_to_local_master
     integer, dimension(master_size*node_size) :: world_to_local

     integer :: master_id, node_id
     integer :: master_data_chunk_size
     integer :: gathered_master_data_chunk_size

     integer :: print_master_from_id
     integer :: print_node_from_id
     integer :: print_master_to_id
     integer :: print_node_to_id
     integer :: print_world_from_id
     integer :: print_world_to_id

     print_master_from_id = 0
     print_node_from_id = 1

     print_master_to_id = 1
     print_node_to_id = 1

     print_world_from_id = 1
     print_world_to_id = 65
   

     !allocate( gather_buff(node_size*(master_size*node_size*10000) ) )
     if ( node_rank == 0 ) then
       write(*,*) "allocating ", node_size*(master_size*node_size*10000)
       allocate( gather_buff(node_size*(master_size*node_size*10000) ) )
       allocate( gather_gather_buff(master_size*(node_size*node_size*10000) ) )
     endif

     ! TODO create look up table from world to node local rank
     local_to_world(node_rank+1) = world_rank
     call mpi_allgather( local_to_world(node_rank+1), &
                        1, &
                        MPI_INTEGER, &
                        local_to_world, &
                        1, &
                        MPI_INTEGER, &
                        node_comm, ierror )

     call mpi_bcast( master_rank, &
                     1, &
                     MPI_INTEGER, &
                     0, &
                     node_comm, & 
                     ierror ) 

     

     ! TODO create look up table from world to node local rank
     world_to_local_master(world_rank+1) = master_rank
     call mpi_allgather( world_to_local_master(world_rank+1), &
                        1, &
                        MPI_INTEGER, &
                        world_to_local_master, &
                        1, &
                        MPI_INTEGER, &
                        world_comm, ierror )

     ! TODO create look up table from world to node local rank
     world_to_local(world_rank+1) = node_rank
     call mpi_allgather( world_to_local(world_rank+1), &
                        1, &
                        MPI_INTEGER, &
                        world_to_local, &
                        1, &
                        MPI_INTEGER, &
                        world_comm, ierror )


      call barrier

    if ( world_rank == 0 ) then
      do idp = 0, world_size-1
       call get_master_and_node_id( master_size, node_size, world_to_local, &
                                    world_to_local_master, idp, master_id, node_id)

        write(*,*) "world_rank:", idp, &
                   "local_master:", master_id, & 
                   "node_rank:", node_id
      enddo
    endif
    ! if ( world_rank == 0 ) then
    !   do i = 0, world_size-1
    !     write(*,*) "world_rank: ", i, " local_master ", world_to_local_master(i+1)
    !   enddo
    ! endif

    ! TODO create look up table

     

     ! redistribute from local processor to local processor
     do i = 1, r%to(iproc)%nn
        to_here(r%from(iproc)%k(i), &
                r%from(iproc)%l(i), &
                r%from(iproc)%m(i)) &
           = from_here(r%to(iproc)%k(i), &
                       r%to(iproc)%l(i), &
                       r%to(iproc)%m(i), &
                       r%to(iproc)%n(i), &
                       r%to(iproc)%o(i))
     end do
      !if ( iproc == 0 ) write(*,*) "running redistribute loop"
      ! redistribute to idpth next processor from idpth preceding processor
      ! or redistribute from idpth preceding processor to idpth next processor
      ! to avoid deadlocks

      call barrier
      do idp = 1, world_size - 1
         ipto = mod(iproc + idp, world_size)
         ipfrom = mod(iproc + world_size - idp, world_size)
         iadp = min(idp, world_size - idp)
         ! send to idpth next processor
         if (r%to(ipto)%nn > 0) then
            call get_master_and_node_id( master_size, node_size, world_to_local, &
                                         world_to_local_master, ipto, master_id, node_id)
            base = master_id*node_size*10000+node_id*10000
            offset = base+1+r%to(ipto)%nn
            !if ( base .gt. size(send_buff) .or. base+1 < 1) then
            !  write(*,*) "master_id ", master_id, " node_id", node_id
            !endif
            do i = 1, r%to(ipto)%nn
               send_buff(base+i) = from_here(r%to(ipto)%k(i), &
                                             r%to(ipto)%l(i), &
                                             r%to(ipto)%m(i), &
                                             r%to(ipto)%n(i), &
                                             r%to(ipto)%o(i))
            end do
            if ( world_rank == print_world_from_id .and. ipto == print_world_to_id ) then
              write(*,*) "master: ", master_id, "node_id", node_id
              write(*,*) "send: ", world_rank, " ", node_rank, " " , send_buff(base+1:base+2)
            endif
         end if
         
      end do


      call barrier
      if ( world_rank == 0 ) then
        write(*,*) "send data prepared"
      endif
      
      master_data_chunk_size=node_size*10000 ! data to transmit to one master process
      gathered_master_data_chunk_size= master_data_chunk_size*node_size ! data gathered on one master from all local processes
      ! multiple gathers
      do i = 0, master_size - 1
        call mpi_gather( send_buff(i*master_data_chunk_size+1:(i+1)*master_data_chunk_size), &
                         master_data_chunk_size, &
                         mpicmplx, &
                         gather_buff(i*gathered_master_data_chunk_size+1:(i+1)*gathered_master_data_chunk_size), & 
                         master_data_chunk_size, &
                         mpicmplx, &
                         0, &
                         node_comm, &
                         ierror ) 

      enddo

      call barrier
      if (world_rank ==0) write(*,*) "gathered all data to node masters"

      if ( world_rank == 0 ) then 
      write(*,*) "data local: ", gather_buff(print_master_to_id*gathered_master_data_chunk_size+ &
                                 print_node_from_id*master_data_chunk_size+ &
                                 print_node_to_id*10000+1 &
                                 :&
                                 print_master_to_id*gathered_master_data_chunk_size+&
                                 print_node_from_id*master_data_chunk_size+&
                                 print_node_to_id*10000+2)
      endif

      if ( node_rank == 0 ) then
        ! TODO complete the target buffer code
        do i = 0, master_size - 1
          call mpi_gather( gather_buff(i*gathered_master_data_chunk_size+1:(i+1)*gathered_master_data_chunk_size), &
                           gathered_master_data_chunk_size, &
                           mpicmplx, &
                           gather_gather_buff, &
                           gathered_master_data_chunk_size, &
                           mpicmplx, &
                           i, &
                           master_comm, &
                           ierror )
        enddo

      endif
      if ( world_rank == 64 ) then 
      write(*,*) "data target master: ", gather_gather_buff(print_master_from_id*gathered_master_data_chunk_size+&
                                         print_node_from_id*master_data_chunk_size+&
                                         print_node_to_id*10000+1 &
                                         :&
                                         print_master_from_id*gathered_master_data_chunk_size+&
                                         print_node_from_id*master_data_chunk_size+&
                                         print_node_to_id*10000+2)
      endif

      call barrier
      if (world_rank ==0) write(*,*) "copied all data to masters"

      do i = 0, master_size - 1
        do j = 0, node_size - 1
          call mpi_scatter( gather_gather_buff(i*gathered_master_data_chunk_size+j*master_data_chunk_size+1 &
                           :&
                           i*gathered_master_data_chunk_size+(j+1)*master_data_chunk_size), &
                           10000, &
                           mpicmplx, &
                           receive_buff(i*master_data_chunk_size+j*10000+1:i*master_data_chunk_size+(j+1)*10000), &
                           10000, &
                           mpicmplx, &
                           0, &
                           node_comm, &
                           ierror )
        enddo
      enddo
      call barrier
      if (world_rank ==0) write(*,*) "copied all data to local processes"
      if ( world_rank == 65 ) then 
        write(*,*) "data target 65: ", receive_buff(print_master_from_id*master_data_chunk_size+print_node_from_id*10000+1:&
                                                    print_master_from_id*master_data_chunk_size+print_node_from_id*10000+2)
      endif

      do idp = 1, world_size - 1
         ipfrom = mod(iproc + world_size - idp, world_size)
         if (r%from(ipfrom)%nn > 0) then
            call get_master_and_node_id( master_size, node_size, world_to_local, &
                                         world_to_local_master, ipfrom, master_id, node_id)
            base = master_id*node_size*10000+node_id*10000
            do i = 1, r%from(ipfrom)%nn
               to_here(r%from(ipfrom)%k(i), &
                       r%from(ipfrom)%l(i), &
                       r%from(ipfrom)%m(i)) &
                  = receive_buff(base+i)

            end do
         end if
      end do
      call barrier
      if (world_rank ==0) write(*,*) "mapped into local arrays"

      !! TODO scatter to all node members
      !do i = 0, node_size - 1
      !  call mpi_scatter( gather_gather_buff,
      !enddo
     
      !do i = 0, master_size - 1

      !
      !if ( world_rank == 0 ) then
      !  master_id = 0
      !  node_id = 0
      !  base = 1*size(send_buff)+master_id*(node_size*10000)+node_id*(10000)
      !  write(*,*) "recv: ", world_rank, " ", node_rank, " ", gather_buff(base+1:base+2)
      !endif

      !call mpi_gather( gather_buff(master_rank*(master_size*node_size*10000)+1:(master_rank+1)*(master_size*node_size*10000), (master_size*node_size*10000), mpicmplx, gather_buff, size(send_buff), mpicmplx, 0, node_comm, ierror) 

      !call barrier

      !call exit(1)


  

  end subroutine parallel_scatter_complex_commsplit_impl


  subroutine parallel_scatter_complex_commsplit(r, from_here, to_here)
     use mpi
     use mp, only: iproc, nproc, send, receive, waitall, waitany, barrier, mp_comm, mpicmplx, &
                   comm_send_complex_array, comm_receive_complex_array
     type(redist_type), intent(in out) :: r

     complex, dimension(r%to_low(1):, &
                        r%to_low(2):, &
                        r%to_low(3):, &
                        r%to_low(4):, &
                        r%to_low(5):), intent(in) :: from_here

     complex, dimension(r%from_low(1):, &
                        r%from_low(2):, &
                        r%from_low(3):), intent(in out) :: to_here

     integer :: ierror
     integer :: world_rank, world_size, world_comm 
     integer :: node_comm
     integer :: node_rank, node_size
     integer :: master_comm
     integer :: master_rank, master_size


     world_comm = mp_comm
     world_rank = iproc
     world_size = nproc

     call mpi_comm_split_type(mp_comm, OMPI_COMM_TYPE_SOCKET, world_rank, MPI_INFO_NULL, node_comm, ierror)
     call mpi_comm_size(node_comm, node_size, ierror)
     call mpi_comm_rank(node_comm, node_rank, ierror)

     call mpi_comm_split(mp_comm, node_rank, world_rank, master_comm, ierror)
     call mpi_comm_size(master_comm, master_size, ierror)
     call mpi_comm_rank(master_comm, master_rank, ierror)

     call parallel_scatter_complex_commsplit_impl( r, from_here, to_here, world_rank, world_size, world_comm, &
                                                   master_rank, master_size, master_comm, node_rank, node_size, node_comm )
 
  end subroutine parallel_scatter_complex_commsplit


   subroutine parallel_scatter_complex(r, from_here, to_here)

      use mpi
      use mp, only: iproc, nproc, send, receive, waitall, waitany, barrier, mp_comm, mpicmplx
      type(redist_type), intent(in out) :: r

      complex, dimension(r%to_low(1):, &
                         r%to_low(2):, &
                         r%to_low(3):, &
                         r%to_low(4):, &
                         r%to_low(5):), intent(in) :: from_here

      complex, dimension(r%from_low(1):, &
                         r%from_low(2):, &
                         r%from_low(3):), intent(in out) :: to_here

      integer :: i, idp, ipto, ipfrom, iadp

      complex, dimension(nproc*10000) :: send_buff
      complex, dimension(nproc*10000) :: receive_buff
      integer, dimension(nproc) :: send_requests
      integer, dimension(nproc) :: receive_requests
      integer :: send_request_idx 
      integer :: receive_request_idx 
      integer, dimension(MPI_STATUS_SIZE) :: statuses
      integer :: idx
      integer :: ierror


      integer :: base 
      integer :: offset

      ! redistribute from local processor to local processor

!      if (iproc == 0) open(unit=70, file = "ipfrom_elements.txt")

!      if (iproc == 0) open(unit=69, file = "ipto_elements.txt")
      do i = 1, r%to(iproc)%nn
         to_here(r%from(iproc)%k(i), &
                 r%from(iproc)%l(i), &
                 r%from(iproc)%m(i)) &
            = from_here(r%to(iproc)%k(i), &
                        r%to(iproc)%l(i), &
                        r%to(iproc)%m(i), &
                        r%to(iproc)%n(i), &
                        r%to(iproc)%o(i))
!                  if (iproc == 0) write(70,*) "ipfrom: ", iproc," i: ", i," k: ", r%from(iproc)%k(i)," l: ", &
!                                               r%from(iproc)%l(i), " m: ",r%from(iproc)%m(i)
!                  if (iproc == 0) write(69,*) "ipto: ", iproc," i: ", i," k: ", r%to(iproc)%k(i)," l: ", r%to(iproc)%l(i), &
!                                              " m: ",r%to(iproc)%m(i)," n: ", r%to(iproc)%n(i)," o: ", r%to(iproc)%o(i)



      end do

      send_request_idx = 0
      receive_request_idx = 0
      !if ( iproc == 0 ) write(*,*) "running redistribute loop"
      ! redistribute to idpth next processor from idpth preceding processor
      ! or redistribute from idpth preceding processor to idpth next processor
      ! to avoid deadlocks

      do idp = 1, nproc - 1
         ipto = mod(iproc + idp, nproc)
         ipfrom = mod(iproc + nproc - idp, nproc)
         iadp = min(idp, nproc - idp)
         ! send to idpth next processor
         if (r%to(ipto)%nn > 0) then
            base = ipto*nproc
            offset = base+1+r%to(ipto)%nn
            do i = 1, r%to(ipto)%nn
               send_buff(base+i) = from_here(r%to(ipto)%k(i), &
                                             r%to(ipto)%l(i), &
                                             r%to(ipto)%m(i), &
                                             r%to(ipto)%n(i), &
                                             r%to(ipto)%o(i))

!                  if (iproc == 0) write(69,*) "ipto: ", ipto," i: ", i," k: ", r%to(ipto)%k(i)," l: ", r%to(ipto)%l(i), &
!                                              " m: ",r%to(ipto)%m(i)," n: ", r%to(ipto)%n(i)," o: ", r%to(ipto)%o(i)
            end do

            ! TODO send to node local rank 0
            

            send_request_idx = send_request_idx + 1
            call send(send_buff(base+1:offset), ipto, iproc*nproc+ipto, send_requests(send_request_idx))
         end if

         ! receive from idpth preceding processor
         if (r%from(ipfrom)%nn > 0) then
            base = ipfrom*nproc
            offset = base+1+r%from(ipfrom)%nn

            ! TODO go through multi block receives and write to corresponding buffer

            receive_request_idx = receive_request_idx + 1
            call receive(receive_buff(base+1:offset), ipfrom, &
                 ipfrom*nproc+iproc, receive_requests(receive_request_idx) )
         end if
         
      end do
!      if (iproc == 0) close(69)

      !if ( iproc == 0 ) write(*,*) "calling waitall"
      !if ( iproc == 0 ) write(*,*) "calling waitall send with ", send_request_idx
      !if ( iproc == 0 ) write(*,*) "calling waitall receive with ", receive_request_idx
      if( send_request_idx > 0) call waitall( send_request_idx, send_requests ) 
      if( receive_request_idx > 0) call waitall( receive_request_idx, receive_requests ) 


      !write(*,*) "status: ", statuses(MPI_SOURCE), " ", statuses(MPI_TAG), " ", statuses(MPI_ERROR)

      do idp = 1, nproc - 1
         ipfrom = mod(iproc + nproc - idp, nproc)
         base = ipfrom*nproc
         offset = base+1+r%from(ipfrom)%nn
         if (r%from(ipfrom)%nn > 0) then
            do i = 1, r%from(ipfrom)%nn
               to_here(r%from(ipfrom)%k(i), &
                       r%from(ipfrom)%l(i), &
                       r%from(ipfrom)%m(i)) &
                  = receive_buff(base+i)
!                  if (iproc == 0) write(70,*) "ipfrom: ", ipfrom," i: ", i," k: ", r%from(ipfrom)%k(i)," l: ", &
!                                               r%from(ipfrom)%l(i), " m: ",r%from(ipfrom)%m(i)

            end do

            !if( ipfrom*nproc+iproc == 120 ) write(*,*) "merge step tag: ", ipfrom*nproc+iproc, "buff: ", receive_buff(base+1)
         end if
      enddo
!      if (iproc == 0) close(70)

      call barrier
      !if ( iproc == 0 ) write(*,*) "done scatter "
   end subroutine parallel_scatter_complex

   subroutine report_map_property(r)

      use mp, only: iproc, nproc, proc0, sum_reduce, max_reduce
      type(redist_type), intent(in) :: r
      type :: redist_prp
         integer :: local_max, local_total
         integer :: comm_max, comm_total
         integer :: elm_max, elm_total
      end type redist_prp
      type(redist_prp) :: prp
      integer :: ip, rank_from, rank_to
      integer, dimension(:), allocatable :: lbd_from, lbd_to

      prp%comm_max = 0
      prp%comm_total = 0
      prp%elm_total = 0

      do ip = 0, nproc - 1
         if (ip == iproc) then
            prp%local_total = r%to(ip)%nn
         else
            if (r%to(ip)%nn > 0) then
               prp%comm_total = prp%comm_total + 1
               prp%elm_total = prp%elm_total + r%to(ip)%nn
            end if
         end if
      end do
      prp%local_max = prp%local_total
      prp%comm_max = prp%comm_total
      prp%elm_max = prp%elm_total

      call max_reduce(prp%local_max, 0)
      call sum_reduce(prp%local_total, 0)
      call max_reduce(prp%comm_max, 0)
      call sum_reduce(prp%comm_total, 0)
      call max_reduce(prp%elm_max, 0)
      call sum_reduce(prp%elm_total, 0)

      if (proc0) then
         rank_from = 1
         if (associated(r%from(0)%l)) rank_from = 2
         if (associated(r%from(0)%m)) rank_from = 3
         if (associated(r%from(0)%n)) rank_from = 4
         rank_to = 1
         if (associated(r%to(0)%l)) rank_to = 2
         if (associated(r%to(0)%m)) rank_to = 3
         if (associated(r%to(0)%n)) rank_to = 4
         allocate (lbd_from(rank_from), lbd_to(rank_to))
         lbd_from = r%from_low(1:rank_from)
         lbd_to = r%to_low(1:rank_to)
         print '(a,i2,a,i2)', 'From rank', rank_from, ' to rank', rank_to
         print '(a,t20,4i10)', 'From lbound (proc0)', r%from_low(1:rank_from)
         print '(a,t20,4i10)', 'To lbound (proc0)', r%to_low(1:rank_to)
         print '(a,t49,a,t64,a)', '--- Redistribution statistics ---', &
            'max', 'avg'
         print '(a,t40,i12,t55,f15.2)', 'Number of local move elements', &
            prp%local_max, real(prp%local_total) / real(nproc)
         print '(a,t40,i12,t60,f10.2)', &
            'Number of inter-processor communications', &
            prp%comm_max, real(prp%comm_total) / real(nproc)
         print '(a,t40,i12,t55,f15.2)', &
            'Number of inter-processor move elements', &
            prp%elm_max, real(prp%elm_total) / real(nproc)
         print *
      end if

   end subroutine report_map_property

end module redistribute
