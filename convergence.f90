!> Routines for estimating convergence of growth rates and frequencies for linear runs

module convergence

   use mpi

   implicit none
   public init_convergence, finish_convergence
   public init_overflow, finish_overflow
   public testing_convergence
   public rescale_g_and_phi
   public convergence_switch
   public overflow_switch

   private

   real :: omega_prec
   real :: overflow_limit
   real :: underflow_limit

   integer :: window_size

   logical :: convergence_switch
   logical :: overflow_switch

   !> Needed for testing convergence of growth rates and frequencies
   complex, dimension(:, :, :), allocatable :: omega_window
   complex, dimension(:, :, :), allocatable :: omega_avg_window

   real, dimension(:, :, :, :, :), allocatable :: g_abs

contains

   subroutine init_convergence

      use mp, only: broadcast, proc0

      !> Read the namelist "convergence_knobs" in the input file
      call read_parameters_convergence
      call allocate_arrays_convergence

      call broadcast(convergence_switch)
      call broadcast(omega_prec)
   end subroutine init_convergence

   subroutine init_overflow

      use mp, only: broadcast, proc0

      !> Read the namelist "overflow_knob" in the input file
      call read_parameters_overflow
      call allocate_arrays_overflow

      call broadcast(overflow_switch)
      call broadcast(overflow_limit)
      call broadcast(underflow_limit)
   end subroutine init_overflow

   !> Allocate the module-level arrays
   subroutine allocate_arrays_convergence

      use kt_grids, only: nakx, naky
      if (.not. allocated(omega_window)) then
         if (convergence_switch) then
            allocate (omega_window(window_size, naky, nakx))
            allocate (omega_avg_window(window_size, naky, nakx))
            omega_window = 0.
            omega_avg_window = 0.
         else
            allocate (omega_window(1, 1, 1))
            allocate (omega_avg_window(1, 1, 1))
            window_size = 1
         end if
      end if

   end subroutine allocate_arrays_convergence

   subroutine allocate_arrays_overflow

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      if (.not. allocated(g_abs)) then
         if (overflow_switch) then
            allocate (g_abs(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            g_abs = 0
         else
            allocate (g_abs(1, 1, 1, 1, 1))
         end if
      end if

   end subroutine allocate_arrays_overflow

   !> Read the module parameters

   subroutine read_parameters_convergence

      use mp, only: proc0
      use file_utils, only: input_unit_exist

      implicit none

      logical :: exist
      integer :: in_file

      namelist /convergence_knobs/ convergence_switch, &
         omega_prec, window_size

      if (proc0) then
         convergence_switch = .false.
         omega_prec = 10e-6
         window_size = 200
         in_file = input_unit_exist("convergence_knobs", exist)
         if (exist) read (unit=in_file, nml=convergence_knobs)

      end if

   end subroutine read_parameters_convergence

   subroutine read_parameters_overflow

      use mp, only: proc0
      use file_utils, only: input_unit_exist

      implicit none

      logical :: exist
      integer :: in_file

      namelist /overflow_knobs/ overflow_switch, &
         overflow_limit, underflow_limit

      if (proc0) then
         overflow_switch = .false.
         overflow_limit = 1e30
         underflow_limit = 1e-8

         in_file = input_unit_exist("overflow_knobs", exist)
         if (exist) read (unit=in_file, nml=overflow_knobs)

      end if

   end subroutine read_parameters_overflow

   !> Subroutine testing if maximum standard deviation of omega calculated for the fixed window_size is smaller than omega_prec

   subroutine testing_convergence(istep, conv_exit)
      use mp, only: proc0, broadcast
      use constants, only: zi
      use fields_arrays, only: phi, apar
      use fields_arrays, only: phi_old
      use stella_time, only: code_dt
      use kt_grids, only: naky, nakx

      use volume_averages, only: fieldline_average

      implicit none

      !> The current timestep
      integer, intent(in) :: istep

      logical, intent(in out) :: conv_exit

      real :: zero
      real :: maxstd
      logical :: reached_convergence

      complex, dimension(:, :), allocatable :: phiavg, phioldavg

      complex, dimension(:, :), allocatable :: omega_avg1
      complex, dimension(:, :), allocatable :: omega_avg2
      complex, dimension(:, :), allocatable :: omega_square

      real, dimension(:, :), allocatable :: omega_std
      real, dimension(:, :), allocatable :: omega_avg_squared
      integer, dimension(:), allocatable ::  posmax
      integer :: posky
      integer :: poskx

      ! calculation of omega requires computation of omega more
      ! frequently than every nwrite time steps
      if (proc0) then
         zero = 100.*epsilon(0.)
         if (istep > 0) then
            allocate (phiavg(naky, nakx))
            allocate (phioldavg(naky, nakx))
            call fieldline_average(phi, phiavg)
            call fieldline_average(phi_old, phioldavg)
            where (abs(phiavg) < zero .or. abs(phioldavg) < zero)
               omega_window(mod(istep, window_size) + 1, :, :) = 0.0
            elsewhere
               omega_window(mod(istep, window_size) + 1, :, :) = log(phiavg / phioldavg) * zi / code_dt
            end where
            deallocate (phiavg, phioldavg)
         end if
      end if

      if (proc0) then
         if (istep > window_size) then
            allocate (omega_avg1(naky, nakx))
            omega_avg1 = sum(omega_window, dim=1)
            omega_avg_window(mod(istep, window_size) + 1, :, :) = omega_avg1 / window_size
            if (istep > 2 * window_size .and. mod(istep, 10) == 0) then
               allocate (omega_square(naky, nakx))
               allocate (omega_avg2(naky, nakx))
               allocate (omega_avg_squared(naky, nakx))
               allocate (posmax(2))
               omega_avg2 = sum(omega_avg_window, dim=1)
               omega_square = sum(omega_avg_window * CONJG(omega_avg_window), dim=1)
               omega_std = omega_square - omega_avg2 * CONJG(omega_avg2) / real(window_size)
               omega_std = omega_std / real(window_size - 1)
               omega_avg_squared = omega_avg2 * CONJG(omega_avg2) / (window_size**2)
               maxstd = sqrt(maxval(omega_std / omega_avg_squared))
               posmax = maxloc(omega_std)
               posky = posmax(1)
               poskx = posmax(2)
               reached_convergence = maxstd < omega_prec
               conv_exit = conv_exit .or. reached_convergence

               !write (*, *)  ' maxstd ', maxstd, ' convergence reached ', reached_convergence, ' exit ', conv_exit, ' posky ', posky, ' poskx ', poskx
               deallocate (posmax)
               deallocate (omega_square)
               deallocate (omega_avg2)
               deallocate (omega_avg_squared)
            end if
            ! do not need omega_avg again this time step
            deallocate (omega_avg1)
         end if
      end if
      !check for NaN in omega_window
      if (proc0 .and. any(isNan(abs(omega_window)))) then
         conv_exit = .true.
         write(*,*) 'NaN in omega_window'
      end if
      call broadcast(conv_exit)

   end subroutine testing_convergence

! Rescale distribution function and phi (not working properly)
   subroutine rescale_g_and_phi(istep, overflow_exit)
      use mp, only: proc0, broadcast, mp_comm, mpireal
      use constants, only: zi
      use dist_fn_arrays, only: gnew, gold
      use fields_arrays, only: phi
      use fields_arrays, only: phi_old
      use kt_grids, only: naky, nakx, ny

      use volume_averages, only: volume_average, fieldline_average

      implicit none

      !> The current timestep
      integer, intent(in) :: istep

      logical, intent(in out) :: overflow_exit

      real :: zero
      real :: locmax, globmax
      real :: locmin, globmin
      integer :: ierror

      locmax = maxval(abs(gnew))
      locmin = minval(abs(gnew))
      zero = 100.*epsilon(0.)
      call mpi_allreduce(locmax, globmax, 1, mpireal, MPI_MAX, mp_comm, ierror)
      call mpi_allreduce(locmin, globmin, 1, mpireal, MPI_MIN, mp_comm, ierror)
      overflow_exit = overflow_exit .or. (globmax /= globmax) ! check for NaN, NaNQ
      !if (proc0.and.mod(istep, 1) == 0) then
      !   write (*, *) 'istep ', istep, ' globmax ', globmax, ' globmin ', globmin
      !end if

      IF (.NOT. overflow_exit .and. globmax > overflow_limit) then
         if (underflow_limit <= 0) then
            gold = gold / globmax * sqrt(1e-6)
            gnew = gnew / globmax * sqrt(1e-6)
            phi = phi / globmax
            phi_old = phi_old / globmax
         else
            !gold = gold/globmax
            gnew = gnew / globmax
            !phi = phi/globmax
            phi_old = phi_old / globmax
         end if
         !if (proc0.and.mod(istep, 1) == 0) write (*, *) 'rescaled'
      END IF
   end subroutine rescale_g_and_phi

   !==============================================
   !======== FINISH CONVERGENCE DIAGNOSTIC ============
   !==============================================
   subroutine finish_convergence

      call deallocate_arrays_convergence

   end subroutine finish_convergence

   subroutine finish_overflow

      call deallocate_arrays_overflow

   end subroutine finish_overflow

   !==============================================
   !============ DEALLCOATE ARRAYS ===============
   !==============================================
   subroutine deallocate_arrays_convergence

      implicit none

      if (allocated(omega_window)) deallocate (omega_window)
      if (allocated(omega_avg_window)) deallocate (omega_avg_window)

   end subroutine deallocate_arrays_convergence

   subroutine deallocate_arrays_overflow

      implicit none

      if (allocated(g_abs)) deallocate (g_abs)

   end subroutine deallocate_arrays_overflow

end module convergence
