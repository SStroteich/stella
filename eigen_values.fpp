!> A module for finding eigenvalues and eigenmodes
!> Requires SLEPC and PETSC
module eigval
   use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
   use petscvec
   use petscmat
   use slepceps

   implicit none

!Allows us to use SLEPC_VERSION_MAJOR etc. to automatically
!disable unsupported features
#include <slepcversion.h>
#include <petscversion.h>

#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>
#elif PETSC_VERSION_LT(3,8,0)
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#else
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#endif

#if SLEPC_VERSION_LT(3,6,0)
#include <finclude/slepcepsdef.h>
#elif SLEPC_VERSION_LT(3,8,0)
#include <slepc/finclude/slepcepsdef.h>
#else
#include <slepc/finclude/slepceps.h>
#endif

   public :: eigval_functional
   public :: init_eigval, finish_eigval, BasicSolve
   public :: read_parameters, wnml_eigval, check_eigval

   public :: eigval_config_type
   public :: set_eigval_config
   public :: get_eigval_config

   private

   !//////////////////
   !Slepc related configuration parameters
   !//////////////////
   !Solver type
   integer, parameter :: &        !Note not all of these may be available
      SolverTypePower = 1, &
      SolverTypeSubspace = 2, &
      SolverTypeArnoldi = 3, &
      SolverTypeLanczos = 4, &
      SolverTypeKrylovSchur = 5, &
      SolverTypeGD = 6, &
      SolverTypeJD = 7, &
      SolverTypeRQCG = 8, &
      SolverTypeCISS = 9, &
      SolverTypeLapack = 10, &
      SolverTypeArpack = 11, &
      SolverTypeBlzpack = 12, &
      SolverTypeTrlan = 13, &
      SolverTypeBlopex = 14, &
      SolverTypePrimme = 15, &
      SolverTypeFeast = 16, &
      SolverTypeNotSpecified = 17 !=>Use slepc default
   integer :: solver_option_switch

   !Extraction type
   integer, parameter :: & !Note not all extractions are supported by every solver
      ExtractionRitz = 1, &
      ExtractionHarmonic = 2, &
      ExtractionHarmonicRelative = 3, &
      ExtractionHarmonicRight = 4, &
      ExtractionHarmonicLargest = 5, &
      ExtractionRefined = 6, &
      ExtractionRefinedHarmonic = 7, &
      ExtractionNotSpecified = 8 !=>Use slepc default
   integer :: extraction_option_switch

   !Which eigenpairs do we search for
   integer, parameter :: &
      WhichLargestMagnitude = 1, &
      WhichSmallestMagnitude = 2, &
      WhichLargestReal = 3, &
      WhichSmallestReal = 4, &
      WhichLargestImaginary = 5, &
      WhichSmallestImaginary = 6, &
      WhichTargetMagnitude = 7, &
      WhichTargetReal = 8, &
      WhichTargetImaginary = 9, & !Complex builds only
      WhichAll = 10, &   !Requires an interval to be set
      WhichUser = 11, &  !Requires a used comparison routine to be defined
      WhichNotSpecified = 12 !=>Use slepc default
   integer :: which_option_switch

   !What sort of spectral transform do we use
   integer, parameter :: &
      TransformShell = 1, &
      TransformShift = 2, &
      TransformInvert = 3, &
      TransformCayley = 4, &
      TransformFold = 5, &
      TransformPrecond = 6, &
      TransformNotSpecified = 7
   integer :: transform_option_switch

   !//////////////////

   !Number of eigenvalues to seek
   integer :: n_eig

   !Maximum number of iterations (calls to advance) allowed
   integer :: max_iter

   !Tolerance to converge to
   real :: tolerance

   !Real and imaginary components of target (roughly where we expect eigvals to be)
   real :: targ_re, targ_im

   !Do we specify the initial condition (through ginit)?
   logical :: use_ginit

   !Number of stella advance steps per SLEPc call to advance_eigen
   !May be useful when looking at marginal modes etc.
   integer :: nadv

   !If true then save restart files for each eigenmode
   logical :: save_restarts

   !Initialisation state
   logical :: initialized = .false.

   !> A custom type to make it easy to encapsulate specific settings
   type EpsSettings
      logical :: use_ginit
      logical :: analyse_ddt_operator
      integer :: n_eig
      integer :: max_iter
      integer :: solver_option_switch
      integer :: extraction_option_switch
      integer :: which_option_switch
      integer :: transform_option_switch
      PetscInt :: local_size, global_size
      real :: tolerance
      real :: targ_re, targ_im
      complex :: target
      PetscScalar :: target_slepc
   end type EpsSettings

   character(len=12), parameter :: nml_name = "eigval_knobs"

   type, extends(abstract_config_type) :: eigval_config_type
      ! namelist : eigval_knobs
      ! indexed : false
      !> Determines which operator SLEPc is finding eigenvalues of.  If
      !> .false. then SLEPc analyses the time advance operator, so
      !> internally works with the eigenvalue `exp(-i*omega*nadv*dt)`.
      !> If .true. then SLEPc analyses a time derivative operator, so
      !> internally works with the eigenvalue `-i*omega`. Note we form
      !> a first order one sided approximation of the time derivative
      !> (e.g. `(g_{n+nadv}-g_{n})/(nadv*code_dt)`) to be analysed.
      logical :: analyse_ddt_operator = .false.
      !> Sets the extraction technique, must be one of:
      !>
      !> - 'default' (use SLEPC default)
      !> - 'slepc_default' (use SLEPC default)
      !> - 'ritz'
      !> - 'harmonic'
      !> - 'harmonic_relative'
      !> - 'harmonic_right'
      !> - 'harmonic_largest'
      !> - 'refined'
      !> - 'refined_harmonic'
      !>
      character(len=20) :: extraction_option = 'default'
      !> Sets the maximum number of SLEPC iterations used.  If not set
      !> (recommended) then let SLEPC decide what to use (varies with
      !> different options).
      integer :: max_iter = PETSC_DECIDE
      !> The number of eigenmodes to search for. The actual number of
      !> modes found may be larger than this.
      integer :: n_eig = 1
      !> How many GS2 timesteps to take each time SLEPc wants to
      !> advance the distribution function. Useful to separate closely
      !> spaced eigenvalues without changing [[knobs:delt]].
      integer :: nadv = 1
      !> If `true` then we save a set of restart files for each eigenmode
      !> found. These are named as the standard restart file
      !> (i.e. influenced by the restart_file input), but have `eig_<id>`
      !> appended near the end, where `<id>` is an integer representing the
      !> eigenmode id. If `save_distfn` of [[gs2_diagnostics_knobs]] is true
      !> then will also save the distribution function files.
      logical :: save_restarts = .false.
      !> Sets the type of solver to use, must be one of:
      !>
      !> - 'default' (Krylov-Schur)
      !> - 'slepc_default' (Krylov-Schur)
      !> - 'power'
      !> - 'subspace'
      !> - 'arnoldi'
      !> - 'lanczos'
      !> - 'krylov'
      !> - 'GD'
      !> - 'JD'
      !> - 'RQCG'
      !> - 'CISS'
      !> - 'lapack'
      !> - 'arpack'
      !> - 'blzpack'
      !> - 'trlan'
      !> - 'blopex'
      !> - 'primme'
      !> - 'feast'
      !>
      !> Not all solver types are compatible with other eigenvalue options,
      !> some options may not be supported in older SLEPC versions and some
      !> may require certain flags to be set when SLEPC is compiled.
      character(len=20) :: solver_option = 'default'
      !> Imaginary part of the eigenvalue target. Often beneficial to
      !> set this fairly large (e.g. 10) when looking for unstable
      !> modes. Used with the `target_*` [[eigval_knobs:which_option]]
      !> mode.
      real :: targ_im = 0.5
      !> Real part of the eigenvalue target. Often beneficial to set
      !> this fairly small (e.g. ~0). Used with the `target_*`
      !> [[eigval_knobs:which_option]] mode.
      real :: targ_re = 0.5
      !> Sets tolerance on SLEPC eigenmode search. Used to determine
      !> when an eigenmode has been found. Note this is the tolerance
      !> based on the SLEPc eigenvalue, which is the time advance
      !> eigenvalue
      !> \(\exp{\left(-i\Omega\textrm{nadv delt}\right)}\)
      !> rather than the GS2 eigenvalue, \(\Omega\).
      real :: tolerance = 1.0d-6
      !> Sets the type of spectral transform to be used. This can help
      !> extract interior eigenvalues. Must be one of
      !>
      !> - 'default' (let SLEPC decide)
      !> - 'slepc_default' (let SLEPC decide)
      !> - 'shell'
      !> - 'shift'
      !> - 'invert'
      !> - 'cayley'
      !> - 'fold'
      !> - 'precond' (not implemented)
      !>
      !> Not all options are available in all versions of the library.
      character(len=20) :: transform_option = 'default'
      !> If `true` then provide an initial guess for the eigenmode
      !> based on using [[init_g]] routines to initialise `g`. Can be
      !> helpful in accelerating initial phase of eigensolver. Can also
      !> be quite useful with `ginit_option='many'` (etc.)  to start an
      !> eigenvalue search from a previously obtained solution,
      !> allowing both eigenmode refinement and sub-dominant mode
      !> detection.
      logical :: use_ginit = .false.
      !> Sets SLEPC mode of operation (i.e. what sort of eigenvalues it
      !>looks for). Note that this refers to the SLEPc eigenvalue,
      !>i.e. the time advance eigenvalue
      !>\(\exp{\left(-i\Omega\textrm{nadv delt}\right)}\). In other
      !>words one should select `'largest_real'` to search for modes
      !>with the largest growth rate. Must be one of
      !>
      !> - 'default' (equivalent to 'target_magnitude')
      !> - 'slepc_default' (let SLEPC decide)
      !> - 'largest_magnitude'
      !> - 'smallest_magnitude'
      !> - 'largest_real'
      !> - 'smallest_real'
      !> - 'largest_imaginary'
      !> - 'smallest_imaginary'
      !> - 'target_magnitude' (complex eigenvalue magnitude closest to magnitude of target)
      !> - 'target_real'
      !> - 'target_imaginary'
      !> - 'all' (only some solver types, e.g. lapack) -- not supported
      !> - 'user' (will use a user specified function to pick between eigenmodes, note not currently implemented)
      !>
      character(len=20) :: which_option = 'default'
   contains
      procedure, public :: read => read_eigval_config
      procedure, public :: write => write_eigval_config
      procedure, public :: reset => reset_eigval_config
      procedure, public :: broadcast => broadcast_eigval_config
      procedure, public, nopass :: get_default_name => get_default_name_eigval_config
      procedure, public, nopass :: get_default_requires_index => get_default_requires_index_eigval_config
   end type eigval_config_type

   type(eigval_config_type) :: eigval_config

contains

! procedures for the eigval_config_type

!> Reads in the eigval_knobs namelist and populates the member variables
   subroutine read_eigval_config(self)
      use file_utils, only: input_unit_exist, get_indexed_namelist_unit
      use mp, only: proc0
      implicit none
      class(eigval_config_type), intent(in out) :: self
      logical :: exist
      integer :: in_file

      ! Note: When this routine is in the module where these variables live
      ! we shadow the module level variables here. This is intentional to provide
      ! isolation and ensure we can move this routine to another module easily.
      character(len=20) :: extraction_option, solver_option, transform_option, which_option
      integer :: max_iter, n_eig, nadv
      logical :: analyse_ddt_operator, save_restarts, use_ginit
      real :: targ_im, targ_re, tolerance

    namelist /eigval_knobs/ analyse_ddt_operator, extraction_option, max_iter, n_eig, nadv, save_restarts, solver_option, targ_im, targ_re, tolerance, transform_option, use_ginit, which_option

      ! Only proc0 reads from file
      if (.not. proc0) return

      ! First set local variables from current values
      analyse_ddt_operator = self%analyse_ddt_operator
      extraction_option = self%extraction_option
      max_iter = self%max_iter
      n_eig = self%n_eig
      nadv = self%nadv
      save_restarts = self%save_restarts
      solver_option = self%solver_option
      targ_im = self%targ_im
      targ_re = self%targ_re
      tolerance = self%tolerance
      transform_option = self%transform_option
      use_ginit = self%use_ginit
      which_option = self%which_option

      ! Now read in the main namelist
      in_file = input_unit_exist(self%get_name(), exist)
      if (exist) read (in_file, nml=eigval_knobs)

      ! Now copy from local variables into type members
      self%analyse_ddt_operator = analyse_ddt_operator
      self%extraction_option = extraction_option
      self%max_iter = max_iter
      self%n_eig = n_eig
      self%nadv = nadv
      self%save_restarts = save_restarts
      self%solver_option = solver_option
      self%targ_im = targ_im
      self%targ_re = targ_re
      self%tolerance = tolerance
      self%transform_option = transform_option
      self%use_ginit = use_ginit
      self%which_option = which_option

      self%exist = exist
   end subroutine read_eigval_config

   !> Writes out a namelist representing the current state of the config object
   subroutine write_eigval_config(self, unit)
      implicit none
      class(eigval_config_type), intent(in) :: self
      integer, intent(in), optional:: unit
      integer :: unit_internal

      unit_internal = 6 ! @todo -- get stdout from file_utils
      if (present(unit)) then
         unit_internal = unit
      end if

      call self%write_namelist_header(unit_internal)
      call self%write_key_val("analyse_ddt_operator", self%analyse_ddt_operator, unit_internal)
      call self%write_key_val("extraction_option", self%extraction_option, unit_internal)
      call self%write_key_val("max_iter", self%max_iter, unit_internal)
      call self%write_key_val("n_eig", self%n_eig, unit_internal)
      call self%write_key_val("nadv", self%nadv, unit_internal)
      call self%write_key_val("save_restarts", self%save_restarts, unit_internal)
      call self%write_key_val("solver_option", self%solver_option, unit_internal)
      call self%write_key_val("targ_im", self%targ_im, unit_internal)
      call self%write_key_val("targ_re", self%targ_re, unit_internal)
      call self%write_key_val("tolerance", self%tolerance, unit_internal)
      call self%write_key_val("transform_option", self%transform_option, unit_internal)
      call self%write_key_val("use_ginit", self%use_ginit, unit_internal)
      call self%write_key_val("which_option", self%which_option, unit_internal)
      call self%write_namelist_footer(unit_internal)
   end subroutine write_eigval_config

   !> Resets the config object to the initial empty state
   subroutine reset_eigval_config(self)
      class(eigval_config_type), intent(in out) :: self
      type(eigval_config_type) :: empty
      select type (self)
      type is (eigval_config_type)
         self = empty
      end select
   end subroutine reset_eigval_config

   !> Broadcasts all config parameters so object is populated identically on
  !! all processors
   subroutine broadcast_eigval_config(self)
      use mp, only: broadcast
      implicit none
      class(eigval_config_type), intent(in out) :: self
      call broadcast(self%analyse_ddt_operator)
      call broadcast(self%extraction_option)
      call broadcast(self%max_iter)
      call broadcast(self%n_eig)
      call broadcast(self%nadv)
      call broadcast(self%save_restarts)
      call broadcast(self%solver_option)
      call broadcast(self%targ_im)
      call broadcast(self%targ_re)
      call broadcast(self%tolerance)
      call broadcast(self%transform_option)
      call broadcast(self%use_ginit)
      call broadcast(self%which_option)

      call broadcast(self%exist)
   end subroutine broadcast_eigval_config

   !> Gets the default name for this namelist
   function get_default_name_eigval_config()
      implicit none
      character(len=CONFIG_MAX_NAME_LEN) :: get_default_name_eigval_config
      get_default_name_eigval_config = "eigval_knobs"
   end function get_default_name_eigval_config

   !> Gets the default requires index for this namelist
   function get_default_requires_index_eigval_config()
      implicit none
      logical :: get_default_requires_index_eigval_config
      get_default_requires_index_eigval_config = .false.
   end function get_default_requires_index_eigval_config

! subroutines for the eigen value calculation

   function eigval_functional()
      logical :: eigval_functional
      eigval_functional = .true.
   end function eigval_functional

   subroutine init_eigval(eigval_config_in)
      use mp, only: mp_comm, proc0
      implicit none
      type(eigval_config_type), intent(in), optional :: eigval_config_in

      PetscErrorCode :: ierr

      if (.not. eigval_functional()) then
         return
      end if
      if (initialized) return
      initialized = .true.
      call read_parameters(eigval_config_in)

      !Set PETSC_COMM_WORLD to be mp_comm so we can use whatever sub-comm
      !needed for list mode to work
      PETSC_COMM_WORLD = mp_comm

      !Initialise slepc
      call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
   end subroutine init_eigval

   subroutine read_parameters(eigval_config_in)
      use file_utils, only: error_unit
      use text_options, only: text_option, get_option_value
      implicit none
      type(eigval_config_type), intent(in), optional :: eigval_config_in

      !NOTE: text_option is case sensitive (currently)
      type(text_option), dimension(18), parameter :: solver_type_opts = &
                                                     (/ &
                                                     text_option('slepc_default', SolverTypeNotSpecified), &
                                                     text_option('default', SolverTypeKrylovSchur), &
                                                     text_option('power', SolverTypePower), &
                                                     text_option('subspace', SolverTypeSubspace), &
                                                     text_option('arnoldi', SolverTypeArnoldi), &
                                                     text_option('lanczos', SolverTypeLanczos), &
                                                     text_option('krylov', SolverTypeKrylovSchur), &
                                                     text_option('GD', SolverTypeGD), &
                                                     text_option('JD', SolverTypeJD), &
                                                     text_option('RQCG', SolverTypeRQCG), &
                                                     text_option('CISS', SolverTypeCISS), &
                                                     text_option('lapack', SolverTypeLapack), &
                                                     text_option('arpack', SolverTypeArpack), &
                                                     text_option('blzpack', SolverTypeBlzpack), &
                                                     text_option('trlan', SolverTypeTrlan), &
                                                     text_option('blopex', SolverTypeBlopex), &
                                                     text_option('primme', SolverTypePrimme), &
                                                     text_option('feast', SolverTypeFeast) &
                                                     /)

      type(text_option), dimension(9), parameter :: extraction_type_opts = &
                                                    (/ &
                                                    text_option('slepc_default', ExtractionNotSpecified), &
                                                    text_option('default', ExtractionNotSpecified), &
                                                    text_option('ritz', ExtractionRitz), &
                                                    text_option('harmonic', ExtractionHarmonic), &
                                                    text_option('harmonic_relative', ExtractionHarmonicRelative), &
                                                    text_option('harmonic_right', ExtractionHarmonicRight), &
                                                    text_option('harmonic_largest', ExtractionHarmonicLargest), &
                                                    text_option('refined', ExtractionRefined), &
                                                    text_option('refined_harmonic', ExtractionRefinedHarmonic) &
                                                    /)

      type(text_option), dimension(13), parameter :: which_type_opts = &
                                                     (/ &
                                                     text_option('slepc_default', WhichNotSpecified), &
                                                     text_option('default', WhichTargetMagnitude), &
                                                     text_option('largest_magnitude', WhichLargestMagnitude), &
                                                     text_option('smallest_magnitude', WhichSmallestMagnitude), &
                                                     text_option('largest_real', WhichLargestReal), &
                                                     text_option('smallest_real', WhichSmallestReal), &
                                                     text_option('largest_imaginary', WhichLargestImaginary), &
                                                     text_option('smallest_imaginary', WhichSmallestImaginary), &
                                                     text_option('target_magnitude', WhichTargetMagnitude), &
                                                     text_option('target_real', WhichTargetReal), &
                                                     text_option('target_imaginary', WhichTargetImaginary), &
                                                     text_option('all', WhichAll), &
                                                     text_option('user', WhichUser) &
                                                     /)

      type(text_option), dimension(8), parameter :: transform_type_opts = &
                                                    (/ &
                                                    text_option('slepc_default', TransformNotSpecified), &
                                                    text_option('default', TransformNotSpecified), &
                                                    text_option('shell', TransformShell), &
                                                    text_option('shift', TransformShift), &
                                                    text_option('invert', TransformInvert), &
                                                    text_option('cayley', TransformCayley), &
                                                    text_option('fold', TransformFold), &
                                                    text_option('precond', TransformPrecond) &
                                                    /)

      character(len=20) :: solver_option, extraction_option, which_option, transform_option

      integer :: ierr
      logical :: exist

      if (present(eigval_config_in)) eigval_config = eigval_config_in

      call eigval_config%init(name='eigval_knobs', requires_index=.false.)

      ! Copy out internal values into module level parameters
      analyse_ddt_operator = eigval_config%analyse_ddt_operator
      extraction_option = eigval_config%extraction_option
      max_iter = eigval_config%max_iter
      n_eig = eigval_config%n_eig
      nadv = eigval_config%nadv
      save_restarts = eigval_config%save_restarts
      solver_option = eigval_config%solver_option
      targ_im = eigval_config%targ_im
      targ_re = eigval_config%targ_re
      tolerance = eigval_config%tolerance
      transform_option = eigval_config%transform_option
      use_ginit = eigval_config%use_ginit
      which_option = eigval_config%which_option

      exist = eigval_config%exist

      !Get error unit for output
      ierr = error_unit()

      !Convert string options to integer flags
      call get_option_value &
         (solver_option, solver_type_opts, solver_option_switch, &
          ierr, "solver_option in "//nml_name, .true.)

      call get_option_value &
         (extraction_option, extraction_type_opts, extraction_option_switch, &
          ierr, "extraction_option in "//nml_name, .true.)

      call get_option_value &
         (which_option, which_type_opts, which_option_switch, &
          ierr, "which_option in "//nml_name, .true.)

      call get_option_value &
         (transform_option, transform_type_opts, transform_option_switch, &
          ierr, "transform_option in "//nml_name, .true.)
   end subroutine read_parameters

   subroutine finish_eigval
      implicit none
      PetscErrorCode :: ierr
      initialized = .false.
      !Finish up slepc
      call SlepcFinalize(ierr)
   end subroutine finish_eigval

   subroutine InitMatrixOperator(mat_operator, Adv_routine, eps_settings)
      implicit none
      Mat, intent(inout) :: mat_operator
      external :: Adv_routine !This returns the result of mat_operator.X (where X is a vector)
      type(EpsSettings), intent(in) :: eps_settings
      PetscErrorCode :: ierr
      !Note, whilst we're defining a matrix (i.e. 2D array) we only specify
      !the local and global length of one side. This is because it is a shell
      !matrix where we only ever deal with vectors which can interact with the
      !matrix and never the matrix itself (hence why Sum(n_loc*n_loc)/=n_glob*n_glob
      !is ok!).

      !Make a shell matrix operator (AdvMat)
      call MatCreateShell(PETSC_COMM_WORLD, eps_settings%local_size, eps_settings%local_size, &
                          eps_settings%global_size, eps_settings%global_size, PETSC_NULL_INTEGER, mat_operator, ierr)

      !Set the shell MATOP_MULT operation, i.e. which routine returns the result
      !of a AdvMat.x
      call MatShellSetOperation(mat_operator, MATOP_MULT, Adv_routine, ierr)
   end subroutine InitMatrixOperator
end module eigval
