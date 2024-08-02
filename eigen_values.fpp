#if WITH_EIG
!> A module for finding eigenvalues and eigenmodes
!> Requires SLEPC and PETSC
module eigen_values
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

   public :: time_eigval

   public :: eigval_functional
   public :: run_eigensolver
   public :: test_eigensolver
   public :: read_parameters, wnml_eigval, check_eigval

   public :: eigval_config_type
   public :: set_eigval_config
   public :: get_eigval_config

   private

   real, dimension(2) :: time_eigval = 0.

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

   !Do we evaluate the time derivative operator or the time advance operator?
   logical :: analyse_ddt_operator

   !If true then save restart files for each eigenmode
   logical :: save_restarts

   !Initialisation state
   logical :: initialized = .false.

   !> A custom type to make it easy to encapsulate specific settings
   type EpsSettings
      logical :: use_ginit
      logical :: analyse_ddt_operator
      PetscInt :: n_eig
      PetscInt :: max_iter
      PetscInt :: solver_option_switch
      PetscInt :: extraction_option_switch
      PetscInt :: which_option_switch
      PetscInt :: transform_option_switch
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
      !> How many stella timesteps to take each time SLEPc wants to
      !> advance the distribution function. Useful to separate closely
      !> spaced eigenvalues without changing [[knobs:delt]].
      integer :: nadv = 1
      !> If `true` then we save a set of restart files for each eigenmode
      !> found. These are named as the standard restart file
      !> (i.e. influenced by the restart_file input), but have `eig_<id>`
      !> appended near the end, where `<id>` is an integer representing the
      !> eigenmode id. If `save_distfn` of [[stella_diagnostics_knobs]] is true
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
      !> rather than the stella eigenvalue, \(\Omega\).
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

   !> Returns true if stella was compiled with WITH_EIG defined
   function eigval_functional()
      logical :: eigval_functional
      eigval_functional = .true.
   end function eigval_functional

   !> Initialise eigenvalue problem
   subroutine init_eigval(eigval_config_in)
      use mp, only: mp_comm, proc0
      use file_utils, only: error_unit
      implicit none
      type(eigval_config_type), intent(in), optional :: eigval_config_in
      PetscErrorCode :: ierr
#ifndef PETSC_USE_COMPLEX
      integer :: unit, ii
      integer, dimension(2) :: units
#endif

      !If we don't have eigenvalue support then ignore any settings in input file
      if (.not. eigval_functional()) then
         return
      end if

      if (initialized) return
      initialized = .true.

      !If petsc has not been compiled with complex support then display
      !a warning that the results are most likely incorrect.
#ifndef PETSC_USE_COMPLEX
      units(1) = error_unit() !Error file
      units(2) = 6            !Screen
      if (proc0) then
         do ii = 1, 2
            unit = units(ii)
            write (unit, '(66("#"))')
            write (unit, '("# ","WARNING : The PETSC library used does not have complex support"," #")')
            write (unit, '("# ","          --> The results are most likely incorrect.          "," #")')
            write (unit, '(66("#"))')
         end do
      end if
#endif

      !Get the input parameters
      call read_parameters(eigval_config_in)

      !Set PETSC_COMM_WORLD to be mp_comm so we can use whatever sub-comm
      !needed for list mode to work
      PETSC_COMM_WORLD = mp_comm

      !Initialise petsc
      !call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      !Initialise slepc
      call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
   end subroutine init_eigval

   !> Clean up eigenvalue problem
   subroutine finish_eigval
      implicit none
      PetscErrorCode :: ierr
      initialized = .false.
      !Finish up slepc
      call SlepcFinalize(ierr)
      !call PetscFinalize(ierr)
   end subroutine finish_eigval

   !> Read the eigval_knobs namelist
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

   !> Write the eigval_knobs namelist
   subroutine wnml_eigval(unit)
      use mp, only: mp_abort
      implicit none
      integer, intent(in) :: unit
      character(len=4) :: inden = "    "
      character(len=30) :: choice
      write (unit, *)
      write (unit, '(" &",A)') nml_name
      !Basic pars
      write (unit, '(A,A,"=",L1)') inden, 'analyse_ddt_operator', analyse_ddt_operator
      write (unit, '(A,A,"=",L1)') inden, 'use_ginit', use_ginit
      write (unit, '(A,A,"=",I0)') inden, 'n_eig', n_eig
      write (unit, '(A,A,"=",I0)') inden, 'max_iter', max_iter
      write (unit, '(A,A,"=",I0)') inden, 'nadv', nadv
      write (unit, '(A,A,"=",F12.5)') inden, 'toleranace', tolerance
      write (unit, '(A,A,"=",F12.5)') inden, 'targ_re', targ_re
      write (unit, '(A,A,"=",F12.5)') inden, 'targ_im', targ_im
      !string pars
      !Solver
      select case (solver_option_switch)
      case (SolverTypePower)
         choice = 'power'
      case (SolverTypeSubspace)
         choice = 'subspace'
      case (SolverTypeArnoldi)
         choice = 'arnoldi'
      case (SolverTypeLanczos)
         choice = 'lanczos'
      case (SolverTypeKrylovSchur)
         choice = 'krylov'
      case (SolverTypeGD)
         choice = 'GD'
      case (SolverTypeJD)
         choice = 'JD'
      case (SolverTypeRQCG)
         choice = 'RQCG'
      case (SolverTypeCISS)
         choice = 'CISS'
      case (SolverTypeLapack)
         choice = 'lapack'
      case (SolverTypeArpack)
         choice = 'arpack'
      case (SolverTypeBlzpack)
         choice = 'blzpack'
      case (SolverTypeTrlan)
         choice = 'trlan'
      case (SolverTypeBlopex)
         choice = 'blopex'
      case (SolverTypePrimme)
         choice = 'primme'
      case (SolverTypeFeast)
         choice = 'feast'
      case (SolverTypeNotSpecified)
         choice = 'slepc_default'
      case default
         !Should never get here
         call mp_abort("Unknown value of solver_option_switch")
      end select
      write (unit, '(A,A,"=",A,A,A)') inden, 'solver_option', '"', choice, '"'

      !Extraction
      select case (extraction_option_switch)
      case (ExtractionRitz)
         choice = 'ritz'
      case (ExtractionHarmonic)
         choice = 'harmonic'
      case (ExtractionHarmonicRelative)
         choice = 'harmonic_relative'
      case (ExtractionHarmonicRight)
         choice = 'harmonic_right'
      case (ExtractionHarmonicLargest)
         choice = 'harmonic_largest'
      case (ExtractionRefined)
         choice = 'refined'
      case (ExtractionRefinedHarmonic)
         choice = 'refined_harmonic'
      case (ExtractionNotSpecified)
         choice = 'slepc_default'
      case default
         !Should never get here
         call mp_abort("Unknown value of extraction_option_switch")
      end select
      write (unit, '(A,A,"=",A,A,A)') inden, 'extraction_option', '"', choice, '"'

      !Which type
      select case (which_option_switch)
      case (WhichLargestMagnitude)
         choice = 'largest_magnitude'
      case (WhichSmallestMagnitude)
         choice = 'smallest_magnitude'
      case (WhichLargestReal)
         choice = 'largest_real'
      case (WhichSmallestReal)
         choice = 'smallest_real'
      case (WhichLargestImaginary)
         choice = 'largest_imaginary'
      case (WhichSmallestImaginary)
         choice = 'smallest_imaginary'
      case (WhichTargetMagnitude)
         choice = 'target_magnitude'
      case (WhichTargetReal)
         choice = 'target_real'
      case (WhichTargetImaginary)
         choice = 'target_imaginary'
      case (WhichAll)
         choice = 'all'
      case (WhichUser)
         choice = 'user'
      case (WhichNotSpecified)
         choice = 'slepc_default'
      case default
         !Should never get here
         call mp_abort("Unknown value of which_option_switch")
      end select
      write (unit, '(A,A,"=",A,A,A)') inden, 'which_option', '"', choice, '"'

      !Transform type
      select case (transform_option_switch)
      case (TransformShell)
         choice = 'shell'
      case (TransformShift)
         choice = 'shift'
      case (TransformInvert)
         choice = 'invert'
      case (TransformCayley)
         choice = 'cayley'
      case (TransformFold)
         choice = 'fold'
      case (TransformPrecond)
         choice = 'precond'
      case (TransformNotSpecified)
         choice = 'slepc_default'
      case default
         !Should never get here
         call mp_abort("Unknown value of transform_option_switch")
      end select
      write (unit, '(A,A,"=",A,A,A)') inden, 'transform_option', '"', choice, '"'

      !Done
      write (unit, '(" /")')

   end subroutine wnml_eigval

   !> Check the eigval settings
   subroutine check_eigval
   end subroutine check_eigval

   !> Create a default settings object
   subroutine InitSettings(eps_settings)
      use stella_time, only: code_dt
      implicit none
      type(EpsSettings), intent(inout) :: eps_settings
      eps_settings%analyse_ddt_operator = analyse_ddt_operator
      eps_settings%use_ginit = use_ginit
      eps_settings%n_eig = n_eig
      eps_settings%max_iter = max_iter
      eps_settings%solver_option_switch = solver_option_switch
      eps_settings%extraction_option_switch = extraction_option_switch
      eps_settings%which_option_switch = which_option_switch
      eps_settings%transform_option_switch = transform_option_switch
      eps_settings%tolerance = tolerance
      eps_settings%targ_re = targ_re
      eps_settings%targ_im = targ_im
      eps_settings%target = cmplx(targ_re, targ_im)
      !Convert stella eigenvalue to slepc eigenvalue
      eps_settings%target_slepc = Eig_StellaToSlepc(eps_settings%target)
      !Initialise the dimension attributes to -1 to indicate unset
      eps_settings%local_size = -1
      eps_settings%global_size = -1
   end subroutine InitSettings

   !> Setup the matrix operator
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

   !> Create a solver and set parameters based on eps_settings
   subroutine InitEPS(eps_solver, mat_operator, eps_settings)
      use init_g, only: ginit
      use mp, only: barrier, proc0
      implicit none
      EPS, intent(inout) :: eps_solver !The solver object to initialise
      Mat, intent(in) :: mat_operator !The matrix to find eigenpairs of
      type(EpsSettings), intent(in) :: eps_settings
      Vec :: init_vec
      PetscErrorCode :: ierr
      logical :: restarted
      integer :: istep0
      PetscInt :: error1

      !Now create the eps solver
      call EPSCreate(PETSC_COMM_WORLD, eps_solver, ierr)

      !Set the matrix operator we want to find eigenvalues of
#if PETSC_VERSION_LT(3,8,0)
      call EPSSetOperators(eps_solver, mat_operator, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetOperators(eps_solver, mat_operator, PETSC_NULL_MAT, ierr)
#endif

      !Set the type of eigenproblem -->Always non-hermittian for us
      call EPSSetProblemType(eps_solver, EPS_NHEP, ierr)

      !Now setup from options
      call SetupEPS(eps_solver, eps_settings)

      !Now initialise if we requested this
      if (eps_settings%use_ginit) then
         !gvmu is initialised in the ginit routine
         ! the layout is later shifted to the kykxz layout

         !Setup vector
#if PETSC_VERSION_LT(3,6,0)
         call MatGetVecs(mat_operator, PETSC_NULL_OBJECT, init_vec, ierr)
#elif PETSC_VERSION_LT(3,8,0)
         call MatCreateVecs(mat_operator, PETSC_NULL_OBJECT, init_vec, ierr)
#else
         call MatCreateVecs(mat_operator, PETSC_NULL_VEC, init_vec, ierr)
#endif

         call GnewToVec(init_vec)

         !Set the initial vector
         error1 = 1
         call EPSSetInitialSpace(eps_solver, error1, init_vec, ierr)

         !Now destroy the vector
         call VecDestroy(init_vec, ierr)
      end if
   end subroutine InitEPS

   !> Setup the passed solver
   subroutine SetupEPS(eps_solver, eps_settings)
      use mp, only: mp_abort
      implicit none
      EPS, intent(inout) :: eps_solver
      type(EpsSettings), intent(in) :: eps_settings
      ST :: st
      PetscErrorCode :: ierr
      EPSType :: SolverType, TransformType
      PetscInt :: ExtractionType, WhichType

      !NOTE: There is no consistency checking here (yet) so it is up to the user
      !to make sure they request a valid combination of parameters

      !First set the number of eigenpairs to find
      call EPSSetDimensions(eps_solver, eps_settings%n_eig, &
                            PETSC_DECIDE, PETSC_DECIDE, ierr)

      !Now set tolerance and iteration limits
      call EPSSetTolerances(eps_solver, eps_settings%tolerance, &
                            eps_settings%max_iter, ierr)

      !Now set what type of eigenpairs we're looking for if we don't ask for slepc_default
      if (eps_settings%which_option_switch /= WhichNotSpecified) then
         select case (eps_settings%which_option_switch)
         case (WhichLargestMagnitude)
            WhichType = EPS_LARGEST_MAGNITUDE
         case (WhichSmallestMagnitude)
            WhichType = EPS_SMALLEST_MAGNITUDE
         case (WhichLargestReal)
            WhichType = EPS_LARGEST_REAL
         case (WhichSmallestReal)
            WhichType = EPS_SMALLEST_REAL
         case (WhichLargestImaginary)
            WhichType = EPS_LARGEST_IMAGINARY
         case (WhichSmallestImaginary)
            WhichType = EPS_SMALLEST_IMAGINARY
         case (WhichTargetMagnitude)
            WhichType = EPS_TARGET_MAGNITUDE
         case (WhichTargetReal)
            WhichType = EPS_TARGET_REAL
         case (WhichTargetImaginary)
            WhichType = EPS_TARGET_IMAGINARY
         case (WhichAll)
            WhichType = EPS_ALL
         case (WhichUser)
            WhichType = EPS_WHICH_USER
         case default
            !Should never get here
            call mp_abort("Unknown value of which_option_switch")
         end select

         call EpsSetWhichEigenpairs(eps_solver, int(WhichType), ierr)
      end if

      !Set the target
      call EpsSetTarget(eps_solver, eps_settings%target_slepc, ierr)

      !Now set the solver type if we don't ask for slepc_default
      if (eps_settings%solver_option_switch /= SolverTypeNotSpecified) then
         select case (eps_settings%solver_option_switch)
         case (SolverTypePower)
            SolverType = EPSPOWER
         case (SolverTypeSubspace)
            SolverType = EPSSUBSPACE
         case (SolverTypeArnoldi)
            SolverType = EPSARNOLDI
         case (SolverTypeLanczos)
            SolverType = EPSLANCZOS
         case (SolverTypeKrylovSchur)
            SolverType = EPSKRYLOVSCHUR
#if SLEPC_VERSION_GE(3,1,0)
         case (SolverTypeGD)
            SolverType = EPSGD
         case (SolverTypeJD)
            SolverType = EPSJD
#endif
#if SLEPC_VERSION_GE(3,3,0)
         case (SolverTypeRQCG)
            SolverType = EPSRQCG
#endif
#if SLEPC_VERSION_GE(3,4,0)
         case (SolverTypeCISS)
            SolverType = EPSCISS
#endif
         case (SolverTypeLapack)
            SolverType = EPSLAPACK
         case (SolverTypeArpack)
            SolverType = EPSARPACK
         case (SolverTypeTrlan)
            SolverType = EPSTRLAN
         case (SolverTypeBlopex)
            SolverType = EPSBLOPEX
         case (SolverTypePrimme)
            SolverType = EPSPRIMME
#if SLEPC_VERSION_GE(3,4,0)
#if SLEPC_VERSION_LE(3,11,0) || SLEPC_VERSION_GE(3,14,0)
         case (SolverTypeFeast)
            SolverType = EPSFEAST
#endif
#endif
         case default
            !Should never get here
            call mp_abort("Unknown value of solver_option_switch")
         end select

         call EPSSetType(eps_solver, SolverType, ierr)
      end if

      !Now set the extraction method
      if (eps_settings%extraction_option_switch /= ExtractionNotSpecified) then
         select case (eps_settings%extraction_option_switch)
         case (ExtractionRitz)
            ExtractionType = EPS_RITZ
         case (ExtractionHarmonic)
            ExtractionType = EPS_HARMONIC
         case (ExtractionHarmonicRelative)
            ExtractionType = EPS_HARMONIC_RELATIVE
         case (ExtractionHarmonicRight)
            ExtractionType = EPS_HARMONIC_RIGHT
         case (ExtractionHarmonicLargest)
            ExtractionType = EPS_HARMONIC_LARGEST
         case (ExtractionRefined)
            ExtractionType = EPS_REFINED
         case (ExtractionRefinedHarmonic)
            ExtractionType = EPS_REFINED_HARMONIC
         case default
            !Should never get here
            call mp_abort("Unknown value of extraction_option_switch")
         end select

         call EPSSetExtraction(eps_solver, int(ExtractionType), ierr)
      end if

      !Now set the spectral transform
      if (eps_settings%transform_option_switch /= TransformNotSpecified) then
         !Get spectral transform object
         call EPSGetSt(eps_solver, st, ierr)

         select case (eps_settings%transform_option_switch)
         case (TransformShell)
            TransformType = STSHELL
         case (TransformShift)
            TransformType = STSHIFT
         case (TransformInvert)
            TransformType = STSINVERT
         case (TransformCayley)
            TransformType = STCAYLEY
#if SLEPC_VERSION_LE(3,4,4)
         case (TransformFold)
            TransformType = STFOLD
#endif
         case (TransformPrecond)
            TransformType = STPRECOND
         case default
            !Should never get here
            call mp_abort("Unknown value of transform_option_switch")
         end select

         !Set the shift
         call STSetShift(st, eps_settings%target_slepc, ierr)

         !Set the type
         call STSetType(st, TransformType, ierr)

         !Set the matrix mode to shell -- This is needed to use a lot of the
         !spectral transforms but is commented out for the moment because a
         !few other settings would need to be set including the KSP and PC
         !types as without implementing other matrix methods (like MATOP_GET_DIAGONAL)
         !the default KSP/PC won't work. To enable this at run time try something like
         ! -st_type sinvert -st_ksp_type gmres -st_pc_type none -st_matmode shell
         !call STSetMatMode(st,ST_MATMODE_SHELL,ierr)
      end if
   end subroutine SetupEPS

   !> Writes the current settings to a file with extension '.eig.solver_settings'
   subroutine ReportSolverSettings(eps_solver, fext)
      use mp, only: proc0
      use file_utils, only: open_output_file, close_output_file
      implicit none
      EPS, intent(in) :: eps_solver
      character(len=*), intent(in), optional :: fext !Optional text to add to file name
      character(len=80) :: extension
      integer :: ounit
      EPSType :: TmpType !String
      EPSExtraction :: Extract !Integer
      PetscInt :: TmpInt
      PetscScalar :: TmpScal !Complex
      PetscReal :: TmpReal
      ST :: st
      PetscErrorCode :: ierr
      integer :: target_type

      !Only proc0 does the writing
      if (.not. proc0) return

      !Set the default extension
      extension = ".eig.solver_settings"

      !Add optional text is passed
      if (present(fext)) then
         extension = trim(fext)//trim(extension)
      end if

      !Open an output file for writing
      call open_output_file(ounit, extension)
      write (ounit, '("Slepc eigensolver settings:")')

      !Now we fetch settings and write them to file
#ifdef PETSC_USE_COMPLEX
      write (ounit, '("   ",A30)') "Compiled with COMPLEX support"
#else
      write (ounit, '("   ",A30)') "Compiled with REAL support"
#endif

      !1. Solver type
      call EPSGetType(eps_solver, TmpType, ierr)
      write (ounit, '("   ",A22,2X,":",2X,A22)') "Type", adjustl(TmpType)

      !2. Number of eigenvalues to look for
      call EPSGetDimensions(eps_solver, TmpInt, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      write (ounit, '("   ",A22,2X,":",2X,I0)') "Number of eigenvalues", TmpInt

      !3. Tolerances
      call EPSGetTolerances(eps_solver, TmpReal, TmpInt, ierr)
      write (ounit, '("   ",A22,2X,":",2X,I0)') "Max iterations", TmpInt
      write (ounit, '("   ",A22,2X,":",1X,E11.4)') "Tolerance", TmpReal

      !4. ExtractionType
      call EPSGetExtraction(eps_solver, Extract, ierr)
      !Could have a select case (in function) to convert integer to string name
      write (ounit, '("   ",A22,2X,":",2X,I0)') "Extraction type", Extract !Note integer, not string

      !5. WhichType
      call EPSGetWhichEigenPairs(eps_solver, target_type, ierr)
      !Could have a select case (in function) to convert integer to string name
      write (ounit, '("   ",A22,2X,":",2X,I0)') "Target type", target_type !Note integer, not string
      !5b. Target value
      call EPSGetTarget(eps_solver, TmpScal, ierr)
#ifdef PETSC_USE_COMPLEX
      write (ounit, '("   ",A22,2X,":",2X,F12.5)') "Real target (slepc)", PetscRealPart(TmpScal)
      write (ounit, '("   ",A22,2X,":",2X,F12.5)') "Imag target (slepc)", PetscImaginaryPart(TmpScal)
#else
      write (ounit, '("   ",A22,2X,":",2X,F12.5)') "Mag. target (slepc)", PetscRealPart(TmpScal)
#endif

      !6. Transform type
      call EPSGetST(eps_solver, st, ierr)
      call STGetType(st, TmpType, ierr)
      write (ounit, '("   ",A22,2X,":",2X,A22)') "Transform type", adjustl(TmpType)

      !7. Analysis of DDT operator

      write (ounit, '("   ",A22,2X,":",2X,L1)') "Analyse DDT operator", analyse_ddt_operator

      !8. Use ginit
      write (ounit, '("   ",A22,2X,":",2X,L1)') "Use ginit", use_ginit

      !9. nadv
      write (ounit, '("   ",A22,2X,":",2X,I0)') "nadv", nadv

      !Close output file
      call close_output_file(ounit)
   end subroutine ReportSolverSettings

   !> Main routine that does the setup, starts the solver and
  !! coordinates reporting of the results
   subroutine BasicSolve
      use mp, only: barrier, proc0
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use stella_layouts, only: vmu_lo

      implicit none
      EPS :: my_solver
      Mat :: my_operator
      PetscInt :: d1_size, d2_size, d3_size, d4_size
      PetscInt :: d5_size_local, d5_size_global
      PetscInt :: n_loc, n_glob
      type(EpsSettings) :: my_settings
      PetscErrorCode :: ierr

      !Define dimensions for problem
      d1_size = naky !Size of ky grid
      d2_size = nakx !Size of kx grid
      d3_size = 2 * nzgrid + 1 !Size of z grid
      d4_size = ntubes !Size of tube grid
      d5_size_local = (1 + vmu_lo%ulim_proc - vmu_lo%llim_proc) !Size of local distributed domain
      if (vmu_lo%ulim_proc < vmu_lo%llim_proc) d5_size_local = 0

      d5_size_global = (1 + vmu_lo%ulim_world - vmu_lo%llim_world) !Size of global distributed domain

      !How big is the "advance operator matrix"
      n_loc = d1_size * d2_size * d3_size * d4_size * d5_size_local
      n_glob = d1_size * d2_size * d3_size * d4_size * d5_size_global
      call barrier
      !Pack the settings type
      if (proc0) write (*, *) "Initialising settings"
      call InitSettings(my_settings)

      !Set the sizes
      my_settings%local_size = n_loc
      my_settings%global_size = n_glob

      !Initialise the matrix operator
      if (proc0) write (*, *) "Initialising matrix operator"
      if (analyse_ddt_operator) then
         call InitMatrixOperator(my_operator, ddt_eigen, my_settings)
      else
         call InitMatrixOperator(my_operator, advance_eigen, my_settings)
      end if

      !Create the solver
      if (proc0) write (*, *) "Initialising solver"
      call InitEPS(my_solver, my_operator, my_settings)

      !Attach a monitor
      !call EPSMonitorSet(my_solver,SimpleMonitor,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)

      !Write settings to file
      if (proc0) write (*, *) "Writing settings to file"
      call ReportSolverSettings(my_solver)

      !Now we're ready to solve
      call barrier
      if (proc0) write (*, *) "Solving"
      call EPSSolve(my_solver, ierr)

      !Now do the reporting and diagnostic output
      if (proc0) write (*, *) "Reporting results"
      call ReportResults(my_solver, my_operator, my_settings)

      call barrier
      !Now destroy objects
      if (proc0) write (*, *) "Destroying objects"
      call EPSDestroy(my_solver, ierr)
      call MatDestroy(my_operator, ierr)
   end subroutine BasicSolve

   !> Routine to write results to screen and netcdf
   !TODO save eigen values to netcdf
   subroutine ReportResults(my_solver, my_operator, my_settings)
      use mp, only: proc0
      use stella_time, only: code_dt
      use run_parameters, only: fphi, fapar, fbpar
      use stella_io, only: EigNetcdfID, init_eigenfunc_file, add_eigenpair_to_file, finish_eigenfunc_file
      use file_utils, only: run_name
      use dist_fn_arrays, only: gnew
      implicit none
      EPS, intent(in) :: my_solver
      Mat, intent(in) :: my_operator
      type(EpsSettings), intent(in) :: my_settings
      PetscErrorCode :: ierr
      PetscInt :: iteration_count, n_converged
      Vec :: eig_vec_r, eig_vec_i
      PetscScalar :: eig_val_r, eig_val_i
      complex :: EigVal
      logical :: all_found
      PetscInt :: ieig
      type(EigNetcdfID) :: io_ids
      character(len=20) :: restart_unique_string
      complex, dimension(:), allocatable :: eigenvalues
      real, dimension(:), allocatable :: local_conv
      integer :: istatus

      !Find out how many iterations were performed
      call EPSGetIterationNumber(my_solver, iteration_count, ierr)

      !Find out how many converged eigenpairs where found
      call EPSGetConverged(my_solver, n_converged, ierr)
      all_found = (n_converged >= my_settings%n_eig)
      if (proc0) write (6, '("Found ",I0," eigenvalues in ",I0," iterations.")') n_converged, iteration_count

      !Initialise the eigenvalues array
      if (.not. allocated(eigenvalues)) then
         allocate (eigenvalues(n_converged))
      end if
      if (.not. allocated(local_conv)) then
         allocate (local_conv(n_converged))
      end if
      if (n_converged > 0) then
         !Initialise the eigenvector arrays
#if PETSC_VERSION_LT(3,6,0)
         call MatGetVecs(my_operator, PETSC_NULL_OBJECT, eig_vec_r, ierr)
         call MatGetVecs(my_operator, PETSC_NULL_OBJECT, eig_vec_i, ierr)
#elif PETSC_VERSION_LT(3,8,0)
         call MatCreateVecs(my_operator, PETSC_NULL_OBJECT, eig_vec_r, ierr)
         call MatCreateVecs(my_operator, PETSC_NULL_OBJECT, eig_vec_i, ierr)
#else
         call MatCreateVecs(my_operator, PETSC_NULL_VEC, eig_vec_r, ierr)
         call MatCreateVecs(my_operator, PETSC_NULL_VEC, eig_vec_i, ierr)
#endif
         if (proc0) call init_eigenfunc_file(trim(run_name)//".eig.out.nc", fphi, fapar, fbpar, io_ids)

         !Now loop over converged values
         do ieig = 0, n_converged - 1
            !Get eigenpair
            call EPSGetEigenpair(my_solver, ieig, eig_val_r, &
                                 eig_val_i, eig_vec_r, eig_vec_i, ierr)

            !NOTE: If petsc has been compiled with complex support then _i variables
            !are empty, else contains imaginary component
#ifdef PETSC_USE_COMPLEX
            EigVal = cmplx(PetscRealPart(eig_val_r), PetscImaginaryPart(eig_val_r))
#else
            EigVal = cmplx(PetscRealPart(eig_val_r), PetscRealPart(eig_val_i))
#endif
            !Convert to stella eigenvalue
            Eigval = Eig_SlepcToStella(EigVal)

            !Get field eigenmodes
            !/First set g
#ifdef PETSC_USE_COMPLEX
            call VecToGnew(eig_vec_r)
#else
            call VecToGnew2(eig_vec_r, eig_vec_i)
#endif

            !Now store the eigenvalue
            eigenvalues(ieig + 1) = EigVal
            local_conv(ieig + 1) = ieig + 1
         end do
         if (proc0) call add_eigenpair_to_file(eigenvalues, io_ids, local_conv)

         if (proc0) call finish_eigenfunc_file(io_ids)

         !Close netcdf file

         call VecDestroy(eig_vec_r, ierr)
         call VecDestroy(eig_vec_i, ierr)
      else
         if (proc0) write (6, '("No converged eigenvalues found.")')
      end if
      if (allocated(eigenvalues)) deallocate (eigenvalues)
      if (allocated(local_conv)) deallocate (local_conv)
   end subroutine ReportResults

   !> Function to convert SLEPc eigenvalue to a stella one
   elemental complex function Eig_SlepcToStella(eig)
      use stella_time, only: code_dt
      implicit none
      complex, intent(in) :: eig
      if (analyse_ddt_operator) then
         Eig_SlepcToStella = eig * cmplx(0.0, 1.0)
      else
         Eig_SlepcToStella = log(eig) * cmplx(0.0, 1.0) / (code_dt * nadv)
      end if
   end function Eig_SlepcToStella

   !> Function to convert a stella eigenvalue to a SLEPc one
   elemental PetscScalar function Eig_StellaToSlepc(eig)
      use stella_time, only: code_dt
      implicit none
      complex, intent(in) :: eig
      !If PETSC doesn't have a complex type then return the
      !magnitude of the eigenvalue.
      if (analyse_ddt_operator) then
#ifdef PETSC_USE_COMPLEX
         Eig_StellaToSlepc = -cmplx(0.0, 1.0) * eig
#else
         Eig_StellaToSlepc = abs(cmplx(0.0, 1.0) * eig)
#endif
      else
#ifdef PETSC_USE_COMPLEX
         Eig_StellaToSlepc = exp(-cmplx(0.0, 1.0) * code_dt * eig * nadv)
#else
         Eig_StellaToSlepc = abs(exp(-cmplx(0.0, 1.0) * code_dt * eig * nadv))
#endif
      end if
   end function Eig_StellaToSlepc

   !> Call back routine to represent linear advance operator
   subroutine advance_eigen(MatOperator, VecIn, Res, ierr)
      use fields, only: init_fields
      use time_advance, only: advance_linear
      use dist_fn_arrays, only: gnew
      use stella_time, only: code_dt
      implicit none
      Mat, intent(in) :: MatOperator
      Vec, intent(inout) :: VecIn, Res
      PetscErrorCode, intent(inout) :: ierr
      integer, parameter :: istep = 2
      integer :: is
      !First unpack input vector into gnew
      call VecToGnew(VecIn)

      !Now set fields to be consistent with gnew
      call init_fields

      !Now do a number of time steps
      do is = 1, nadv
         !Note by using a fixed istep we
         !force the explicit terms to be excluded
         !except for the very first call.
         call advance_linear(istep)
      end do

      !Now pack gnew into output
      call GnewToVec(Res)
   end subroutine advance_eigen

   !> Call back routine to represent time derivative operator
   subroutine ddt_eigen(MatOperator, VecIn, Res, ierr)
      use fields, only: init_fields
      use time_advance, only: advance_linear
      use dist_fn_arrays, only: gnew
      use stella_time, only: code_dt
      use zgrid, only: nzgrid
      implicit none
      Mat, intent(in) :: MatOperator
      Vec, intent(inout) :: VecIn, Res
      PetscErrorCode, intent(inout) :: ierr
      integer, parameter :: istep = 2
      integer :: is
      complex, dimension(:, :, :, :, :), allocatable :: ginit
      !First unpack input vector into gnew
      call VecToGnew(VecIn)

      ! Store a copy of the initial state
      ginit = gnew

      !Now set fields to be consistent with gnew
      call init_fields

      !Now do a number of time steps
      do is = 1, nadv
         !Note by using a fixed istep we
         !force the explicit terms to be excluded
         !except for the very first call.
         call advance_linear(istep)
      end do

      ! Form a crude estimate of the time derivative
      gnew = (gnew - ginit) / (code_dt * nadv)

      !Now pack gnew into output
      call GnewToVec(Res)
   end subroutine ddt_eigen

   !> Unpacks a SLEPc vector into gnew
   subroutine VecToGnew(VecIn)
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use stella_layouts, only: vmu_lo
      use dist_fn_arrays, only: gnew
      implicit none
      Vec, intent(inout) :: VecIn
      PetscScalar, pointer :: VecInArr(:)
      PetscInt :: d1_size, d2_size, d3_size, d4_size, d5_size_local, d5_size_global
      PetscErrorCode :: ierr
      integer :: itube, iz, ikx, iky, ivmu, local_index

      !No local data
      if (vmu_lo%ulim_proc < vmu_lo%llim_proc) return

      !Define dimensions
      d1_size = naky !Size of ky grid
      d2_size = nakx !Size of kx grid
      d3_size = 2 * nzgrid + 1 !Size of z grid
      d4_size = ntubes !Size of tube grid
      d5_size_local = (1 + vmu_lo%ulim_proc - vmu_lo%llim_proc) !Size of local distributed domain
      d5_size_global = (1 + vmu_lo%ulim_world - vmu_lo%llim_world) !Size of global distributed domain

      !Get a pointer to the data
      call VecGetArrayReadF90(VecIn, VecInArr, ierr)

      !Extract
      local_index = 0
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
         do itube = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     local_index = 1 + local_index
                     gnew(iky, ikx, iz, itube, ivmu) = VecInArr(local_index)
                  end do
               end do
            end do
         end do
      end do

      !Get rid of pointer
      call VecRestoreArrayReadF90(VecIn, VecInArr, ierr)
   end subroutine VecToGnew

#ifndef PETSC_USE_COMPLEX
   subroutine VecToGnew2(VecIn, VecInI)
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use stella_layouts, only: vmu_lo
      use dist_fn_arrays, only: gnew
      implicit none
      Vec, intent(inout) :: VecIn, VecInI
      PetscScalar, pointer :: VecInArr(:), VecInIArr(:)
      PetscInt :: d1_size, d2_size, d3_size, d4_size, d5_size_local, d5_size_global
      PetscErrorCode :: ierr
      integer :: itube, iz, ikx, iky, ivmu, local_index

      !No local data
      if (vmu_lo%ulim_proc < vmu_lo%llim_proc) return

      !Define dimensions
      d1_size = naky !Size of ky grid
      d2_size = nakx !Size of kx grid
      d3_size = 2 * nzgrid + 1 !Size of z grid
      d4_size = ntubes !Size of tube grid
      d5_size_local = (1 + vmu_lo%ulim_proc - vmu_lo%llim_proc) !Size of local distributed domain
      d5_size_global = (1 + vmu_lo%ulim_world - vmu_lo%llim_world) !Size of global distributed domain

      !Get a pointer to the data
      call VecGetArrayReadF90(VecIn, VecInArr, ierr)
      call VecGetArrayReadF90(VecInI, VecInIArr, ierr)

      !Extract
      local_index = 0
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
         do itube = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     !Form local index (note we could just having a running counter which we
                     !increment on each loop
                     local_index = 1 + local_index
                     gnew(iky, ikx, iz, itube, ivmu) = cmplx(PetscRealPart(VecInArr(local_index)), PetscRealPart(VecInIArr(local_index)))
                  end do
               end do
            end do
         end do

         !Get rid of pointer
         call VecRestoreArrayReadF90(VecIn, VecInArr, ierr)
         call VecRestoreArrayReadF90(VecInI, VecInIArr, ierr)
         end subroutine VecToGnew2
#endif

         !> Packs gnew into a SLEPc vector
         subroutine GnewToVec(VecIn)
            use zgrid, only: nzgrid, ntubes
            use kt_grids, only: naky, nakx
            use stella_layouts, only: vmu_lo
            use dist_fn_arrays, only: gnew
            use mp, only: barrier, proc0, nproc, iproc, mp_abort
            implicit none
            Vec, intent(inout) :: VecIn
            PetscScalar, pointer :: VecInArr(:)
            PetscInt :: d1_size, d2_size, d3_size, d4_size, d5_size_local, d5_size_global
            PetscErrorCode :: ierr
            integer :: itube, iz, ikx, iky, ivmu, local_index

            !No local data
            if (vmu_lo%ulim_proc < vmu_lo%llim_proc) return

            !Define dimensions
            d1_size = naky !Size of ky grid
            d2_size = nakx !Size of kx grid
            d3_size = 2 * nzgrid + 1 !Size of z grid
            d4_size = ntubes !Size of tube grid
            d5_size_local = (1 + vmu_lo%ulim_proc - vmu_lo%llim_proc) !Size of local distributed domain
            d5_size_global = (1 + vmu_lo%ulim_world - vmu_lo%llim_world) !Size of global distributed domain
            !Get a pointer to the data
            call VecGetArrayF90(VecIn, VecInArr, ierr)
            !Extract
            local_index = 0
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
               do itube = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     do ikx = 1, nakx
                        do iky = 1, naky
                           !local_index= (ivmu-vmu_lo%llim_proc)*d1_size*d2_size*d3_size*d4_size + &
                           !             (itube-1)*d1_size*d2_size*d3_size + &
                           !             (iz+nzgrid)*d1_size*d2_size + &
                           !             (ikx-1)*d1_size + &
                           !             iky
                           local_index = local_index + 1
                           VecInArr(local_index) = gnew(iky, ikx, iz, itube, ivmu)
                        end do
                     end do
                  end do
               end do
            end do
            !Get rid of pointer
            call VecRestoreArrayF90(VecIn, VecInArr, ierr)
         end subroutine GnewToVec

         !> An example of a custom SLEPc monitor -- may not work well
         PetscErrorCode function SimpleMonitor(solver, iters, nconv, eig_r, eig_i, err_est, n_est, mm)
            use mp, only: proc0
            use stella_time, only: code_dt
            implicit none
            EPS, intent(inout) :: solver
            PetscInt, intent(inout) :: iters, nconv, n_est
            PetscScalar, dimension(:), intent(inout) :: eig_r, eig_i
            PetscReal, dimension(:), intent(inout) :: err_est
            integer, pointer, intent(inout), optional :: mm
            integer, save :: nconv_prev = 0, iter_prev = 0
            integer :: new_conv, new_iter

            !Set return value
            SimpleMonitor = 0

            !If no change in nconv then return
            if (nconv_prev == nconv) then
               return
            end if

            !Calculate how many extra modes found and how many iterations it took
            new_conv = nconv - nconv_prev
            new_iter = iters - iter_prev
            if (proc0) write (6, '("Found ",I0," eigenvalues in ",I0," iterations:")') new_conv, new_iter

            !Update previous value
            nconv_prev = nconv
            iter_prev = iters
         end function SimpleMonitor

         !> Set the module level config type
         !> Will abort if the module has already been initialised to avoid
         !> inconsistencies.
         subroutine set_eigval_config(eigval_config_in)
            use mp, only: mp_abort
            type(eigval_config_type), intent(in), optional :: eigval_config_in
            if (initialized) then
               call mp_abort("Trying to set eigval_config when already initialized.")
            end if
            if (present(eigval_config_in)) then
               eigval_config = eigval_config_in
            end if
         end subroutine set_eigval_config

         !> Get the module level config instance
         function get_eigval_config()
            type(eigval_config_type) :: get_eigval_config
            get_eigval_config = eigval_config
         end function get_eigval_config

         !---------------------------------------
         ! Following is for the config_type
         !---------------------------------------

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

         subroutine run_eigensolver
            use job_manage, only: time_message
            use mp, only: mp_abort, proc0, barrier

            !Start timer
            call time_message(.false., time_eigval, ' Eigensolver')

            call barrier
            if (proc0) write (6, *) 'Initialising eigensolver'
            !Initialise slepc and the default/input eigensolver parameters
            call init_eigval

            call barrier
            if (proc0) write (6, *) 'Running eigensolver'
            !Create a solver based on input paramters, use it to solve and
            !then destroy it.
            call BasicSolve

            call barrier
            if (proc0) write (6, *) 'Finishing eigensolver results'
            !Tidy up
            call finish_eigval

            !Stop timer
            call time_message(.false., time_eigval, ' Eigensolver')

         end subroutine run_eigensolver

         subroutine test_eigensolver
            use mp, only: mp_comm, proc0, nproc, barrier, iproc
            implicit none
            PetscErrorCode :: ierr
            Mat :: my_operator
            PetscInt :: n, Istart, Iend, i, index
            PetscInt :: n_converged, iteration_count
            EPS :: my_solver
            PetscScalar :: eig_val_r, eig_val_i
            real :: start_time, end_time

            !Test the eigensolver
            PETSC_COMM_WORLD = mp_comm
            n = 10000

            !Initialise slepc
            call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
            if (proc0) write (*, *) "Slepc initialised"
            call cpu_time(start_time)
            call MatCreate(PETSC_COMM_WORLD, my_operator, ierr)
            call MatSetSizes(my_operator, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
            call MatSetFromOptions(my_operator, ierr)
            call MatSetUp(my_operator, ierr)
            if (proc0) write (*, *) "Matrix created"
            call MatGetOwnershipRange(my_operator, Istart, Iend, ierr)
            write (*, *) "Rank: ", iproc, " Istart: ", Istart, " Iend: ", Iend
            call barrier
            do i = Istart, Iend
               index = i - 1
               call MatSetValue(my_operator, index, index, complex(index * 1.0, (n - index - 1) * (1.0)), INSERT_VALUES, ierr)
            end do
            call MatAssemblyBegin(my_operator, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(my_operator, MAT_FINAL_ASSEMBLY, ierr)

            if (proc0) write (*, *) "Matrix assembled"
            call EPSCreate(PETSC_COMM_WORLD, my_solver, ierr)
            call EPSSetOperators(my_solver, my_operator, PETSC_NULL_MAT, ierr)
            !call EpsSetWhichEigenpairs(my_solver, EPS_LARGEST_REAL, ierr)
            call EPSSetFromOptions(my_solver, ierr)
            call ReportSolverSettings(my_solver)
            if (proc0) write (*, *) "Solver created"
            call EPSSolve(my_solver, ierr)
            if (proc0) write (*, *) "Solved"
            call EPSGetConverged(my_solver, n_converged, ierr)
            call EPSGetIterationNumber(my_solver, iteration_count, ierr)
            if (proc0) write (*, *) "Converged: ", n_converged, " after ", iteration_count, " iterations"

            do i = 0, n_converged - 1
               call EPSGetEigenvalue(my_solver, i, eig_val_r, eig_val_i, ierr)
               if (proc0) write (*, *) "Eigenvalue: ", eig_val_r, " with magnitude: ", abs(eig_val_r)
            end do
            call EPSDestroy(my_solver, ierr)
            call MatDestroy(my_operator, ierr)
            if (proc0) write (*, *) "Matrix destroyed"
            call cpu_time(end_time)
            if (proc0) write (*, *) "Time taken: ", end_time - start_time
            if (proc0) write (*, *) "Total CPU Time taken: ", (end_time - start_time) * nproc
            call SlepcFinalize(ierr)
            if (proc0) write (*, *) "Slepc finalised"

         end subroutine test_eigensolver

         subroutine advance_test(MatOperator, VecIn, Res, ierr)
            use dist_fn_arrays, only: gnew
            implicit none
            Mat, intent(in) :: MatOperator
            Vec, intent(inout) :: VecIn, Res
            PetscErrorCode, intent(inout) :: ierr
            integer, parameter :: istep = 1
            integer :: is
            !First unpack input vector into gnew
            call VecToGnew(VecIn)

            !Now set fields to be consistent with gnew

            !Now do a number of time steps
            do is = 1, nadv
               !Note by using a fixed istep we
               !force the explicit terms to be excluded
               !except for the very first call.
               gnew = gnew * complex(0, 1)
            end do

            !Now pack gnew into output
            call GnewToVec(Res)
         end subroutine advance_test

         end module eigen_values

#else

!###############################################################
!THE SECTION BELOW PROVIDES THE EIGENVAL MODULE IN THE ABSENCE OF
!SLEPC/PETSC (i.e. WITH_EIG not defined).

!> A stub module representing the eigenvalue solver
         module eigen_values
            implicit none
            private
            public :: eigval_functional
            public :: time_eigval
            public :: run_eigensolver
            public :: test_eigensolver
            real, dimension(2) :: time_eigval = 0.

         contains
            !> Returns true if stella was compiled with WITH_EIG defined -- here forced to .false.
  !! as compiled without SLEPc support
            function eigval_functional()
               logical :: eigval_functional
               eigval_functional = .false.
            end function eigval_functional

            subroutine run_eigensolver

               use job_manage, only: time_message
               use mp, only: mp_abort, proc0
               if (proc0) write (*, *) "Stella was compiled without Eigenvalue solver support."
               call mp_abort("Require slepc/petsc")

            end subroutine run_eigensolver

            subroutine test_eigensolver

               use job_manage, only: time_message
               use mp, only: mp_abort, proc0
               if (proc0) write (*, *) "Stella was compiled without Eigenvalue solver support."
               call mp_abort("Require slepc/petsc")

            end subroutine test_eigensolver
         end module eigen_values

#endif
