!> Energy diagnostics data container and lifecycle helpers
module energy_diagnostic

   implicit none

   public :: energy_diagnostics_type, energy_diag
   public :: init_energy_diagnostic, finish_energy_diagnostic
   public :: get_free_energy

   type :: energy_diagnostics_type
      real, dimension(:), allocatable :: energy_total, dedt_total
      real, dimension(:), allocatable :: drive_term, diss_perp, diss_zed, diss_vpa
      real, dimension(:), allocatable :: drifts_term, streaming_term, nonlinear_term, mirror_term
      real, dimension(:), allocatable :: weights_energy, factor_spec

      real, dimension(:), allocatable :: energy_total_vmu

      complex, dimension(:, :, :, :, :), allocatable :: velocity_integral1
      complex, dimension(:, :, :), allocatable :: spatial_integral1

      real, dimension(:, :, :, :, :), allocatable :: free_energy_kxkyz, dedt_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: diss_perp_kxkyz, diss_zed_kxkyz, diss_vpa_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: drive_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: drifts_kxkyz, streaming_kxkyz, nonlinear_kxkyz, mirror_kxkyz

      real, dimension(:, :, :), allocatable :: free_energy_vmu, dedt_vmu
      real, dimension(:, :, :), allocatable :: diss_perp_vmu, diss_zed_vmu, diss_vpa_vmu
      real, dimension(:, :, :), allocatable :: drive_vmu
      real, dimension(:, :, :), allocatable :: drifts_vmu, streaming_vmu, nonlinear_vmu, mirror_vmu
   end type energy_diagnostics_type

   type(energy_diagnostics_type), save :: energy_diag

   logical, parameter :: debug = .false.

contains

   subroutine init_energy_diagnostic(write_energy_kxkyz, write_energy_vmu)

      use species, only: nspec
      use kt_grids, only: nakx, naky
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nmu, nvpa

      implicit none

      logical, intent(in) :: write_energy_kxkyz
      logical, intent(in) :: write_energy_vmu

      if (.not. allocated(energy_diag%energy_total)) allocate (energy_diag%energy_total(nspec)); energy_diag%energy_total = 0.
      if (.not. allocated(energy_diag%dedt_total)) allocate (energy_diag%dedt_total(nspec)); energy_diag%dedt_total = 0.
      if (.not. allocated(energy_diag%drive_term)) allocate (energy_diag%drive_term(nspec)); energy_diag%drive_term = 0.
      if (.not. allocated(energy_diag%diss_perp)) allocate (energy_diag%diss_perp(nspec)); energy_diag%diss_perp = 0.
      if (.not. allocated(energy_diag%diss_zed)) allocate (energy_diag%diss_zed(nspec)); energy_diag%diss_zed = 0.
      if (.not. allocated(energy_diag%diss_vpa)) allocate (energy_diag%diss_vpa(nspec)); energy_diag%diss_vpa = 0.
      if (.not. allocated(energy_diag%drifts_term)) allocate (energy_diag%drifts_term(nspec)); energy_diag%drifts_term = 0.
      if (.not. allocated(energy_diag%streaming_term)) allocate (energy_diag%streaming_term(nspec)); energy_diag%streaming_term = 0.
      if (.not. allocated(energy_diag%nonlinear_term)) allocate (energy_diag%nonlinear_term(nspec)); energy_diag%nonlinear_term = 0.
      if (.not. allocated(energy_diag%mirror_term)) allocate (energy_diag%mirror_term(nspec)); energy_diag%mirror_term = 0.

      if (.not. allocated(energy_diag%weights_energy)) allocate (energy_diag%weights_energy(nspec)); energy_diag%weights_energy = 1.
      if (.not. allocated(energy_diag%factor_spec)) allocate (energy_diag%factor_spec(nspec)); energy_diag%factor_spec = 0.

      if (.not. allocated(energy_diag%energy_total_vmu)) allocate (energy_diag%energy_total_vmu(nspec)); energy_diag%energy_total_vmu = 0.

      if (.not. allocated(energy_diag%velocity_integral1)) then
         allocate (energy_diag%velocity_integral1(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         energy_diag%velocity_integral1 = 0.
      end if
      if (.not. allocated(energy_diag%spatial_integral1)) then
         allocate (energy_diag%spatial_integral1(nvpa, nmu, nspec))
         energy_diag%spatial_integral1 = 0.
      end if

      if (write_energy_kxkyz) then
         if (.not. allocated(energy_diag%free_energy_kxkyz)) allocate (energy_diag%free_energy_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%dedt_kxkyz)) allocate (energy_diag%dedt_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%diss_perp_kxkyz)) allocate (energy_diag%diss_perp_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%diss_zed_kxkyz)) allocate (energy_diag%diss_zed_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%diss_vpa_kxkyz)) allocate (energy_diag%diss_vpa_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%drive_kxkyz)) allocate (energy_diag%drive_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%drifts_kxkyz)) allocate (energy_diag%drifts_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%streaming_kxkyz)) allocate (energy_diag%streaming_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%nonlinear_kxkyz)) allocate (energy_diag%nonlinear_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         if (.not. allocated(energy_diag%mirror_kxkyz)) allocate (energy_diag%mirror_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      end if

      if (write_energy_vmu) then
         if (.not. allocated(energy_diag%free_energy_vmu)) allocate (energy_diag%free_energy_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%dedt_vmu)) allocate (energy_diag%dedt_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%diss_perp_vmu)) allocate (energy_diag%diss_perp_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%diss_zed_vmu)) allocate (energy_diag%diss_zed_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%diss_vpa_vmu)) allocate (energy_diag%diss_vpa_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%drive_vmu)) allocate (energy_diag%drive_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%drifts_vmu)) allocate (energy_diag%drifts_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%streaming_vmu)) allocate (energy_diag%streaming_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%nonlinear_vmu)) allocate (energy_diag%nonlinear_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%mirror_vmu)) allocate (energy_diag%mirror_vmu(nvpa, nmu, nspec))
      end if

   end subroutine init_energy_diagnostic

   !the subroutine takes input values of g, phi, factor_spec and returns sum_spec, sum_total and the array term_kxkyz
   !
   subroutine get_one_energy_term_kxkyz(h, term, factor_spec, sum_spec, sum_total, term_kxkyz)
      use mp, only: proc0, sum_allreduce

      use dist_fn_arrays, only: g0
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: nspec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: wgts_mu, wgts_vpa

      use vpamu_grids, only: integrate_vmu
      use volume_averages, only: mode_fac
      use kt_grids, only: naky, nakx
      use stella_geometry, only: dVolume
      use volume_averages, only: volume_total

      implicit none
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h, term

      real, dimension(nspec), intent(in) :: factor_spec
      real, dimension(nspec), intent(out) :: sum_spec
      real, intent(out) :: sum_total
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: term_kxkyz

      integer :: ivmu, imu, iv, iz, it, is, ia, ikx, iky

      energy_diag%weights_energy = 1.
      sum_spec = 0.
      term_kxkyz = 0.
      energy_diag%velocity_integral1 = 0.
      sum_total = 0.
      g0 = 0.

      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do iz = -nzgrid, nzgrid
            g0(:, :, iz, :, ivmu) = term(:, :, iz, :, ivmu) * conjg(h(:, :, iz, :, ivmu)) &
                                    / (maxwell_fac(is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is))
            energy_diag%velocity_integral1(:, :, iz, :, is) = energy_diag%velocity_integral1(:, :, iz, :, is) + &
                                                          wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g0(:, :, iz, :, ivmu) * energy_diag%weights_energy(is)
         end do
      end do

      call sum_allreduce(energy_diag%velocity_integral1)

      if (proc0) then
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     do iky = 1, naky
               term_kxkyz(iky, ikx, iz, it, is) = 0.5 * mode_fac(iky) * (real(factor_spec(is) * energy_diag%velocity_integral1(iky, ikx, iz, it, is)))
                        sum_spec(is) = sum_spec(is) + term_kxkyz(iky, ikx, iz, it, is) * dVolume(ia, ikx, iz)
                     end do
                  end do
               end do
            end do
            sum_spec(is) = sum_spec(is) / volume_total
            sum_total = sum_total + sum_spec(is)
         end do
      end if
      g0 = 0.

   end subroutine get_one_energy_term_kxkyz

   subroutine get_one_energy_term_vmu(h, term, factor_spec, term_vmu)
      use mp, only: proc0, sum_allreduce
      use dist_fn_arrays, only: g0, gvmu0
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: is_idx, ikx_idx, iky_idx, iz_idx, it_idx
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: nspec
      use zgrid, only: ntubes, nzgrid
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: nvpa, nmu

      use volume_averages, only: mode_fac
      use stella_geometry, only: dVolume, bmag
      use volume_averages, only: volume_total

      use redistribute, only: scatter
      use dist_redistribute, only: kxkyz2vmu

      implicit none
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h, term

      real, dimension(nspec), intent(in) :: factor_spec

      real, dimension(:, :, :), intent(out) :: term_vmu

      integer :: ivmu, imu, iv, is
      integer :: ikxkyz, iz, it, ia, ikx, iky

      energy_diag%weights_energy = 1.

      term_vmu = 0.

      energy_diag%spatial_integral1 = 0.
      g0 = 0.

      ia = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0(:, :, iz, it, ivmu) = term(:, :, iz, it, ivmu) * conjg(h(:, :, iz, it, ivmu))
            end do
         end do
      end do

      call scatter(kxkyz2vmu, g0, gvmu0)

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         is = is_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            do iv = 1, nvpa    
               energy_diag%spatial_integral1(iv, imu, is) = energy_diag%spatial_integral1(iv, imu, is) + mode_fac(iky) * bmag(ia,iz) * &
                     factor_spec(is) * gvmu0(iv, imu, ikxkyz) * dVolume(ia, ikx, iz) / (maxwell_fac(is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is))
            end do
         end do
      end do

      energy_diag%spatial_integral1 = energy_diag%spatial_integral1 / volume_total
      call sum_allreduce(energy_diag%spatial_integral1)

      if (proc0) then
         do is = 1, nspec
            do imu = 1, nmu
               do iv = 1, nvpa
                  term_vmu(iv, imu, is) = real(energy_diag%spatial_integral1(iv, imu, is))
               end do
            end do
         end do
      end if

      g0 = 0.
      gvmu0 = 0.

   end subroutine get_one_energy_term_vmu

! Here sum_spec and sum_total are integrated over vpa and mu which is not the case in the subroutine above.
! This subroutine is kept for reference since the total integrals are coinciding with the one from get_one_energy_term_kxkyz
!
   subroutine get_one_energy_term_vmu_ref(h, term, factor_spec, sum_spec, sum_total, term_vmu)
      use mp, only: proc0, sum_allreduce
      use dist_fn_arrays, only: g0, gvmu0
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: is_idx, ikx_idx, iky_idx, iz_idx, it_idx
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: nspec
      use zgrid, only: ntubes, nzgrid
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: wgts_mu, wgts_vpa, wgts_mu_bare
      use volume_averages, only: mode_fac
      use stella_geometry, only: dVolume, bmag
      use volume_averages, only: volume_total

      use redistribute, only: gather, scatter
      use dist_redistribute, only: kxkyz2vmu

      implicit none
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h, term

      real, dimension(nspec), intent(out) :: sum_spec
      real, dimension(nspec), intent(in) :: factor_spec
      real, intent(out) :: sum_total

      real, dimension(:, :, :), intent(out) :: term_vmu

      integer :: ivmu, imu, iv, is
      integer :: ikxkyz, iz, it, ia, ikx, iky

      energy_diag%weights_energy = 1.
      sum_spec = 0.
      term_vmu = 0.
      sum_total = 0.
      energy_diag%spatial_integral1 = 0.
      g0 = 0.

      ia = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0(:, :, iz, it, ivmu) = term(:, :, iz, it, ivmu) * conjg(h(:, :, iz, it, ivmu))
            end do
         end do
      end do

      call scatter(kxkyz2vmu, g0, gvmu0)

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         is = is_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            do iv = 1, nvpa
               energy_diag%spatial_integral1(iv, imu, is) = energy_diag%spatial_integral1(iv, imu, is) + 0.5 * mode_fac(iky) * bmag(ia, iz) * &
               factor_spec(is) * gvmu0(iv, imu, ikxkyz) * dVolume(ia, ikx, iz) / (maxwell_fac(is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is))
            end do
         end do
      end do

      energy_diag%spatial_integral1 = energy_diag%spatial_integral1 / volume_total
      call sum_allreduce(energy_diag%spatial_integral1)

      if (proc0) then
         do is = 1, nspec
            do imu = 1, nmu
               do iv = 1, nvpa
                  term_vmu(iv, imu, is) = real(energy_diag%spatial_integral1(iv, imu, is))
                  sum_spec(is) = sum_spec(is) + 2 * wgts_mu_bare(imu) * wgts_vpa(iv) * term_vmu(iv, imu, is) * energy_diag%weights_energy(is)
               end do
            end do
            sum_total = sum_total + sum_spec(is)
         end do
      end if

      g0 = 0.
      gvmu0 = 0.

   end subroutine get_one_energy_term_vmu_ref

   !> Calculate free energy, the drive term and the dissipation
   !>
   subroutine get_free_energy(h, g, phi, istep, energy_unit, write_energy_vmu)

      use mp, only: proc0
      use dist_fn_arrays, only: g1, kperp2, gold2
      use fields_arrays, only: phi_zero
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: spec, nspec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: mu, vpa, nmu, nvpa
      use run_parameters, only: fphi
      use kt_grids, only: naky, nakx
      use stella_time, only: code_time, code_dt
      use stella_geometry, only: b_dot_grad_z, dbdzed
      use physics_flags, only: nonlinear
      use hyper, only: D_hyper, k2max, advance_hyper_zed, advance_hyper_vpa, hyp_vpa, hyp_zed
      use redistribute, only: gather, scatter
      use dist_redistribute, only: kxkyz2vmu
      use stella_layouts, only: kxyz_lo, kxkyz_lo, vmu_lo
      use time_advance, only: advance_wdriftx_explicit, advance_wdrifty_explicit, advance_ExB_nonlinearity
      use time_advance, only: advance_wstar_explicit
      use parallel_streaming, only: advance_parallel_streaming_explicit
      use mirror_terms, only: advance_mirror_explicit
      use g_tofrom_h, only: g_to_h

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: g, h
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi

      integer, intent(in) :: istep
      integer, intent(in) :: energy_unit
      logical, intent(in) :: write_energy_vmu

      integer :: ivmu, imu, iv, iz, it, is, ia, ikx, iky, ikxkyz
      real :: energy_sum, energy_sum_vmu
      real :: dedt_sum

      real :: diss_perp_sum
      real :: diss_zed_sum
      real :: diss_vpa_sum
      real :: drive_sum
      real :: drifts_sum
      real :: streaming_sum
      real :: nonlinear_sum
      real :: mirror_sum

      real :: total_sum

      logical :: restart_time_step
      restart_time_step = .false.

      energy_diag%free_energy_kxkyz = 0.
      energy_diag%dedt_kxkyz = 0.
      energy_diag%diss_perp_kxkyz = 0.
      energy_diag%diss_zed_kxkyz = 0.
      energy_diag%diss_vpa_kxkyz = 0.
      energy_diag%drive_kxkyz = 0.
      energy_diag%drifts_kxkyz = 0.
      energy_diag%streaming_kxkyz = 0.
      energy_diag%mirror_kxkyz = 0.
      energy_diag%nonlinear_kxkyz = 0.

      energy_sum = 0.
      energy_sum_vmu = 0.
      dedt_sum = 0.

      diss_perp_sum = 0.
      diss_zed_sum = 0.
      diss_vpa_sum = 0.

      drive_sum = 0.
      drifts_sum = 0.
      streaming_sum = 0.
      nonlinear_sum = 0.
      mirror_sum = 0.

      total_sum = 0.

      ia = 1

      phi_zero = 0.
      g1 = 0.

      ! FLAG - electrostatic for now
      ! get electrostatic contributions to energy terms
      if (fphi > epsilon(0.0)) then

         ! Calculate free energy
         ! This is g * h_conj

         energy_diag%factor_spec = spec%dens * spec%temp
         call get_one_energy_term_kxkyz(h, g, energy_diag%factor_spec, energy_diag%energy_total, energy_sum, energy_diag%free_energy_kxkyz)

         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, g, energy_diag%factor_spec, energy_diag%free_energy_vmu)
         end if

         if (debug) then
            if (write_energy_vmu) then
            call get_one_energy_term_vmu_ref(h, g, energy_diag%factor_spec, energy_diag%energy_total_vmu, energy_sum_vmu, energy_diag%free_energy_vmu)
               if (proc0) then
                  write (*, *) 'Free energy from kxkyz: ', energy_sum
                  if (write_energy_vmu) then
                     write (*, *) 'Free energy from vmu: ', energy_sum_vmu
                     write (*, *) 'Check: ratio = ', energy_sum / energy_sum_vmu
                  end if
               end if
            end if
         end if

         ! Calculate dE/dt
         ! This is (g - g_old) / dt * h_conj
         energy_diag%factor_spec = spec%dens * spec%temp
         g1 = (g - gold2) / code_dt

         call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%dedt_total, dedt_sum, energy_diag%dedt_kxkyz)

         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%dedt_vmu)
         end if
         ! Calculate dissipation perpendicular
         ! This is - D_hyper * k_perp^4 * g * h_conj
         energy_diag%factor_spec = -D_hyper * spec%dens * spec%temp
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g1(:, :, iz, it, ivmu) = g(:, :, iz, it, ivmu) * (kperp2(:, :, ia, iz) / k2max)**2
               end do
            end do
         end do
         call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%diss_perp, diss_perp_sum, energy_diag%diss_perp_kxkyz)
         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%diss_perp_vmu)
         end if
         ! Calculate numerical dissipation in the zed direction
         ! This is - code_dt * D_zed * delzed(0)**4 / 16 * dgdz * h_conj
         if (hyp_zed) then
            g1 = 0
            call advance_hyper_zed(g, g1)
            g1 = g1 / code_dt
            energy_diag%factor_spec = spec%dens * spec%temp

            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%diss_zed, diss_zed_sum, energy_diag%diss_zed_kxkyz)
            if (write_energy_vmu) then
               call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%diss_zed_vmu)
            end if
         end if
         !Calculate numerical dissipation in the parallel velocity
         ! This is - code_dt * D_vpa * delvpa(0)**4 / 16 * dgvpa * h_conj
         if (hyp_vpa) then

            g1 = 0
            call advance_hyper_vpa(g, g1)
            g1 = g1 / code_dt
            energy_diag%factor_spec = spec%dens * spec%temp

            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%diss_vpa, diss_vpa_sum, energy_diag%diss_vpa_kxkyz)
            if (write_energy_vmu) then
               call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%diss_vpa_vmu)
            end if
         end if

         ! Calculate drive

         g1 = 0.
         call advance_wstar_explicit(phi, g1)
         g1 = g1 / code_dt
         energy_diag%factor_spec = spec%dens * spec%temp

         call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%drive_term, drive_sum, energy_diag%drive_kxkyz)
         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%drive_vmu)
         end if

         ! Calculate drifts

         g1 = 0
         call advance_wdriftx_explicit(g, phi, g1)
         call advance_wdrifty_explicit(g, phi, g1)
         g1 = g1 / code_dt
         energy_diag%factor_spec = spec%dens * spec%temp

         call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%drifts_term, drifts_sum, energy_diag%drifts_kxkyz)
         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%drifts_vmu)
         end if

         ! Calculate streaming

         g1 = 0
         call advance_parallel_streaming_explicit(g, phi, g1)
         g1 = g1 / code_dt
         energy_diag%factor_spec = spec%dens * spec%temp

         call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%streaming_term, streaming_sum, energy_diag%streaming_kxkyz)
         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%streaming_vmu)
         end if

         ! Calculate mirror

         g1 = 0
         call advance_mirror_explicit(g, g1)
         g1 = g1 / code_dt
         energy_diag%factor_spec = spec%dens * spec%temp

         call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%mirror_term, mirror_sum, energy_diag%mirror_kxkyz)
         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%mirror_vmu)
         end if

         !Calculate nonlinearity

         if (nonlinear) then
            g1 = 0
            call advance_ExB_nonlinearity(g, g1, restart_time_step, istep)
            g1 = g1 / code_dt
            energy_diag%factor_spec = spec%dens * spec%temp

            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%nonlinear_term, nonlinear_sum, energy_diag%nonlinear_kxkyz)
            if (write_energy_vmu) then
               call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%nonlinear_vmu)
            end if
         end if
      end if

      total_sum = diss_perp_sum + diss_zed_sum + diss_vpa_sum + drive_sum + drifts_sum + streaming_sum + mirror_sum + nonlinear_sum

      if (proc0) then
      write (energy_unit, '(12e20.8E3)') code_time, energy_sum, dedt_sum, diss_perp_sum, diss_zed_sum, diss_vpa_sum, drive_sum, drifts_sum, streaming_sum, mirror_sum &
            , nonlinear_sum, total_sum
         call flush (energy_unit)
      end if

   end subroutine get_free_energy

   subroutine finish_energy_diagnostic
      implicit none

      if (allocated(energy_diag%energy_total)) deallocate (energy_diag%energy_total)
      if (allocated(energy_diag%dedt_total)) deallocate (energy_diag%dedt_total)
      if (allocated(energy_diag%drive_term)) deallocate (energy_diag%drive_term)
      if (allocated(energy_diag%diss_perp)) deallocate (energy_diag%diss_perp)
      if (allocated(energy_diag%diss_zed)) deallocate (energy_diag%diss_zed)
      if (allocated(energy_diag%diss_vpa)) deallocate (energy_diag%diss_vpa)
      if (allocated(energy_diag%drifts_term)) deallocate (energy_diag%drifts_term)
      if (allocated(energy_diag%streaming_term)) deallocate (energy_diag%streaming_term)
      if (allocated(energy_diag%nonlinear_term)) deallocate (energy_diag%nonlinear_term)
      if (allocated(energy_diag%mirror_term)) deallocate (energy_diag%mirror_term)
      if (allocated(energy_diag%weights_energy)) deallocate (energy_diag%weights_energy)
      if (allocated(energy_diag%factor_spec)) deallocate (energy_diag%factor_spec)

      if (allocated(energy_diag%free_energy_kxkyz)) deallocate (energy_diag%free_energy_kxkyz)
      if (allocated(energy_diag%dedt_kxkyz)) deallocate (energy_diag%dedt_kxkyz)
      if (allocated(energy_diag%drive_kxkyz)) deallocate (energy_diag%drive_kxkyz)
      if (allocated(energy_diag%diss_perp_kxkyz)) deallocate (energy_diag%diss_perp_kxkyz)
      if (allocated(energy_diag%diss_zed_kxkyz)) deallocate (energy_diag%diss_zed_kxkyz)
      if (allocated(energy_diag%diss_vpa_kxkyz)) deallocate (energy_diag%diss_vpa_kxkyz)
      if (allocated(energy_diag%drifts_kxkyz)) deallocate (energy_diag%drifts_kxkyz)
      if (allocated(energy_diag%streaming_kxkyz)) deallocate (energy_diag%streaming_kxkyz)
      if (allocated(energy_diag%nonlinear_kxkyz)) deallocate (energy_diag%nonlinear_kxkyz)
      if (allocated(energy_diag%mirror_kxkyz)) deallocate (energy_diag%mirror_kxkyz)

      if (allocated(energy_diag%energy_total_vmu)) deallocate (energy_diag%energy_total_vmu)

      if (allocated(energy_diag%free_energy_vmu)) deallocate (energy_diag%free_energy_vmu)
      if (allocated(energy_diag%dedt_vmu)) deallocate (energy_diag%dedt_vmu)
      if (allocated(energy_diag%drive_vmu)) deallocate (energy_diag%drive_vmu)
      if (allocated(energy_diag%diss_perp_vmu)) deallocate (energy_diag%diss_perp_vmu)
      if (allocated(energy_diag%diss_zed_vmu)) deallocate (energy_diag%diss_zed_vmu)
      if (allocated(energy_diag%diss_vpa_vmu)) deallocate (energy_diag%diss_vpa_vmu)
      if (allocated(energy_diag%drifts_vmu)) deallocate (energy_diag%drifts_vmu)
      if (allocated(energy_diag%streaming_vmu)) deallocate (energy_diag%streaming_vmu)
      if (allocated(energy_diag%nonlinear_vmu)) deallocate (energy_diag%nonlinear_vmu)
      if (allocated(energy_diag%mirror_vmu)) deallocate (energy_diag%mirror_vmu)
      if (allocated(energy_diag%velocity_integral1)) deallocate (energy_diag%velocity_integral1)

   end subroutine finish_energy_diagnostic

end module energy_diagnostic
