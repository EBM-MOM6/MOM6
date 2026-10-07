! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the MOM6 version of the estuary box model parameterization,
!! which is based on the algorithm described in Sun, Q., Whitney, M. M., Bryan, F. O.,
!! & Tseng, Y. (2017). A box model for representing estuarine physical processes in
!! Earth system models. Ocean Modelling, 112, 139-153.
!! https://doi.org/10.1016/j.ocemod.2017.03.004
module MOM_EBM

use MOM_diag_mediator,         only : post_data, register_static_field, diag_ctrl
use MOM_diag_mediator,         only : register_diag_field, time_type
use MOM_error_handler,         only : MOM_error, WARNING, FATAL, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_io,                    only : file_exists, field_exists, field_size, open_file_to_read
use MOM_io,                    only : close_file_to_read, read_variable, stdout, slasher
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, remapping_core_h
use MOM_remapping,             only : remappingSchemesDoc, remappingDefaultScheme
use MOM_verticalGrid,          only : verticalGrid_type

#include <MOM_memory.h>

implicit none ; private

public EBM_init, calculate_EBM, EBM_is_used, post_EBM_diagnostics
public EBM_unit_tests

!> Control structure including parameters for the estuary box model.
type, public :: EBM_cs ; private

  real :: tide_amp !< Averaged tidal amplitude at estuary mouth [m]
  real, allocatable, dimension(:,:) :: W_h !< Estuary head width [m]
  real :: H        !< Estuary averaged depth [m]
  real, allocatable, dimension(:,:) :: a2  !< A constant of tidal diffusion [nondim]
  real, allocatable, dimension(:,:) :: a1  !< A constant of estuarine mixing length [nondim]
  real, allocatable, dimension(:,:) :: h0  !< A constant of ratio of geometry: h_l/H [nondim]
  real :: g        !< Gravitational acceleration [m s-2]
  real :: rho_ref    !< Reference density for linear equation of state [kg m-3]
  real :: beta_S   !< Saline contraction coefficient [ppt-1]
  real :: Sc       !< Schmidt number [nondim]
  logical :: do_exchange !< If true, apply the EBM exchange flow [nondim]
  real, allocatable, dimension(:,:) :: H_est !< Estuary averaged depth [m]
  real, allocatable, dimension(:,:) :: H_U  !< Thickness of EBM discharge into MOM6 (upper) [m]
  real, allocatable, dimension(:,:) :: H_L  !< Thickness of MOM6 discharging to EBM (lower) [m]
  real, allocatable, dimension(:,:,:) :: wf_U !< Weighting function for the upper layer [m-1]
  real, allocatable, dimension(:,:,:) :: wf_L !< Weighting function for the lower layer [m-1]
  integer :: deg              !< Degree of polynomial reconstruction [nondim]
  real    :: H_subroundoff   !< A thickness that is so small that it can be added to a thickness of
                             !! Angstrom or larger without changing it at the bit level [H ~> m or kg m-2].
  type(remapping_CS)                :: remap_CS !< Control structure to hold remapping configuration.
  type(diag_ctrl), pointer :: diag => NULL() !< Structure used to regulate diagnostic output
  integer :: id_ebm_S_upper = -1 !< Diagnostic ID for EBM upper layer salinity
  integer :: id_ebm_S_lower = -1 !< Diagnostic ID for EBM lower layer salinity
  integer :: id_ebm_Q_u = -1     !< Diagnostic ID for EBM upper layer volume flux
  integer :: id_ebm_Q_l = -1     !< Diagnostic ID for EBM lower layer volume flux
  integer :: id_ebm_S_u = -1     !< Diagnostic ID for EBM upper layer outflow salinity
  real, allocatable, dimension(:,:) :: ebm_S_upper !< EBM upper layer salinity [S ~> ppt]
  real, allocatable, dimension(:,:) :: ebm_S_lower !< EBM lower layer salinity [S ~> ppt]
  real, allocatable, dimension(:,:) :: ebm_Q_u     !< EBM upper layer volume flux [m3 s-1]
  real, allocatable, dimension(:,:) :: ebm_Q_l     !< EBM lower layer volume flux [m3 s-1]
  real, allocatable, dimension(:,:) :: ebm_S_u     !< EBM upper layer outflow salinity [S ~> ppt]
  real, allocatable, dimension(:) :: wf_U_work    !< Work array for remapped upper weighting function [nondim]
  real, allocatable, dimension(:) :: wf_L_work    !< Work array for remapped lower weighting function [nondim]

end type EBM_cs

character(len=40) :: mdl = "MOM_EBM"  !< This module's name.

contains

!> Initializes the estuary box model parameterization.
!! Returns .true. if the parameterization is enabled.
logical function EBM_init(Time, param_file, G, GV, diag, CS)

  type(time_type), target, intent(in)    :: Time      !< The current model time
  type(param_file_type),  intent(in)    :: param_file !< Run-time parameter file handle
  type(ocean_grid_type),  intent(in)    :: G          !< The ocean's grid structure
  type(verticalGrid_type), intent(in)   :: GV         !< The ocean's vertical grid structure
  type(diag_ctrl), target, intent(inout) :: diag      !< Structure used to regulate diagnostic output
  type(EBM_cs),           intent(inout) :: CS         !< EBM control structure

  ! local variables
  character(len=80)  :: string              ! Temporary strings
  character(len=200) :: ebm_edits_file      ! Name of the EBM depth edits file
  character(len=200) :: inputdir            ! Path to the input directory
  logical            :: boundary_extrap     ! Controls if boundary extrapolation is used
  logical            :: om4_remap_via_sub_cells ! Use the OM4-era remap_via_sub_cells
  real :: W_h_dflt                           ! Default estuary head width [m]
  real :: a1_dflt                            ! Default mixing length constant [nondim]
  real :: a2_dflt                            ! Default tidal diffusion constant [nondim]
  real :: h0_dflt                            ! Default geometry ratio h_l/H [nondim]
  real, dimension(:), allocatable :: new_depth ! The new values of estuary depth [m]
  real, dimension(:), allocatable :: new_W_h   ! Override values for W_h [m]
  real, dimension(:), allocatable :: new_a1    ! Override values for a1 [nondim]
  real, dimension(:), allocatable :: new_a2    ! Override values for a2 [nondim]
  real, dimension(:), allocatable :: new_h0    ! Override values for h0 [nondim]
  integer, dimension(:), allocatable :: ig, jg ! The global indices of points to modify
  integer :: id                              ! Diagnostic id for static fields
  integer :: i, j, n, ncid, n_edits, i_file, j_file, ndims, sizes(8)

  ! This include declares and sets the variable "version".
# include "version_variable.h"

  call get_param(param_file, mdl, "USE_EBM", EBM_init, default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
       "Estuary box model parameterization", all_default=.not.EBM_init)
  call get_param(param_file, mdl, "USE_EBM", EBM_init, &
                 "If true, enables the estuary box model (EBM) parameterization. ", &
                 default=.false.)

  if (.not. EBM_init) return

  if (.not.GV%Boussinesq) call MOM_error(FATAL, trim(mdl)// &
       ": the EBM parameterization currently requires Boussinesq mode.")

  call get_param(param_file, mdl, "EBM_TIDAL_AMPLITUDE", CS%tide_amp, &
                 "Averaged tidal amplitude at the estuary mouth.", &
                 units="m", default=1.0)

  call get_param(param_file, mdl, "EBM_HEAD_WIDTH", W_h_dflt, &
                 "Estuary head width.", &
                 units="m", default=2000.0)

  call get_param(param_file, mdl, "EBM_DEPTH", CS%H, &
                 "Estuary averaged depth.", &
                 units="m", default=10.0)

  call get_param(param_file, mdl, "EBM_A1", a1_dflt, &
                 "A constant of estuarine mixing length.", &
                 units="nondim", default=0.876)

  call get_param(param_file, mdl, "EBM_A2", a2_dflt, &
                 "A constant of tidal diffusion.", &
                 units="nondim", default=0.0)

  call get_param(param_file, mdl, "EBM_H0", h0_dflt, &
                 "A constant of ratio of geometry: h_l/H.", &
                 units="nondim", default=0.5)

  call get_param(param_file, mdl, "EBM_G", CS%g, &
                 "Gravitational acceleration used in the EBM.", &
                 units="m s-2", default=9.8)

  call get_param(param_file, mdl, "EBM_RHO_REF", CS%rho_ref, &
                 "Reference density for the linear EoS used in the EBM.", &
                 units="kg m-3", default=1000.0)

  call get_param(param_file, mdl, "EBM_BETA_S", CS%beta_S, &
                 "Saline contraction coefficient used in the EBM.", &
                 units="ppt-1", default=7.7e-4)

  call get_param(param_file, mdl, "EBM_SC", CS%Sc, &
                 "Schmidt number used in the EBM.", &
                 units="nondim", default=2.2)

  call get_param(param_file, mdl, "EBM_EXCHANGE", CS%do_exchange, &
                 "If true, apply the EBM exchange flow between the upper "//&
                 "and lower estuary layers.", &
                 default=.true.)

  ! Remapping initialization
  CS%H_subroundoff = GV%H_subroundoff
  call get_param(param_file, mdl, "EBM_BOUNDARY_EXTRAP", boundary_extrap, &
                 "Use boundary extrapolation in EBM remapping.", &
                 default=.false.)
  call get_param(param_file, mdl, "EBM_REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 do_not_log=.true., default=.true.)
  call get_param(param_file, mdl, "EBM_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for the EBM. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
  call initialize_remapping(CS%remap_CS, string, boundary_extrapolation=boundary_extrap, &
                            om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                            check_reconstruction=.false., check_remapping=.false., &
                            h_neglect=CS%H_subroundoff, h_neglect_edge=CS%H_subroundoff)
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)

  ! Allocate and initialize arrays
  allocate(CS%H_est(SZI_(G),SZJ_(G)))
  allocate(CS%H_U(SZI_(G),SZJ_(G)))
  allocate(CS%H_L(SZI_(G),SZJ_(G)))
  allocate(CS%wf_U(SZI_(G),SZJ_(G),3))
  allocate(CS%wf_L(SZI_(G),SZJ_(G),3))
  allocate(CS%W_h(SZI_(G),SZJ_(G)))
  allocate(CS%a1(SZI_(G),SZJ_(G)))
  allocate(CS%a2(SZI_(G),SZJ_(G)))
  allocate(CS%h0(SZI_(G),SZJ_(G)))
  allocate(CS%wf_U_work(GV%ke), source=0.0)
  allocate(CS%wf_L_work(GV%ke), source=0.0)
  CS%H_est(:,:) = CS%H
  CS%W_h(:,:) = W_h_dflt
  CS%a1(:,:) = a1_dflt
  CS%a2(:,:) = a2_dflt
  CS%h0(:,:) = h0_dflt
  CS%H_U(:,:) = 0.0
  CS%H_L(:,:) = 0.0
  CS%wf_U(:,:,:) = 0.0
  CS%wf_L(:,:,:) = 0.0

  ! Apply optional point-by-point overrides of EBM parameters from a file
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "EBM_EDITS_FILE", ebm_edits_file, &
                 "The file from which to read per-column EBM parameter overrides. "//&
                 "The file must be in NetCDF format with scalar variables 'ni' and "//&
                 "'nj' matching the global grid dimensions, and 1-D arrays 'iEdit' "//&
                 "and 'jEdit' (of size nEdits) giving the global i-index and j-index "//&
                 "(both using Python 0-based indexing). Optional 1-D arrays of size "//&
                 "nEdits: 'zEdit' (estuary depth [m]), 'W_hEdit' (head width [m]), "//&
                 "'a1Edit' (mixing length constant), 'a2Edit' (tidal diffusion "//&
                 "constant), 'h0Edit' (geometry ratio h_l/H).", &
                 default="")

  if (len_trim(ebm_edits_file) > 0) then
    ebm_edits_file = trim(inputdir)//trim(ebm_edits_file)
    if (is_root_pe()) then
      if (.not.file_exists(ebm_edits_file, G%Domain)) &
        call MOM_error(FATAL, trim(mdl)//': Unable to find file '//trim(ebm_edits_file))
      call open_file_to_read(ebm_edits_file, ncid)
    else
      ncid = -1
    endif

    ! Read and check the values of ni and nj in the file for consistency with this configuration.
    call read_variable(ebm_edits_file, 'ni', i_file, ncid_in=ncid)
    call read_variable(ebm_edits_file, 'nj', j_file, ncid_in=ncid)
    if (i_file /= G%ieg) call MOM_error(FATAL, trim(mdl)//': Incompatible i-dimension of grid in '//&
                                        trim(ebm_edits_file))
    if (j_file /= G%jeg) call MOM_error(FATAL, trim(mdl)//': Incompatible j-dimension of grid in '//&
                                        trim(ebm_edits_file))

    ! Get nEdits from iEdit
    call field_size(ebm_edits_file, 'iEdit', sizes, ndims=ndims, ncid_in=ncid)
    if (ndims /= 1) call MOM_error(FATAL, "The variable iEdit has an "//&
              "unexpected number of dimensions in "//trim(ebm_edits_file))
    n_edits = sizes(1)
    allocate(ig(n_edits), jg(n_edits))

    ! Read iEdit and jEdit (required)
    call read_variable(ebm_edits_file, 'iEdit', ig, ncid_in=ncid)
    call read_variable(ebm_edits_file, 'jEdit', jg, ncid_in=ncid)

    ! Read optional override arrays
    if (field_exists(ebm_edits_file, 'zEdit')) then
      allocate(new_depth(n_edits))
      call read_variable(ebm_edits_file, 'zEdit', new_depth, ncid_in=ncid)
    endif
    if (field_exists(ebm_edits_file, 'W_hEdit')) then
      allocate(new_W_h(n_edits))
      call read_variable(ebm_edits_file, 'W_hEdit', new_W_h, ncid_in=ncid)
    endif
    if (field_exists(ebm_edits_file, 'a1Edit')) then
      allocate(new_a1(n_edits))
      call read_variable(ebm_edits_file, 'a1Edit', new_a1, ncid_in=ncid)
    endif
    if (field_exists(ebm_edits_file, 'a2Edit')) then
      allocate(new_a2(n_edits))
      call read_variable(ebm_edits_file, 'a2Edit', new_a2, ncid_in=ncid)
    endif
    if (field_exists(ebm_edits_file, 'h0Edit')) then
      allocate(new_h0(n_edits))
      call read_variable(ebm_edits_file, 'h0Edit', new_h0, ncid_in=ncid)
    endif
    call close_file_to_read(ncid, ebm_edits_file)

    do n = 1, n_edits
      i = ig(n) - G%idg_offset + 1 ! +1 for python indexing
      j = jg(n) - G%jdg_offset + 1
      if (i>=G%isc .and. i<=G%iec .and. j>=G%jsc .and. j<=G%jec) then
        if (allocated(new_depth)) then
          write(stdout,'(a,3i5,f8.2,a,f8.2,2i4)') &
            'EBM depth edit: ', n, ig(n), jg(n), CS%H_est(i,j), '->', abs(new_depth(n)), i, j
          CS%H_est(i,j) = abs(new_depth(n))
        endif
        if (allocated(new_W_h)) then
          write(stdout,'(a,3i5,ES12.4,a,ES12.4)') &
            'EBM W_h edit: ', n, ig(n), jg(n), CS%W_h(i,j), '->', new_W_h(n)
          CS%W_h(i,j) = new_W_h(n)
        endif
        if (allocated(new_a1)) then
          write(stdout,'(a,3i5,ES12.4,a,ES12.4)') &
            'EBM a1 edit: ', n, ig(n), jg(n), CS%a1(i,j), '->', new_a1(n)
          CS%a1(i,j) = new_a1(n)
        endif
        if (allocated(new_a2)) then
          write(stdout,'(a,3i5,ES12.4,a,ES12.4)') &
            'EBM a2 edit: ', n, ig(n), jg(n), CS%a2(i,j), '->', new_a2(n)
          CS%a2(i,j) = new_a2(n)
        endif
        if (allocated(new_h0)) then
          write(stdout,'(a,3i5,ES12.4,a,ES12.4)') &
            'EBM h0 edit: ', n, ig(n), jg(n), CS%h0(i,j), '->', new_h0(n)
          CS%h0(i,j) = new_h0(n)
        endif
      endif
    enddo

    deallocate(ig, jg)
    if (allocated(new_depth)) deallocate(new_depth)
    if (allocated(new_W_h)) deallocate(new_W_h)
    if (allocated(new_a1)) deallocate(new_a1)
    if (allocated(new_a2)) deallocate(new_a2)
    if (allocated(new_h0)) deallocate(new_h0)
  endif

  ! TODO: revisit H_U and H_L later with limiters?
  do j=G%jsc,G%jec
    do i=G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        if (CS%H_est(i,j) <= G%bathyT(i,j)) then
          CS%H_U(i,j) = CS%H_est(i,j) * 0.5
          CS%H_L(i,j) = CS%H_est(i,j) * 0.5
        else
          CS%H_U(i,j) = G%bathyT(i,j) * 0.5
          CS%H_L(i,j) = G%bathyT(i,j) * 0.5
        endif
        CS%wf_U(i,j,1) = 1.0/CS%H_U(i,j)
        CS%wf_L(i,j,2) = 1.0/CS%H_L(i,j)
        CS%wf_U(i,j,3) = 0.0
        CS%wf_L(i,j,3) = 0.0
      endif
    enddo
  enddo

  ! Post the static estuary depth fields
  id = register_static_field('ocean_model', 'ebm_depth', diag%axesT1, &
        'Estuary averaged depth used by the EBM', 'm')
  if (id > 0) call post_data(id, CS%H_est, diag, .true.)

  id = register_static_field('ocean_model', 'ebm_H_U', diag%axesT1, &
        'EBM upper layer thickness', 'm')
  if (id > 0) call post_data(id, CS%H_U, diag, .true.)

  id = register_static_field('ocean_model', 'ebm_H_L', diag%axesT1, &
        'EBM lower layer thickness', 'm')
  if (id > 0) call post_data(id, CS%H_L, diag, .true.)

  id = register_static_field('ocean_model', 'ebm_W_h', diag%axesT1, &
        'Estuary head width used by the EBM', 'm')
  if (id > 0) call post_data(id, CS%W_h, diag, .true.)

  id = register_static_field('ocean_model', 'ebm_a1', diag%axesT1, &
        'Estuarine mixing length constant used by the EBM', 'nondim')
  if (id > 0) call post_data(id, CS%a1, diag, .true.)

  id = register_static_field('ocean_model', 'ebm_a2', diag%axesT1, &
        'Tidal diffusion constant used by the EBM', 'nondim')
  if (id > 0) call post_data(id, CS%a2, diag, .true.)

  id = register_static_field('ocean_model', 'ebm_h0', diag%axesT1, &
        'Geometry ratio h_l/H used by the EBM', 'nondim')
  if (id > 0) call post_data(id, CS%h0, diag, .true.)

  ! Register time-dependent EBM diagnostics
  CS%diag => diag
  CS%id_ebm_S_upper = register_diag_field('ocean_model', 'ebm_S_upper', diag%axesT1, &
        Time, 'EBM upper layer salinity', 'ppt')
  if (CS%id_ebm_S_upper > 0) then
    allocate(CS%ebm_S_upper(SZI_(G),SZJ_(G)), source=0.0)
  endif
  CS%id_ebm_S_lower = register_diag_field('ocean_model', 'ebm_S_lower', diag%axesT1, &
        Time, 'EBM lower layer salinity', 'ppt')
  if (CS%id_ebm_S_lower > 0) then
    allocate(CS%ebm_S_lower(SZI_(G),SZJ_(G)), source=0.0)
  endif
  CS%id_ebm_Q_u = register_diag_field('ocean_model', 'ebm_Q_u', diag%axesT1, &
        Time, 'EBM upper layer volume flux', 'm3 s-1')
  if (CS%id_ebm_Q_u > 0) then
    allocate(CS%ebm_Q_u(SZI_(G),SZJ_(G)), source=0.0)
  endif
  CS%id_ebm_Q_l = register_diag_field('ocean_model', 'ebm_Q_l', diag%axesT1, &
        Time, 'EBM lower layer volume flux', 'm3 s-1')
  if (CS%id_ebm_Q_l > 0) then
    allocate(CS%ebm_Q_l(SZI_(G),SZJ_(G)), source=0.0)
  endif
  CS%id_ebm_S_u = register_diag_field('ocean_model', 'ebm_S_u', diag%axesT1, &
        Time, 'EBM upper layer outflow salinity', 'ppt')
  if (CS%id_ebm_S_u > 0) then
    allocate(CS%ebm_S_u(SZI_(G),SZJ_(G)), source=0.0)
  endif

end function EBM_init

!> Calculates estuary box model exchange and distributes river runoff using
!! weighting functions, updating layer thickness, temperature, salinity, and net mass flux.
!! The estuary exchange fluxes (Q_u, Q_l, S_u) are computed by estuary_box_model
!! using the EBM parameters in CS together with the provided river discharge and
!! lower-layer salinity.
subroutine calculate_EBM(CS, G, i, j, Idt, lrunoff, EnthalpyConst, netMassIn, T2d_col, S_col, h2d_col)

  type(EBM_cs), intent(inout)   :: CS            !< EBM control structure
  type(ocean_grid_type), intent(in) :: G          !< The ocean's grid structure
  integer,      intent(in)      :: i             !< i-index of the current column [nondim]
  integer,      intent(in)      :: j             !< j-index of the current column [nondim]
  real,         intent(in)      :: Idt           !< The inverse of the timestep [T-1 ~> s-1]
  real,         intent(in)      :: lrunoff       !< River runoff for this column [H ~> m or kg m-2]
  real,         intent(in)    :: EnthalpyConst !< Enthalpy constant [nondim]
  real,         intent(inout) :: netMassIn     !< Net mass entering this column [H ~> m or kg m-2]
  real,         intent(inout) :: T2d_col(:)    !< Temperature in the column [C ~> degC]
  real,         intent(inout) :: S_col(:)      !< Salinity in the column [S ~> ppt]
  real,         intent(inout) :: h2d_col(:)    !< Layer thickness in the column [H ~> m or kg m-2]

  ! local variables
  real :: Q_u        !< EBM upper layer volume flux [m3 s-1]
  real :: Q_l        !< EBM lower layer volume flux [m3 s-1]
  real :: S_u        !< EBM upper layer salinity [S ~> ppt]
  real :: dThickness !< Change in layer thickness [H ~> m or kg m-2]
  real :: dTemp      !< Integrated change in layer temperature [C H ~> degC m or degC kg m-2]
  real :: dSalt      !< Integrated change in layer salinity [S H ~> ppt m or ppt kg m-2]
  real :: sum_dThickness !< Accumulated thickness change over the column [H ~> m or kg m-2]
  real :: sum_dTemp      !< Accumulated temperature content change over the column [C H ~> degC m]
  real :: sum_dSalt      !< Accumulated salt content change over the column [S H ~> ppt m]
  real :: Temp_in    !< Temperature of the incoming mass flux [C ~> degC]
  real :: hOld       !< Original layer thickness before update [H ~> m or kg m-2]
  real :: Ithickness !< Inverse of the updated layer thickness [H-1 ~> m-1 or m2 kg-1]
  integer, parameter :: nz_ebm = 3  !< Number of layers in the EBM grid + 1  [nondim]
                                    !! The extra layer is needed to avoid propagating the lower layer value.
  integer :: nz      !< Number of layers in the native grid  [nondim]
  integer :: k       !< Layer index [nondim]

  real, dimension(nz_ebm) :: dz_ebm !< EBM layer thicknesses [m]
  real, dimension(nz_ebm) :: ebm_wf_U  !< EBM weighting function upper layer (pointwise) [m-1]
  real, dimension(nz_ebm) :: ebm_wf_L  !< EBM weighting function lower layer (pointwise) [m-1]
  real, dimension(nz_ebm) :: ebm_S     !<  Salinity in the column on the EBM grid [S ~> ppt]
  real, dimension(nz_ebm) :: ebm_T     !<  Temperature in the column on the EBM grid [C ~> degC]
  real :: integral_U  !< Sum of wf_U over the column, used for normalization [nondim]
  real :: integral_L  !< Sum of wf_L over the column, used for normalization [nondim]

  ! Vertical grid EBM
  dz_ebm(1) = CS%H_U(i,j)
  dz_ebm(2) = CS%H_L(i,j)
  dz_ebm(3) = 1.0E4 ! avoid propagating the lower layer value.

  ! Weighting functions EBM
  ebm_wf_U(1) = CS%wf_U(i,j,1)
  ebm_wf_U(2) = CS%wf_U(i,j,2)
  ebm_wf_L(1) = CS%wf_L(i,j,1)
  ebm_wf_L(2) = CS%wf_L(i,j,2)
  ebm_wf_U(3) = CS%wf_U(i,j,3)
  ebm_wf_L(3) = CS%wf_L(i,j,3)
  ! Weighting functions native

  ! Number of layers in the native grid
  nz = size(T2d_col)

  ! Initialize work arrays
  CS%wf_U_work(:) = 0.0
  CS%wf_L_work(:) = 0.0

  ! 1) Distribute river runoff over upper layer

  ! 1a) remap weighting functions
  call remapping_core_h(CS%remap_cs, nz_ebm, dz_ebm(:), ebm_wf_U(:), nz, h2d_col(:), CS%wf_U_work(:))
  call remapping_core_h(CS%remap_cs, nz_ebm, dz_ebm(:), ebm_wf_L(:), nz, h2d_col(:), CS%wf_L_work(:))

  ! Multiply weights by layer thickness to make them nondim
  do k = 1, nz
    CS%wf_U_work(k) = CS%wf_U_work(k) * h2d_col(k)
    CS%wf_L_work(k) = CS%wf_L_work(k) * h2d_col(k)
  enddo

  ! check that the SUM(wf_U) is 1.
  call check_wf_integral(G, i, j, nz_ebm, dz_ebm, ebm_wf_U, nz, h2d_col, CS%wf_U_work, 'wf_U')
  call check_wf_integral(G, i, j, nz_ebm, dz_ebm, ebm_wf_L, nz, h2d_col, CS%wf_L_work, 'wf_L')

  ! Normalize weighting functions to sum exactly to 1 to prevent
  ! floating-point drift from accumulating into conservation errors.
  integral_U = 0.0
  integral_L = 0.0
  do k = 1, nz
    integral_U = integral_U + CS%wf_U_work(k)
    integral_L = integral_L + CS%wf_L_work(k)
  enddo
  do k = 1, nz
    CS%wf_U_work(k) = CS%wf_U_work(k) / integral_U
    CS%wf_L_work(k) = CS%wf_L_work(k) / integral_L
  enddo

  ! Compute salinity on EBM grid [S ~> ppt]
  call remapping_core_h(CS%remap_cs, nz, h2d_col(:), S_col(:), nz_ebm, dz_ebm(:), ebm_S(:))
  call remapping_core_h(CS%remap_cs, nz, h2d_col(:), T2d_col(:), nz_ebm, dz_ebm(:), ebm_T(:))

  ! Store EBM layer salinities for diagnostics
  if (CS%id_ebm_S_upper > 0) CS%ebm_S_upper(i,j) = ebm_S(1)
  if (CS%id_ebm_S_lower > 0) CS%ebm_S_lower(i,j) = ebm_S(2)

  ! initialize accumulators
  sum_dThickness = 0.0

  ! Distribute river runoff over upper layer
  do k = 1, nz
    dThickness = lrunoff * CS%wf_U_work(k)
    sum_dThickness = sum_dThickness + dThickness
    dTemp = 0.
    dSalt = 0.

    netMassIn = netMassIn - dThickness
    Temp_in  = T2d_col(k)

    ! GMM, TODO:
    ! We will need to change how enthalpy is done with the EBM.
    ! EnthalpyConst should be 1 and we need to remove the
    ! heat_content_lrunoff from the coupler elsewhere.

    dTemp = dTemp + dThickness * Temp_in * EnthalpyConst

    hOld = h2d_col(k)
    h2d_col(k) = h2d_col(k) + dThickness
    if (h2d_col(k) > 0.0) then
      Ithickness = 1.0 / h2d_col(k)
      if (dThickness /= 0. .or. dTemp /= 0.) T2d_col(k) = (hOld * T2d_col(k) + dTemp) * Ithickness
      if (dThickness /= 0. .or. dSalt /= 0.) S_col(k)   = (hOld * S_col(k)   + dSalt) * Ithickness
    end if
  end do

  ! Check conservation of river input distribution
  if (abs((sum_dThickness) - lrunoff) > 1.0e-10) then
    if (is_root_pe()) then
      write(stdout,'(A)') "=== MOM_EBM river input conservation FAILURE ==="
      write(stdout,'(A,I6,A,I6)') "  i          : ", i, "  j : ", j
      write(stdout,'(A,F12.4,A,F12.4)') "  lon        : ", G%geoLonT(i,j), &
                                         "  lat : ", G%geoLatT(i,j)
      write(stdout,'(A,ES15.8)') "  sum_dThickness - lrunoff : ", (sum_dThickness-lrunoff)
      write(stdout,'(A,ES15.8)') "  sum_dThickness           : ", sum_dThickness
      write(stdout,'(A,ES15.8)') "  lrunoff           : ", lrunoff
      do k = 1, nz
        write(stdout,'(A,I4,A,2ES15.8)') "    k=", k, " : ", CS%wf_U_work(k), CS%wf_L_work(k)
      enddo
      call MOM_error(FATAL, "MOM_EBM calculate_EBM: river input distribution is not conservative.")
    endif
  endif

  ! 2) Exchange
  if (CS%do_exchange) then

    ! GMM, todo
    ! EBM parameterization currently requires Boussinesq mode
    ! Runoff is in m ==> converting to m3 s-1
    call estuary_box_model(CS, G, i, j, lrunoff*G%areaT(i,j)*Idt, ebm_S(2), Q_u, Q_l, S_u)
    ! GMM, TODO: S_u is masked out if lrunoff=0.

    ! Store EBM exchange diagnostics
    if (CS%id_ebm_Q_u > 0) CS%ebm_Q_u(i,j) = Q_u
    if (CS%id_ebm_Q_l > 0) CS%ebm_Q_l(i,j) = Q_l
    if (CS%id_ebm_S_u > 0) CS%ebm_S_u(i,j) = S_u

    ! convert Q_U and Q_L back to m
    Q_u = Q_u/(G%areaT(i,j)*Idt)
    Q_l = Q_l/(G%areaT(i,j)*Idt)

    ! Apply exchange flow

    ! initialize accumulators
    sum_dThickness = 0.0
    sum_dTemp = 0.0
    sum_dSalt = 0.0

    ! 2a) Upper and Lower layer combined
    do k = 1, nz
      dThickness = Q_l * (CS%wf_L_work(k) - CS%wf_U_work(k))
      dTemp = dThickness * ebm_T(2)
      dSalt = dThickness * ebm_S(2)
      sum_dThickness = sum_dThickness + dThickness
      sum_dTemp = sum_dTemp + dTemp
      sum_dSalt = sum_dSalt + dSalt
      hOld = h2d_col(k)
      h2d_col(k) = h2d_col(k) + dThickness
      if (h2d_col(k) > 0.0) then
        Ithickness = 1.0 / h2d_col(k)
        if (dThickness /= 0. .or. dTemp /= 0.) T2d_col(k) = (hOld * T2d_col(k) + dTemp) * Ithickness
        if (dThickness /= 0. .or. dSalt /= 0.) S_col(k)   = (hOld * S_col(k)   + dSalt) * Ithickness
      end if
    end do

    ! Check conservation of exchange flow
    if (abs(sum_dThickness) > 1.0e-10 .or. abs(sum_dTemp) > 1.0e-10 .or. abs(sum_dSalt) > 1.0e-10) then
      if (is_root_pe()) then
        write(stdout,'(A)') "=== MOM_EBM exchange conservation FAILURE ==="
        write(stdout,'(A,I6,A,I6)') "  i          : ", i, "  j : ", j
        write(stdout,'(A,F12.4,A,F12.4)') "  lon        : ", G%geoLonT(i,j), &
                                           "  lat : ", G%geoLatT(i,j)
        write(stdout,'(A,ES15.8)') "  sum_dThickness : ", sum_dThickness
        write(stdout,'(A,ES15.8)') "  sum_dTemp      : ", sum_dTemp
        write(stdout,'(A,ES15.8)') "  sum_dSalt      : ", sum_dSalt
        write(stdout,'(A,ES15.8)') "  Q_l            : ", Q_l
        do k = 1, nz
          write(stdout,'(A,I4,A,2ES15.8)') "    k=", k, " : ", CS%wf_U_work(k), CS%wf_L_work(k)
        enddo
        call MOM_error(FATAL, "MOM_EBM calculate_EBM: exchange flow is not conservative.")
      endif
    endif

  endif ! do_exchange

end subroutine calculate_EBM

!> Check that SUM(wf(k), k=1..nz) = 1.
!! After remapping and multiplication by layer thickness, the weighting function
!! must sum to 1 to ensure conservation. If not, detailed debug info is printed
!! on root PE and a FATAL error is issued.
subroutine check_wf_integral(G, i, j, nz_ebm, dz_ebm, ebm_wf, nz, h, wf, varname)
  type(ocean_grid_type), intent(in) :: G       !< The ocean's grid structure
  integer,          intent(in) :: i            !< i-index of the current column [nondim]
  integer,          intent(in) :: j            !< j-index of the current column [nondim]
  integer,          intent(in) :: nz_ebm       !< Number of EBM layers [nondim]
  real,             intent(in) :: dz_ebm(:)    !< EBM layer thicknesses [m]
  real,             intent(in) :: ebm_wf(:)    !< Weighting function on EBM grid [m-1]
  integer,          intent(in) :: nz           !< Number of native layers [nondim]
  real,             intent(in) :: h(:)         !< Native layer thicknesses [m]
  real,             intent(in) :: wf(:)        !< Remapped weighting function on native grid [nondim]
  character(len=*), intent(in) :: varname      !< Name of the variable being checked

  ! local variables
  real, parameter :: tol = 1.0e-10 ! Tolerance for the integral check [nondim]
  real    :: integral ! Integral of wf over the column [nondim]
  integer :: k        ! Layer index [nondim]
  character(len=256) :: mesg ! Error message string

  integral = 0.0
  do k = 1, nz
    integral = integral + wf(k)
  enddo

  if (abs(integral - 1.0) > tol) then
    if (is_root_pe()) then
      write(stdout,'(A)') "=== MOM_EBM check_wf_integral FAILURE ==="
      write(stdout,'(A,A)')       "  variable   : ", trim(varname)
      write(stdout,'(A,I6,A,I6)') "  i          : ", i, "  j : ", j
      write(stdout,'(A,F12.4,A,F12.4)') "  lon        : ", G%geoLonT(i,j), &
                                         "  lat : ", G%geoLatT(i,j)
      write(stdout,'(A,I4)')      "  nz_ebm     : ", nz_ebm
      write(stdout,'(A)')         "  dz_ebm     : "
      do k = 1, nz_ebm
        write(stdout,'(A,I4,A,ES15.8)') "    k=", k, " : ", dz_ebm(k)
      enddo
      write(stdout,'(A)')         "  ebm_wf     : "
      do k = 1, nz_ebm
        write(stdout,'(A,I4,A,ES15.8)') "    k=", k, " : ", ebm_wf(k)
      enddo
      write(stdout,'(A,I4)')      "  nz         : ", nz
      write(stdout,'(A)')         "  h          : "
      do k = 1, nz
        write(stdout,'(A,I4,A,ES15.8)') "    k=", k, " : ", h(k)
      enddo
      write(stdout,'(A)')         "  wf         : "
      do k = 1, nz
        write(stdout,'(A,I4,A,ES15.8)') "    k=", k, " : ", wf(k)
      enddo
      write(mesg, '("SUM(wf) = ",ES15.8," differs from 1 by ",ES15.8," (tol = ",ES15.8,").")') &
           integral, abs(integral - 1.0), tol
      call MOM_error(FATAL, "MOM_EBM check_wf_integral: variable '"//trim(varname)//"': "//&
           trim(mesg)//" The remapped weighting function is not conservative.")
    endif
  endif

end subroutine check_wf_integral


!> Post time-dependent EBM diagnostics. This should be called after the
!! column loop that calls calculate_EBM has completed.
subroutine post_EBM_diagnostics(CS)
  type(EBM_cs), intent(in) :: CS !< EBM control structure

  if (CS%id_ebm_S_upper > 0) call post_data(CS%id_ebm_S_upper, CS%ebm_S_upper, CS%diag)
  if (CS%id_ebm_S_lower > 0) call post_data(CS%id_ebm_S_lower, CS%ebm_S_lower, CS%diag)
  if (CS%id_ebm_Q_u > 0) call post_data(CS%id_ebm_Q_u, CS%ebm_Q_u, CS%diag)
  if (CS%id_ebm_Q_l > 0) call post_data(CS%id_ebm_Q_l, CS%ebm_Q_l, CS%diag)
  if (CS%id_ebm_S_u > 0) call post_data(CS%id_ebm_S_u, CS%ebm_S_u, CS%diag)

end subroutine post_EBM_diagnostics

!> Reads the parameter "USE_EBM" and returns state.
!! This function allows other modules to know whether this parameterization will
!! be used without needing to duplicate the log entry.
logical function EBM_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  call get_param(param_file, mdl, "USE_EBM", EBM_is_used, &
                 default=.false., do_not_log=.true.)

end function EBM_is_used

!> Calculate the estuary box model (EBM) exchange.
!!
!! The EBM is assumed steady state, flat bottom and flat surface with a straight
!! channel and rectangular cross-section. It is built on three global conservation
!! laws: water mass, water volume, and potential energy conservation.
!!
!! Estuary Box geometry:
!!              ^ z
!!              |____________________         _______
!!              |                    |               |
!!  Q_u,S_u <--+--  upper layer   <-+-- Q_r         |
!!              |--------------------|        -+-    | H
!!     Q_l,S_l -+->  lower layer     |         | h_l |
!!           ___|____________________|        _|_____|
!!          x   0                   -LE
!!
!! Note: the negative lower layer volume flux Q_l leaves the ocean; the positive
!! upper layer volume flux Q_u flows into the ocean.
!!
!! Reference: Sun, Q., Whitney, M. M., Bryan, F. O., & Tseng, Y. (2017).
!! A box model for representing estuarine physical processes in
!! Earth system models. Ocean Modelling, 112, 139-153.
!! https://doi.org/10.1016/j.ocemod.2017.03.004
subroutine estuary_box_model(CS, G, ig, jg, Q_r, S_l, Q_u, Q_l, S_u)

  type(EBM_cs), intent(in)  :: CS   !< EBM control structure
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure
  integer,      intent(in)  :: ig   !< i-index of the current column [nondim]
  integer,      intent(in)  :: jg   !< j-index of the current column [nondim]
  real,         intent(in)  :: Q_r  !< River discharge [m3 s-1]
  real,         intent(in)  :: S_l  !< Salinity at estuary lower layer [ppt]
  real,         intent(out) :: Q_u  !< Upper layer volume flux [m3 s-1]
  real,         intent(out) :: Q_l  !< Lower layer volume flux [m3 s-1]
  real,         intent(out) :: S_u  !< Salinity at estuary upper layer [ppt]

  ! local variables
  real :: ERR_EBM     ! Closure error of the EBM potential energy budget [kg m s-3]
  real :: rho_r       ! River water density [kg m-3]
  real :: rho_l       ! Lower layer inflow density [kg m-3]
  real :: rho_u       ! Upper layer outflow density [kg m-3]
  real :: u_t         ! Tidal current amplitude near bottom [m s-1]
  real :: u_r         ! Riverine velocity at head of EBM [m s-1]
  real :: u_l         ! Estuarine inflow velocity at mouth of EBM [m s-1]
  real :: u_u         ! Estuarine outflow velocity at mouth of EBM [m s-1]
  real :: u_bar       ! Layer-averaged net velocity in EBM [m s-1]
  real :: c_wave      ! Densimetric wave phase speed [m s-1]
  real :: ur0         ! Densimetric riverine Froude number [nondim]
  real :: ut0         ! Densimetric tidal current Froude number [nondim]
  real :: ul0         ! Dimensionless lower layer inflow Froude number [nondim]
  real :: uu0         ! Dimensionless upper layer outflow Froude number [nondim]
  real :: R0          ! Layer densimetric riverine Froude number [nondim]
  real :: T0          ! Layer densimetric tidal Froude number [nondim]
  real :: h_l         ! Lower layer water depth of EBM [m]
  real :: a, b, c, d  ! Normalized coefficients of the cubic equation for ul0 [nondim]
  real :: AD          ! Advective potential energy flux term for EBM closure check [kg m s-3]
  real :: HD          ! Horizontal diffusive potential energy flux term [kg m s-3]
  real :: VD          ! Vertical diffusive potential energy flux term [kg m s-3]
  real :: LF          ! Lateral friction potential energy flux term [kg m s-3]
  real, dimension(3,2) :: roots   ! Roots of the 3rd-order polynomial [nondim];
                                  !! roots(:,1) = real part, roots(:,2) = imaginary part
  integer, dimension(3) :: mask   ! Selection mask for physically valid roots [nondim]
  integer :: i, n                 ! Loop and counter indices [nondim]

  real, parameter :: PI = 4.0 * atan(1.0)  !< Ratio of circumference to diameter [nondim]

  ! Skip computations if S_l <= 0, using limit as S_l->0
  if (S_l <= 0.0) then
    Q_u = Q_r
    Q_l = 0.0
    S_u = 0.0
    return
  end if

  ! Water densities
  rho_r = CS%rho_ref
  rho_l = CS%rho_ref * (1.0 + CS%beta_S * S_l)

  ! River and tidal velocities
  u_t    = -CS%tide_amp * sqrt(CS%g / CS%H_est(ig,jg))        ! Tidal velocity (toward river)
  u_r    = Q_r / (CS%W_h(ig,jg) * CS%H_est(ig,jg) * (1.0 - CS%h0(ig,jg)))  ! Riverine velocity at head of upper layer
  c_wave = sqrt(CS%beta_S * S_l * CS%g * CS%H_est(ig,jg))     ! Densimetric wave phase speed

  ! Dimensionless parameters
  ur0 = u_r / c_wave
  ut0 = u_t / c_wave
  R0  = ur0 * (1.0 - CS%h0(ig,jg))
  T0  = ut0 * (1.0 - CS%h0(ig,jg)) / PI

  ! Coefficients of the cubic equation for dimensionless lower layer inflow (ul0)
  a = -CS%h0(ig,jg)**3.0

  b = 2.0 * CS%h0(ig,jg)**2.0 * ((2.0 - CS%h0(ig,jg)) * R0 - CS%a2(ig,jg) * T0)

  c = 0.096 * CS%a1(ig,jg) * CS%h0(ig,jg) * (CS%Sc**2.0 * R0)**(-1.0/3.0) * R0 &
    - CS%h0(ig,jg) * ((2.0 - CS%h0(ig,jg)) * R0 * (R0 - 2.0 * CS%a2(ig,jg) * T0) &
    + CS%a2(ig,jg)**2.0 * T0**2.0)

  d = -0.048 * CS%a1(ig,jg) * (CS%Sc**2.0 * R0)**(-1.0/3.0) &
    * R0 * (R0 - 2.0 * CS%a2(ig,jg) * T0)

  call cubsolve(b/a, c/a, d/a, roots)

  ! Select the physically valid root: real, negative lower layer inflow
  mask = 0
  n    = 0
  do i = 1, 3
    if (roots(i,1) < 0.0 .and. roots(i,2) == 0.0) then
      mask(i) = 1
      n = n + 1
    end if
  end do

  if (n == 0) then
    if (is_root_pe()) then
      write(stdout,'(A)') "=== MOM_EBM estuary_box_model: no valid EBM solution found ==="
      write(stdout,'(A,I6,A,I6)') "  i          : ", ig, "  j : ", jg
      write(stdout,'(A,F12.4,A,F12.4)') "  lon        : ", G%geoLonT(ig,jg), &
                                         "  lat : ", G%geoLatT(ig,jg)
      write(stdout,'(A,ES15.8)') "  Q_r        : ", Q_r
      write(stdout,'(A,ES15.8)') "  S_l        : ", S_l
      write(stdout,'(A,ES15.8)') "  tide_amp   : ", CS%tide_amp
      write(stdout,'(A,ES15.8)') "  W_h        : ", CS%W_h(ig,jg)
      write(stdout,'(A,ES15.8)') "  H_est      : ", CS%H_est(ig,jg)
      write(stdout,'(A,ES15.8)') "  a1         : ", CS%a1(ig,jg)
      write(stdout,'(A,ES15.8)') "  a2         : ", CS%a2(ig,jg)
      write(stdout,'(A,ES15.8)') "  h0         : ", CS%h0(ig,jg)
      write(stdout,'(A,ES15.8)') "  g          : ", CS%g
      write(stdout,'(A,ES15.8)') "  rho_ref    : ", CS%rho_ref
      write(stdout,'(A,ES15.8)') "  beta_S     : ", CS%beta_S
      write(stdout,'(A,ES15.8)') "  Sc         : ", CS%Sc
      write(stdout,'(A,ES15.8)') "  R0         : ", R0
      write(stdout,'(A,ES15.8)') "  T0         : ", T0
      do i = 1, 3
        write(stdout,'(A,I2,A,ES15.8,A,ES15.8)') "  root(", i, ") : ", roots(i,1), " + i*", roots(i,2)
      enddo
      call MOM_error(WARNING, "MOM_EBM estuary_box_model: no valid EBM solution found.")
    endif
    ul0 = 0.0
  else if (n == 1) then
    ul0 = sum(roots(1:3,1) * real(mask))
  else
    if (is_root_pe()) then
      write(stdout,'(A)') "=== MOM_EBM estuary_box_model: multiple valid EBM solutions found ==="
      write(stdout,'(A,I6,A,I6)') "  i          : ", ig, "  j : ", jg
      write(stdout,'(A,F12.4,A,F12.4)') "  lon        : ", G%geoLonT(ig,jg), &
                                         "  lat : ", G%geoLatT(ig,jg)
      write(stdout,'(A,ES15.8)') "  Q_r        : ", Q_r
      write(stdout,'(A,ES15.8)') "  S_l        : ", S_l
      write(stdout,'(A,ES15.8)') "  tide_amp   : ", CS%tide_amp
      write(stdout,'(A,ES15.8)') "  W_h        : ", CS%W_h(ig,jg)
      write(stdout,'(A,ES15.8)') "  H_est      : ", CS%H_est(ig,jg)
      write(stdout,'(A,ES15.8)') "  a1         : ", CS%a1(ig,jg)
      write(stdout,'(A,ES15.8)') "  a2         : ", CS%a2(ig,jg)
      write(stdout,'(A,ES15.8)') "  h0         : ", CS%h0(ig,jg)
      write(stdout,'(A,ES15.8)') "  g          : ", CS%g
      write(stdout,'(A,ES15.8)') "  rho_ref    : ", CS%rho_ref
      write(stdout,'(A,ES15.8)') "  beta_S     : ", CS%beta_S
      write(stdout,'(A,ES15.8)') "  Sc         : ", CS%Sc
      write(stdout,'(A,ES15.8)') "  R0         : ", R0
      write(stdout,'(A,ES15.8)') "  T0         : ", T0
      do i = 1, 3
        write(stdout,'(A,I2,A,ES15.8,A,ES15.8)') "  root(", i, ") : ", roots(i,1), " + i*", roots(i,2)
      enddo
      call MOM_error(WARNING, "MOM_EBM estuary_box_model: multiple valid EBM solutions found.")
    endif
    ul0 = 0.0
  end if

  ! Upper layer salinity and volume fluxes at EBM mouth
  uu0 = R0 / (1.0 - CS%h0(ig,jg)) - CS%h0(ig,jg) / (1.0 - CS%h0(ig,jg)) * ul0
  S_u = (-S_l * ul0 * CS%h0(ig,jg) - S_l * CS%a2(ig,jg) * T0) / &
        (R0 - ul0 * CS%h0(ig,jg) - CS%a2(ig,jg) * T0)
  Q_l = ul0 * CS%h0(ig,jg) * CS%H_est(ig,jg) * CS%W_h(ig,jg) * c_wave
  Q_u = uu0 * (1.0 - CS%h0(ig,jg)) * CS%H_est(ig,jg) * CS%W_h(ig,jg) * c_wave

  ! Verify closure of EBM potential energy budget
  u_l   = ul0 * c_wave
  u_u   = uu0 * c_wave
  u_bar = Q_r / (CS%W_h(ig,jg) * CS%H_est(ig,jg))
  h_l   = CS%H_est(ig,jg) * CS%h0(ig,jg)
  rho_u = CS%rho_ref * (1.0 + CS%beta_S * S_u)

  AD = 0.5 * CS%g * rho_l * u_l * h_l**2.0 &
     + 0.5 * CS%g * (rho_u * u_u - rho_r * u_r) * (CS%H_est(ig,jg)**2.0 - h_l**2.0)

  HD = -0.5 * CS%a2(ig,jg) * CS%g * (rho_l - rho_u) * (CS%H_est(ig,jg)**2.0 - h_l**2.0) * u_t / PI

  VD = -0.5 * CS%g * (rho_u - rho_l) &
     * (rho_l + rho_u - 2.0 * rho_r) / (rho_l - rho_r) &
     * 0.024 * CS%a1(ig,jg) * CS%H_est(ig,jg)**2.0 &
     * (c_wave**4.0 / (u_bar * CS%Sc**2.0))**(1.0/3.0)

  LF = -0.25 * CS%g * Q_l / CS%W_h(ig,jg) &
     * ( (rho_u**2.0 + 2.0 * rho_l * rho_r - 2.0 * rho_u * rho_r) &
         * (CS%H_est(ig,jg) - h_l) &
       - rho_r**2.0 * CS%H_est(ig,jg) + rho_l**2.0 * h_l) / (rho_l - rho_r)

  ! GMM
  ! TODO: make this a 2D field and add option to save as a diagnostic?
  ERR_EBM = AD - HD - VD - LF

  ! Effective upper layer salinity for MOM6; Q_l is negative
  S_u = -Q_l * S_l / Q_u

end subroutine estuary_box_model

!> Solves the depressed cubic equation x^3 + ax^2 + bx + c = 0 analytically.
!! Returns all three roots; roots(:,1) are the real parts and roots(:,2) the
!! imaginary parts.
subroutine cubsolve(a, b, c, roots)

  real, intent(in)  :: a     !< Coefficient of x^2 in x^3 + ax^2 + bx + c = 0 [nondim]
  real, intent(in)  :: b     !< Coefficient of x in x^3 + ax^2 + bx + c = 0 [nondim]
  real, intent(in)  :: c     !< Constant term in x^3 + ax^2 + bx + c = 0 [nondim]
  real, dimension(3,2), intent(out) :: roots !< Roots [nondim]: column 1 = real, column 2 = imaginary

  ! local variables
  real :: Q     ! Intermediate cubic parameter [nondim]
  real :: R     ! Intermediate cubic parameter [nondim]
  real :: Rsqu  ! R squared [nondim]
  real :: Qcub  ! Q cubed [nondim]
  real :: SQ    ! Square root of Q [nondim]
  real :: theta ! Angle for the three-root case [nondim]
  real :: X     ! Intermediate root quantity [nondim]
  real :: Y     ! Intermediate root quantity [nondim]
  real :: XY    ! Sum X + Y [nondim]

  real, parameter :: PI = 4.0 * atan(1.0) !< Ratio of circumference to diameter [nondim]

  Q    = (a**2.0 - 3.0 * b) / 9.0
  R    = (2.0 * a**3.0 - 9.0 * a * b + 27.0 * c) / 54.0
  Rsqu = R**2.0
  Qcub = Q**3.0

  if (Rsqu < Qcub) then  ! Three distinct real roots
    theta      = acos(R / sqrt(Qcub))
    SQ         = sqrt(Q)
    roots(1,1) = -2.0 * SQ * cos(theta / 3.0) - a / 3.0
    roots(2,1) = -2.0 * SQ * cos((theta + 2.0 * PI) / 3.0) - a / 3.0
    roots(3,1) = -2.0 * SQ * cos((theta - 2.0 * PI) / 3.0) - a / 3.0
    roots(1,2) = 0.0
    roots(2,2) = 0.0
    roots(3,2) = 0.0
    return
  end if

  ! One real root and two conjugate complex roots
  X = -(abs(R) + sqrt(Rsqu - Qcub))**(1.0/3.0)
  if (R < 0.0) X = -X

  if (X == 0.0) then
    Y = 0.0
  else
    Y = Q / X
  end if

  XY         = X + Y
  roots(1,1) = XY - a / 3.0
  roots(1,2) = 0.0
  roots(2,1) = -0.5 * XY - a / 3.0
  roots(3,1) = -0.5 * XY - a / 3.0
  roots(2,2) =  sqrt(3.0) * (X - Y) / 2.0
  roots(3,2) = -sqrt(3.0) * (X - Y) / 2.0

end subroutine cubsolve

!> Unit tests for the EBM module. Returns true if any test fails.
logical function EBM_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, output additional information for debugging unit tests

  ! local variables
  real, dimension(3,2) :: roots   ! Roots from cubsolve [nondim]
  real :: tol                     ! Tolerance for floating-point comparisons [nondim]
  character(len=80) :: test_name  ! Brief description of the unit test
  type(EBM_cs) :: test_CS         ! Minimal CS for estuary_box_model tests
  type(ocean_grid_type) :: test_G ! Minimal G for estuary_box_model tests
  real :: Q_r_test                ! Test river discharge [m3 s-1]
  real :: S_l_test                ! Test lower layer salinity [ppt]
  real :: Q_u_test                ! Computed upper layer volume flux [m3 s-1]
  real :: Q_l_test                ! Computed lower layer volume flux [m3 s-1]
  real :: S_u_test                ! Computed upper layer salinity [ppt]
  real :: Q_u_expected            ! Expected upper layer volume flux [m3 s-1]
  real :: Q_l_expected            ! Expected lower layer volume flux [m3 s-1]
  real :: S_u_expected            ! Expected upper layer salinity [ppt]

  EBM_unit_tests = .false.
  tol = 1.0e-5
  write(stdout,*) '==== MOM_EBM ======================================='

  ! --- cubsolve tests ---

  ! Test 1: Three distinct real roots
  ! (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6 => a=-6, b=11, c=-6
  test_name = 'cubsolve: three real roots (x-1)(x-2)(x-3)'
  call cubsolve(-6.0, 11.0, -6.0, roots)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_cubsolve_real(roots, (/3.0, 1.0, 2.0/), tol, test_name, verbose)

  ! Test 2: One real root and two complex conjugate roots
  ! (x-1)(x^2+1) = x^3 - x^2 + x - 1 => a=-1, b=1, c=-1
  test_name = 'cubsolve: one real + two complex roots'
  call cubsolve(-1.0, 1.0, -1.0, roots)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_cubsolve_single_real(roots, 1.0, tol, test_name, verbose)

  ! Test 3: Triple root at x=2
  ! (x-2)^3 = x^3 - 6x^2 + 12x - 8 => a=-6, b=12, c=-8
  test_name = 'cubsolve: triple root at x=2'
  call cubsolve(-6.0, 12.0, -8.0, roots)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_cubsolve_single_real(roots, 2.0, tol, test_name, verbose)

  ! Test 4: Roots include zero
  ! x(x-1)(x+1) = x^3 - x => a=0, b=-1, c=0
  test_name = 'cubsolve: roots at -1, 0, 1'
  call cubsolve(0.0, -1.0, 0.0, roots)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_cubsolve_real(roots, (/1.0, -1.0, 0.0/), tol, test_name, verbose)

  ! Test 5: Negative real roots (relevant to EBM: ul0 must be real and negative)
  ! (x+1)(x+2)(x-3) = x^3 + 0x^2 - 7x - 6 => a=0, b=-7, c=-6 (not depressed but general form)
  ! Actually: x^3 - 7x - 6, a=0, b=-7, c=-6
  test_name = 'cubsolve: roots with two negative values'
  call cubsolve(0.0, -7.0, -6.0, roots)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_cubsolve_real(roots, (/3.0, -1.0, -2.0/), tol, test_name, verbose)

  if (.not. EBM_unit_tests) write(stdout,*) 'Passed cubsolve tests'

  ! --- estuary_box_model tests ---
  ! Infrastructure: set up a minimal CS and G for testing.
  call setup_EBM_test_CS(test_CS, test_G)

  ! Test: Amazon-like river with Sun et al. Table 2 parameters
  ! test_name = 'estuary_box_model: Amazon-like'
  ! call estuary_box_model(test_CS, test_G, 1, 1, Q_r_test, S_l_test, Q_u_test, Q_l_test, S_u_test)
  ! EBM_unit_tests = EBM_unit_tests .or. &
  !   test_ebm_scalar(Q_u_test, Q_u_expected, tol, 'Q_u', test_name, verbose)
  ! EBM_unit_tests = EBM_unit_tests .or. &
  !   test_ebm_scalar(Q_l_test, Q_l_expected, tol, 'Q_l', test_name, verbose)
  ! EBM_unit_tests = EBM_unit_tests .or. &
  !   test_ebm_scalar(S_u_test, S_u_expected, tol, 'S_u', test_name, verbose)

  ! Test: Columbia river, Sun et al. Fig 3.
  test_name = 'estuary_box_model: Columbia river'
  Q_r_test = 7.5E3 ! m3/s
  S_l_test = 32.   ! psu
  test_CS%H_est(1,1) = 10.93  ! m
  test_CS%W_h(1,1)   = 3.67E3 ! m
  test_CS%a1(1,1)     = 1.2
  test_CS%a2(1,1)     = 0.93
  Q_u_expected = 6000.
  Q_l_expected = Q_u_expected - Q_r_test
  S_u_expected = 15.
  tol = 1000.
  call estuary_box_model(test_CS, test_G, 1, 1, Q_r_test, S_l_test, Q_u_test, Q_l_test, S_u_test)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_ebm_scalar(Q_u_test, Q_u_expected, tol, 'Q_u', test_name, verbose)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_ebm_scalar(Q_l_test, Q_l_expected, tol, 'Q_l', test_name, verbose)
  EBM_unit_tests = EBM_unit_tests .or. &
    test_ebm_scalar(S_u_test, S_u_expected, tol, 'S_u', test_name, verbose)
  if (.not. EBM_unit_tests) write(stdout,*) 'Passed estuary_box_model tests'

  call teardown_EBM_test_CS(test_CS, test_G)

end function EBM_unit_tests

!> Returns true if any real root from cubsolve does not match the expected values.
!! All three roots are expected to be real (imaginary parts = 0).
logical function test_cubsolve_real(roots, expected, tol, test_name, verbose)
  real, dimension(3,2), intent(in) :: roots    !< Computed roots [nondim]
  real, dimension(3),   intent(in) :: expected !< Expected real roots [nondim]
  real,                 intent(in) :: tol      !< Tolerance for comparison [nondim]
  character(len=*),     intent(in) :: test_name !< Brief description of the unit test
  logical,              intent(in) :: verbose  !< If true, always print output

  ! local variables
  real, dimension(3) :: calc_sorted, exp_sorted ! Sorted arrays for comparison [nondim]
  integer :: i

  test_cubsolve_real = .false.

  ! Check that all imaginary parts are zero
  do i = 1, 3
    if (abs(roots(i,2)) > tol) then
      test_cubsolve_real = .true.
      write(stdout,*) "MOM_EBM UNIT TEST FAILED: ", trim(test_name)
      write(stdout,'(A,I2,A,ES15.8)') "  Root ", i, " has nonzero imaginary part: ", roots(i,2)
      return
    endif
  enddo

  ! Sort both arrays for comparison (simple bubble sort for 3 elements)
  calc_sorted = roots(:,1)
  exp_sorted = expected
  call sort3(calc_sorted)
  call sort3(exp_sorted)

  do i = 1, 3
    if (abs(calc_sorted(i) - exp_sorted(i)) > tol) then
      test_cubsolve_real = .true.
      write(stdout,*) "MOM_EBM UNIT TEST FAILED: ", trim(test_name)
      write(stdout,'(A,I2,A,ES15.8,A,ES15.8)') "  Sorted root ", i, &
        ": calc=", calc_sorted(i), " expected=", exp_sorted(i)
    endif
  enddo

  if (verbose .and. .not. test_cubsolve_real) then
    write(stdout,*) "  PASSED: ", trim(test_name)
    do i = 1, 3
      write(stdout,'(A,I2,A,ES15.8,A,ES15.8)') "    root ", i, &
        ": calc=", calc_sorted(i), " expected=", exp_sorted(i)
    enddo
  endif

end function test_cubsolve_real

!> Returns true if the single real root from cubsolve does not match the expected value.
!! Only the root(s) with zero imaginary part are checked.
logical function test_cubsolve_single_real(roots, expected, tol, test_name, verbose)
  real, dimension(3,2), intent(in) :: roots     !< Computed roots [nondim]
  real,                 intent(in) :: expected  !< Expected real root value [nondim]
  real,                 intent(in) :: tol       !< Tolerance for comparison [nondim]
  character(len=*),     intent(in) :: test_name !< Brief description of the unit test
  logical,              intent(in) :: verbose   !< If true, always print output

  ! local variables
  integer :: i
  logical :: found

  test_cubsolve_single_real = .false.
  found = .false.

  do i = 1, 3
    if (abs(roots(i,2)) < tol) then
      found = .true.
      if (abs(roots(i,1) - expected) > tol) then
        test_cubsolve_single_real = .true.
        write(stdout,*) "MOM_EBM UNIT TEST FAILED: ", trim(test_name)
        write(stdout,'(A,ES15.8,A,ES15.8)') "  Real root: calc=", roots(i,1), " expected=", expected
      elseif (verbose) then
        write(stdout,*) "  PASSED: ", trim(test_name)
        write(stdout,'(A,ES15.8,A,ES15.8)') "    Real root: calc=", roots(i,1), " expected=", expected
      endif
    endif
  enddo

  if (.not. found) then
    test_cubsolve_single_real = .true.
    write(stdout,*) "MOM_EBM UNIT TEST FAILED: ", trim(test_name)
    write(stdout,*) "  No real root found"
  endif

end function test_cubsolve_single_real

!> Set up a minimal EBM control structure and ocean grid for unit testing.
!! Allocates 1x1 2D arrays in CS and G with Amazon-like parameters from
!! Sun et al. (2017) Table 2.
subroutine setup_EBM_test_CS(CS, G)
  type(EBM_cs),        intent(inout) :: CS !< EBM control structure to initialize
  type(ocean_grid_type), intent(inout) :: G  !< Ocean grid structure to initialize

  ! Scalar EBM parameters
  CS%tide_amp = 1.0      ! [m]
  CS%g        = 9.8      ! [m s-2]
  CS%rho_ref  = 1000.0   ! [kg m-3]
  CS%beta_S   = 7.7e-4   ! [ppt-1]
  CS%Sc       = 2.2      ! [nondim]

  ! Per-column arrays (1x1)
  allocate(CS%H_est(1,1)) ; CS%H_est(1,1) = 21.8    ! Amazon depth [m]
  allocate(CS%W_h(1,1))   ; CS%W_h(1,1)   = 50000.0 ! Amazon head width [m]
  allocate(CS%a1(1,1))    ; CS%a1(1,1)     = 1.0     ! Amazon a1 [nondim]
  allocate(CS%a2(1,1))    ; CS%a2(1,1)     = 0.0     ! Amazon a2 [nondim]
  allocate(CS%h0(1,1))    ; CS%h0(1,1)     = 0.5     ! Amazon h0 [nondim]

  ! Minimal ocean grid (only geoLonT, geoLatT are accessed by estuary_box_model
  ! in warning messages)
  allocate(G%geoLonT(1,1)) ; G%geoLonT(1,1) = 310.67  ! Amazon lon
  allocate(G%geoLatT(1,1)) ; G%geoLatT(1,1) = 0.0     ! Amazon lat

end subroutine setup_EBM_test_CS

!> Deallocate arrays created by setup_EBM_test_CS.
subroutine teardown_EBM_test_CS(CS, G)
  type(EBM_cs),        intent(inout) :: CS !< EBM control structure to clean up
  type(ocean_grid_type), intent(inout) :: G  !< Ocean grid structure to clean up

  deallocate(CS%H_est, CS%W_h, CS%a1, CS%a2, CS%h0)
  deallocate(G%geoLonT, G%geoLatT)

end subroutine teardown_EBM_test_CS

!> Returns true if a computed scalar does not match the expected value within tolerance.
logical function test_ebm_scalar(calc, expected, tol, varname, test_name, verbose)
  real,             intent(in) :: calc      !< Computed value [arbitrary]
  real,             intent(in) :: expected  !< Expected value [arbitrary]
  real,             intent(in) :: tol       !< Tolerance for comparison [nondim]
  character(len=*), intent(in) :: varname   !< Name of the variable being checked
  character(len=*), intent(in) :: test_name !< Brief description of the unit test
  logical,          intent(in) :: verbose   !< If true, always print output

  real :: rel_err ! Relative error [nondim]

  test_ebm_scalar = .false.
  if (expected /= 0.0) then
    rel_err = abs(calc - expected) / abs(expected)
  else
    rel_err = abs(calc - expected)
  endif

  if (rel_err > tol) then
    test_ebm_scalar = .true.
    write(stdout,*) "MOM_EBM UNIT TEST FAILED: ", trim(test_name)
    write(stdout,'(A,A,A,ES15.8,A,ES15.8,A,ES10.3)') "  ", trim(varname), &
      ": calc=", calc, " expected=", expected, " rel_err=", rel_err
  elseif (verbose) then
    write(stdout,'(A,A,A,ES15.8,A,ES15.8)') "  PASSED ", trim(varname), &
      ": calc=", calc, " expected=", expected
  endif

end function test_ebm_scalar

!> Sort a 3-element array in ascending order (simple swap sort).
subroutine sort3(a)
  real, dimension(3), intent(inout) :: a !< Array to sort [arbitrary]
  real :: tmp

  if (a(1) > a(2)) then ; tmp = a(1) ; a(1) = a(2) ; a(2) = tmp ; endif
  if (a(2) > a(3)) then ; tmp = a(2) ; a(2) = a(3) ; a(3) = tmp ; endif
  if (a(1) > a(2)) then ; tmp = a(1) ; a(1) = a(2) ; a(2) = tmp ; endif

end subroutine sort3

end module MOM_EBM
