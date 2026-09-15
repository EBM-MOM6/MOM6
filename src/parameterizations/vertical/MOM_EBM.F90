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
use MOM_error_handler,         only : MOM_error, WARNING, FATAL, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_io,                    only : file_exists, field_size, open_file_to_read
use MOM_io,                    only : close_file_to_read, read_variable, stdout, slasher
use MOM_remapping,             only : remapping_CS, initialize_remapping, reintegrate_column
use MOM_remapping,             only : extract_member_remapping_CS, remapping_core_h
use MOM_remapping,             only : remappingSchemesDoc, remappingDefaultScheme
use MOM_verticalGrid,          only : verticalGrid_type

#include <MOM_memory.h>

implicit none ; private

public EBM_init, calculate_EBM, EBM_is_used

!> Control structure including parameters for the estuary box model.
type, public :: EBM_cs ; private

  real :: tide_amp !< Averaged tidal amplitude at estuary mouth [m]
  real :: W_h      !< Estuary head width [m]
  real :: H        !< Estuary averaged depth [m]
  real :: a2       !< A constant of tidal diffusion [nondim]
  real :: a1       !< A constant of estuarine mixing length [nondim]
  real :: h0       !< A constant of ratio of geometry: h_l/H [nondim]
  real :: g        !< Gravitational acceleration [m s-2]
  real :: rho_ref    !< Reference density for linear equation of state [kg m-3]
  real :: beta_S   !< Saline contraction coefficient [ppt-1]
  real :: Sc       !< Schmidt number [nondim]
  real, allocatable, dimension(:,:) :: H_est !< Estuary averaged depth [m]
  real, allocatable, dimension(:,:) :: H_U  !< Thickness of EBM discharge into MOM6 (upper) [m]
  real, allocatable, dimension(:,:) :: H_L  !< Thickness of MOM6 discharging to EBM (lower) [m]
  real, allocatable, dimension(:,:,:) :: wf_U !< Weighting function for the upper layer [nondim]
  real, allocatable, dimension(:,:,:) :: wf_L !< Weighting function for the lower layer [nondim]
  integer :: deg              !< Degree of polynomial reconstruction [nondim]
  real    :: H_subroundoff   !< A thickness that is so small that it can be added to a thickness of
                             !! Angstrom or larger without changing it at the bit level [H ~> m or kg m-2].
  type(remapping_CS)                :: remap_CS !< Control structure to hold remapping configuration.

end type EBM_cs

character(len=40) :: mdl = "MOM_EBM"  !< This module's name.

contains

!> Initializes the estuary box model parameterization.
!! Returns .true. if the parameterization is enabled.
logical function EBM_init(param_file, G, GV, diag, CS)

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
  real, dimension(:), allocatable :: new_depth ! The new values of estuary depth [m]
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

  call get_param(param_file, mdl, "EBM_TIDAL_AMPLITUDE", CS%tide_amp, &
                 "Averaged tidal amplitude at the estuary mouth.", &
                 units="m", default=1.0)

  call get_param(param_file, mdl, "EBM_HEAD_WIDTH", CS%W_h, &
                 "Estuary head width.", &
                 units="m", default=2000.0)

  call get_param(param_file, mdl, "EBM_DEPTH", CS%H, &
                 "Estuary averaged depth.", &
                 units="m", default=10.0)

  call get_param(param_file, mdl, "EBM_A1", CS%a1, &
                 "A constant of estuarine mixing length.", &
                 units="nondim", default=0.876)

  call get_param(param_file, mdl, "EBM_A2", CS%a2, &
                 "A constant of tidal diffusion.", &
                 units="nondim", default=0.0)

  call get_param(param_file, mdl, "EBM_H0", CS%h0, &
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
  allocate(CS%wf_U(SZI_(G),SZJ_(G),2))
  allocate(CS%wf_L(SZI_(G),SZJ_(G),2))
  CS%H_est(:,:) = CS%H
  CS%H_U(:,:) = 0.0
  CS%H_L(:,:) = 0.0
  CS%wf_U(:,:,:) = 0.0
  CS%wf_L(:,:,:) = 0.0

  ! Apply optional point-by-point overrides of H_est from a file
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "EBM_DEPTH_EDITS_FILE", ebm_edits_file, &
                 "The file from which to read a list of i,j,depth overrides for the "//&
                 "estuary averaged depth (H_est). The file must be in NetCDF format "//&
                 "with scalar variables 'ni' and 'nj' matching the global grid "//&
                 "dimensions, and 1-D arrays 'iEdit', 'jEdit', and 'zEdit' (all of "//&
                 "size nEdits) giving the global i-index, global j-index (both using "//&
                 "Python 0-based indexing), and new estuary depth [m] for each edit.", &
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

    ! Get nEdits
    call field_size(ebm_edits_file, 'zEdit', sizes, ndims=ndims, ncid_in=ncid)
    if (ndims /= 1) call MOM_error(FATAL, "The variable zEdit has an "//&
              "unexpected number of dimensions in "//trim(ebm_edits_file))
    n_edits = sizes(1)
    allocate(ig(n_edits), jg(n_edits), new_depth(n_edits))

    ! Read iEdit, jEdit and zEdit
    call read_variable(ebm_edits_file, 'iEdit', ig, ncid_in=ncid)
    call read_variable(ebm_edits_file, 'jEdit', jg, ncid_in=ncid)
    call read_variable(ebm_edits_file, 'zEdit', new_depth, ncid_in=ncid)
    call close_file_to_read(ncid, ebm_edits_file)

    do n = 1, n_edits
      i = ig(n) - G%idg_offset + 1 ! +1 for python indexing
      j = jg(n) - G%jdg_offset + 1
      if (i>=G%isc .and. i<=G%iec .and. j>=G%jsc .and. j<=G%jec) then
        write(stdout,'(a,3i5,f8.2,a,f8.2,2i4)') &
          'EBM depth edit: ', n, ig(n), jg(n), CS%H_est(i,j), '->', abs(new_depth(n)), i, j
        CS%H_est(i,j) = abs(new_depth(n))
      endif
    enddo

    deallocate(ig, jg, new_depth)
  endif

  ! TODO: revisit H_U and H_L later with limiters?
  do j=G%jsc,G%jec
    do i=G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        CS%H_U(i,j) = CS%H_est(i,j) * 0.5
        CS%H_L(i,j) = CS%H_est(i,j) * 0.5
        CS%wf_U(i,j,1) = 1.0/CS%H_U(i,j)
        CS%wf_L(i,j,2) = 1.0/CS%H_L(i,j)
      endif
    enddo
  enddo

  ! Post the static estuary depth field
  id = register_static_field('ocean_model', 'ebm_depth', diag%axesT1, &
        'Estuary averaged depth used by the EBM', 'm')
  if (id > 0) call post_data(id, CS%H_est, diag, .true.)

end function EBM_init

!> Calculates estuary box model exchange and distributes river runoff over the
!! first 4 layers, updating layer thickness, temperature, salinity, and net mass flux.
!! The estuary exchange fluxes (Q_u, Q_l, S_u) are computed by estuary_box_model
!! using the EBM parameters in CS together with the provided river discharge and
!! lower-layer salinity.
subroutine calculate_EBM(CS, i, j, lrunoff, EnthalpyConst, netMassIn, T2d_col, S_col, h2d_col)

  type(EBM_cs), intent(in)    :: CS            !< EBM control structure
  integer,      intent(in)    :: i             !< i-index of the current column [nondim]
  integer,      intent(in)    :: j             !< j-index of the current column [nondim]
  real,         intent(in)    :: lrunoff       !< River runoff for this column [H ~> m or kg m-2]
  real,         intent(in)    :: EnthalpyConst !< Enthalpy constant [nondim]
  real,         intent(inout) :: netMassIn     !< Net mass entering this column [H ~> m or kg m-2]
  real,         intent(inout) :: T2d_col(:)    !< Temperature in first 4 layers [C ~> degC]
  real,         intent(inout) :: S_col(:)      !< Salinity in first 4 layers [S ~> ppt]
  real,         intent(inout) :: h2d_col(:)    !< Layer thickness in first 4 layers [H ~> m or kg m-2]

  ! local variables
  real :: Q_u        ! EBM upper layer volume flux [m3 s-1]
  real :: Q_l        ! EBM lower layer volume flux [m3 s-1]
  real :: S_l        ! EBM lower layer salinity    [S ~> ppt]
  real :: S_u        ! EBM upper layer salinity [S ~> ppt]
  real :: dThickness ! Change in layer thickness [H ~> m or kg m-2]
  real :: dTemp      ! Integrated change in layer temperature [C H ~> degC m or degC kg m-2]
  real :: dSalt      ! Integrated change in layer salinity [S H ~> ppt m or ppt kg m-2]
  real :: Temp_in    ! Temperature of the incoming mass flux [C ~> degC]
  real :: Salin_in   ! Salinity of the incoming mass flux [S ~> ppt]
  real :: hOld       ! Original layer thickness before update [H ~> m or kg m-2]
  real :: Ithickness ! Inverse of the updated layer thickness [H-1 ~> m-1 or m2 kg-1]
  integer :: k       ! Layer index [nondim]
  real, dimension(2) :: dz_ebm ! EBM layer thicknesses [m]

  ! EBM grid
  dz_ebm(1) = CS%H_U(i,j)
  dz_ebm(2) = CS%H_L(i,j)

  ! GMM, TODO: need to specify Hu, instead of hard code thickness...

  ! Compute salinity of lower layer [S ~> ppt]
  ! s_l = get_lower_layer_salinity(S_col, h2d_col)

  ! Distribute river runoff over upper layer
  do k = 1, size(T2d_col)
    dThickness = lrunoff * 0.25  ! Each of the 4 layers receives 1/4 of the runoff
    dTemp = 0.
    dSalt = 0.

    netMassIn = netMassIn - dThickness
    Temp_in  = T2d_col(k)
    dTemp = dTemp + dThickness * Temp_in * EnthalpyConst
    ! GMM, this is not used. Delete?
    Salin_in = 0.0

    hOld = h2d_col(k)
    h2d_col(k) = h2d_col(k) + dThickness
    if (h2d_col(k) > 0.0) then
      Ithickness = 1.0 / h2d_col(k)
      if (dThickness /= 0. .or. dTemp /= 0.) T2d_col(k) = (hOld * T2d_col(k) + dTemp) * Ithickness
      if (dThickness /= 0. .or. dSalt /= 0.) S_col(k)   = (hOld * S_col(k)   + dSalt) * Ithickness
    end if
  end do

  ! GMM, todo
  !call estuary_box_model(CS, lrunoff, S_l, Q_u, Q_l, S_u)

end subroutine calculate_EBM

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
subroutine estuary_box_model(CS, Q_r, S_l, Q_u, Q_l, S_u)

  type(EBM_cs), intent(in)  :: CS   !< EBM control structure
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
  u_t    = -CS%tide_amp * sqrt(CS%g / CS%H)        ! Tidal velocity (toward river)
  u_r    = Q_r / (CS%W_h * CS%H * (1.0 - CS%h0))  ! Riverine velocity at head of upper layer
  c_wave = sqrt(CS%beta_S * S_l * CS%g * CS%H)     ! Densimetric wave phase speed

  ! Dimensionless parameters
  ur0 = u_r / c_wave
  ut0 = u_t / c_wave
  R0  = ur0 * (1.0 - CS%h0)
  T0  = ut0 * (1.0 - CS%h0) / PI

  ! Coefficients of the cubic equation for dimensionless lower layer inflow (ul0)
  a = -CS%h0**3.0

  b = 2.0 * CS%h0**2.0 * ((2.0 - CS%h0) * R0 - CS%a2 * T0)

  c = 0.096 * CS%a1 * CS%h0 * (CS%Sc**2.0 * R0)**(-1.0/3.0) * R0 &
    - CS%h0 * ((2.0 - CS%h0) * R0 * (R0 - 2.0 * CS%a2 * T0) + CS%a2**2.0 * T0**2.0)

  d = -0.048 * CS%a1 * (CS%Sc**2.0 * R0)**(-1.0/3.0) &
    * R0 * (R0 - 2.0 * CS%a2 * T0)

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
    call MOM_error(WARNING, "MOM_EBM estuary_box_model: no valid EBM solution found.")
    ul0 = 0.0
  else if (n == 1) then
    ul0 = sum(roots(1:3,1) * real(mask))
  else
    call MOM_error(WARNING, "MOM_EBM estuary_box_model: multiple valid EBM solutions found.")
    ul0 = 0.0
  end if

  ! Upper layer salinity and volume fluxes at EBM mouth
  uu0 = R0 / (1.0 - CS%h0) - CS%h0 / (1.0 - CS%h0) * ul0
  S_u = (-S_l * ul0 * CS%h0 - S_l * CS%a2 * T0) / (R0 - ul0 * CS%h0 - CS%a2 * T0)
  Q_l = ul0 * CS%h0 * CS%H * CS%W_h * c_wave
  Q_u = uu0 * (1.0 - CS%h0) * CS%H * CS%W_h * c_wave

  ! Verify closure of EBM potential energy budget
  u_l   = ul0 * c_wave
  u_u   = uu0 * c_wave
  u_bar = Q_r / (CS%W_h * CS%H)
  h_l   = CS%H * CS%h0
  rho_u = CS%rho_ref * (1.0 + CS%beta_S * S_u)

  AD = 0.5 * CS%g * rho_l * u_l * h_l**2.0 &
     + 0.5 * CS%g * (rho_u * u_u - rho_r * u_r) * (CS%H**2.0 - h_l**2.0)

  HD = -0.5 * CS%a2 * CS%g * (rho_l - rho_u) * (CS%H**2.0 - h_l**2.0) * u_t / PI

  VD = -0.5 * CS%g * (rho_u - rho_l) &
     * (rho_l + rho_u - 2.0 * rho_r) / (rho_l - rho_r) &
     * 0.024 * CS%a1 * CS%H**2.0 &
     * (c_wave**4.0 / (u_bar * CS%Sc**2.0))**(1.0/3.0)

  LF = -0.25 * CS%g * Q_l / CS%W_h &
     * ( (rho_u**2.0 + 2.0 * rho_l * rho_r - 2.0 * rho_u * rho_r) &
         * (CS%H - h_l) &
       - rho_r**2.0 * CS%H + rho_l**2.0 * h_l) / (rho_l - rho_r)

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

end module MOM_EBM
