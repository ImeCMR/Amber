module xray_bulk_model_module
  
  use xray_pure_utils, only : real_kind
  use xray_contracts_module
  use xray_unit_cell_module
  
  implicit none
  private
  
  public :: add_bulk_contribution_and_rescale
  public :: get_f_scale
  public :: finalize
  public :: init
  
  integer :: model_id
  integer, parameter :: none_id = 0
  integer, parameter :: afonine_2013_id = 1
  integer, parameter :: simple_id = 2
  
  character(len = *), parameter :: none_name = "none"
  character(len = *), parameter :: afonine_2013_name = "afonine-2013" ! Named after https://doi.org/10.1107/S0907444913000462
  character(len = *), parameter :: simple_name = "simple"

contains
  
  subroutine init(model_name, mask_update_period, scale_update_period, &
      & resolution_high, hkl, unit_cell, atm_atomicnumber, k_sol, b_sol, &
      & solvent_mask_adjustment, solvent_mask_probe_radius )
    use xray_bulk_model_afonine_2013_module, only : init_afonine => init
    use xray_bulk_model_none_module, only : init_none => init
    use xray_bulk_model_simple_module, only : init_simple => init
    implicit none
    integer, intent(in) :: mask_update_period
    integer, intent(in) :: scale_update_period
    character(len = *), intent(in) :: model_name
    double precision, intent(in) :: resolution_high
    class(unit_cell_t), intent(in) :: unit_cell
    integer, intent(in) :: hkl(:, :)
    integer, intent(in) :: atm_atomicnumber(:)
    real(real_kind), intent(in) :: k_sol
    real(real_kind), intent(in) :: b_sol
    real(real_kind), intent(in) :: solvent_mask_adjustment
    real(real_kind), intent(in) :: solvent_mask_probe_radius
    
    model_id = model_name_to_id(model_name)
    
    select case (model_id)
    case (none_id)
      call init_none(scale_update_period)
    case (afonine_2013_id)
      call init_afonine(mask_update_period, scale_update_period, resolution_high, hkl, unit_cell, atm_atomicnumber, &
          & solvent_mask_adjustment, solvent_mask_probe_radius)
    case (simple_id)
      call init_simple(k_sol, b_sol, mask_update_period, scale_update_period, resolution_high, hkl, unit_cell, atm_atomicnumber, &
          & solvent_mask_adjustment, solvent_mask_probe_radius)
    case default
      call check_requirement(.FALSE., "Bad model id")
    end select
  end subroutine init
  
  subroutine finalize()
    use xray_bulk_model_afonine_2013_module, only : finalize_afonine => finalize
    use xray_bulk_model_none_module, only : finalize_none => finalize
    use xray_bulk_model_simple_module, only : finalize_simple => finalize
    implicit none
    
    select case (model_id)
    case (none_id)
      call finalize_none()
    case (afonine_2013_id)
      call finalize_afonine()
    case (simple_id)
      call finalize_simple()
    case default
      call check_requirement(.FALSE., "Bad model id")
    end select
  end subroutine finalize
  
  subroutine add_bulk_contribution_and_rescale(frac, current_step, absFobs, Fcalc, mSS4, hkl)
    use xray_bulk_model_afonine_2013_module, only : afonine_f => add_bulk_contribution_and_rescale
    use xray_bulk_model_none_module, only : none_f => add_bulk_contribution_and_rescale
    use xray_bulk_model_simple_module, only : simple_f => add_bulk_contribution_and_rescale
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(inout) :: Fcalc(size(absFobs)) !< input: Fcalc=Fprot, output Fcalc=Fcalc
    real(real_kind), intent(in) :: mSS4(:)
    integer, intent(in) :: hkl(:, :)
    logical, save :: first_call = .TRUE.

    call check_requirement(.not. first_call .or. current_step == 0, &
        & "First call of `xray_bulk_model_module::add_bulk_contribution_and_rescale(...)` &
        & must be made with current_step=0")
    first_call = .FALSE.

    select case (model_id)
    case (none_id)
      call none_f(current_step, absFobs, Fcalc)
    case (afonine_2013_id)
      call afonine_f(frac, current_step, absFobs, Fcalc, mSS4, hkl)
    case (simple_id)
      call simple_f(frac, current_step, absFobs, Fcalc, mSS4)
    case default
      call check_requirement(.FALSE., "Bad model id")
    end select
  end subroutine add_bulk_contribution_and_rescale
  
  function get_f_scale(n_hkl) result(result)
    use xray_bulk_model_afonine_2013_module, only : afonine_f => get_f_scale
    use xray_bulk_model_none_module, only : none_f => get_f_scale
    use xray_bulk_model_simple_module, only : simple_f => get_f_scale
    implicit none
    integer, intent(in) :: n_hkl
    real(real_kind) :: result(n_hkl)
    
    select case (model_id)
    case (none_id)
      result = none_f(n_hkl)
    case (afonine_2013_id)
      result = afonine_f(n_hkl)
    case (simple_id)
      result = simple_f(n_hkl)
    case default
      call check_requirement(.FALSE., "Bad model id")
    end select
  end function get_f_scale
  
  function model_name_to_id(name) result(result)
    implicit none
    character(len = *), intent(in) :: name
    integer :: result
    
    select case (trim(name))
    case (none_name)
      result = none_id
    case (afonine_2013_name)
      result = afonine_2013_id
    case (simple_name)
      result = simple_id
    case default
      result = -1 ! to suppress warning
      call check_requirement(.FALSE., "Unknown bulk solvent model name: '" // trim(name) // "'")
    end select
  end function model_name_to_id

end module xray_bulk_model_module