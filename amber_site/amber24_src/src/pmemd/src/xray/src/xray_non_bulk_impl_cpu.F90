module xray_non_bulk_impl_cpu_module
  
  use xray_pure_utils, only : real_kind
  use xray_contracts_module
  use xray_non_bulk_data_module
  
  implicit none
  
  private
  
  public :: init
  public :: finalize
  public :: calc_f_non_bulk
  public :: get_f_non_bulk

  real(real_kind), parameter :: m_twopi = 2 * 3.1415926535897932384626433832795d0 ! FIXME: use constants module

contains
  
  subroutine init(hkl_, mSS4_, b_factor_, scatter_type_index_, occupancy_)
    implicit none
    integer, intent(in), target :: hkl_(:, :)
    real(real_kind), intent(in), target :: mSS4_(:)
    real(real_kind), intent(in) :: b_factor_(:)
    integer, intent(in) :: scatter_type_index_(:)
    real(real_kind), intent(in) :: occupancy_(:)
    
    call check_precondition(size(hkl_, 1) == 3)
    call check_precondition(size(hkl_, 2) == size(mSS4_))
    call check_precondition(size(b_factor_) == size(scatter_type_index_))
    call check_precondition(size(b_factor_) == size(occupancy_))

    
    hkl => hkl_
    mSS4 => mSS4_
    b_factor = b_factor_
    scatter_type_index = scatter_type_index_
    
    occupancy = occupancy_

    allocate(F_non_bulk(size(mSS4_)))
    allocate(f(size(b_factor_)))
    allocate(angle(size(b_factor_)))
  
  end subroutine init
  
  subroutine finalize()
    implicit none
    hkl => null()
    mSS4 => null()
    xyz => null()
    deallocate(b_factor)
    deallocate(scatter_type_index)
    if (allocated(occupancy)) then
      deallocate(occupancy)
    endif
    deallocate(F_non_bulk)
    deallocate(f)
    deallocate(angle)
  end subroutine finalize
  
  !-------------------------------------------------------------------
  ! Caution: Some literature uses S to represent S^2
  !
  ! (for d*, see p. 93 of Glusker, Lewis, Rossi)
  ! S == d* = sqrt(sum(HKL * orth_to_frac)^2) = sqrt(-4*mSS4)
  ! mSS4 = -S*S/4
  ! mSS4 is more relevant to the formulas used than S.
  
  function get_f_non_bulk() result(result)
    complex(real_kind) :: result(size(F_non_bulk))
    result(:) = F_non_bulk(:)
  end function get_f_non_bulk
  
  subroutine calc_f_non_bulk(frac)
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    ! locals
    integer :: ihkl
    
    call check_precondition(associated(hkl))
    call check_precondition(associated(mSS4))
    call check_precondition(allocated(atomic_scatter_factor))
    call check_precondition(allocated(b_factor))
    call check_precondition(allocated(scatter_type_index))
    
    call check_precondition(size(frac, 1) == 3)
    call check_precondition(size(frac, 2) == size(b_factor))
    call check_precondition(size(frac, 2) == size(scatter_type_index))
    call check_precondition(size(hkl, 2) == size(atomic_scatter_factor, 1))
    
    !$omp parallel do private(ihkl,i,f,angle)
    do ihkl = 1, size(hkl, 2)
      
      ! Fhkl = SUM( fj * exp(2 * M_PI * i * (h * xj + k * yj + l * zj)) ),
      !      j = 1,num_selected_atoms
      ! where:
      !    The sum is versus j, over all selected atoms
      !    fj is the atomic scatter for atom j:   atomic_scatter_factor(j)
      !    h,k,l are from the hkl list for index ihkl:   hkl(1:3,ihkl)
      !    x,y,z are coordinates for atom j:   xyz(1:3,j)
      !        xyz(:) may be a reduced list.
      !
      ! Rather than using a complex exponential where the real part is
      ! always zero, this is optimized to calculate sin and cosine parts,
      ! then convert to a complex number
      ! after the A and B components are summed over all selected atoms.
      ! This can be written as:
      !
      ! Ahkl = SUM( fj * cos(2 * M_PI * (h * xj + k * yj + l * zj)) ),
      ! Bhkl = SUM( fj * sin(2 * M_PI * (h * xj + k * yj + l * zj)) ),
      !    j = 1,num_selected_atoms
      
      f(:) = exp(mSS4(ihkl) * b_factor(:)) &
         & * atomic_scatter_factor(ihkl, scatter_type_index(:)) &
         & * occupancy(:)

      angle(:) = matmul(M_TWOPI * hkl(1:3, ihkl), frac(1:3, :))
      
      F_non_bulk(ihkl) = cmplx(sum(f(:) * cos(angle(:))), &
          sum(f(:) * sin(angle(:))), real_kind)
    
    end do
    !$omp end parallel do
  
  end subroutine calc_f_non_bulk


end module xray_non_bulk_impl_cpu_module