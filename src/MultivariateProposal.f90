module MultivariateProposal
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use precision, only: fp_kind
  use RNG
  use PriorProposal
  implicit none
contains

  !------------------------------------------------------------------
  ! Subroutine: CholeskyDecomposition
  ! Purpose: Compute the Cholesky decomposition of an n x n matrix A.
  !          On output, L is lower triangular such that A = L * L^T.
  !------------------------------------------------------------------
  subroutine CholeskyDecomposition(n, A, L, info)
    implicit none
    integer, intent(in) :: n
    real(fp_kind), intent(in) :: A(n,n)
    real(fp_kind), intent(out) :: L(n,n)
    integer, intent(out) :: info
    integer :: i, j, k
    real(fp_kind) :: sum

    L = 0.0d0
    info = 0
    do i = 1, n
       do j = 1, i
          sum = A(i,j)
          do k = 1, j-1
             sum = sum - L(i,k)*L(j,k)
          end do
          if (i == j) then
             if (sum <= 0.0d0) then
                info = i
                return
             end if
             L(i,j) = sqrt(sum)
          else
             L(i,j) = sum / L(j,j)
          end if
       end do
    end do
  end subroutine CholeskyDecomposition

  !------------------------------------------------------------------
  ! Subroutine: MultivariateNormalSample
  ! Purpose: Sample from a multivariate normal distribution 
  !          N(0, Sigma) using the Cholesky factorization.
  !------------------------------------------------------------------
  subroutine MultivariateNormalSample(Sigma, z)
    implicit none
    real(fp_kind), intent(in) :: Sigma(:,:)
    real(fp_kind), intent(out), allocatable :: z(:)
    integer :: n, info, i
    real(fp_kind), allocatable :: L(:,:), r(:)

    n = size(Sigma, 1)
    allocate(L(n,n))
    allocate(r(n))
    call CholeskyDecomposition(n, Sigma, L, info)
    if (info /= 0) then
       print *, "Error: Cholesky decomposition failed at row", info
       stop
    end if
    do i = 1, n
       call NormalRandom(r(i))
    end do
    allocate(z(n))
    z = matmul(L, r)
    deallocate(L, r)
  end subroutine MultivariateNormalSample

  !------------------------------------------------------------------
  ! Function: GetProposalCovariance
  ! Purpose: Construct a diagonal covariance matrix for the proposal,
  !          using the prior sigma vector. (This can be replaced by a full
  !          covariance if available.)
  !------------------------------------------------------------------
  function GetProposalCovariance(prior_sigma,l) result(Sigma)
    implicit none
    real(fp_kind), intent(in) :: prior_sigma(:)
    integer(int32), intent(in) :: l
    real(fp_kind), allocatable :: Sigma(:,:)
    integer :: n, i
    n = size(prior_sigma)
    allocate(Sigma(n,n))
    Sigma = 0.0d0
    do i = 1, n
       Sigma(i,i) = GetProposalSigma(prior_sigma(i),l)**2
    end do
  end function GetProposalCovariance

  !------------------------------------------------------------------
  ! Subroutine: ProposeNewMV
  ! Purpose: Generate a new candidate for a block using a multivariate 
  !          normal proposal: proposed_block = current_block + z,
  !          where z ~ N(0, Sigma). The covariance matrix Sigma is obtained
  !          from the prior sigma vector for the block.
  ! Inputs:
  !   current_block : current parameter values in the block.
  !   prior_sigma   : prior sigma values for each parameter in the block.
  ! Output:
  !   proposed_block: newly proposed parameter vector.
  !------------------------------------------------------------------
  subroutine ProposeNewMV(current_block, prior_sigma, proposed_block,l)
    implicit none
    real(fp_kind), intent(in) :: current_block(:)
    real(fp_kind), intent(in) :: prior_sigma(:)
    integer(int32), intent(in) :: l
    real(fp_kind), intent(out), allocatable :: proposed_block(:)
    real(fp_kind), allocatable :: Sigma(:,:), z(:)
    integer :: block_size

    block_size = size(current_block)
    if (size(prior_sigma) /= block_size) then
       print *, "Error: Size mismatch in ProposeNewMV."
       stop
    end if

    ! Build the proposal covariance matrix.
    Sigma = GetProposalCovariance(prior_sigma,l)
    call MultivariateNormalSample(Sigma, z)
    allocate(proposed_block(block_size))
    proposed_block = current_block + z
    deallocate(Sigma, z)
  end subroutine ProposeNewMV

end module MultivariateProposal
