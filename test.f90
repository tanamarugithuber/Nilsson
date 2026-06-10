program test_basis
    use iso_fortran_env, only: dp => real64
    use nilsson_basis_mod
    implicit none

    integer :: Nmax
    integer :: omega2_target
    integer :: nbasis
    type(basis_state), allocatable :: basis(:)

    Nmax = 6

    ! Omega = 1/2 のブロック
    omega2_target = 1

    call make_basis(Nmax, basis, nbasis)
    call print_basis(basis, nbasis)

end program test_basis