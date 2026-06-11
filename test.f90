program test_basis
    use iso_fortran_env, only: dp => real64
    use nilsson_basis_mod
    use nilsson_matrix_mod
    implicit none

    integer :: Nmax
    integer :: omega2_target
    integer :: nbasis
    type(basis_state), allocatable :: basis(:)
    type(nuclear) :: nucleus
    type(nilsson_parameters) :: params
    integer :: lwork
    integer :: info,i,k
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: H(:,:), eig(:), j2(:, :), expj2(:), j_(:)
    real(dp) :: j_1, j_2, m_1, m_2, J, M
    real(dp) :: CG
    

    Nmax = 7
    nucleus%A = 100
    ! Omega = 1/2 のブロック
    omega2_target = 1
    params%eta = 0.0
    
    ! CG =  compute_exp_Y20(0,0,2,0)
    ! print *, "Spherical harmonic coefficient: ", CG

    call make_basis(Nmax, basis, nbasis)
    call print_basis(basis, nbasis)

    print *, "Number of basis states: ", nbasis
    allocate(H(nbasis,nbasis))
    allocate(eig(nbasis))
    allocate(j2(nbasis, nbasis))
    allocate(expj2(nbasis))
    allocate(j_(nbasis))
    print *, "Computing matrix elements..."
    call compute_nilsson_parameters(nucleus, params)
    call compute_matrix_element(params, basis, nbasis, H)
    print *, "Computing j2 matrix elements..."
    call compute_j2_matrix_element(params, basis, nbasis, j2)

     ! output the Hamiltonian matrix to file
    
    lwork = max(1,3*nbasis)
    allocate(work(lwork))

    call dsyev( &
        'V', &
        'U', &
        nbasis, &
        H, &
        nbasis, &
        eig, &
        work, &
        lwork, &
        info )

    eig(:) = eig(:) / params%hbar_omega
    expj2(:) = 0.0_dp
        do i = 1, nbasis
            expj2(i) = dot_product( &
             H(:,i), &
             matmul(j2,H(:,i)))
             j_(i) = 0.5_dp * (sqrt(1.0_dp + 4.0_dp * expj2(i)) - 1.0_dp)
             write(*,*) "State ", i, ": Eigenvalue = ", eig(i), " j2 = ", expj2(i), " j = ", j_(i)

             ! output the eigenvalues and corresponding basis states to file
            open(unit=20, file="eigenvalues.txt", status="replace")
            write(20, '(A)') " i    Eigen_energy   j^2   j"
            write(20, '(A)') "-----------------------------------------"
            do k = 1, nbasis
            write(20, '(I5, F15.6, 2F15.6)') k, eig(k), expj2(k), j_(k)
            end do
            close(20)
        end do


    

    !  ! print eigenvalues and corresponding basis states
    !  print *, "Eigen_energies (in units of hbar*omega):"
    ! do i = 1, nbasis
    !     print *, "  ", i, eig(i), basis(i)%N, basis(i)%l, basis(i)%Delta, basis(i)%Sigma2, basis(i)%omega2
    ! end do
    !     ! print eigenvalues and corresponding basis states to file
    !     open(unit=20, file="eigenvalues.txt", status="replace")
    !     write(20, '(A)') " i    Eigen_energy   N   l   Delta   Sigma2   omega2"
    !     write(20, '(A)') "---------------------------------------------------"
    !     do i = 1, nbasis
    !         write(20, '(I5, F15.6, 5I8)') i, eig(i), basis(i)%N, basis(i)%l, basis(i)%Delta, basis(i)%Sigma2, basis(i)%omega2
    !     end do
    !     close(20)
end program test_basis