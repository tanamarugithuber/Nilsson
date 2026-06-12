program test_basis
    use iso_fortran_env, only: dp => real64
    use nilsson_basis_mod
    use nilsson_matrix_mod
    implicit none

    integer :: Nmax
    integer :: omega2_target
    integer :: nbasis
    integer :: vec_unit
    type(basis_state), allocatable :: basis(:)
    type(nuclear) :: nucleus
    type(nilsson_parameters) :: params
    integer :: lwork
    integer :: info,i,k,l,N_label
    real(dp) :: delta
    integer  :: state_index
    integer  :: basis_index
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: H(:,:), eig(:), j2(:, :), expj2(:), j_(:)
    real(dp) :: j_1, j_2, m_1, m_2, j_exp, M,N_exp, j2_exp, omega_exp, sigma_exp
    real(dp) :: CG
    character(len=20) :: delta_str
    

    Nmax = 7
    nucleus%A = 100
    ! ! Omega = 1/2 のブロック
    ! omega2_target = 1
    params%delta = 0.0_dp
    
    ! CG =  compute_exp_Y20(0,0,2,0)
    ! print *, "Spherical harmonic coefficient: ", CG

    call make_basis(Nmax, basis, nbasis)
    call print_basis(basis, nbasis)
    print *, "Number of basis states: ", nbasis
    allocate(H(nbasis,nbasis))
    allocate(eig(nbasis))
    allocate(j2(nbasis, nbasis))
    ! allocate(expj2(nbasis))
    ! allocate(j_(nbasis))
    lwork = max(1,3*nbasis)
    allocate(work(lwork))

    do l = 1, 61
        params%delta = -0.3_dp + 0.01_dp * real(l-1, dp)
    

        print *, "Computing matrix elements..."
        call compute_nilsson_parameters(nucleus, params)
        call compute_matrix_element(params, basis, nbasis, H)
        print *, "Computing j2 matrix elements..."
        call compute_j2_matrix_element(params, basis, nbasis, j2)

        ! output the Hamiltonian matrix to file
        


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

        if (info /= 0) then
            print *, "DSYEV failed. info = ", info
            stop
        end if
        
        write(delta_str, '(SP,F5.2)') params%delta
        delta_str = adjustl(delta_str)

        open(newunit=vec_unit, file='eigenvector_'//trim(delta_str)//'.dat', status="replace")

        write(vec_unit,'(A)') "# delta state basis N l Lambda Sigma2 coeff"

        do k = 1, nbasis          ! k: 固有状態番号
            do i = 1, nbasis      ! i: 基底番号

                write(vec_unit,'(F8.3,2I6,4I6,F15.8)') &
                    params%delta, &
                    k, &
                    i, &
                    basis(i)%N, &
                    basis(i)%l, &
                    basis(i)%Delta, &
                    basis(i)%Sigma2, &
                    H(i,k)

            end do
        end do

        close(vec_unit)


        eig(:) = eig(:) / params%hbar_omega
        j2_exp = 0.0_dp
        

        open(unit=20, file="eigenvalues"//trim(delta_str)//".dat", status="replace")
                write(20, '(A)') " i    Eigen_energy   N  j^2   j sigma omega"
                write(20, '(A)') "-----------------------------------------"
            do i = 1, nbasis
                j2_exp = dot_product( &
                H(:,i), &
                matmul(j2,H(:,i)))

                j_exp = 0.5_dp * (sqrt(1.0_dp + 4.0_dp * j2_exp) - 1.0_dp)
                j_exp = 0.5_dp * nint(2.0_dp * j_exp)


                N_exp = 0.0_dp
                omega_exp = 0.0_dp
                sigma_exp = 0.0_dp
                do k = 1, nbasis
                    N_exp = N_exp + H(k,i)**2 * real(basis(k)%N, dp)
                    omega_exp = omega_exp + H(k,i)**2 * (  real(basis(k)%delta, dp) &
                    + 0.5_dp * real(basis(k)%Sigma2, dp) )
                    sigma_exp = sigma_exp + H(k,i)**2 * real(basis(k)%Sigma2, dp) * 0.5_dp
                end do
                N_label = nint(N_exp)
                ! sigma_exp = 0.5_dp * nint(sigma_exp*2.0_dp)
                ! omega_exp = 0.5_dp * nint(omega_exp*2.0_dp)
                write(*,*) "State ", i, ": Eigenvalue = ", eig(i), " N = ", N_exp, " j2 = ", j2_exp, " j = ", j_exp

                ! output the eigenvalues and corresponding basis states to file
                
                write(20, '(I5, F15.6, I5, F15.6, F15.6, F15.6, F15.6)') i, eig(i), N_label, j2_exp, j_exp, sigma_exp, omega_exp
                
            
            end do
            close(20)
            H(:,:) = 0.0_dp
            j2(:,:) = 0.0_dp
            work(:) = 0.0_dp
            j2_exp = 0.0_dp
            j_exp = 0.0_dp


    end do

end program test_basis