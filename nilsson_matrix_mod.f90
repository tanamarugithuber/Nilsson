module nilsson_matrix_mod
    use iso_fortran_env, only: dp => real64
    use nilsson_basis_mod
    implicit none
    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
    public :: pi
    type nuclear
        integer :: A 
        integer :: Z
        integer :: N
        real(dp) :: hbar_omega0 
    end type nuclear

    type nilsson_parameters
        real(dp) :: mu
        real(dp) :: eta
        real(dp) :: delta
        real(dp) :: kappa
        real(dp) :: hbar_omega0
        real(dp) :: hbar_omega
    end type nilsson_parameters

    contains
        subroutine compute_nilsson_parameters(nucleus, params)
            implicit none
            type(nuclear), intent(in) :: nucleus
            ! type(basis_state), intent(in) :: basis(:)
            ! integer, intent(in) :: size(basis)
            type(nilsson_parameters), intent(out) :: params

            real(dp) :: A, Z, N
            ! real(dp) :: hbar_omega0

            A = real(nucleus%A, dp)
            Z = real(nucleus%Z, dp)
            N = real(nucleus%N, dp)
            params%hbar_omega0 = 41.0_dp / (A**(1.0_dp/3.0_dp))


            

            params%kappa = 0.05
            params%eta = params%delta*(1.0_dp - 4.0_dp/3.0_dp * params%delta**2 &
            - 16.0_dp/27.0_dp * params%delta**3)**(-1.0_dp/6.0_dp)/params%kappa
            params%hbar_omega = params%hbar_omega0 *(1.0_dp - 4.0_dp/3.0_dp &
            * params%delta**2 - 16.0_dp/27.0_dp * params%delta**3)**(-1.0_dp/6.0_dp)

        end subroutine compute_nilsson_parameters

        subroutine compute_matrix_element(params, basis, nbasis, matrix)
            implicit none
            type(nilsson_parameters), intent(in) :: params
            real(dp) :: mu, eta, delta, kappa, hbar_omega, hbar_omega0
            type(basis_state), intent(in) :: basis(:)
            integer, intent(in) :: nbasis
            real(dp), intent(out) :: matrix(nbasis, nbasis)

            integer :: i, j
            real(dp) :: l2avg, r2y20

            mu = params%mu
            eta = params%eta
            delta = params%delta
            kappa = params%kappa
            hbar_omega = params%hbar_omega
            hbar_omega0 = params%hbar_omega0

            matrix(:,:) = 0.0_dp

            do j = 1, nbasis
                do i = 1, nbasis
                    select case (basis(j)%N)
                    case (0:2)
                        mu = 0.0_dp
                    case (3)
                        mu = 0.35_dp
                    case (4)
                        mu = 0.45_dp
                    case (5:6)
                        mu = 0.45_dp
                    case (7)
                        mu = 0.40_dp
                    end select

                    ! oscillator term
                    if (same_basis(basis(i), basis(j))) then
                        matrix(i,j) = matrix(i,j) &
                            + hbar_omega * (real(basis(j)%N, dp) + 1.5_dp)
                    end if

                    ! l dot s
                    if (basis(i)%N == basis(j)%N .and. basis(i)%l == basis(j)%l) then

                        ! diagonal l_z s_z
                        if (basis(i)%Delta == basis(j)%Delta .and. &
                            basis(i)%Sigma2 == basis(j)%Sigma2) then

                            matrix(i,j) = matrix(i,j) &
                                - 2.0_dp * kappa * hbar_omega0 &
                                * real(basis(j)%Delta, dp) &
                                * 0.5_dp * real(basis(j)%Sigma2, dp)

                        end if

                        ! l_+ s_-
                        if (basis(i)%Delta == basis(j)%Delta + 1 .and. &
                            basis(i)%Sigma2 == basis(j)%Sigma2 - 2) then

                            matrix(i,j) = matrix(i,j) &
                                - kappa * hbar_omega0 &
                                * sqrt(real((basis(j)%l - basis(j)%Delta) &
                                        * (basis(j)%l + basis(j)%Delta + 1), dp))
                        end if

                        ! l_- s_+
                        if (basis(i)%Delta == basis(j)%Delta - 1 .and. &
                            basis(i)%Sigma2 == basis(j)%Sigma2 + 2) then

                            matrix(i,j) = matrix(i,j) &
                                - kappa * hbar_omega0 &
                                * sqrt(real((basis(j)%l + basis(j)%Delta) &
                                        * (basis(j)%l - basis(j)%Delta + 1), dp))
                        end if

                    end if

                    ! l^2 term
                    if (same_basis(basis(i), basis(j))) then
                        ! l2avg = 0.5_dp * real(basis(j)%N * (basis(j)%N + 3), dp)

                        matrix(i,j) = matrix(i,j) &
                            - kappa * mu * hbar_omega0 &
                            * real(basis(j)%l * (basis(j)%l + 1), dp) 
                    end if

                    ! deformation term
                    if (basis(i)%Sigma2 == basis(j)%Sigma2 .and. &
                        basis(i)%Delta  == basis(j)%Delta) then

                        r2y20 = exp_r_2(basis(i)%N, basis(i)%l, &
                                        basis(j)%N, basis(j)%l) &
                            * compute_exp_Y20(basis(i)%l, basis(i)%Delta, &
                                                basis(j)%l, basis(j)%Delta)

                        matrix(i,j) = matrix(i,j) &
                            - delta * hbar_omega0 &
                            * (4.0_dp/3.0_dp) * sqrt(pi/5.0_dp) &
                            * r2y20
                    end if

                end do
            end do

        end subroutine compute_matrix_element

        logical function same_basis(state1, state2)
            implicit none
            type(basis_state), intent(in) :: state1, state2

            same_basis = (state1%N == state2%N) .and. &
                         (state1%l == state2%l) .and. &
                         (state1%Delta == state2%Delta) .and. &
                         (state1%Sigma2 == state2%Sigma2)

        end function same_basis

        ! function compute_r2_Y20(N1, l1, Delta1, Sigma2_1, N2, l2, Delta2, Sigma2_2) result(value)
        !     implicit none
        !     integer :: N1, l1, Delta1, Sigma2_1, N2, l2, Delta2, Sigma2_2
        !     real(dp) :: value

        !     ! This is a placeholder function. The actual implementation would depend on the specific form of the r^2 Y20 operator and the basis states.
        !     value = exp_r_2()
        ! end function compute_r2_Y20

        function compute_exp_Y20(l1, Delta1, l2, Delta2) result(value)
            implicit none
            integer :: l1, Delta1, l2, Delta2
            real(dp) :: value

            ! This is a placeholder function. The actual implementation would depend on the specific form of the Y20 operator and the basis states.
            value = sqrt(5.0_dp/4.0_dp/pi) *sqrt(real(2* l2 + 1, dp)&
            /real(2*l1 + 1, dp)) &
            * compute_Clebsch_Gordan(real(l2, dp), 2.0_dp, real(Delta2, dp), &
            0.0_dp, real(l1, dp), real(Delta1, dp)) &
            * compute_Clebsch_Gordan(real(l2, dp), 2.0_dp,&
             0.0_dp, 0.0_dp, real(l1, dp), 0.0_dp) 
        end function compute_exp_Y20

        function compute_Clebsch_Gordan(j_1, j_2, m_1, m_2, J, M) result(CG)
            implicit none
            real(dp), intent(in) :: j_1, j_2, m_1, m_2, J, M
            real(dp) :: CG, sum_over_k, pref1, pref2
            integer :: k, k_min, k_max
            integer :: a1, a2, a3, a4, a5, a6

            CG = 0.0_dp

            ! Selection rules
            if (abs((m_1 + m_2) - M) > 1.0e-12_dp) return
            if (J < abs(j_1 - j_2) .or. J > j_1 + j_2) return
            if (abs(m_1) > j_1 .or. abs(m_2) > j_2 .or. abs(M) > J) return

            ! The following Racah formula is written with factorials.
            ! Since factorial(n) = gamma(n+1), every gamma argument below is +1.
            k_min = max(0, &
                        nint(j_2 - J - m_1), &
                        nint(j_1 - J + m_2))

            k_max = min(nint(j_1 + j_2 - J), &
                        nint(j_1 - m_1), &
                        nint(j_2 + m_2))

            if (k_min > k_max) return

            sum_over_k = 0.0_dp

            do k = k_min, k_max
                a1 = k
                a2 = nint(j_1 + j_2 - J - k)
                a3 = nint(j_1 - m_1 - k)
                a4 = nint(j_2 + m_2 - k)
                a5 = nint(J - j_2 + m_1 + k)
                a6 = nint(J - j_1 - m_2 + k)

                if (min(a1,a2,a3,a4,a5,a6) < 0) cycle

                sum_over_k = sum_over_k + (-1.0_dp)**k &
                    / ( gamma(real(a1 + 1, dp)) &
                      * gamma(real(a2 + 1, dp)) &
                      * gamma(real(a3 + 1, dp)) &
                      * gamma(real(a4 + 1, dp)) &
                      * gamma(real(a5 + 1, dp)) &
                      * gamma(real(a6 + 1, dp)) )
            end do

            pref1 = sqrt( (2.0_dp*J + 1.0_dp) &
                * gamma(J + j_1 - j_2 + 1.0_dp) &
                * gamma(J - j_1 + j_2 + 1.0_dp) &
                * gamma(j_1 + j_2 - J + 1.0_dp) &
                / gamma(j_1 + j_2 + J + 2.0_dp) )

            pref2 = sqrt( gamma(J + M + 1.0_dp) &
                * gamma(J - M + 1.0_dp) &
                * gamma(j_1 - m_1 + 1.0_dp) &
                * gamma(j_1 + m_1 + 1.0_dp) &
                * gamma(j_2 - m_2 + 1.0_dp) &
                * gamma(j_2 + m_2 + 1.0_dp) )

            CG = pref1 * pref2 * sum_over_k

        end function compute_Clebsch_Gordan

        function exp_r_2(N1,l1,N2,l2) result(value)
            implicit none
            integer :: N1, l1, N2, l2,n_r2
            real(dp) :: value

            ! This is a placeholder function. The actual implementation would depend on the specific form of the radial wavefunctions and the operator exp(r^2).
            ! n_r1 = (N1 - l1) / 2
            n_r2 = (N2 - l2) / 2
            value = 0.0_dp
            if (N1 == N2 .and. l1 == l2) then
                value = N2 + 1.5_dp
            else if (N1 == N2 .and. l1 == l2 - 2) then
                value = sqrt(real((n_r2 + 1) * (n_r2 + l2 + 0.5_dp), dp))*2.0_dp
            else if (N1 == N2-2 .and. l1 == l2 ) then
                value = sqrt(real(n_r2 * (n_r2 + l2 + 0.5_dp), dp))
            else if (N1 == N2-2 .and. l1 == l2 - 2) then
                value = sqrt(real((n_r2 + 1 - 0.5_dp) * (n_r2 + l2 + 0.5_dp), dp))
            else if (N1 == N2-2 .and. l1 == l2 + 2) then
                value = sqrt(real(n_r2 *(n_r2 - 1), dp))
            end if
        end function exp_r_2
        
        subroutine diagonalize_matrix(matrix, nbasis, eigenvalues, eigenvectors)
            implicit none
            real(dp), intent(inout) :: matrix(nbasis, nbasis)
            integer, intent(in) :: nbasis
            real(dp), intent(out) :: eigenvalues(nbasis)
            real(dp), intent(out) :: eigenvectors(nbasis, nbasis)

            ! This is a placeholder subroutine. The actual implementation would depend on the specific diagonalization method used (e.g., LAPACK routines).
            ! For example, you could use the DSYEV routine from LAPACK to diagonalize the matrix and obtain the eigenvalues and eigenvectors.


        end subroutine diagonalize_matrix

        subroutine compute_j2_matrix_element(params, basis, nbasis, matrix)
            implicit none

            type(nilsson_parameters), intent(in) :: params
            type(basis_state), intent(in) :: basis(:)
            integer, intent(in) :: nbasis
            real(dp), intent(out) :: matrix(nbasis, nbasis)

            integer :: i, j
            real(dp) :: l, lambda, sigma
            real(dp) :: ls_elem

            matrix(:,:) = 0.0_dp

            do j = 1, nbasis
                do i = 1, nbasis

                    ls_elem = 0.0_dp

                    if (basis(i)%N == basis(j)%N .and. &
                        basis(i)%l == basis(j)%l) then

                        l = real(basis(j)%l, dp)
                        lambda = real(basis(j)%Delta, dp)
                        sigma = 0.5_dp * real(basis(j)%Sigma2, dp)

                        ! diagonal part: L_z S_z
                        if (basis(i)%Delta == basis(j)%Delta .and. &
                            basis(i)%Sigma2 == basis(j)%Sigma2) then

                            ls_elem = ls_elem + lambda * sigma
                        end if

                        ! 1/2 L_+ S_-
                        if (basis(i)%Delta == basis(j)%Delta + 1 .and. &
                            basis(i)%Sigma2 == basis(j)%Sigma2 - 2) then

                            ls_elem = ls_elem + 0.5_dp * &
                                sqrt(real((basis(j)%l - basis(j)%Delta) * &
                                        (basis(j)%l + basis(j)%Delta + 1), dp))
                        end if

                        ! 1/2 L_- S_+
                        if (basis(i)%Delta == basis(j)%Delta - 1 .and. &
                            basis(i)%Sigma2 == basis(j)%Sigma2 + 2) then

                            ls_elem = ls_elem + 0.5_dp * &
                                sqrt(real((basis(j)%l + basis(j)%Delta) * &
                                        (basis(j)%l - basis(j)%Delta + 1), dp))
                        end if

                    end if

                    ! J^2 = L^2 + S^2 + 2 L.S
                    matrix(i,j) = matrix(i,j) + 2.0_dp * ls_elem

                    if (same_basis(basis(i), basis(j))) then
                        matrix(i,j) = matrix(i,j) + &
                            real(basis(j)%l * (basis(j)%l + 1), dp) + 0.75_dp
                    end if

                end do
            end do
        end subroutine compute_j2_matrix_element

end module nilsson_matrix_mod