module nilsson_basis_mod
    use iso_fortran_env, only: dp => real64
    implicit none

    type basis_state
        integer :: N
        integer :: l
        integer :: Delta ! Delta = lz
        integer :: Sigma2   ! Sigma2 = 2*Sigma = +1 or -1
        integer :: omega2 ! omega2 = 2*Omega = 2*Delta + Sigma2
    end type basis_state

contains

    subroutine make_basis(Nmax, basis, nbasis)
        implicit none

        integer, intent(in) :: Nmax
        ! integer, intent(in) :: omega2_target
        type(basis_state), allocatable, intent(out) :: basis(:)
        integer, intent(out) :: nbasis

        type(basis_state), allocatable :: tmp(:)
        integer :: N, l, Delta, Sigma2
        integer :: count
        integer :: omega2_target, omega2_target_max

        allocate(tmp(10000))
        count = 0
        ! omega2_target = Nmax*2 + 1

        do N = 0, Nmax

            ! harmonic oscillator condition: N = 2*n_r + l
            ! therefore l = N, N-2, N-4, ...
            omega2_target_max = N*2 + 1
            do omega2_target = 1, omega2_target_max, 2
                do l = N, 0, -2

                    do Delta = -l, l

                        do Sigma2 = -1, 1, 2

                            if (2*Delta + Sigma2 == omega2_target) then
                                count = count + 1

                                tmp(count)%N = N
                                tmp(count)%l = l
                                tmp(count)%Delta = Delta
                                tmp(count)%Sigma2 = Sigma2
                                tmp(count)%omega2 = omega2_target
                            end if

                        end do
                    end do
                end do
            end do
        end do

        nbasis = count
        allocate(basis(nbasis))
        basis(:) = tmp(1:nbasis)

        deallocate(tmp)

    end subroutine make_basis


    subroutine print_basis(basis, nbasis)
        implicit none

        type(basis_state), intent(in) :: basis(:)
        integer, intent(in) :: nbasis
        integer :: i

        print *, " i       N       l      Delta     Sigma2     omega2"
        print *, "---------------------------------------------"

        do i = 1, nbasis
            print '(6I8)', i, basis(i)%N, basis(i)%l, basis(i)%Delta, &
                          basis(i)%Sigma2, basis(i)%omega2
        end do

        ! output to file
        open(unit=10, file="basis_states.txt", status="replace")
        write(10, '(A)') " i       N       l      Delta     Sigma2     omega2"
        write(10, '(A)') "---------------------------------------------"
        do i = 1, nbasis
            write(10, '(6I8)') i, basis(i)%N, basis(i)%l, basis(i)%Delta, basis(i)%Sigma2, basis(i)%omega2
        end do
        close(10)
            

    end subroutine print_basis

end module nilsson_basis_mod