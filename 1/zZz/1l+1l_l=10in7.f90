program row_add_cpu
    implicit none

    integer, parameter :: l = 10000000
    real(8), allocatable :: A(:), B(:), C(:)
    integer :: i
    real(8) :: t1, t2

    allocate(A(l), B(l), C(l))

    ! Инициализация
    do i = 1, l
        A(i) = 1.0d0
        B(i) = 2.0d0
    end do

    ! Таймер старт
    call cpu_time(t1)

    ! Сложение
    do i = 1, l
        C(i) = A(i) + B(i)
    end do

    ! Таймер стоп
    call cpu_time(t2)

    print *, "Time (s): ", t2 - t1

    ! GFLOPS
    print *, "GFLOPS: ", dble(l) / (t2 - t1) / 1.0d9

    deallocate(A, B, C)

end program row_add_cpu