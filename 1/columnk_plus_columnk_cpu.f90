program columnk_plus_columnk_cpu
    implicit none

    integer :: k, i
    real(8), allocatable :: A(:), B(:), C(:)
    real(8) :: t1, t2
    character(len=100) :: filename
    integer :: unit

    ! ===== ввод k =====
    print *, "Enter vector length k:"
    read *, k

    if (k <= 0) then
        print *, "Invalid k"
        stop
    end if

    ! ===== аллокация =====
    allocate(A(k), B(k), C(k))

    ! ===== генерация случайных чисел =====
    call random_seed()

    do i = 1, k
        call random_number(A(i))
        call random_number(B(i))
    end do

    ! ===== таймер старт =====
    call cpu_time(t1)

    ! ===== сложение =====
    do i = 1, k
        C(i) = A(i) + B(i)
    end do

    ! ===== таймер стоп =====
    call cpu_time(t2)

    ! ===== output файл =====
    filename = "columnk_plus_columnk_cpu.out"
    open(newunit=unit, file=filename, status="replace", action="write")

    write(unit,*) "k =", k

    write(unit,*) "Vector A:"
    do i = 1, k
        write(unit,*) A(i)
    end do

    write(unit,*) "Vector B:"
    do i = 1, k
        write(unit,*) B(i)
    end do

    write(unit,*) "Vector C = A + B:"
    do i = 1, k
        write(unit,*) C(i)
    end do

    write(unit,*) "Time (s): ", t2 - t1
    write(unit,*) "GFLOPS: ", dble(k) / (t2 - t1) / 1.0d9

    close(unit)

    print *, "Done. Output written to ", filename

    deallocate(A, B, C)

end program columnk_plus_columnk_cpu
