program rowl_plus_rowl_cpu
    implicit none

    integer :: l, i
    real(8), allocatable :: A(:), B(:), C(:)
    real(8) :: t1, t2
    character(len=100) :: filename
    integer :: unit

    ! ===== ввод l =====
    print *, "Enter vector length l:"
    read *, l

    if (l <= 0) then
        print *, "Invalid l"
        stop
    end if

    ! ===== аллокация =====
    allocate(A(l), B(l), C(l))

    ! ===== генерация случайных чисел =====
    call random_seed()

    do i = 1, l
        call random_number(A(i))
        call random_number(B(i))
    end do

    ! ===== таймер старт =====
    call cpu_time(t1)

    ! ===== сложение =====
    do i = 1, l
        C(i) = A(i) + B(i)
    end do

    ! ===== таймер стоп =====
    call cpu_time(t2)

    ! ===== output файл =====
    filename = "rowl_plus_rowl_cpu.out"
    open(newunit=unit, file=filename, status="replace", action="write")

    write(unit,*) "l =", l

    write(unit,*) "Vector A:"
    do i = 1, l
        write(unit,*) A(i)
    end do

    write(unit,*) "Vector B:"
    do i = 1, l
        write(unit,*) B(i)
    end do

    write(unit,*) "Vector C = A + B:"
    do i = 1, l
        write(unit,*) C(i)
    end do

    write(unit,*) "Time (s): ", t2 - t1
    write(unit,*) "GFLOPS: ", dble(l) / (t2 - t1) / 1.0d9

    close(unit)

    print *, "Done. Output written to ", filename

    deallocate(A, B, C)

end program rowl_plus_rowl_cpu
