program matrixkl_plus_matrixkl_cpu
    implicit none

    integer :: k, l, i, j
    real(8), allocatable :: A(:,:), B(:,:), C(:,:)
    real(8) :: t1, t2
    character(len=100) :: filename
    integer :: unit

    ! ===== ввод k и l =====
    print *, "Enter matrix dimensions k (rows) and l (cols):"
    read *, k, l

    if (k <= 0 .or. l <= 0) then
        print *, "Invalid k or l"
        stop
    end if

    ! ===== аллокация =====
    allocate(A(k,l), B(k,l), C(k,l))

    ! ===== генерация случайных чисел =====
    call random_seed()

    do j = 1, l
        do i = 1, k
            call random_number(A(i,j))
            call random_number(B(i,j))
        end do
    end do

    ! ===== таймер старт =====
    call cpu_time(t1)

    ! ===== сложение =====
    do j = 1, l
        do i = 1, k
            C(i,j) = A(i,j) + B(i,j)
        end do
    end do

    ! ===== таймер стоп =====
    call cpu_time(t2)

    ! ===== output файл =====
    filename = "matrixkl_plus_matrixkl_cpu.out"
    open(newunit=unit, file=filename, status="replace", action="write")

    write(unit,*) "k =", k
    write(unit,*) "l =", l

    write(unit,*) "Matrix A:"
    do i = 1, k
        write(unit,*) (A(i,j), j=1,l)
    end do

    write(unit,*) "Matrix B:"
    do i = 1, k
        write(unit,*) (B(i,j), j=1,l)
    end do

    write(unit,*) "Matrix C = A + B:"
    do i = 1, k
        write(unit,*) (C(i,j), j=1,l)
    end do

    write(unit,*) "Time (s): ", t2 - t1
    write(unit,*) "GFLOPS: ", dble(k*l) / (t2 - t1) / 1.0d9

    close(unit)

    print *, "Done. Output written to ", filename

    deallocate(A, B, C)

end program matrixkl_plus_matrixkl_cpu
