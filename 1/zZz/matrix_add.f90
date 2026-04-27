cat > matrix_add.f90 << 'EOF' && nvfortran -cuda -o matrix_add matrix_add.f90 && ./matrix_add
module matrix_kernel_module
    use cudafor
    implicit none
contains
    attributes(global) subroutine add_matrix_kernel(a, b, c, n, m)
        real, device :: a(:,:), b(:,:), c(:,:)
        integer, value :: n, m
        integer :: i, j
        i = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        j = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        if (i <= n .and. j <= m) then
            c(i, j) = a(i, j) + b(i, j)
        end if
    end subroutine add_matrix_kernel
end module matrix_kernel_module
program matrix_add_interactive
    use cudafor
    use matrix_kernel_module
    implicit none
    integer :: n, m, i, j, istat
    real, allocatable :: A(:,:), B(:,:), C_cpu(:,:), C_gpu(:,:)
    real, device, allocatable :: A_d(:,:), B_d(:,:), C_d(:,:)
    real :: max_diff
    print *, '============================================='
    print *, 'СЛОЖЕНИЕ МАТРИЦ НА CPU И GPU'
    print *, '============================================='
    print *, ''
    print *, 'Введите количество СТРОК матрицы:'
    read *, n
    print *, 'Введите количество СТОЛБЦОВ матрицы:'
    read *, m
    print *, ''
    print *, 'Размер матрицы:', n, 'x', m
    print *, 'Всего элементов:', n * m
    print *, ''
    allocate(A(n,m), B(n,m), C_cpu(n,m), C_gpu(n,m), stat=istat)
    if (istat /= 0) stop 'Ошибка памяти на CPU'
    allocate(A_d(n,m), B_d(n,m), C_d(n,m), stat=istat)
    if (istat /= 0) stop 'Ошибка памяти на GPU'
    call random_seed()
    call random_number(A)
    call random_number(B)
    
    ! Вывод матрицы A
    print *, '=== МАТРИЦА A ==='
    do i = 1, n
        write(*, '(100f10.4)') (A(i, j), j = 1, m)
    end do
    print *, ''
    
    ! Вывод матрицы B
    print *, '=== МАТРИЦА B ==='
    do i = 1, n
        write(*, '(100f10.4)') (B(i, j), j = 1, m)
    end do
    print *, ''
    
    ! CPU сложение
    do i = 1, n
        do j = 1, m
            C_cpu(i,j) = A(i,j) + B(i,j)
        end do
    end do
    
    ! GPU сложение
    A_d = A
    B_d = B
    
    call add_matrix_kernel<<<dim3((m+15)/16, (n+15)/16, 1), dim3(16,16,1)>>>(A_d, B_d, C_d, n, m)
    
    istat = cudaDeviceSynchronize()
    C_gpu = C_d
    
    ! Сравнение
    max_diff = 0.0
    do i = 1, n
        do j = 1, m
            if (abs(C_cpu(i,j) - C_gpu(i,j)) > max_diff) then
                max_diff = abs(C_cpu(i,j) - C_gpu(i,j))
            end if
        end do
    end do
    
    print *, '=== РЕЗУЛЬТАТЫ ==='
    print *, 'Максимальное отклонение CPU vs GPU:', max_diff
    
    if (max_diff < 1e-5) then
        print *, '✓ РЕЗУЛЬТАТЫ СОВПАДАЮТ!'
    else
        print *, '✗ ОШИБКА! Результаты не совпадают!'
    end if
    
    ! Вывод результирующей матрицы (GPU)
    print *, ''
    print *, '=== РЕЗУЛЬТИРУЮЩАЯ МАТРИЦА C = A + B (GPU) ==='
    do i = 1, n
        write(*, '(100f10.4)') (C_gpu(i, j), j = 1, m)
    end do
    
    deallocate(A, B, C_cpu, C_gpu, A_d, B_d, C_d)
    
    print *, ''
    print *, '============================================='
    print *, 'ПРОГРАММА ЗАВЕРШЕНА!'
    print *, '============================================='
end program matrix_add_interactive
EOF
