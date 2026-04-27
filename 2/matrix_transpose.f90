module matrix_kernel_module
    use cudafor
    implicit none
contains

    ! Наивное ядро транспонирования (каждый поток обрабатывает один элемент)
    attributes(global) subroutine mat_transpose_naive_kernel(a, b, rows, cols)
        real, device, intent(in)  :: a(:,:)   ! rows x cols
        real, device, intent(out) :: b(:,:)   ! cols x rows
        integer, value, intent(in) :: rows, cols
        integer :: row, col
        
        row = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        col = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        
        if (row <= rows .and. col <= cols) then
            b(col, row) = a(row, col)
        end if
    end subroutine mat_transpose_naive_kernel

end module matrix_kernel_module

program matrix_transpose_main
    use cudafor
    use matrix_kernel_module
    implicit none
    
    integer :: rows, cols, i, j, istat
    real, allocatable :: A(:,:), B_cpu(:,:), B_gpu(:,:)
    real, device, allocatable :: A_d(:,:), B_d(:,:)
    real :: max_diff
    type(dim3) :: blocks, threads
    
    print *, '============================================='
    print *, 'НАИВНОЕ ТРАНСПОНИРОВАНИЕ МАТРИЦЫ НА GPU'
    print *, '============================================='
    print *, ''
    
    print *, 'Введите количество СТРОК матрицы:'
    read *, rows
    print *, 'Введите количество СТОЛБЦОВ матрицы:'
    read *, cols
    print *, ''
    print *, 'Размер матрицы:', rows, 'x', cols
    print *, 'Всего элементов:', rows * cols
    print *, ''
    
    ! Выделение памяти на CPU
    allocate(A(rows, cols), B_cpu(cols, rows), B_gpu(cols, rows), stat=istat)
    if (istat /= 0) stop 'Ошибка памяти на CPU'
    
    ! Выделение памяти на GPU
    allocate(A_d(rows, cols), B_d(cols, rows), stat=istat)
    if (istat /= 0) stop 'Ошибка памяти на GPU'
    
    ! Заполнение матрицы A случайными числами
    call random_seed()
    call random_number(A)
    
    ! ====================================================================
    ! CPU транспонирование
    ! ====================================================================
    do i = 1, rows
        do j = 1, cols
            B_cpu(j, i) = A(i, j)
        end do
    end do
    
    ! ====================================================================
    ! GPU наивное транспонирование
    ! ====================================================================
    A_d = A
    
    ! Настройка сетки: один поток на элемент
    threads = dim3(16, 16, 1)
    blocks = dim3((rows + 15) / 16, (cols + 15) / 16, 1)
    
    call mat_transpose_naive_kernel<<<blocks, threads>>>(A_d, B_d, rows, cols)
    
    istat = cudaDeviceSynchronize()
    B_gpu = B_d
    
    ! ====================================================================
    ! Сравнение результатов CPU и GPU
    ! ====================================================================
    max_diff = 0.0
    do i = 1, cols
        do j = 1, rows
            if (abs(B_cpu(i, j) - B_gpu(i, j)) > max_diff) then
                max_diff = abs(B_cpu(i, j) - B_gpu(i, j))
            end if
        end do
    end do
    
    ! ====================================================================
    ! Вывод результатов
    ! ====================================================================
    print *, '=== ИСХОДНАЯ МАТРИЦА A (', rows, 'x', cols, ') ==='
    do i = 1, rows
        write(*, '(100f8.4)') (A(i, j), j = 1, cols)
    end do
    print *, ''
    
    print *, '=== ТРАНСПОНИРОВАННАЯ МАТРИЦА B (', cols, 'x', rows, ') ==='
    do i = 1, cols
        write(*, '(100f8.4)') (B_gpu(i, j), j = 1, rows)
    end do
    print *, ''
    
    print *, '=== РЕЗУЛЬТАТЫ ==='
    print *, 'Максимальное отклонение CPU vs GPU:', max_diff
    
    if (max_diff < 1e-5) then
        print *, '✓ РЕЗУЛЬТАТЫ СОВПАДАЮТ!'
    else
        print *, '✗ ОШИБКА! Результаты не совпадают!'
    end if
    
    ! ====================================================================
    ! Очистка памяти
    ! ====================================================================
    deallocate(A, B_cpu, B_gpu, A_d, B_d)
    
    print *, ''
    print *, '============================================='
    print *, 'ПРОГРАММА ЗАВЕРШЕНА!'
    print *, '============================================='
    
end program matrix_transpose_main
