module matrix_kernel_module
    use cudafor
    implicit none
    
    integer, parameter :: TILE_DIM = 32
    integer, parameter :: BLOCK_ROWS = 8
    
contains

    ! ============================================================
    ! Наивное ядро (без shared memory)
    ! ============================================================
    attributes(global) subroutine mat_transpose_naive_kernel(a, b, rows, cols)
        real, device, intent(in)  :: a(:,:)
        real, device, intent(out) :: b(:,:)
        integer, value, intent(in) :: rows, cols
        integer :: row, col
        
        row = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        col = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        
        if (row <= rows .and. col <= cols) then
            b(col, row) = a(row, col)
        end if
    end subroutine mat_transpose_naive_kernel

    ! ============================================================
    ! Тайловое транспонирование с shared memory (оптимальное)
    ! ============================================================
    attributes(global) subroutine mat_transpose_tiled_kernel(a, b, rows, cols)
        real, device, intent(in)  :: a(:,:)
        real, device, intent(out) :: b(:,:)
        integer, value, intent(in) :: rows, cols
        real, shared :: tile(TILE_DIM, TILE_DIM)
        integer :: tx, ty, bx, by, row, col
        integer :: i
        
        tx = threadIdx%x
        ty = threadIdx%y
        bx = blockIdx%x
        by = blockIdx%y
        
        row = (bx - 1) * TILE_DIM + tx
        col = (by - 1) * TILE_DIM + ty
        
        do i = 0, TILE_DIM - 1, BLOCK_ROWS
            if (row <= rows .and. col + i <= cols) then
                tile(tx, ty + i) = a(row, col + i)
            end if
        end do
        call syncthreads()
        
        col = (by - 1) * TILE_DIM + tx
        row = (bx - 1) * TILE_DIM + ty
        
        do i = 0, TILE_DIM - 1, BLOCK_ROWS
            if (col <= cols .and. row + i <= rows) then
                b(col, row + i) = tile(ty + i, tx)
            end if
        end do
    end subroutine mat_transpose_tiled_kernel

    ! ============================================================
    ! Эталонное копирование через shared memory (для сравнения)
    ! ============================================================
    attributes(global) subroutine copy_shared_kernel(odata, idata, rows, cols)
        real, device, intent(in)  :: idata(:,:)
        real, device, intent(out) :: odata(:,:)
        integer, value, intent(in) :: rows, cols
        real, shared :: tile(TILE_DIM, TILE_DIM)
        integer :: tx, ty, bx, by, row, col
        integer :: i
        
        tx = threadIdx%x
        ty = threadIdx%y
        bx = blockIdx%x
        by = blockIdx%y
        
        row = (bx - 1) * TILE_DIM + tx
        col = (by - 1) * TILE_DIM + ty
        
        do i = 0, TILE_DIM - 1, BLOCK_ROWS
            if (row <= rows .and. col + i <= cols) then
                tile(tx, ty + i) = idata(row, col + i)
            end if
        end do
        call syncthreads()
        
        do i = 0, TILE_DIM - 1, BLOCK_ROWS
            if (row <= rows .and. col + i <= cols) then
                odata(row, col + i) = tile(tx, ty + i)
            end if
        end do
    end subroutine copy_shared_kernel

end module matrix_kernel_module

program matrix_transpose_main
    use cudafor
    use matrix_kernel_module
    implicit none
    
    integer :: rows, cols, i, j, istat
    real, allocatable :: A(:,:), B_naive(:,:), B_tiled(:,:)
    real, device, allocatable :: A_d(:,:), B_naive_d(:,:), B_tiled_d(:,:)
    real :: max_diff_tiled
    type(dim3) :: blocks_naive, blocks_tiled, threads_naive, threads_tiled
    real :: start_time, end_time
    real :: time_naive, time_tiled
    
    print *, '============================================================'
    print *, 'ОПТИМАЛЬНОЕ ТРАНСПОНИРОВАНИЕ МАТРИЦЫ С SHARED MEMORY'
    print *, '============================================================'
    print *, ''
    
    print *, 'Введите количество СТРОК матрицы:'
    read *, rows
    print *, 'Введите количество СТОЛБЦОВ матрицы:'
    read *, cols
    print *, ''
    print *, 'Размер матрицы:', rows, 'x', cols
    print *, 'Всего элементов:', rows * cols
    print *, ''
    
    allocate(A(rows, cols), B_naive(cols, rows), B_tiled(cols, rows), stat=istat)
    if (istat /= 0) stop 'Ошибка памяти на CPU'
    
    allocate(A_d(rows, cols), B_naive_d(cols, rows), B_tiled_d(cols, rows), stat=istat)
    if (istat /= 0) stop 'Ошибка памяти на GPU'
    
    call random_seed()
    call random_number(A)
    
    ! CPU транспонирование
    do i = 1, rows
        do j = 1, cols
            B_naive(j, i) = A(i, j)
        end do
    end do
    
    A_d = A
    
    ! ====================================================================
    ! GPU: Наивное транспонирование
    ! ====================================================================
    threads_naive = dim3(16, 16, 1)
    blocks_naive = dim3((rows + 15) / 16, (cols + 15) / 16, 1)
    
    call cpu_time(start_time)
    call mat_transpose_naive_kernel<<<blocks_naive, threads_naive>>>(A_d, B_naive_d, rows, cols)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_naive = end_time - start_time
    B_naive = B_naive_d
    
    ! ====================================================================
    ! GPU: Тайловое транспонирование
    ! ====================================================================
    threads_tiled = dim3(TILE_DIM, BLOCK_ROWS, 1)
    blocks_tiled = dim3((rows + TILE_DIM - 1) / TILE_DIM, (cols + TILE_DIM - 1) / TILE_DIM, 1)
    
    call cpu_time(start_time)
    call mat_transpose_tiled_kernel<<<blocks_tiled, threads_tiled>>>(A_d, B_tiled_d, rows, cols)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_tiled = end_time - start_time
    B_tiled = B_tiled_d
    
    ! ====================================================================
    ! Сравнение результатов
    ! ====================================================================
    max_diff_tiled = 0.0
    do i = 1, cols
        do j = 1, rows
            if (abs(B_naive(i, j) - B_tiled(i, j)) > max_diff_tiled) then
                max_diff_tiled = abs(B_naive(i, j) - B_tiled(i, j))
            end if
        end do
    end do
    
    ! ====================================================================
    ! Вывод
    ! ====================================================================
    print *, '=== ИСХОДНАЯ МАТРИЦА A (', rows, 'x', cols, ') ==='
    do i = 1, min(rows, 10)
        write(*, '(100f8.4)') (A(i, j), j = 1, min(cols, 10))
    end do
    if (rows > 10 .or. cols > 10) print *, '... (показаны первые 10x10)'
    print *, ''
    
    print *, '=== ТРАНСПОНИРОВАННАЯ МАТРИЦА (тайловое ядро) ==='
    do i = 1, min(cols, 10)
        write(*, '(100f8.4)') (B_tiled(i, j), j = 1, min(rows, 10))
    end do
    if (cols > 10 .or. rows > 10) print *, '... (показаны первые 10x10)'
    print *, ''
    
    print *, '============================================================'
    print *, 'РЕЗУЛЬТАТЫ СРАВНЕНИЯ'
    print *, '============================================================'
    print '(a, f10.4, a)', ' Наивное ядро:          ', time_naive * 1000, ' мс'
    print '(a, f10.4, a)', ' Тайловое ядро (32x8):   ', time_tiled * 1000, ' мс'
    print *, ''
    
    if (time_tiled < time_naive) then
        print '(a, f5.1, a)', ' ✓ Тайловое ядро БЫСТРЕЕ наивного в ', time_naive / time_tiled, ' раз(а)'
    else
        print '(a, f5.1, a)', ' ✗ Тайловое ядро МЕДЛЕННЕЕ наивного в ', time_tiled / time_naive, ' раз(а)'
    end if
    print *, ''
    
    print *, 'Максимальное отклонение (наивное vs тайловое):', max_diff_tiled
    if (max_diff_tiled < 1e-5) then
        print *, '✓ РЕЗУЛЬТАТЫ СОВПАДАЮТ между ядрами!'
    else
        print *, '✗ ОШИБКА! Результаты не совпадают!'
    end if
    
    deallocate(A, B_naive, B_tiled, A_d, B_naive_d, B_tiled_d)
    
    print *, ''
    print *, '============================================================'
    print *, 'ПРОГРАММА ЗАВЕРШЕНА!'
    print *, '============================================================'
    
end program matrix_transpose_main
