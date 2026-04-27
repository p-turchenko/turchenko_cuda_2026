module matrix_kernel_module
    use cudafor
    implicit none
    
    integer, parameter :: TILE = 16
    integer, parameter :: DP = kind(1.d0)
    
contains

    ! ============================================================
    ! ПЛОХОЕ НАИВНОЕ ЯДРО (bad naive)
    ! i = f(ty), j = f(tx) → некоалесцированный доступ
    ! ============================================================
    attributes(global) subroutine gemm_bad_naive_kernel(a, b, c, m, n, k)
        real(DP), device, intent(in)  :: a(:,:), b(:,:)
        real(DP), device, intent(out) :: c(:,:)
        integer, value, intent(in) :: m, n, k
        integer :: i, j, p
        real(DP) :: acc
        
        i = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        j = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        
        if (i <= m .and. j <= n) then
            acc = 0.d0
            do p = 1, k
                acc = acc + a(i, p) * b(p, j)
            end do
            c(i, j) = acc
        end if
    end subroutine gemm_bad_naive_kernel

    ! ============================================================
    ! ХОРОШЕЕ НАИВНОЕ ЯДРО (good naive)
    ! i = f(tx), j = f(ty) → коалесцированный доступ
    ! ============================================================
    attributes(global) subroutine gemm_good_naive_kernel(a, b, c, m, n, k)
        real(DP), device, intent(in)  :: a(:,:), b(:,:)
        real(DP), device, intent(out) :: c(:,:)
        integer, value, intent(in) :: m, n, k
        integer :: i, j, p
        real(DP) :: acc
        
        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        
        if (i <= m .and. j <= n) then
            acc = 0.d0
            do p = 1, k
                acc = acc + a(i, p) * b(p, j)
            end do
            c(i, j) = acc
        end if
    end subroutine gemm_good_naive_kernel

    ! ============================================================
    ! ТАЙЛОВОЕ ЯДРО (tiled with shared memory)
    ! ============================================================
    attributes(global) subroutine gemm_tiled_kernel(a, b, c, m, n, k)
        real(DP), device, intent(in)  :: a(:,:), b(:,:)
        real(DP), device, intent(out) :: c(:,:)
        integer, value, intent(in) :: m, n, k
        real(DP), shared :: as(TILE, TILE), bs(TILE, TILE)
        integer :: tx, ty, i, j, p, t, kk
        real(DP) :: acc
        
        tx = threadIdx%x
        ty = threadIdx%y
        i = (blockIdx%x - 1) * TILE + tx
        j = (blockIdx%y - 1) * TILE + ty
        acc = 0.d0
        
        do t = 1, k, TILE
            ! Загрузка тайла A
            kk = t + ty - 1
            if (i <= m .and. kk <= k) then
                as(tx, ty) = a(i, kk)
            else
                as(tx, ty) = 0.d0
            end if
            
            ! Загрузка тайла B
            kk = t + tx - 1
            if (kk <= k .and. j <= n) then
                bs(tx, ty) = b(kk, j)
            else
                bs(tx, ty) = 0.d0
            end if
            
            call syncthreads()
            
            ! Перемножение тайлов
            do p = 1, TILE
                acc = acc + as(tx, p) * bs(p, ty)
            end do
            
            call syncthreads()
        end do
        
        if (i <= m .and. j <= n) then
            c(i, j) = acc
        end if
    end subroutine gemm_tiled_kernel

end module matrix_kernel_module

program gemm_main
    use cudafor
    use matrix_kernel_module
    implicit none
    
    integer, parameter :: N = 1024
    integer, parameter :: M = N, K = N, N_mat = N
    integer :: i, j, istat
    real(DP), allocatable :: A(:,:), B(:,:), C_cpu(:,:)
    real(DP), allocatable :: C_bad_gpu(:,:), C_good_gpu(:,:), C_tiled_gpu(:,:)
    real(DP), device, allocatable :: A_d(:,:), B_d(:,:), C_d(:,:)
    real(DP) :: max_diff_bad, max_diff_good, max_diff_tiled
    type(dim3) :: grid, threads
    real :: start_time, end_time
    real :: time_cpu, time_bad, time_good, time_tiled
    real(DP) :: gflops_cpu, gflops_bad, gflops_good, gflops_tiled
    
    print *, '============================================================'
    print *, 'GEMM: ПЕРЕМНОЖЕНИЕ МАТРИЦ (DGEMM)'
    print *, '============================================================'
    print '(a, i0, a, i0, a)', ' Матрица: ', N, ' x ', N, ' (double)'
    print '(a, i0, a)', ' Всего элементов: ', N * N, ' в каждой матрице'
    print '(a, f0.1, a)', ' Память на матрицу: ', N*N*8.0/1024.0/1024.0, ' MB'
    print *, ''
    
    allocate(A(M, K), B(K, N_mat), C_cpu(M, N_mat))
    allocate(C_bad_gpu(M, N_mat), C_good_gpu(M, N_mat), C_tiled_gpu(M, N_mat))
    allocate(A_d(M, K), B_d(K, N_mat), C_d(M, N_mat))
    
    call random_seed()
    call random_number(A)
    call random_number(B)
    
    ! ====================================================================
    ! CPU: Наивное перемножение
    ! ====================================================================
    call cpu_time(start_time)
    do i = 1, M
        do j = 1, N_mat
            C_cpu(i, j) = sum(A(i, :) * B(:, j))
        end do
    end do
    call cpu_time(end_time)
    time_cpu = end_time - start_time
    gflops_cpu = 2.d0 * M * N_mat * K / time_cpu / 1e9
    
    ! Копирование на GPU
    A_d = A
    B_d = B
    
    ! ====================================================================
    ! GPU: Плохое наивное ядро (bad access)
    ! ====================================================================
    threads = dim3(16, 16, 1)
    grid = dim3((N_mat + 15) / 16, (M + 15) / 16, 1)
    
    call cpu_time(start_time)
    call gemm_bad_naive_kernel<<<grid, threads>>>(A_d, B_d, C_d, M, N_mat, K)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_bad = end_time - start_time
    gflops_bad = 2.d0 * M * N_mat * K / time_bad / 1e9
    C_bad_gpu = C_d
    
    ! ====================================================================
    ! GPU: Хорошее наивное ядро (good access)
    ! ====================================================================
    call cpu_time(start_time)
    call gemm_good_naive_kernel<<<grid, threads>>>(A_d, B_d, C_d, M, N_mat, K)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_good = end_time - start_time
    gflops_good = 2.d0 * M * N_mat * K / time_good / 1e9
    C_good_gpu = C_d
    
    ! ====================================================================
    ! GPU: Тайловое ядро (tiled)
    ! ====================================================================
    threads = dim3(TILE, TILE, 1)
    grid = dim3((M + TILE - 1) / TILE, (N_mat + TILE - 1) / TILE, 1)
    
    call cpu_time(start_time)
    call gemm_tiled_kernel<<<grid, threads>>>(A_d, B_d, C_d, M, N_mat, K)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_tiled = end_time - start_time
    gflops_tiled = 2.d0 * M * N_mat * K / time_tiled / 1e9
    C_tiled_gpu = C_d
    
    ! ====================================================================
    ! Сравнение результатов
    ! ====================================================================
    max_diff_bad = maxval(abs(C_cpu - C_bad_gpu))
    max_diff_good = maxval(abs(C_cpu - C_good_gpu))
    max_diff_tiled = maxval(abs(C_cpu - C_tiled_gpu))
    
    ! ====================================================================
    ! Вывод результатов
    ! ====================================================================
    print *, '============================================================'
    print *, 'РЕЗУЛЬТАТЫ'
    print *, '============================================================'
    print '(a, f10.3, a, f10.2, a)', ' CPU naive:    ', time_cpu * 1000, ' ms, ', gflops_cpu, ' GFLOPS'
    print *, ''
    print '(a, f10.3, a, f10.2, a)', ' GPU bad:      ', time_bad * 1000, ' ms, ', gflops_bad, ' GFLOPS'
    print '(a, f10.3, a, f10.2, a)', ' GPU good:     ', time_good * 1000, ' ms, ', gflops_good, ' GFLOPS'
    print '(a, f10.3, a, f10.2, a)', ' GPU tiled:    ', time_tiled * 1000, ' ms, ', gflops_tiled, ' GFLOPS'
    print *, ''
    
    print '(a, f5.1, a)', ' Ускорение good vs bad:    ', time_bad / time_good, 'x'
    print '(a, f5.1, a)', ' Ускорение tiled vs good:  ', time_good / time_tiled, 'x'
    print *, ''
    
    print '(a, e10.2)', ' Макс. отклонение (bad vs CPU):   ', max_diff_bad
    print '(a, e10.2)', ' Макс. отклонение (good vs CPU):  ', max_diff_good
    print '(a, e10.2)', ' Макс. отклонение (tiled vs CPU): ', max_diff_tiled
    
    if (max(max_diff_bad, max_diff_good, max_diff_tiled) < 1e-8) then
        print *, '✓ ВСЕ РЕЗУЛЬТАТЫ СОВПАДАЮТ!'
    else
        print *, '✗ ОШИБКА! Результаты не совпадают!'
    end if
    
    print *, ''
    print *, '============================================================'
    print *, 'ПРИМЕЧАНИЕ:'
    print *, '  - bad vs good: разница только в порядке индексов'
    print *, '  - tiled использует shared memory'
    print *, '============================================================'
    
    deallocate(A, B, C_cpu, C_bad_gpu, C_good_gpu, C_tiled_gpu, A_d, B_d, C_d)
    
end program gemm_main
