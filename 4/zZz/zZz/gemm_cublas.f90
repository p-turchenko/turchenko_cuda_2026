module matrix_kernel_module
    use cudafor
    implicit none
    
    integer, parameter :: TILE = 16
    integer, parameter :: DP = kind(1.d0)
    
contains

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
            kk = t + ty - 1
            if (i <= m .and. kk <= k) then
                as(tx, ty) = a(i, kk)
            else
                as(tx, ty) = 0.d0
            end if
            
            kk = t + tx - 1
            if (kk <= k .and. j <= n) then
                bs(tx, ty) = b(kk, j)
            else
                bs(tx, ty) = 0.d0
            end if
            
            call syncthreads()
            
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
    
    ! Внешние функции cuBLAS
    interface
        function cublasCreate(handle) bind(C, name='cublasCreate_v2')
            use iso_c_binding
            integer(c_int) :: cublasCreate
            integer(c_intptr_t) :: handle
        end function cublasCreate
        
        function cublasDestroy(handle) bind(C, name='cublasDestroy_v2')
            use iso_c_binding
            integer(c_int) :: cublasDestroy
            integer(c_intptr_t) :: handle
        end function cublasDestroy
        
        function cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
            bind(C, name='cublasDgemm_v2')
            use iso_c_binding
            integer(c_int) :: cublasDgemm
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: transa, transb
            integer(c_int), value :: m, n, k
            real(c_double), value :: alpha
            type(c_ptr), value :: A
            integer(c_int), value :: lda
            type(c_ptr), value :: B
            integer(c_int), value :: ldb
            real(c_double), value :: beta
            type(c_ptr), value :: C
            integer(c_int), value :: ldc
        end function cublasDgemm
    end interface
    
    integer, parameter :: N = 1024
    integer, parameter :: M = N, K = N, N_mat = N
    integer :: i, j, istat
    real(DP), allocatable, target :: A(:,:), B(:,:), C_cpu(:,:)
    real(DP), allocatable :: C_bad_gpu(:,:), C_good_gpu(:,:), C_tiled_gpu(:,:), C_cublas_gpu(:,:)
    real(DP), device, allocatable, target :: A_d(:,:), B_d(:,:), C_d(:,:)
    real(DP) :: alpha, beta
    real(DP) :: max_diff_bad, max_diff_good, max_diff_tiled, max_diff_cublas
    type(dim3) :: grid, threads
    real :: start_time, end_time
    real :: time_cpu, time_bad, time_good, time_tiled, time_cublas
    real(DP) :: gflops_cpu, gflops_bad, gflops_good, gflops_tiled, gflops_cublas
    integer(c_intptr_t) :: handle
    integer(c_int) :: status
    integer(c_int), parameter :: CUBLAS_OP_N = 0
    
    print *, '============================================================'
    print *, 'GEMM: ПЕРЕМНОЖЕНИЕ МАТРИЦ (DGEMM) с cuBLAS'
    print *, '============================================================'
    print '(a, i0, a, i0, a)', ' Матрица: ', N, ' x ', N, ' (double)'
    print '(a, i0, a)', ' Всего элементов: ', N * N, ' в каждой матрице'
    print '(a, f0.1, a)', ' Память на матрицу: ', N*N*8.0/1024.0/1024.0, ' MB'
    print *, ''
    
    allocate(A(M, K), B(K, N_mat), C_cpu(M, N_mat))
    allocate(C_bad_gpu(M, N_mat), C_good_gpu(M, N_mat), C_tiled_gpu(M, N_mat), C_cublas_gpu(M, N_mat))
    allocate(A_d(M, K), B_d(K, N_mat), C_d(M, N_mat))
    
    call random_seed()
    call random_number(A)
    call random_number(B)
    
    ! CPU
    call cpu_time(start_time)
    do i = 1, M
        do j = 1, N_mat
            C_cpu(i, j) = sum(A(i, :) * B(:, j))
        end do
    end do
    call cpu_time(end_time)
    time_cpu = end_time - start_time
    gflops_cpu = 2.d0 * M * N_mat * K / time_cpu / 1e9
    
    A_d = A
    B_d = B
    
    ! GPU bad
    threads = dim3(16, 16, 1)
    grid = dim3((N_mat + 15) / 16, (M + 15) / 16, 1)
    
    call cpu_time(start_time)
    call gemm_bad_naive_kernel<<<grid, threads>>>(A_d, B_d, C_d, M, N_mat, K)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_bad = end_time - start_time
    gflops_bad = 2.d0 * M * N_mat * K / time_bad / 1e9
    C_bad_gpu = C_d
    
    ! GPU good
    call cpu_time(start_time)
    call gemm_good_naive_kernel<<<grid, threads>>>(A_d, B_d, C_d, M, N_mat, K)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_good = end_time - start_time
    gflops_good = 2.d0 * M * N_mat * K / time_good / 1e9
    C_good_gpu = C_d
    
    ! GPU tiled
    threads = dim3(TILE, TILE, 1)
    grid = dim3((M + TILE - 1) / TILE, (N_mat + TILE - 1) / TILE, 1)
    
    call cpu_time(start_time)
    call gemm_tiled_kernel<<<grid, threads>>>(A_d, B_d, C_d, M, N_mat, K)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_tiled = end_time - start_time
    gflops_tiled = 2.d0 * M * N_mat * K / time_tiled / 1e9
    C_tiled_gpu = C_d
    
    ! GPU cuBLAS
    status = cublasCreate(handle)
    alpha = 1.d0
    beta = 0.d0
    
    call cpu_time(start_time)
    status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                         M, N_mat, K, alpha, c_loc(A_d), M, &
                         c_loc(B_d), K, beta, c_loc(C_d), M)
    istat = cudaDeviceSynchronize()
    call cpu_time(end_time)
    time_cublas = end_time - start_time
    gflops_cublas = 2.d0 * M * N_mat * K / time_cublas / 1e9
    C_cublas_gpu = C_d
    
    status = cublasDestroy(handle)
    
    ! Сравнение
    max_diff_bad = maxval(abs(C_cpu - C_bad_gpu))
    max_diff_good = maxval(abs(C_cpu - C_good_gpu))
    max_diff_tiled = maxval(abs(C_cpu - C_tiled_gpu))
    max_diff_cublas = maxval(abs(C_cpu - C_cublas_gpu))
    
    print *, '============================================================'
    print *, 'РЕЗУЛЬТАТЫ'
    print *, '============================================================'
    print '(a, f10.3, a, f10.2, a)', ' CPU naive:      ', time_cpu * 1000, ' ms, ', gflops_cpu, ' GFLOPS'
    print *, ''
    print '(a, f10.3, a, f10.2, a)', ' GPU bad:        ', time_bad * 1000, ' ms, ', gflops_bad, ' GFLOPS'
    print '(a, f10.3, a, f10.2, a)', ' GPU good:       ', time_good * 1000, ' ms, ', gflops_good, ' GFLOPS'
    print '(a, f10.3, a, f10.2, a)', ' GPU tiled:      ', time_tiled * 1000, ' ms, ', gflops_tiled, ' GFLOPS'
    print '(a, f10.3, a, f10.2, a)', ' GPU cuBLAS:     ', time_cublas * 1000, ' ms, ', gflops_cublas, ' GFLOPS'
    print *, ''
    
    print '(a, f5.1, a)', ' good vs bad:      ', time_bad / time_good, 'x'
    print '(a, f5.1, a)', ' tiled vs good:    ', time_good / time_tiled, 'x'
    if (time_cublas > 0) then
        print '(a, f5.1, a)', ' cuBLAS vs tiled:  ', time_tiled / time_cublas, 'x'
    end if
    print *, ''
    
    print '(a, e10.2)', ' diff bad:    ', max_diff_bad
    print '(a, e10.2)', ' diff good:   ', max_diff_good
    print '(a, e10.2)', ' diff tiled:  ', max_diff_tiled
    print '(a, e10.2)', ' diff cuBLAS: ', max_diff_cublas
    
    if (max(max(max(max_diff_bad, max_diff_good), max_diff_tiled), max_diff_cublas) < 1e-8) then
        print *, '✓ ВСЕ РЕЗУЛЬТАТЫ СОВПАДАЮТ!'
    else
        print *, '✗ ОШИБКА! Результаты не совпадают!'
    end if
    
    print *, ''
    print *, '============================================================'
    print *, 'ПРИМЕЧАНИЕ:'
    print *, '  - bad vs good: разница только в порядке индексов (2x)'
    print *, '  - tiled использует shared memory'
    print *, '  - cuBLAS — оптимизированная библиотека NVIDIA'
    print *, '============================================================'
    
    deallocate(A, B, C_cpu, C_bad_gpu, C_good_gpu, C_tiled_gpu, C_cublas_gpu, A_d, B_d, C_d)
    
end program gemm_main
