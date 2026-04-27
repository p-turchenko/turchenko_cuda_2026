    call gemm_bad_naive_kernel<<<dim3((M+15)/16,(N_mat+15)/16,1), dim3(16,16,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_bad, startEvent, stopEvent)
    time_bad = time_bad / 1000.0
    gflops_bad = 2.d0*M*N_mat*K/time_bad/1e9
    C_bad = C_d
    
    ! good kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_good_naive_kernel<<<dim3((M+15)/16,(N_mat+15)/16,1), dim3(16,16,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_good, startEvent, stopEvent)
    time_good = time_good / 1000.0
    gflops_good = 2.d0*M*N_mat*K/time_good/1e9
    C_good = C_d
    
    ! tiled kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_tiled_kernel<<<dim3((M+TILE-1)/TILE,(N_mat+TILE-1)/TILE,1), dim3(TILE,TILE,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_tiled, startEvent, stopEvent)
    time_tiled = time_tiled / 1000.0
    gflops_tiled = 2.d0*M*N_mat*K/time_tiled/1e9
    C_tiled = C_d
    
    ! cuBLAS - optimized version
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    istat = cublasDgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                           M, N_mat, K, &
                           alpha, A_d, M, B_d, K, beta, C_d, M)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_cublas, startEvent, stopEvent)
    time_cublas = time_cublas / 1000.0
    gflops_cublas = 2.d0*M*N_mat*K/time_cublas/1e9
    C_cublas = C_d
    
    print *, 'CPU:    ', time_cpu, 'sec, ', gflops_cpu, 'GFLOPS'
    print *, 'bad:    ', time_bad, 'sec, ', gflops_bad, 'GFLOPS'
    print *, 'good:   ', time_good, 'sec, ', gflops_good, 'GFLOPS'
    print *, 'tiled:  ', time_tiled, 'sec, ', gflops_tiled, 'GFLOPS'
    print *, 'cublas: ', time_cublas, 'sec, ', gflops_cublas, 'GFLOPS'
    
    ! Check results
    print *, 'Max error bad vs CPU:   ', maxval(abs(C_bad - C_cpu))
    print *, 'Max error good vs CPU:  ', maxval(abs(C_good - C_cpu))
    print *, 'Max error tiled vs CPU: ', maxval(abs(C_tiled - C_cpu))
    print *, 'Max error cublas vs CPU:', maxval(abs(C_cublas - C_cpu))
    
    ! Cleanup
    istat = cublasDestroy(handle)
    istat = cudaEventDestroy(startEvent)
    istat = cudaEventDestroy(stopEvent)
    
end program
EOF

mc
./gemm_no_cublas
cat > gemm_cublas.cuf << 'EOF' && nvfortran -O3 -o gemm_cublas gemm_cublas.cuf -I/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/lib64 -cudalib=cublas -lcudart && ./gemm_cublasmodule matrix_kernel_module
    use cudafor
    implicit none
    
    integer, parameter :: TILE = 16
    integer, parameter :: DP = kind(1.d0)
    
contains

    attributes(global) subroutine gemm_bad_naive_kernel(a, b, c, m, n, k)
        real(DP), device :: a(:,:), b(:,:), c(:,:)
        integer, value :: m, n, k
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
    end subroutine

    attributes(global) subroutine gemm_good_naive_kernel(a, b, c, m, n, k)
        real(DP), device :: a(:,:), b(:,:), c(:,:)
        integer, value :: m, n, k
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
    end subroutine

    attributes(global) subroutine gemm_tiled_kernel(a, b, c, m, n, k)
        real(DP), device :: a(:,:), b(:,:), c(:,:)
        integer, value :: m, n, k
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
    end subroutine

end module

program gemm_cublas
    use cudafor
    use cublas_v2
    use matrix_kernel_module
    implicit none
    
    integer, parameter :: N = 1024
    integer, parameter :: M = N, K = N, N_mat = N
    integer :: i, j, istat
    real(DP), allocatable :: A(:,:), B(:,:), C_cpu(:,:)
    real(DP), allocatable :: C_bad(:,:), C_good(:,:), C_tiled(:,:), C_cublas(:,:)
    real(DP), device, allocatable :: A_d(:,:), B_d(:,:), C_d(:,:)
    
    type(cudaEvent) :: startEvent, stopEvent
    type(cublasHandle) :: handle
    real :: time_bad, time_good, time_tiled, time_cpu, time_cublas
    real(DP) :: gflops_bad, gflops_good, gflops_tiled, gflops_cpu, gflops_cublas
    real(DP), parameter :: alpha = 1.d0, beta = 0.d0
    
    allocate(A(M,K), B(K,N_mat), C_cpu(M,N_mat))
    allocate(C_bad(M,N_mat), C_good(M,N_mat), C_tiled(M,N_mat), C_cublas(M,N_mat))
    allocate(A_d(M,K), B_d(K,N_mat), C_d(M,N_mat))
    
    call random_number(A)
    call random_number(B)
    
    ! Initialize CUDA events
    istat = cudaEventCreate(startEvent)
    istat = cudaEventCreate(stopEvent)
    
    ! Initialize cuBLAS handle once
    istat = cublasCreate(handle)
    if (istat /= CUBLAS_STATUS_SUCCESS) then
        write(*,*) 'cublasCreate failed'
        stop
    end if
    
    ! CPU computation
    istat = cudaEventRecord(startEvent, 0)
    do j = 1, N_mat
        do i = 1, M
            C_cpu(i,j) = sum(A(i,:) * B(:,j))
        end do
    end do
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_cpu, startEvent, stopEvent)
    time_cpu = time_cpu / 1000.0  ! Convert ms to seconds
    gflops_cpu = 2.d0*M*N_mat*K/time_cpu/1e9
    
    A_d = A
    B_d = B
    
    ! bad kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_bad_naive_kernel<<<dim3((M+15)/16,(N_mat+15)/16,1), dim3(16,16,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_bad, startEvent, stopEvent)
    time_bad = time_bad / 1000.0
    gflops_bad = 2.d0*M*N_mat*K/time_bad/1e9
    C_bad = C_d
    
    ! good kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_good_naive_kernel<<<dim3((M+15)/16,(N_mat+15)/16,1), dim3(16,16,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_good, startEvent, stopEvent)
    time_good = time_good / 1000.0
    gflops_good = 2.d0*M*N_mat*K/time_good/1e9
    C_good = C_d
    
    ! tiled kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_tiled_kernel<<<dim3((M+TILE-1)/TILE,(N_mat+TILE-1)/TILE,1), dim3(TILE,TILE,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_tiled, startEvent, stopEvent)
    time_tiled = time_tiled / 1000.0
    gflops_tiled = 2.d0*M*N_mat*K/time_tiled/1e9
    C_tiled = C_d
    
    ! cuBLAS - using pre-created handle
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    istat = cublasDgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                           M, N_mat, K, &
                           alpha, A_d, M, B_d, K, beta, C_d, M)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_cublas, startEvent, stopEvent)
    time_cublas = time_cublas / 1000.0
    gflops_cublas = 2.d0*M*N_mat*K/time_cublas/1e9
    C_cublas = C_d
    
    print *, 'Results for matrix multiplication ', M, 'x', N_mat, 'x', K
    print *, '--------------------------------------------------'
    print *, 'CPU:    ', time_cpu, 'sec, ', gflops_cpu, 'GFLOPS'
    print *, 'bad:    ', time_bad, 'sec, ', gflops_bad, 'GFLOPS'
    print *, 'good:   ', time_good, 'sec, ', gflops_good, 'GFLOPS'
    print *, 'tiled:  ', time_tiled, 'sec, ', gflops_tiled, 'GFLOPS'
    print *, 'cublas: ', time_cublas, 'sec, ', gflops_cublas, 'GFLOPS'
    
    ! Cleanup
    istat = cublasDestroy(handle)
    istat = cudaEventDestroy(startEvent)
    istat = cudaEventDestroy(stopEvent)
    
    deallocate(A, B, C_cpu)
    deallocate(C_bad, C_good, C_tiled, C_cublas)
    deallocate(A_d, B_d, C_d)
    
end program gemm_cublas
EOF

./gemm_cublas
mc
cat > gemm_cublas.cuf << 'EOF' && nvfortran -O3 -o gemm_cublas gemm_cublas.cuf -I/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/lib64 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/cuda/lib64 -cudalib=cublas -lcudart && ./gemm_cublasmodule matrix_kernel_module
    use cudafor
    implicit none
    
    integer, parameter :: TILE = 16
    integer, parameter :: DP = kind(1.d0)
    
contains

    attributes(global) subroutine gemm_bad_naive_kernel(a, b, c, m, n, k)
        real(DP), device :: a(:,:), b(:,:), c(:,:)
        integer, value :: m, n, k
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
    end subroutine

    attributes(global) subroutine gemm_good_naive_kernel(a, b, c, m, n, k)
        real(DP), device :: a(:,:), b(:,:), c(:,:)
        integer, value :: m, n, k
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
    end subroutine

    attributes(global) subroutine gemm_tiled_kernel(a, b, c, m, n, k)
        real(DP), device :: a(:,:), b(:,:), c(:,:)
        integer, value :: m, n, k
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
    end subroutine

end module

program gemm_cublas
    use cudafor
    use cublas_v2
    use matrix_kernel_module
    implicit none
    
    integer, parameter :: N = 1024
    integer, parameter :: M = N, K = N, N_mat = N
    integer :: i, j, istat
    real(DP), allocatable :: A(:,:), B(:,:), C_cpu(:,:)
    real(DP), allocatable :: C_bad(:,:), C_good(:,:), C_tiled(:,:), C_cublas(:,:)
    real(DP), device, allocatable :: A_d(:,:), B_d(:,:), C_d(:,:)
    
    type(cudaEvent) :: startEvent, stopEvent
    type(cublasHandle) :: handle
    real :: time_bad, time_good, time_tiled, time_cpu, time_cublas
    real(DP) :: gflops_bad, gflops_good, gflops_tiled, gflops_cpu, gflops_cublas
    real(DP), parameter :: alpha = 1.d0, beta = 0.d0
    
    allocate(A(M,K), B(K,N_mat), C_cpu(M,N_mat))
    allocate(C_bad(M,N_mat), C_good(M,N_mat), C_tiled(M,N_mat), C_cublas(M,N_mat))
    allocate(A_d(M,K), B_d(K,N_mat), C_d(M,N_mat))
    
    call random_number(A)
    call random_number(B)
    
    ! Initialize CUDA events
    istat = cudaEventCreate(startEvent)
    istat = cudaEventCreate(stopEvent)
    
    ! Initialize cuBLAS handle once
    istat = cublasCreate(handle)
    if (istat /= CUBLAS_STATUS_SUCCESS) then
        write(*,*) 'cublasCreate failed'
        stop
    end if
    
    A_d = A
    B_d = B
    C_d = 0.d0
    
    ! Warmup - run cuBLAS once to initialize GPU and avoid first-run overhead
    istat = cublasDgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                           M, N_mat, K, &
                           alpha, A_d, M, B_d, K, beta, C_d, M)
    istat = cudaDeviceSynchronize()
    
    ! CPU computation
    istat = cudaEventRecord(startEvent, 0)
    do j = 1, N_mat
        do i = 1, M
            C_cpu(i,j) = sum(A(i,:) * B(:,j))
        end do
    end do
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_cpu, startEvent, stopEvent)
    time_cpu = time_cpu / 1000.0  ! Convert ms to seconds
    gflops_cpu = 2.d0*M*N_mat*K/time_cpu/1e9
    
    ! bad kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_bad_naive_kernel<<<dim3((M+15)/16,(N_mat+15)/16,1), dim3(16,16,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_bad, startEvent, stopEvent)
    time_bad = time_bad / 1000.0
    gflops_bad = 2.d0*M*N_mat*K/time_bad/1e9
    C_bad = C_d
    
    ! good kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_good_naive_kernel<<<dim3((M+15)/16,(N_mat+15)/16,1), dim3(16,16,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_good, startEvent, stopEvent)
    time_good = time_good / 1000.0
    gflops_good = 2.d0*M*N_mat*K/time_good/1e9
    C_good = C_d
    
    ! tiled kernel
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    call gemm_tiled_kernel<<<dim3((M+TILE-1)/TILE,(N_mat+TILE-1)/TILE,1), dim3(TILE,TILE,1)>>>(A_d,B_d,C_d,M,N_mat,K)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_tiled, startEvent, stopEvent)
    time_tiled = time_tiled / 1000.0
    gflops_tiled = 2.d0*M*N_mat*K/time_tiled/1e9
    C_tiled = C_d
    
    ! cuBLAS - using pre-created handle
    C_d = 0.d0
    istat = cudaEventRecord(startEvent, 0)
    istat = cublasDgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                           M, N_mat, K, &
                           alpha, A_d, M, B_d, K, beta, C_d, M)
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)
    istat = cudaEventElapsedTime(time_cublas, startEvent, stopEvent)
    time_cublas = time_cublas / 1000.0
    gflops_cublas = 2.d0*M*N_mat*K/time_cublas/1e9
    C_cublas = C_d
    
    print *, 'Results for matrix multiplication ', M, 'x', N_mat, 'x', K
    print *, '--------------------------------------------------'
    print *, 'CPU:    ', time_cpu, 'sec, ', gflops_cpu, 'GFLOPS'
    print *, 'bad:    ', time_bad, 'sec, ', gflops_bad, 'GFLOPS'
    print *, 'good:   ', time_good, 'sec, ', gflops_good, 'GFLOPS'
    print *, 'tiled:  ', time_tiled, 'sec, ', gflops_tiled, 'GFLOPS'
    print *, 'cublas: ', time_cublas, 'sec, ', gflops_cublas, 'GFLOPS'
    
    ! Cleanup
    istat = cublasDestroy(handle)
    istat = cudaEventDestroy(startEvent)
    istat = cudaEventDestroy(stopEvent)
    
    deallocate(A, B, C_cpu)
    deallocate(C_bad, C_good, C_tiled, C_cublas)
    deallocate(A_d, B_d, C_d)
    
end program gemm_cublas
EOF

./gemm_no_cublas
./gemm_cublas
mc
mc
./gemm_cublas
mc
./gemm_cublas
nvprof ./gemm_cublas
./concurrent_gemm.exe
mc
