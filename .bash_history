    print *, '========================================================================'
    
    deallocate(h_A, h_B, h_C_serial, h_C_concurrent)
    deallocate(d_A, d_B_all, d_C_all)
    
end program gemm_streams_main
EOF

cd ~/5 && cat > gemm_streams.cuf << 'EOF' && nvfortran -cuda -o gemm_streams gemm_streams.cuf -lcublas -lcudart && ./gemm_streams
program gemm_streams_main
    use cudafor
    implicit none

    integer, parameter :: DP = kind(1.d0)
    
    ! Размеры матриц (константы)
    integer, parameter :: M = 512
    integer, parameter :: N = 512
    integer, parameter :: K_LARGE = 4096
    integer, parameter :: NSTREAMS = 4
    integer, parameter :: SZPART_MAX = (K_LARGE + NSTREAMS - 1) / NSTREAMS
    
    integer :: i, j, istat, istream, szpart, shift
    real(DP), allocatable, target :: h_A(:,:), h_B(:,:)
    real(DP), allocatable, target :: h_C_serial(:,:), h_C_concurrent(:,:)
    real(DP), device, allocatable, target :: d_A(:,:)
    real(DP), device, allocatable, target :: d_B_all(:,:), d_C_all(:,:)
    real(DP) :: alpha, beta, max_diff
    real :: start_time, end_time
    real :: time_serial, time_concurrent
    real(DP) :: gflops_serial, gflops_concurrent
    
    ! CUDA Streams и cuBLAS handles
    integer(kind=cuda_stream_kind) :: stream(NSTREAMS)
    integer(kind=c_intptr_t) :: handle(NSTREAMS)
    
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
            integer(c_intptr_t), value :: handle
        end function cublasDestroy
        
        function cublasSetStream(handle, stream) bind(C, name='cublasSetStream_v2')
            use iso_c_binding
            integer(c_int) :: cublasSetStream
            integer(c_intptr_t), value :: handle
            integer(cuda_stream_kind), value :: stream
        end function cublasSetStream
        
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
    
    integer(c_int), parameter :: CUBLAS_OP_N = 0
    
    print *, '========================================================================'
    print *, 'CUDA STREAMS: SERIAL vs CONCURRENT GEMM'
    print *, '========================================================================'
    print '(a, i0, a, i0, a)', ' A: ', M, ' x ', N
    print '(a, i0, a, i0, a)', ' B: ', N, ' x ', K_LARGE
    print '(a, i0, a, i0, a)', ' C: ', M, ' x ', K_LARGE
    print '(a, i0)',            ' Количество streams: ', NSTREAMS
    print '(a, i0, a)',         ' Максимальный размер части: ', SZPART_MAX, ' столбцов'
    print *, ''
    
    ! ====================================================================
    ! Выделение памяти на CPU
    ! ====================================================================
    allocate(h_A(M, N), h_B(N, K_LARGE))
    allocate(h_C_serial(M, K_LARGE), h_C_concurrent(M, K_LARGE))
    
    ! ====================================================================
    ! Выделение памяти на GPU (с константными размерами)
    ! ====================================================================
    allocate(d_A(M, N))
    allocate(d_B_all(N, SZPART_MAX * NSTREAMS))
    allocate(d_C_all(M, SZPART_MAX * NSTREAMS))
    
    ! ====================================================================
    ! Инициализация случайными числами
    ! ====================================================================
    call random_seed()
    call random_number(h_A)
    call random_number(h_B)
    
    ! ====================================================================
    ! SERIAL MODEL
    ! ====================================================================
    print *, '=== SERIAL MODEL (синхронная обработка) ==='
    print *, ''
    
    d_A = h_A
    alpha = 1.d0
    beta  = 0.d0
    
    call cpu_time(start_time)
    
    do istream = 0, NSTREAMS - 1
        szpart = K_LARGE / NSTREAMS
        if (istream == NSTREAMS - 1) szpart = szpart + mod(K_LARGE, NSTREAMS)
        shift = (K_LARGE / NSTREAMS) * istream
        
        print '(a, i0, a, i0, a, i0)', '  Часть ', istream, ': columns ', shift+1, '-', shift+szpart
        
        ! Синхронное копирование части B
        d_B_all(1:N, istream*SZPART_MAX+1 : istream*SZPART_MAX+szpart) = &
            h_B(1:N, shift+1 : shift+szpart)
        
        ! GEMM
        call cublasDgemm(0_c_intptr_t, CUBLAS_OP_N, CUBLAS_OP_N, &
                         M, szpart, N, alpha, &
                         c_loc(d_A), M, &
                         c_loc(d_B_all(1, istream*SZPART_MAX+1)), N, &
                         beta, &
                         c_loc(d_C_all(1, istream*SZPART_MAX+1)), M)
        
        istat = cudaDeviceSynchronize()
        
        ! Синхронное копирование результата
        h_C_serial(1:M, shift+1 : shift+szpart) = &
            d_C_all(1:M, istream*SZPART_MAX+1 : istream*SZPART_MAX+szpart)
    end do
    
    call cpu_time(end_time)
    time_serial = end_time - start_time
    gflops_serial = 2.d0 * M * N * K_LARGE / time_serial / 1e9
    
    print '(a, f10.3, a)', '  Время: ', time_serial * 1000, ' ms'
    print '(a, f10.2, a)', '  GFLOPS: ', gflops_serial, ''
    print *, ''
    
    ! ====================================================================
    ! CONCURRENT MODEL
    ! ====================================================================
    print *, '=== CONCURRENT MODEL (асинхронная обработка) ==='
    print *, ''
    
    do istream = 1, NSTREAMS
        istat = cudaStreamCreate(stream(istream))
        istat = cublasCreate(handle(istream))
        istat = cublasSetStream(handle(istream), stream(istream))
    end do
    
    d_A = h_A
    alpha = 1.d0
    beta  = 0.d0
    
    call cpu_time(start_time)
    
    do istream = 0, NSTREAMS - 1
        szpart = K_LARGE / NSTREAMS
        if (istream == NSTREAMS - 1) szpart = szpart + mod(K_LARGE, NSTREAMS)
        shift = (K_LARGE / NSTREAMS) * istream
        
        print '(a, i0, a, i0, a, i0)', '  Часть ', istream, ': columns ', shift+1, '-', shift+szpart
        
        ! Асинхронное копирование H2D
        istat = cudaMemcpyAsync(d_B_all(1, istream*SZPART_MAX+1), &
                                c_loc(h_B(1, shift+1)), &
                                N * szpart * 8, &
                                cudaMemcpyHostToDevice, &
                                stream(istream+1))
        
        ! GEMM в том же stream
        istat = cublasDgemm(handle(istream+1), CUBLAS_OP_N, CUBLAS_OP_N, &
                            M, szpart, N, alpha, &
                            c_loc(d_A), M, &
                            c_loc(d_B_all(1, istream*SZPART_MAX+1)), N, &
                            beta, &
                            c_loc(d_C_all(1, istream*SZPART_MAX+1)), M)
        
        ! Асинхронное копирование D2H
        istat = cudaMemcpyAsync(c_loc(h_C_concurrent(1, shift+1)), &
                                d_C_all(1, istream*SZPART_MAX+1), &
                                M * szpart * 8, &
                                cudaMemcpyDeviceToHost, &
                                stream(istream+1))
    end do
    
    istat = cudaDeviceSynchronize()
    
    call cpu_time(end_time)
    time_concurrent = end_time - start_time
    gflops_concurrent = 2.d0 * M * N * K_LARGE / time_concurrent / 1e9
    
    print '(a, f10.3, a)', '  Время: ', time_concurrent * 1000, ' ms'
    print '(a, f10.2, a)', '  GFLOPS: ', gflops_concurrent, ''
    print *, ''
    
    do istream = 1, NSTREAMS
        istat = cublasDestroy(handle(istream))
        istat = cudaStreamDestroy(stream(istream))
    end do
    
    ! ====================================================================
    ! Проверка
    ! ====================================================================
    print *, '=== ПРОВЕРКА ==='
    max_diff = 0.d0
    do i = 1, M
        do j = 1, K_LARGE
            if (abs(h_C_serial(i, j) - h_C_concurrent(i, j)) > max_diff) then
                max_diff = abs(h_C_serial(i, j) - h_C_concurrent(i, j))
            end if
        end do
    end do
    
    print '(a, e10.2)', '  Макс. отклонение: ', max_diff
    if (max_diff < 1e-8) then
        print *, '  ✓ РЕЗУЛЬТАТЫ СОВПАДАЮТ'
    else
        print *, '  ✗ ОШИБКА'
    end if
    print *, ''
    
    ! ====================================================================
    ! Итоги
    ! ====================================================================
    print *, '========================================================================'
    print *, 'РЕЗУЛЬТАТЫ'
    print *, '========================================================================'
    print '(a, f10.3, a, f12.2, a)', ' SERIAL:     ', time_serial * 1000, ' ms, ', gflops_serial, ' GFLOPS'
    print '(a, f10.3, a, f12.2, a)', ' CONCURRENT: ', time_concurrent * 1000, ' ms, ', gflops_concurrent, ' GFLOPS'
    print *, ''
    print '(a, f5.2, a)', ' УСКОРЕНИЕ: ', time_serial / time_concurrent, 'x'
    print *, ''
    print *, '========================================================================'
    
    deallocate(h_A, h_B, h_C_serial, h_C_concurrent)
    deallocate(d_A, d_B_all, d_C_all)
    
end program gemm_streams_main
EOF

cd ~/5 && cat > gemm_streams.cuf << 'EOF' && nvfortran -cuda -o gemm_streams gemm_streams.cuf -lcublas -lcudart && ./gemm_streams
program gemm_streams_main
    use cudafor
    implicit none

    integer, parameter :: DP = kind(1.d0)
    
    ! Размеры матриц
    integer, parameter :: M = 512
    integer, parameter :: N = 512
    integer, parameter :: K_LARGE = 4096
    integer, parameter :: NSTREAMS = 4
    integer, parameter :: SZPART_MAX = (K_LARGE + NSTREAMS - 1) / NSTREAMS
    
    integer :: i, j, istat, istream, szpart, shift
    real(DP), allocatable, target :: h_A(:,:), h_B(:,:)
    real(DP), allocatable, target :: h_C_serial(:,:), h_C_concurrent(:,:)
    real(DP), device, allocatable, target :: d_A(:,:)
    real(DP), device, allocatable, target :: d_B(:,:), d_C(:,:)
    real(DP) :: alpha, beta, max_diff
    real :: start_time, end_time
    real :: time_serial, time_concurrent
    real(DP) :: gflops_serial, gflops_concurrent
    
    ! CUDA Streams и cuBLAS handles
    integer(kind=cuda_stream_kind) :: stream(NSTREAMS)
    integer(kind=c_intptr_t) :: handle(NSTREAMS)
    
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
            integer(c_intptr_t), value :: handle
        end function cublasDestroy
        
        function cublasSetStream(handle, stream) bind(C, name='cublasSetStream_v2')
            use iso_c_binding
            integer(c_int) :: cublasSetStream
            integer(c_intptr_t), value :: handle
            integer(cuda_stream_kind), value :: stream
        end function cublasSetStream
        
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
    
    integer(c_int), parameter :: CUBLAS_OP_N = 0
    
    print *, '========================================================================'
    print *, 'CUDA STREAMS: SERIAL vs CONCURRENT GEMM'
    print *, '========================================================================'
    print '(a, i0, a, i0, a)', ' A: ', M, ' x ', N
    print '(a, i0, a, i0, a)', ' B: ', N, ' x ', K_LARGE
    print '(a, i0, a, i0, a)', ' C: ', M, ' x ', K_LARGE
    print '(a, i0)',            ' Количество streams: ', NSTREAMS
    print '(a, i0, a)',         ' Максимальный размер части: ', SZPART_MAX, ' столбцов'
    print *, ''
    
    ! ====================================================================
    ! Выделение памяти на CPU
    ! ====================================================================
    allocate(h_A(M, N), h_B(N, K_LARGE))
    allocate(h_C_serial(M, K_LARGE), h_C_concurrent(M, K_LARGE))
    
    ! ====================================================================
    ! Выделение памяти на GPU (отдельные буферы для каждой части)
    ! ====================================================================
    allocate(d_A(M, N))
    allocate(d_B(N, SZPART_MAX), d_C(M, SZPART_MAX))
    
    ! ====================================================================
    ! Инициализация случайными числами
    ! ====================================================================
    call random_seed()
    call random_number(h_A)
    call random_number(h_B)
    
    ! ====================================================================
    ! SERIAL MODEL
    ! ====================================================================
    print *, '=== SERIAL MODEL (синхронная обработка) ==='
    print *, ''
    
    d_A = h_A
    alpha = 1.d0
    beta  = 0.d0
    
    call cpu_time(start_time)
    
    do istream = 0, NSTREAMS - 1
        szpart = K_LARGE / NSTREAMS
        if (istream == NSTREAMS - 1) szpart = szpart + mod(K_LARGE, NSTREAMS)
        shift = (K_LARGE / NSTREAMS) * istream
        
        print '(a, i0, a, i0, a, i0)', '  Часть ', istream, ': columns ', shift+1, '-', shift+szpart
        
        ! Синхронное копирование части B
        d_B(:, 1:szpart) = h_B(:, shift+1 : shift+szpart)
        
        ! GEMM
        call cublasDgemm(0_c_intptr_t, CUBLAS_OP_N, CUBLAS_OP_N, &
                         M, szpart, N, alpha, &
                         c_loc(d_A), M, &
                         c_loc(d_B), N, &
                         beta, &
                         c_loc(d_C), M)
        
        istat = cudaDeviceSynchronize()
        
        ! Синхронное копирование результата
        h_C_serial(:, shift+1 : shift+szpart) = d_C(:, 1:szpart)
    end do
    
    call cpu_time(end_time)
    time_serial = end_time - start_time
    gflops_serial = 2.d0 * M * N * K_LARGE / time_serial / 1e9
    
    print '(a, f10.3, a)', '  Время: ', time_serial * 1000, ' ms'
    print '(a, f10.2, a)', '  GFLOPS: ', gflops_serial, ''
    print *, ''
    
    ! ====================================================================
    ! CONCURRENT MODEL
    ! ====================================================================
    print *, '=== CONCURRENT MODEL (асинхронная обработка) ==='
    print *, ''
    
    do istream = 1, NSTREAMS
        istat = cudaStreamCreate(stream(istream))
        istat = cublasCreate(handle(istream))
        istat = cublasSetStream(handle(istream), stream(istream))
    end do
    
    d_A = h_A
    alpha = 1.d0
    beta  = 0.d0
    
    call cpu_time(start_time)
    
    do istream = 0, NSTREAMS - 1
        szpart = K_LARGE / NSTREAMS
        if (istream == NSTREAMS - 1) szpart = szpart + mod(K_LARGE, NSTREAMS)
        shift = (K_LARGE / NSTREAMS) * istream
        
        print '(a, i0, a, i0, a, i0)', '  Часть ', istream, ': columns ', shift+1, '-', shift+szpart
        
        ! Асинхронное копирование H2D
        istat = cudaMemcpyAsync(d_B, c_loc(h_B(1, shift+1)), &
                                N * szpart * 8, &
                                cudaMemcpyHostToDevice, &
                                stream(istream+1))
        
        ! GEMM в том же stream
        istat = cublasDgemm(handle(istream+1), CUBLAS_OP_N, CUBLAS_OP_N, &
                            M, szpart, N, alpha, &
                            c_loc(d_A), M, &
                            c_loc(d_B), N, &
                            beta, &
                            c_loc(d_C), M)
        
        ! Асинхронное копирование D2H
        istat = cudaMemcpyAsync(c_loc(h_C_concurrent(1, shift+1)), &
                                d_C, &
                                M * szpart * 8, &
                                cudaMemcpyDeviceToHost, &
                                stream(istream+1))
    end do
    
    istat = cudaDeviceSynchronize()
    
    call cpu_time(end_time)
    time_concurrent = end_time - start_time
    gflops_concurrent = 2.d0 * M * N * K_LARGE / time_concurrent / 1e9
    
    print '(a, f10.3, a)', '  Время: ', time_concurrent * 1000, ' ms'
    print '(a, f10.2, a)', '  GFLOPS: ', gflops_concurrent, ''
    print *, ''
    
    do istream = 1, NSTREAMS
        istat = cublasDestroy(handle(istream))
        istat = cudaStreamDestroy(stream(istream))
    end do
    
    ! ====================================================================
    ! Проверка
    ! ====================================================================
    print *, '=== ПРОВЕРКА ==='
    max_diff = 0.d0
    do i = 1, M
        do j = 1, K_LARGE
            if (abs(h_C_serial(i, j) - h_C_concurrent(i, j)) > max_diff) then
                max_diff = abs(h_C_serial(i, j) - h_C_concurrent(i, j))
            end if
        end do
    end do
    
    print '(a, e10.2)', '  Макс. отклонение: ', max_diff
    if (max_diff < 1e-8) then
        print *, '  ✓ РЕЗУЛЬТАТЫ СОВПАДАЮТ'
    else
        print *, '  ✗ ОШИБКА'
    end if
    print *, ''
    
    ! ====================================================================
    ! Итоги
    ! ====================================================================
    print *, '========================================================================'
    print *, 'РЕЗУЛЬТАТЫ'
    print *, '========================================================================'
    print '(a, f10.3, a, f12.2, a)', ' SERIAL:     ', time_serial * 1000, ' ms, ', gflops_serial, ' GFLOPS'
    print '(a, f10.3, a, f12.2, a)', ' CONCURRENT: ', time_concurrent * 1000, ' ms, ', gflops_concurrent, ' GFLOPS'
    print *, ''
    print '(a, f5.2, a)', ' УСКОРЕНИЕ: ', time_serial / time_concurrent, 'x'
    print *, ''
    print *, '========================================================================'
    
    deallocate(h_A, h_B, h_C_serial, h_C_concurrent)
    deallocate(d_A, d_B, d_C)
    
end program gemm_streams_main
EOF

