module stream_kernel_module
    use cudafor
    implicit none
    
    integer, parameter :: DP = kind(1.d0)
    
contains

    ! Простое ядро для проверки (не используется в GEMM, но оставлено для совместимости)
    attributes(global) subroutine dummy_kernel(a, b, c, n, m)
        real(DP), device, intent(in)  :: a(:,:), b(:,:)
        real(DP), device, intent(out) :: c(:,:)
        integer, value, intent(in) :: n, m
        integer :: i, j
        
        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        
        if (i <= n .and. j <= m) then
            c(i, j) = a(i, j) + b(i, j)
        end if
    end subroutine dummy_kernel

end module stream_kernel_module

program gemm_streams_main
    use cudafor
    use cublas
    use stream_kernel_module
    implicit none
    
    integer, parameter :: M = 512
    integer, parameter :: N = 512
    integer, parameter :: K_LARGE = 4096
    integer, parameter :: NSTREAMS = 4
    
    integer :: i, istat, istream, szpart, shift, szpart_max
    integer :: i0, j0, p
    real(DP), allocatable, target :: h_A(:,:), h_B(:,:), h_C_serial(:,:), h_C_concurrent(:,:)
    real(DP), device, allocatable, target :: d_A(:,:), d_B_all(:,:), d_C_all(:,:)
    real(DP) :: alpha, beta
    real(DP) :: max_diff
    real :: start_time, end_time
    real :: time_serial, time_concurrent
    real(DP) :: gflops_serial, gflops_concurrent
    
    ! CUDA Streams и cuBLAS handles
    integer(cuda_stream_kind), allocatable :: stream(:)
    type(cublasHandle),        allocatable :: handle(:)
    
    ! Указатели на слоты в больших буферах
    real(DP), device, pointer :: d_B_ptr(:,:), d_C_ptr(:,:)
    
    print *, '========================================================================'
    print *, 'CUDA STREAMS: SERIAL vs CONCURRENT GEMM'
    print *, '========================================================================'
    print '(a, i0, a, i0, a)', ' A: ', M, ' x ', N
    print '(a, i0, a, i0, a)', ' B: ', N, ' x ', K_LARGE
    print '(a, i0)',            ' C: ', M, ' x ', K_LARGE
    print '(a, i0)',            ' Количество streams: ', NSTREAMS
    print *, ''
    print '(a, f0.1, a)', ' Память на матрицу A: ', M*N*8.0/1024.0/1024.0, ' MB'
    print '(a, f0.1, a)', ' Память на матрицу B: ', N*K_LARGE*8.0/1024.0/1024.0, ' MB'
    print '(a, f0.1, a)', ' Память на матрицу C: ', M*K_LARGE*8.0/1024.0/1024.0, ' MB'
    print *, ''
    
    ! ====================================================================
    ! Выделение памяти на CPU
    ! ====================================================================
    allocate(h_A(M, N), h_B(N, K_LARGE), h_C_serial(M, K_LARGE), h_C_concurrent(M, K_LARGE))
    
    ! ====================================================================
    ! Выделение памяти на GPU
    ! ====================================================================
    allocate(d_A(M, N))
    
    ! Максимальный размер части (округление вверх для безопасности)
    szpart_max = (K_LARGE + NSTREAMS - 1) / NSTREAMS
    print '(a, i0, a)', ' Максимальный размер части: ', szpart_max, ' столбцов'
    print *, ''
    
    allocate(d_B_all(N, szpart_max * NSTREAMS))
    allocate(d_C_all(M, szpart_max * NSTREAMS))
    
    ! ====================================================================
    ! Инициализация случайными числами
    ! ====================================================================
    call random_seed()
    call random_number(h_A)
    call random_number(h_B)
    
    ! ====================================================================
    ! Копирование A на GPU (один раз для всех частей)
    ! ====================================================================
    d_A = h_A
    
    ! ====================================================================
    ! SERIAL MODEL: все части обрабатываются последовательно
    ! ====================================================================
    print *, '=== SERIAL MODEL (синхронная обработка) ==='
    print *, ''
    
    alpha = 1.d0
    beta  = 0.d0
    
    call cpu_time(start_time)
    
    do istream = 0, NSTREAMS - 1
        ! Вычисление размера и смещения части
        szpart = K_LARGE / NSTREAMS
        if (istream < mod(K_LARGE, NSTREAMS)) szpart = szpart + 1
        
        shift = (K_LARGE / NSTREAMS) * istream
        if (istream >= mod(K_LARGE, NSTREAMS)) shift = shift + mod(K_LARGE, NSTREAMS)
        
        print '(a, i0, a, i0, a, i0)', '  Stream ', istream, ': szpart=', szpart, ' shift=', shift
        
        ! Синхронное копирование H2D
        d_B_all(1:N, istream*szpart_max+1 : istream*szpart_max+szpart) = &
            h_B(1:N, shift+1 : shift+szpart)
        
        ! Указатель на слот
        d_B_ptr => d_B_all(1:N, istream*szpart_max+1 : istream*szpart_max+szpart)
        d_C_ptr => d_C_all(1:M, istream*szpart_max+1 : istream*szpart_max+szpart)
        
        ! GEMM через cuBLAS (используем стандартный handle)
        call cublasDgemm('N', 'N', M, szpart, N, alpha, d_A, M, d_B_ptr, N, beta, d_C_ptr, M)
        
        ! Синхронизация перед D2H
        istat = cudaDeviceSynchronize()
        
        ! Синхронное копирование D2H
        h_C_serial(1:M, shift+1 : shift+szpart) = d_C_ptr(1:M, 1:szpart)
    end do
    
    call cpu_time(end_time)
    time_serial = end_time - start_time
    gflops_serial = 2.d0 * M * N * K_LARGE / time_serial / 1e9
    
    print '(a, f10.3, a)', '  Полное время: ', time_serial * 1000, ' ms'
    print '(a, f10.2, a)', '  Производительность: ', gflops_serial, ' GFLOPS'
    print *, ''
    
    ! ====================================================================
    ! CONCURRENT MODEL: асинхронная обработка с независимыми streams
    ! ====================================================================
    print *, '=== CONCURRENT MODEL (асинхронная обработка) ==='
    print *, ''
    
    ! Выделение памяти под streams и handles
    allocate(stream(NSTREAMS), handle(NSTREAMS))
    
    ! Создание streams и привязка cuBLAS handles
    do istream = 0, NSTREAMS - 1
        istat = cudaStreamCreate(stream(istream))
        istat = cublasCreate(handle(istream))
        istat = cublasSetStream(handle(istream), stream(istream))
        
        if (istat /= 0) then
            print *, 'Ошибка создания stream/handle для', istream
            stop
        end if
    end do
    
    call cpu_time(start_time)
    
    do istream = 0, NSTREAMS - 1
        ! Вычисление размера и смещения части
        szpart = K_LARGE / NSTREAMS
        if (istream < mod(K_LARGE, NSTREAMS)) szpart = szpart + 1
        
        shift = (K_LARGE / NSTREAMS) * istream
        if (istream >= mod(K_LARGE, NSTREAMS)) shift = shift + mod(K_LARGE, NSTREAMS)
        
        print '(a, i0, a, i0, a, i0)', '  Stream ', istream, ': szpart=', szpart, ' shift=', shift
        
        ! Указатель на слот в большом буфере
        d_B_ptr => d_B_all(1:N, istream*szpart_max+1 : istream*szpart_max+szpart)
        d_C_ptr => d_C_all(1:M, istream*szpart_max+1 : istream*szpart_max+szpart)
        
        ! Асинхронное копирование H2D (не блокирует CPU)
        istat = cudaMemcpyAsync(d_B_ptr, c_loc(h_B(1, shift+1)), &
                                N * szpart * 8, cudaMemcpyHostToDevice, stream(istream))
        
        ! GEMM в том же stream (начнётся после завершения H2D)
        istat = cublasDgemm(handle(istream), CUBLAS_OP_N, CUBLAS_OP_N, &
                            M, szpart, N, alpha, d_A, M, d_B_ptr, N, beta, d_C_ptr, M)
        
        ! Асинхронное копирование D2H (начнётся после GEMM)
        istat = cudaMemcpyAsync(c_loc(h_C_concurrent(1, shift+1)), d_C_ptr, &
                                M * szpart * 8, cudaMemcpyDeviceToHost, stream(istream))
    end do
    
    ! Одна синхронизация в конце — ждём завершения ВСЕХ streams
    istat = cudaDeviceSynchronize()
    
    call cpu_time(end_time)
    time_concurrent = end_time - start_time
    gflops_concurrent = 2.d0 * M * N * K_LARGE / time_concurrent / 1e9
    
    print '(a, f10.3, a)', '  Полное время: ', time_concurrent * 1000, ' ms'
    print '(a, f10.2, a)', '  Производительность: ', gflops_concurrent, ' GFLOPS'
    print *, ''
    
    ! ====================================================================
    ! Уничтожение streams и handles
    ! ====================================================================
    do istream = 0, NSTREAMS - 1
        istat = cublasDestroy(handle(istream))
        istat = cudaStreamDestroy(stream(istream))
    end do
    
    deallocate(stream, handle)
    
    ! ====================================================================
    ! Сравнение результатов serial и concurrent
    ! ====================================================================
    print *, '=== ПРОВЕРКА КОРРЕКТНОСТИ ==='
    print *, ''
    
    max_diff = 0.d0
    do i = 1, M
        do j = 1, K_LARGE
            if (abs(h_C_serial(i, j) - h_C_concurrent(i, j)) > max_diff) then
                max_diff = abs(h_C_serial(i, j) - h_C_concurrent(i, j))
            end if
        end do
    end do
    
    print '(a, e10.2)', '  Максимальное отклонение (serial vs concurrent): ', max_diff
    
    if (max_diff < 1e-8) then
        print *, '  ✓ РЕЗУЛЬТАТЫ СОВПАДАЮТ!'
    else
        print *, '  ✗ ОШИБКА! Результаты не совпадают!'
    end if
    print *, ''
    
    ! ====================================================================
    ! Вывод итоговых результатов
    ! ====================================================================
    print *, '========================================================================'
    print *, 'РЕЗУЛЬТАТЫ'
    print *, '========================================================================'
    print '(a, f10.3, a, f12.2, a)', ' SERIAL:     ', time_serial * 1000, ' ms, ', gflops_serial, ' GFLOPS'
    print '(a, f10.3, a, f12.2, a)', ' CONCURRENT: ', time_concurrent * 1000, ' ms, ', gflops_concurrent, ' GFLOPS'
    print *, ''
    print '(a, f5.2, a)', ' Ускорение (serial/concurrent): ', time_serial / time_concurrent, 'x'
    print '(a, f5.2, a)', ' Эффективность (speedup/NSTREAMS): ', (time_serial / time_concurrent) / NSTREAMS, 'x'
    print *, ''
    
    if (time_serial / time_concurrent > 3.5) then
        print *, ' ✓ Отличный scalability! Почти линейное ускорение.'
    else if (time_serial / time_concurrent > 2.5) then
        print *, ' ✓ Хороший scalability.'
    else
        print *, ' ~ Умеренный scalability.'
    end if
    
    print *, ''
    print *, '========================================================================'
    print *, 'КЛЮЧЕВЫЕ УРОКИ:'
    print *, '  1. Асинхронный H2D/GEMM/D2H с разными streams перекрываются'
    print *, '  2. Отдельный cublas handle на каждую часть с привязкой к своему stream'
    print *, '  3. Одна синхронизация в конце вместо NSTREAMS'
    print *, '  4. Буферы для каждой части не должны пересекаться'
    print *, '========================================================================'
    
    deallocate(h_A, h_B, h_C_serial, h_C_concurrent)
    deallocate(d_A, d_B_all, d_C_all)
    
end program gemm_streams_main
