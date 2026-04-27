./rowl_plus_rowl_cpu
git push -f origin main
git commit -m "Сверхпервый пуш"
git add .
git commit -m "Сверхпервый пуш"
git push -f origin main
mc
mc
gfortran columnk_plus_columnk_cpu.f90 -O3 -o columnk_plus_columnk_cpu
./columnk_plus_columnk_cpu
git add .
git commit -m "Сложение столбцов высоты k"
git push -f origin main
mc
gfortran matrixkl_plus_matrixkl_cpu.f90 -O3 -o matrixkl_plus_matrixkl_cpu_cpu
./matrixkl_plus_matrixkl_cpu_cpu
mc
gfortran matrixkl_plus_matrixkl_cpu.f90 -O3 -o matrixkl_plus_matrixkl_cpu
./matrixkl_plus_matrixkl_cpu
git add .
git commit -m "Сложение матриц k на l на cpu"
git push -f origin main
mc
mc
rowl_plus_rowl_gpu.f90 -O3 -o rowl_plus_rowl_gpu
rowl_plus_rowl_gpu.cuf -O3 -o rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
mc
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu.exe

~
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
git add .
git commit -m "Сложение строк длины l на gpu"
git push -f origin main
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 columnk_plus_columnk_gpu.cuf -o columnk_plus_columnk_gpu
./columnk_plus_columnk_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
./rowl_plus_rowl_cpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu_2.cuf -o rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu_2.cuf -o rowl_plus_rowl_gpu_2
mc
mc
nvfortran -O3 -o rowl_plus_rowl rowl_plus_rowl.cuf
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 vec_add_module.cuf -o rowl_plus_rowl_gpu
mc
gfortran makefile.f90 -O3 -o 1l+1l_cpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 vec_add_module.cuf -o rowl_plus_rowl_gpu
nvfortran vec_add.cuf -o vec_add
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 vec_add_module.cuf -o rowl_plus_rowl_gpu
mc
nvfortran -O3 -cuda -cudalib=cublas -gpu=cc80 rowl_plus_rowl_gpu.cuf -o rowl_plus_rowl_gpu
./rowl_plus_rowl_gpu
make
./vec_mat_add
ьс
mc
nvfortran -O3 -cuda -gpu=cc80 rowl_plus_rowl_gpu_2.cuf -o rowl_plus_rowl_gpu_2
./rowl_plus_rowl_gpu_2
nvfortran -O3 -cuda -gpu=cc80 rowl_plus_rowl_gpu_2.cuf -o rowl_plus_rowl_gpu_2
./rowl_plus_rowl_gpu_2
mc
# Компиляция
nvfortran -c utils_module.f90
nvfortran -c vec_add_cpu.f90
nvfortran -c vec_add_module.cuf
nvfortran -o vec_add main.f90 utils_module.o vec_add_cpu.o vec_add_module.o
# Запуск
mc
nvfortran vec_add.f90 -o vec_add
mc
cd ~/task1 && gfortran test_1.f90 -o test_1 && ./test_1 && gfortran test_2.f90 -o test_2 && ./test_2 && nvfortran -cuda cuda_matrix_add_random.f90 -o cuda_matrix_add_random && ./cuda_matrix_add_random && nvfortran -cuda matrix_add_interactive.f90 -o matrix_add_interactive && ./matrix_add_interactive
cd ~/1 -o test_1 && ./test_1 && gfortran test_2.f90 -o test_2 && ./test_2 && nvfortran -cuda cuda_matrix_add_random.f90 -o cuda_matrix_add_random && ./cuda_matrix_add_random && nvfortran -cuda matrix_add_interactive.f90 -o matrix_add_interactive && ./matrix_add_interactive
mc
nvfortran -cuda -Mfixed -o matrix_add matrix_add.f90 && ./matrix_add
nvfortran -cuda -o matrix_add matrix_add.f90
./matrix_add.f90
./matrix_add
mc
! ================================================================
! ЗАДАНИЕ 1: СЛОЖЕНИЕ МАТРИЦ ПРОИЗВОЛЬНОГО РАЗМЕРА (CPU + GPU)
! ================================================================
module matrix_kernel_module
    use cudafor
    implicit none
    
contains
    attributes(global) subroutine add_matrix_kernel(a, b, c, n, m)
        real, device, intent(in)  :: a(:,:)
        real, device, intent(in)  :: b(:,:)
        real, device, intent(out) :: c(:,:)
        integer, value, intent(in) :: n, m
        
        integer :: i, j
        
        i = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        j = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        
        if (i <= n .and. j <= m) then             c(i, j) = a(i, j) + b(i, j)
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
    
    if (max_diff < 1e-5) then         print *, '✓ РЕЗУЛЬТАТЫ СОВПАДАЮТ!';     else         print *, '✗ ОШИБКА! Результаты не совпадают!';     end if         ! Вывод первых 3x3 элементов;     print *, '';     print *, 'Пример (первые 3x3):';     print *, '   i   j   A(i,j)   +   B(i,j)   =   C(i,j) (GPU)'         do i = 1, min(3,n)
        do j = 1, min(3,m)
            print '(i4,i4,3f10.4)', i, j, A(i,j), B(i,j), C_gpu(i,j)
        end do
    end do
    
    deallocate(A, B, C_cpu, C_gpu, A_d, B_d, C_d)
    
    print *, ''
    print *, '============================================='
    print *, 'ПРОГРАММА ЗАВЕРШЕНА!'
    print *, '============================================='
end program matrix_add_interactivewqw
! ЗАДАНИЕ 1: СЛОЖЕНИЕ МАТРИЦ ПРОИЗВОЛЬНОГО РАЗМЕРА (CPU + GPU)
! ================================================================
module matrix_kernel_moduleo
mc
cat > matrix_add.f90 << 'EOF' && nvfortran -cuda -o matrix_add matrix_add.f90 && echo -e "3\n3" | ./matrix_add
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
    do i = 1, n
        do j = 1, m
            C_cpu(i,j) = A(i,j) + B(i,j)
        end do
    end do
    A_d = A
    B_d = B
    call add_matrix_kernel<<<dim3((m+15)/16, (n+15)/16, 1), dim3(16,16,1)>>>(A_d, B_d, C_d, n, m)
    istat = cudaDeviceSynchronize()
    C_gpu = C_d
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
    print *, ''
    print *, 'Пример (первые 3x3):'
    print *, '   i   j   A(i,j)   +   B(i,j)   =   C(i,j) (GPU)'
    do i = 1, min(3,n)
        do j = 1, min(3,m)
            print '(i4,i4,3f10.4)', i, j, A(i,j), B(i,j), C_gpu(i,j)
        end do
    end do
    deallocate(A, B, C_cpu, C_gpu, A_d, B_d, C_d)
    print *, ''
    print *, '============================================='
    print *, 'ПРОГРАММА ЗАВЕРШЕНА!'
    print *, '============================================='
end program matrix_add_interactive
EOF

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
    do i = 1, n
        do j = 1, m
            C_cpu(i,j) = A(i,j) + B(i,j)
        end do
    end do
    A_d = A
    B_d = B
    call add_matrix_kernel<<<dim3((m+15)/16, (n+15)/16, 1), dim3(16,16,1)>>>(A_d, B_d, C_d, n, m)
    istat = cudaDeviceSynchronize()
    C_gpu = C_d
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
    print *, ''
    print *, 'Пример (первые 3x3):'
    print *, '   i   j   A(i,j)   +   B(i,j)   =   C(i,j) (GPU)'
    do i = 1, min(3,n)
        do j = 1, min(3,m)
            print '(i4,i4,3f10.4)', i, j, A(i,j), B(i,j), C_gpu(i,j)
        end do
    end do
    deallocate(A, B, C_cpu, C_gpu, A_d, B_d, C_d)
    print *, ''
    print *, '============================================='
    print *, 'ПРОГРАММА ЗАВЕРШЕНА!'
    print *, '============================================='
end program matrix_add_interactive
EOF

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

mc
