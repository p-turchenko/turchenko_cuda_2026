program main

    use cudafor
    use utils_module
    use vec_add_module      ! GPU-версия: vec_add_gpu
    use vec_add_cpu_module  ! CPU-версия: vec_add_cpu
    use mat_add_module      ! GPU-версия: mat_add_gpu
    use mat_add_cpu_module  ! CPU-версия: mat_add_cpu

    implicit none

    integer, parameter :: N_VEC = 1024*1024
    integer, parameter :: N_MAT = 1024        
    integer, parameter :: THREADS_PER_BLOCK = 256
    integer, parameter :: BLOCK_DIM = 16

    real, allocatable :: va(:), vb(:), vc_gpu(:), vc_cpu(:)
    real, allocatable :: ma(:,:), mb(:,:), mc_gpu(:,:), mc_cpu(:,:)
    integer :: istat

   !==================================================================
    write(*,'(a)') '=== Vector addition ==='
    allocate(va(N_VEC), vb(N_VEC), vc_gpu(N_VEC), vc_cpu(N_VEC), stat=istat)
    if (istat /= 0) then
        write(*,'(a)') 'Memory allocation error in main (vector arrays)'
        stop 1
    end if

    call rnd_fill(va)
    call rnd_fill(vb)

    call vec_add_gpu(va, vb, vc_gpu, THREADS_PER_BLOCK)
    call vec_add_cpu(va, vb, vc_cpu)

    if (check_equal_vec(vc_gpu, vc_cpu)) then
        print *, 'Vectors are equal'
    else
        print *, 'Vectors are not equal'
    end if

    deallocate(va, vb, vc_gpu, vc_cpu)
   !==================================================================

   !==================================================================
    write(*,'(a)') '=== Matrix addition ==='
    allocate(ma(N_MAT,N_MAT), mb(N_MAT,N_MAT), mc_gpu(N_MAT,N_MAT), mc_cpu(N_MAT,N_MAT), stat=istat)
    if (istat /= 0) then
        write(*,'(a)') 'Memory allocation error in main (matrix arrays)'
        stop 1
    end if

    call rnd_fill(ma)
    call rnd_fill(mb)

    call mat_add_gpu(ma, mb, mc_gpu, BLOCK_DIM)
    call mat_add_cpu(ma, mb, mc_cpu)

    if (check_equal_mat(mc_gpu, mc_cpu)) then
        print *, 'Arrays are equal'
    else
        print *, 'Arrays are not equal'
    end if

    deallocate(ma, mb, mc_gpu, mc_cpu)
   !==================================================================
end program main
