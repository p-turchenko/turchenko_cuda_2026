module mat_add_cpu_module
    implicit none
    private
    public :: mat_add_cpu
contains

    ! Эталонная CPU-реализация: сложение матриц c = a + b.
    subroutine mat_add_cpu(a, b, c)
        real, intent(in)  :: a(:,:), b(:,:)
        real, intent(out) :: c(:,:)
        c = a + b
    end subroutine mat_add_cpu

end module mat_add_cpu_module
