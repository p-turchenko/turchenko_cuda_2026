module vec_add_cpu_module
    implicit none
    private
    public :: vec_add_cpu
contains

    ! Эталонная CPU-реализация: простое сложение векторов c = a + b.
    subroutine vec_add_cpu(a, b, c)
        real, intent(in)  :: a(:), b(:)
        real, intent(out) :: c(:)
        c = a + b
    end subroutine vec_add_cpu

end module vec_add_cpu_module
