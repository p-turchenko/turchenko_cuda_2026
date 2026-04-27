module utils_module
    implicit none
    private
    public :: rnd_fill, check_equal_vec, check_equal_mat

    interface rnd_fill
        module procedure rnd_fill_1d, rnd_fill_2d
    end interface rnd_fill

contains
    ! Заполнение вектора случайными значениями в диапазоне [0, 100].
    subroutine rnd_fill_1d(a)
        real, intent(out) :: a(:)
        call random_number(a)
        a = 100.0 * a
    end subroutine rnd_fill_1d

    ! Заполнение матрицы случайными значениями в диапазоне [0, 100].
    subroutine rnd_fill_2d(a)
        real, intent(out) :: a(:,:)
        call random_number(a)
        a = 100.0 * a
    end subroutine rnd_fill_2d

    ! Поэлементное сравнение векторов с допуском на погрешность вычислений.
    logical function check_equal_vec(a, b) result(equal)
        real, intent(in) :: a(:), b(:)
        equal = all(abs(a - b) <= 1.0e-4)
    end function check_equal_vec

    ! Поэлементное сравнение матриц с допуском на погрешность вычислений.
    logical function check_equal_mat(a, b) result(equal)
        real, intent(in) :: a(:,:), b(:,:)
        equal = all(abs(a - b) <= 1.0e-4)
    end function check_equal_mat

end module utils_module
