use mke_elements
use make_system
use solve_system
    ! Как-то хорошо бы распараллелить процессы — хотя бы построение картинок и собственно расчеты
    ! Думать над названиями
    double precision :: tolschina = 3.0, rho = 10.0
    call Razbienie(100, tolschina)          ! Разбиение области на конечные элементы Nx,l
    call MakeSys(0.99*tolschina, rho)      ! Расчет матрицы системы и столбца свободных членов l_gr, rho
    call SolveSys()                         ! Решение полученной СЛАУ

end