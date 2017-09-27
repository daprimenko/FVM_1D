module solve_system
    use make_system
    use mke_elements
    use ifport
    
    double precision, allocatable :: a(:)
    double precision, allocatable :: u(:,:)
    contains
    
    subroutine SolveSys()
        ! Решает систему линейных уравнений по найденным раньше матрице K и столбцу свободных членов p
    
        ! Пока что будет прогонка, потом использовать разложение Холецкого или что-то в этом роде.
        ! Думать, как в общем случае использовать симметрию матрицы
        
        double precision, allocatable :: ap(:), bp(:)   !Массивы прогоночных коэффициентов
        integer i
        
        allocate(ap(Nx-1), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva <<ap>> № ', ierr
            stop
        endif
        
        allocate(bp(Nx-1), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva <<bp>> № ', ierr
            stop
        endif
        
        allocate(a(0:Nx), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva <<a>> № ', ierr
            stop
        endif
        
        !Прямой ход прогонки
        
        !Первые прогоночные коэффициенты
        
        ! Граничные условия первого рода, поэтому пока что нумерацию начинаем с единицы.
        !### Это все будет меняться ###
        ap(1) = -K(1,2)/K(1,1)
        bp(1) = p(1)/K(1,1)
                      
        do i = 2, Nx-2
        
            !Прогоночные коэффициенты
            ap(i) = -K(i,i+1) / ( K(i-1,i)*ap(i-1) + K(i,i) )
            bp(i) = ( p(i) - K(i-1,i)*bp(i-1) ) / ( K(i-1,i)*ap(i-1) + K(i,i) )
            
        end do 
        
        !Обратный ход прогонки
        
        a(Nx-1) = ( p(Nx-1) - K(Nx-1,Nx-2)*bp(Nx-2) ) / ( K(Nx-1,Nx-2)*ap(Nx-2) + K(Nx-1,Nx-1) )
                
        do i = 2, Nx-1
            m = Nx-i
            a(m)=ap(m)*a(m+1)+bp(m)
        end do
        
        !### Вытаскиваем крайние значения из граничных условий ###
        a(0) = ul
        a(Nx) = ur
        
        call KoeffToFunc()
        call VisFunc()
        
    end subroutine SolveSys
    
    subroutine KoeffToFunc()
        ! Суммирует базисные функции с найденными коэффициентами a_i. В результате получаем искомую функцию u.
        
        integer i
        ! Сколько точек брать для функции?
        ! Так как интерполяция линейная, смысла брать больше, чем Nx, нет — 
        ! графопостроитель сделает такую же линейную интерполяцию
        
        allocate(u(0:Nx, 2), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva <<u>> № ', ierr
            stop
        endif
        
        u(0:Nx,1) = uzly
        u(0:Nx,2) = a
        
        deallocate(a, uzly)
                
    end subroutine KoeffToFunc
    
    subroutine VisFunc()
        ! Визуализирует рассчитанyю функцию u.
        ! Делать через Python
        !open(2, file = 'inf.dat', form = 'binary')
        !write(2) Nx
        !close(2)
        logical(4) res
        
        open(2, file = 'u.dat', form = 'binary')
        write(2) u(0:Nx,2)
        close(2)
        deallocate(u)
        
        res = systemqq('python visualize_function_1D.py')
        
    end subroutine VisFunc
    
end module solve_system    