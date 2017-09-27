! Непонятно, в каком количестве делать проверки. Может быть, напихать везд try-except-ов, а может быть, убрать вообще все защиты

module mke_elements
    use ifport
    
    double precision, allocatable :: uzly(:)    ! Массив с координатами узловых точек, интегралы от плотности заряда и производной базисной функции
    double precision, allocatable :: hx(:)      ! Длины конечных элементов
    double precision, parameter :: MaxAdjacentRatio = 2.0
    
    contains
    
    subroutine Razbienie(Nx, l)  ! Nx — число отрезков, l — толщина области
    
        ! Выполняет разбиение области на элементы
        ! Координаты узлов записываются в массив uzly
        ! Сделать запись номеров элементов
        integer Nx, ierr
        double precision l
        
        !Размещение массива с координатами узлов
        allocate(uzly(0:Nx), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva № ', ierr
            stop
        endif
        
        !Строится массив, содержащий координаты узловых точек
        h = l/Nx
        allocate(hx(Nx), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva № ', ierr
            stop
        endif
        hx = h      !### Здесь можно менять на произвольный шаг ###
        hx(1) = 0.5*h
        hx(2) = 1.5*h
        uzly(0) = 0
        do i = 1,Nx
            uzly(i) = uzly(i-1) + hx(i)
        enddo
        
        call VisUzly()
        
        !Массив не удаляем — он еще пригодится
    end subroutine Razbienie
    
    
    subroutine VisUzly()
        ! Вызывать python с помощью 
        ! use msflib — runqq или лучше systemqq 
        ! есть еще какой-то просто system, но вроде он на Linux
        ! или еще use kernel32 + какая-то команда WinAPI
    
        logical(4) res
        open(2, file='partition_info.dat', form = 'binary')
        write(2) size(uzly)
        close(2)
        
        open(2, file = 'partition.dat', form = 'binary')
        write(2) uzly
        close(2)
        
        !res = systemqq('python visualize_points_1D.py')
        
    end subroutine VisUzly
    
    subroutine SegmentPartition(x0, xN, h0, hN)
        double precision x0, xN, h0, hN !Начало и конец фрагмента, начальный и конечный шаг
        
        integer N                   ! Количество отрезков разбиения        
        double precision h, q, l    ! Шаг, коэффициент изменения шага и длина фрагмента
        double precision, allocatable :: x(:)
        
        ! Для таких данных существует единственная пара значений q, N. При этом N может получиться дробным.
        ! Думать, что с этим делать
        ! Пока что округляем N и делаем все отрезки, кроме последнего. Последний отрезок — как получится.
        ! Чтобы получалось адекватно, нужно давать достаточно адекватные входные данные
        
        ! Еще вариант: не давать правой границы, вместо этого — коэффициент изменения шага и количество отрезков
        ! Тогда правая граница будет некрасивая, и в итоге на самом краю счетной области получится непонятный отрезок
        
        ! Еще можно итеративный процесс: округлили N — пересчитали q — снова округлили N, и так пока N не перестанет меняться
        
        l = xN - x0
        q = (l - h0) / (l - hN)        
        if(q > 2) then 
            print*, 'Слишком резкое изменение шага!'
        endif
        N = 1 + idnint( log(hN/h0) / log(q) )
        
        allocate(x(N))
        x(1) = x0
        h = h0
        do i = 2, N-1
            x(i) = x(i-1) + h
            h = h*q
        enddo
        x(N) = xN
        
        if((abs(x(N) - x(N-1)) > (2*hN)) .or. (abs(x(N) - x(N-1)) < 0.5*hN)) then
            print*, 'Слишком резкое изменение шага на границе!'
        endif
        
            
    end subroutine SegmentPartition
        
end module mke_elements