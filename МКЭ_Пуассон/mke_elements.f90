! Непонятно, в каком количестве делать проверки. Может быть, напихать везд try-except-ов, а может быть, убрать вообще все защиты
! Придумать оптимальный формат описания сетки — править GetPartitionInfo
    
module mke_elements
    use ifport
    
    double precision, allocatable :: uzly(:)    ! Массив с координатами узловых точек, интегралы от плотности заряда и производной базисной функции
    double precision, allocatable :: Points(:)  ! Массив с координатами узловых точек, интегралы от плотности заряда и производной базисной функции
    double precision, allocatable :: hx(:)      ! Длины конечных элементов
    double precision, parameter :: MaxAdjacentRatio = 2.0
    integer NumberOfFragment
    character(len=*), parameter :: PartitionFileName = "partition_info.txt"
    
    contains
    
    subroutine Partition()
    
        ! Выполняет разбиение области на элементы
     
        double precision, allocatable :: CurrentFragmentPoints(:), info(:,:)
        integer i
        
        call GetPartitionInfo(PartitionFileName, NumberOfFragments, info)
        
        do i = 1, NumberOfFragments
            call FragmentPartition(CurrentFragmentPoints, info(:,i))
            call AddArrayToArray(CurrentFragmentPoints, Points)
            deallocate(CurrentFragmentPoints)
        enddo
        
        call VisUzly()
        
        !Массив не удаляем — он еще пригодится
    end subroutine Partition
     
     
    subroutine VisUzly()
        ! Вызывать python с помощью 
        ! use msflib — runqq или лучше systemqq 
        ! есть еще какой-то просто system, но вроде он на Linux
        ! или еще use kernel32 + какая-то команда WinAPI
    
        logical(4) res
        open(2, file='number_of_points.dat', form = 'binary')
        write(2) size(Points)
        close(2)
        
        open(2, file = 'partition.dat', form = 'binary')
        write(2) Points
        close(2)
        
        res = systemqq('python visualize_points_1D.py')
        
    end subroutine VisUzly
    
    
    subroutine GetPartitionInfo(fname, N, inf)
        character(len=*) fname
        integer N, ios
        double precision, allocatable :: inf(:,:)
        double precision curinfo(4)
        
        open(2, file=fname)
        !read(2,*, IOSTAT=ios) N
        !allocate(inf(4, N))
        !do i = 1, N
        !    read(2,*, IOSTAT=ios) inf(:, i)
        !enddo
        
        ios = 0
        N = 0
        do while(ios == 0)
            read(2, *, IOSTAT=ios) curinfo
            call AddColumnToMatrix(curinfo, inf)
            N = N+1
        enddo
        N = N - 1
        
        close(2)
    end subroutine GetPartitionInfo
    
    
    subroutine FragmentPartition(x, parameters)
        double precision x0, xN, h0, hN, parameters(4)        !Начало и конец фрагмента, начальный и конечный шаг
        double precision, allocatable :: x(:)   ! Временный массив, в который будут писаться координаты
        
        integer N                   ! Количество отрезков разбиения        
        double precision h, q, l    ! Шаг, коэффициент изменения шага и длина фрагмента
        
        ! Для таких данных существует единственная пара значений q, N. При этом N может получиться дробным.
        ! Думать, что с этим делать
        ! Пока что округляем N и делаем все отрезки, кроме последнего. Последний отрезок — как получится.
        ! Чтобы получалось адекватно, нужно давать достаточно адекватные входные данные
        
        ! Еще вариант: не давать правой границы, вместо этого — коэффициент изменения шага и количество отрезков
        ! Тогда правая граница будет некрасивая, и в итоге на самом краю счетной области получится непонятный отрезок
        
        ! Еще можно итеративный процесс: округлили N — пересчитали q — снова округлили N, и так пока N не перестанет меняться
        
        x0 = parameters(1)
        xN = parameters(2)
        h0 = parameters(3)
        hN = parameters(4)
        
        !! Сделать отдельный блок для постоянного шага
        l = xN - x0
        q = (l - h0) / (l - hN)        
        if(q > MaxAdjacentRatio) then 
            print*, 'Слишком резкое изменение шага!'
        endif
        
        if (q == 1) then
            N = 1 + l/h0
        else
            N = 2 + idnint( log(hN/h0) / log(q) )
        endif
        
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
    end subroutine FragmentPartition
    
    
    subroutine AddArrayToArray(subarray, array)
        double precision, allocatable :: subarray(:), array(:)
        double precision, allocatable :: temparray(:)
        integer oldsize, newsize
        
        if(.not.allocated(subarray)) then
            print*, 'Добавляемый массив не размещен!'
            stop
        endif

        if(allocated(array)) then
            oldsize = size(array)
            newsize = oldsize + size(subarray)
            allocate(temparray(newsize))
            temparray(1:oldsize) = array
            temparray(oldsize+1:newsize) = subarray

            deallocate(array)
            call move_alloc(temparray, array)
        else
            allocate(array(size(subarray)))
            array = subarray
            !!! сделать просто call move_alloc(temparray, array)?
        end if
    end subroutine AddArrayToArray
    
    
    subroutine AddColumnToMatrix(column, matrix)
        double precision, allocatable :: matrix(:,:), newmatrix(:,:)
        double precision column(:)
        integer m, n
        
        m = size(column)       
        
        if(.not. allocated(matrix)) then
            allocate(matrix(m, 1))
            matrix(:,1) = column
        else
            if( size(matrix, 1) /= m) then
                print*, "Добавляемый столбец не согласован с матрицей по размеру"
                stop
            endif
            
            n = size(matrix, 2)
            allocate(newmatrix(m, n+1))
            newmatrix(:, 1:n) = matrix
            newmatrix(:, n+1) = column
            deallocate(matrix)
            call move_alloc(newmatrix, matrix)            
        endif
        
    end subroutine AddColumnToMatrix

    
end module mke_elements