 module make_system
    use mke_elements
    
    double precision, parameter :: ul = 2.0, ur = 0.0
    double precision, allocatable :: p(:), K(:,:)
    integer Nx
    
    contains
    
    subroutine MakeSys(l_gr, rho)
        double precision l_gr, rho
        
        Nx = size(uzly)-1
        call IntDphiDx()             ! Расчет матрицы системы
        call IntRhoPhi(l_gr, rho)    ! Расчет столбца свободных членов
        call GranUsl()               ! Учет граничных условий
        
        deallocate(hx)
        
    end subroutine MakeSys
    
    subroutine IntDphiDx()
        ! Рассчитывает коэффициенты K_ij — интегралы от произведения dphi_i/dx*dphi_j/dx
        ! Значения коэффициентов записываются в массив K, который представляет собой матрицу системы
        integer ierr, Nx
            
        Nx = size(uzly)-1
        allocate(K((Nx-1),0:Nx), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva № ', ierr
            stop
        endif
        
        K = 0
        !!! Разобраться здесь с нумерацией
        do i = 1, Nx-1
            K(i, i-1) = -1/hx(i)
            K(i, i) = 1/hx(i) + 1/hx(i+1)
            K(i, i+1) = -1/hx(i+1)
        enddo
        
    end subroutine IntDphiDx
    
    
    subroutine IntRhoPhi(l_gr, rho)   ! Координата границы, плотность заряда (деленная на \epsilon_0)
        ! Рассчитывает p_i — интегралы от произведения плотности на базисные функции
        ! Значения интегралов записываются в массив p — столбец свободных членов
        integer ierr, n_gr, Nx
        double precision l_gr, rho
        
        ! Проверка того, что посчитан массив с узлами
        if (.not. allocated(uzly)) then
            print*, 'massiv s koordinatami uzlov ne opredelen'
            stop
        endif
        
        !! Размещаем массив для p
        Nx = size(uzly)-1   ! Чтобы не мучиться с передачей числа узлов через параметры процедуры
        allocate(p(Nx-1), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva № ', ierr
            stop
        endif
        
        !! Присваиваем p значения.
        ! Код только для постоянного p, для сложной зависимости нужно переписывать
        
        n_gr = nint(l_gr*Nx/uzly(Nx)) ! Возможно, лучше все же передавать в явном виде Nx и l — это решать по мере развития программы
        
        p = 0
        do i = 1, (n_gr-1)
            p(i) = 0.5*rho*(hx(i)+hx(i+1))
        enddo
        p(n_gr) = 0.5*rho*hx(n_gr)
        
    end subroutine IntRhoPhi
    
    subroutine GranUsl()
        !### Учет граничных условий ###
        p(1) = p(1) - K(1,0)*ul                 ! меняем свободные члены в первом и последнем уравнении
        p(Nx-1) = p(Nx-1) - K(Nx-1,Nx)*ur
        
        ! В матрице K остались лишние столбцы — первый и последний.
        !       Можно делать матрицу меньшего размера, но тогда будет сложнее код: 
        !       придется отдельно считать \alpha_1 и \gamma_{N-1} и отдельно просчитывать первую и последнюю строку 
        !       
        !       Можно просто учитывать нумерацию в SolveSyS
        !       ПОКА ЧТО ТАК И ДЕЛАЕМ
        !
    end subroutine GranUsl
    
end module make_system