 module make_system
    use mke_elements
    
    double precision, parameter :: ul = 2.0, ur = 0.0
    double precision, allocatable :: p(:), K(:,:)
    integer Nx
    
    contains
    
    subroutine MakeSys(l_gr, rho)
        double precision l_gr, rho
        
        Nx = size(uzly)-1
        call IntDphiDx()             ! ������ ������� �������
        call IntRhoPhi(l_gr, rho)    ! ������ ������� ��������� ������
        call GranUsl()               ! ���� ��������� �������
        
        deallocate(hx)
        
    end subroutine MakeSys
    
    subroutine IntDphiDx()
        ! ������������ ������������ K_ij � ��������� �� ������������ dphi_i/dx*dphi_j/dx
        ! �������� ������������� ������������ � ������ K, ������� ������������ ����� ������� �������
        integer ierr, Nx
            
        Nx = size(uzly)-1
        allocate(K((Nx-1),0:Nx), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva � ', ierr
            stop
        endif
        
        K = 0
        !!! ����������� ����� � ����������
        do i = 1, Nx-1
            K(i, i-1) = -1/hx(i)
            K(i, i) = 1/hx(i) + 1/hx(i+1)
            K(i, i+1) = -1/hx(i+1)
        enddo
        
    end subroutine IntDphiDx
    
    
    subroutine IntRhoPhi(l_gr, rho)   ! ���������� �������, ��������� ������ (�������� �� \epsilon_0)
        ! ������������ p_i � ��������� �� ������������ ��������� �� �������� �������
        ! �������� ���������� ������������ � ������ p � ������� ��������� ������
        integer ierr, n_gr, Nx
        double precision l_gr, rho
        
        ! �������� ����, ��� �������� ������ � ������
        if (.not. allocated(uzly)) then
            print*, 'massiv s koordinatami uzlov ne opredelen'
            stop
        endif
        
        !! ��������� ������ ��� p
        Nx = size(uzly)-1   ! ����� �� �������� � ��������� ����� ����� ����� ��������� ���������
        allocate(p(Nx-1), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva � ', ierr
            stop
        endif
        
        !! ����������� p ��������.
        ! ��� ������ ��� ����������� p, ��� ������� ����������� ����� ������������
        
        n_gr = nint(l_gr*Nx/uzly(Nx)) ! ��������, ����� ��� �� ���������� � ����� ���� Nx � l � ��� ������ �� ���� �������� ���������
        
        p = 0
        do i = 1, (n_gr-1)
            p(i) = 0.5*rho*(hx(i)+hx(i+1))
        enddo
        p(n_gr) = 0.5*rho*hx(n_gr)
        
    end subroutine IntRhoPhi
    
    subroutine GranUsl()
        !### ���� ��������� ������� ###
        p(1) = p(1) - K(1,0)*ul                 ! ������ ��������� ����� � ������ � ��������� ���������
        p(Nx-1) = p(Nx-1) - K(Nx-1,Nx)*ur
        
        ! � ������� K �������� ������ ������� � ������ � ���������.
        !       ����� ������ ������� �������� �������, �� ����� ����� ������� ���: 
        !       �������� �������� ������� \alpha_1 � \gamma_{N-1} � �������� ������������ ������ � ��������� ������ 
        !       
        !       ����� ������ ��������� ��������� � SolveSyS
        !       ���� ��� ��� � ������
        !
    end subroutine GranUsl
    
end module make_system