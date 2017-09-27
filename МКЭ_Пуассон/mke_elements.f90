! ���������, � ����� ���������� ������ ��������. ����� ����, �������� ���� try-except-��, � ����� ����, ������ ������ ��� ������

module mke_elements
    use ifport
    
    double precision, allocatable :: uzly(:)    ! ������ � ������������ ������� �����, ��������� �� ��������� ������ � ����������� �������� �������
    double precision, allocatable :: hx(:)      ! ����� �������� ���������
    double precision, parameter :: MaxAdjacentRatio = 2.0
    
    contains
    
    subroutine Razbienie(Nx, l)  ! Nx � ����� ��������, l � ������� �������
    
        ! ��������� ��������� ������� �� ��������
        ! ���������� ����� ������������ � ������ uzly
        ! ������� ������ ������� ���������
        integer Nx, ierr
        double precision l
        
        !���������� ������� � ������������ �����
        allocate(uzly(0:Nx), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva � ', ierr
            stop
        endif
        
        !�������� ������, ���������� ���������� ������� �����
        h = l/Nx
        allocate(hx(Nx), stat = ierr)
        if (ierr .ne. 0) then
            print*, 'Oshibka razmescheniya massiva � ', ierr
            stop
        endif
        hx = h      !### ����� ����� ������ �� ������������ ��� ###
        hx(1) = 0.5*h
        hx(2) = 1.5*h
        uzly(0) = 0
        do i = 1,Nx
            uzly(i) = uzly(i-1) + hx(i)
        enddo
        
        call VisUzly()
        
        !������ �� ������� � �� ��� ����������
    end subroutine Razbienie
    
    
    subroutine VisUzly()
        ! �������� python � ������� 
        ! use msflib � runqq ��� ����� systemqq 
        ! ���� ��� �����-�� ������ system, �� ����� �� �� Linux
        ! ��� ��� use kernel32 + �����-�� ������� WinAPI
    
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
        double precision x0, xN, h0, hN !������ � ����� ���������, ��������� � �������� ���
        
        integer N                   ! ���������� �������� ���������        
        double precision h, q, l    ! ���, ����������� ��������� ���� � ����� ���������
        double precision, allocatable :: x(:)
        
        ! ��� ����� ������ ���������� ������������ ���� �������� q, N. ��� ���� N ����� ���������� �������.
        ! ������, ��� � ���� ������
        ! ���� ��� ��������� N � ������ ��� �������, ����� ����������. ��������� ������� � ��� ���������.
        ! ����� ���������� ���������, ����� ������ ���������� ���������� ������� ������
        
        ! ��� �������: �� ������ ������ �������, ������ ����� � ����������� ��������� ���� � ���������� ��������
        ! ����� ������ ������� ����� ����������, � � ����� �� ����� ���� ������� ������� ��������� ���������� �������
        
        ! ��� ����� ����������� �������: ��������� N � ����������� q � ����� ��������� N, � ��� ���� N �� ���������� ��������
        
        l = xN - x0
        q = (l - h0) / (l - hN)        
        if(q > 2) then 
            print*, '������� ������ ��������� ����!'
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
            print*, '������� ������ ��������� ���� �� �������!'
        endif
        
            
    end subroutine SegmentPartition
        
end module mke_elements