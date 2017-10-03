! ���������, � ����� ���������� ������ ��������. ����� ����, �������� ���� try-except-��, � ����� ����, ������ ������ ��� ������
! ��������� ����������� ������ �������� ����� � ������� GetPartitionInfo
    
module mke_elements
    use ifport
    
    double precision, allocatable :: uzly(:)    ! ������ � ������������ ������� �����, ��������� �� ��������� ������ � ����������� �������� �������
    double precision, allocatable :: Points(:)  ! ������ � ������������ ������� �����, ��������� �� ��������� ������ � ����������� �������� �������
    double precision, allocatable :: hx(:)      ! ����� �������� ���������
    double precision, parameter :: MaxAdjacentRatio = 2.0
    integer NumberOfFragment
    character(len=*), parameter :: PartitionFileName = "partition_info.txt"
    
    contains
    
    subroutine Partition()
    
        ! ��������� ��������� ������� �� ��������
     
        double precision, allocatable :: CurrentFragmentPoints(:), info(:,:)
        integer i
        
        call GetPartitionInfo(PartitionFileName, NumberOfFragments, info)
        
        do i = 1, NumberOfFragments
            call FragmentPartition(CurrentFragmentPoints, info(:,i))
            call AddArrayToArray(CurrentFragmentPoints, Points)
            deallocate(CurrentFragmentPoints)
        enddo
        
        call VisUzly()
        
        !������ �� ������� � �� ��� ����������
    end subroutine Partition
     
     
    subroutine VisUzly()
        ! �������� python � ������� 
        ! use msflib � runqq ��� ����� systemqq 
        ! ���� ��� �����-�� ������ system, �� ����� �� �� Linux
        ! ��� ��� use kernel32 + �����-�� ������� WinAPI
    
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
        double precision x0, xN, h0, hN, parameters(4)        !������ � ����� ���������, ��������� � �������� ���
        double precision, allocatable :: x(:)   ! ��������� ������, � ������� ����� �������� ����������
        
        integer N                   ! ���������� �������� ���������        
        double precision h, q, l    ! ���, ����������� ��������� ���� � ����� ���������
        
        ! ��� ����� ������ ���������� ������������ ���� �������� q, N. ��� ���� N ����� ���������� �������.
        ! ������, ��� � ���� ������
        ! ���� ��� ��������� N � ������ ��� �������, ����� ����������. ��������� ������� � ��� ���������.
        ! ����� ���������� ���������, ����� ������ ���������� ���������� ������� ������
        
        ! ��� �������: �� ������ ������ �������, ������ ����� � ����������� ��������� ���� � ���������� ��������
        ! ����� ������ ������� ����� ����������, � � ����� �� ����� ���� ������� ������� ��������� ���������� �������
        
        ! ��� ����� ����������� �������: ��������� N � ����������� q � ����� ��������� N, � ��� ���� N �� ���������� ��������
        
        x0 = parameters(1)
        xN = parameters(2)
        h0 = parameters(3)
        hN = parameters(4)
        
        !! ������� ��������� ���� ��� ����������� ����
        l = xN - x0
        q = (l - h0) / (l - hN)        
        if(q > MaxAdjacentRatio) then 
            print*, '������� ������ ��������� ����!'
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
            print*, '������� ������ ��������� ���� �� �������!'
        endif
    end subroutine FragmentPartition
    
    
    subroutine AddArrayToArray(subarray, array)
        double precision, allocatable :: subarray(:), array(:)
        double precision, allocatable :: temparray(:)
        integer oldsize, newsize
        
        if(.not.allocated(subarray)) then
            print*, '����������� ������ �� ��������!'
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
            !!! ������� ������ call move_alloc(temparray, array)?
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
                print*, "����������� ������� �� ���������� � �������� �� �������"
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