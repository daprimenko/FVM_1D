use mke_elements
use make_system
use solve_system
    ! ���-�� ������ �� �������������� �������� � ���� �� ���������� �������� � ���������� �������
    ! ������ ��� ����������
    double precision :: tolschina = 3.0, rho = 10.0
    call Partition()          ! ��������� ������� �� �������� �������� Nx,l
    !call MakeSys(0.99*tolschina, rho)      ! ������ ������� ������� � ������� ��������� ������ l_gr, rho
    !call SolveSys()                         ! ������� ���������� ����

end