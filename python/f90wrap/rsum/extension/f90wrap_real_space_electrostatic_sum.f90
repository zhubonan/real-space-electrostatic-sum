! Module real_space_electrostatic_sum defined in file real_space_electrostatic_sum.f90

subroutine f90wrap_energy(a1, a2, a3, n, rx, ry, rz, z, rc, rd, e, n0, n1, n2, n3)
    use real_space_electrostatic_sum, only: energy
    implicit none
    
    real(8), intent(in), dimension(3) :: a1
    real(8), intent(in), dimension(3) :: a2
    real(8), intent(in), dimension(3) :: a3
    integer, intent(in) :: n
    real(8), intent(in), dimension(n0) :: rx
    real(8), intent(in), dimension(n1) :: ry
    real(8), intent(in), dimension(n2) :: rz
    real(8), intent(in), dimension(n3) :: z
    real(8), intent(in) :: rc
    real(8), intent(in) :: rd
    real(8), intent(out) :: e
    integer :: n0
    !f2py intent(hide), depend(rx) :: n0 = shape(rx,0)
    integer :: n1
    !f2py intent(hide), depend(ry) :: n1 = shape(ry,0)
    integer :: n2
    !f2py intent(hide), depend(rz) :: n2 = shape(rz,0)
    integer :: n3
    !f2py intent(hide), depend(z) :: n3 = shape(z,0)
    call energy(a1=a1, a2=a2, a3=a3, n=n, rx=rx, ry=ry, rz=rz, z=z, rc=rc, rd=rd, e=e)
end subroutine f90wrap_energy

subroutine f90wrap_add_energy_corr_term(zi, qi, rho, rd, corr_term)
    use real_space_electrostatic_sum, only: add_energy_corr_term
    implicit none
    
    real(8), intent(in) :: zi
    real(8), intent(in) :: qi
    real(8), intent(in) :: rho
    real(8), intent(in) :: rd
    real(8), intent(out) :: corr_term
    call add_energy_corr_term(zi=zi, qi=qi, rho=rho, rd=rd, corr_term=corr_term)
end subroutine f90wrap_add_energy_corr_term

subroutine f90wrap_force(a1, a2, a3, n, rx, ry, rz, z, rc, rd, fx, fy, fz, n0, n1, n2, n3, n4, n5, n6)
    use real_space_electrostatic_sum, only: force
    implicit none
    
    real(8), intent(in), dimension(3) :: a1
    real(8), intent(in), dimension(3) :: a2
    real(8), intent(in), dimension(3) :: a3
    integer, intent(in) :: n
    real(8), intent(in), dimension(n0) :: rx
    real(8), intent(in), dimension(n1) :: ry
    real(8), intent(in), dimension(n2) :: rz
    real(8), intent(in), dimension(n3) :: z
    real(8), intent(in) :: rc
    real(8), intent(in) :: rd
    real(8), intent(inout), dimension(n4) :: fx
    real(8), intent(inout), dimension(n5) :: fy
    real(8), intent(inout), dimension(n6) :: fz
    integer :: n0
    !f2py intent(hide), depend(rx) :: n0 = shape(rx,0)
    integer :: n1
    !f2py intent(hide), depend(ry) :: n1 = shape(ry,0)
    integer :: n2
    !f2py intent(hide), depend(rz) :: n2 = shape(rz,0)
    integer :: n3
    !f2py intent(hide), depend(z) :: n3 = shape(z,0)
    integer :: n4
    !f2py intent(hide), depend(fx) :: n4 = shape(fx,0)
    integer :: n5
    !f2py intent(hide), depend(fy) :: n5 = shape(fy,0)
    integer :: n6
    !f2py intent(hide), depend(fz) :: n6 = shape(fz,0)
    call force(a1=a1, a2=a2, a3=a3, n=n, rx=rx, ry=ry, rz=rz, z=z, rc=rc, rd=rd, fx=fx, fy=fy, fz=fz)
end subroutine f90wrap_force

subroutine f90wrap_stress(a1, a2, a3, n, rx, ry, rz, z, rc, rd, s, n0, n1, n2, n3)
    use real_space_electrostatic_sum, only: stress
    implicit none
    
    real(8), intent(in), dimension(3) :: a1
    real(8), intent(in), dimension(3) :: a2
    real(8), intent(in), dimension(3) :: a3
    integer, intent(in) :: n
    real(8), intent(in), dimension(n0) :: rx
    real(8), intent(in), dimension(n1) :: ry
    real(8), intent(in), dimension(n2) :: rz
    real(8), intent(in), dimension(n3) :: z
    real(8), intent(in) :: rc
    real(8), intent(in) :: rd
    real(8), dimension(6), intent(inout) :: s
    integer :: n0
    !f2py intent(hide), depend(rx) :: n0 = shape(rx,0)
    integer :: n1
    !f2py intent(hide), depend(ry) :: n1 = shape(ry,0)
    integer :: n2
    !f2py intent(hide), depend(rz) :: n2 = shape(rz,0)
    integer :: n3
    !f2py intent(hide), depend(z) :: n3 = shape(z,0)
    call stress(a1=a1, a2=a2, a3=a3, n=n, rx=rx, ry=ry, rz=rz, z=z, rc=rc, rd=rd, s=s)
end subroutine f90wrap_stress

subroutine f90wrap_invert_3x3(a, b)
    use real_space_electrostatic_sum, only: invert_3x3
    implicit none
    
    real(8), intent(in), dimension(3,3) :: a
    real(8), dimension(3,3), intent(inout) :: b
    call invert_3x3(a=a, b=b)
end subroutine f90wrap_invert_3x3

subroutine f90wrap_real_space_electrostatic_sum__get__dp(f90wrap_dp)
    use real_space_electrostatic_sum, only: real_space_electrostatic_sum_dp => dp
    implicit none
    integer, intent(out) :: f90wrap_dp
    
    f90wrap_dp = real_space_electrostatic_sum_dp
end subroutine f90wrap_real_space_electrostatic_sum__get__dp

subroutine f90wrap_real_space_electrostatic_sum__get__pi(f90wrap_pi)
    use real_space_electrostatic_sum, only: real_space_electrostatic_sum_pi => pi
    implicit none
    real(8), intent(out) :: f90wrap_pi
    
    f90wrap_pi = real_space_electrostatic_sum_pi
end subroutine f90wrap_real_space_electrostatic_sum__get__pi

subroutine f90wrap_real_space_electrostatic_sum__get__sqrt_pi(f90wrap_sqrt_pi)
    use real_space_electrostatic_sum, only: real_space_electrostatic_sum_sqrt_pi => sqrt_pi
    implicit none
    real(8), intent(out) :: f90wrap_sqrt_pi
    
    f90wrap_sqrt_pi = real_space_electrostatic_sum_sqrt_pi
end subroutine f90wrap_real_space_electrostatic_sum__get__sqrt_pi

subroutine f90wrap_real_space_electrostatic_sum__get__one_third(f90wrap_one_third)
    use real_space_electrostatic_sum, only: real_space_electrostatic_sum_one_third => one_third
    implicit none
    real(8), intent(out) :: f90wrap_one_third
    
    f90wrap_one_third = real_space_electrostatic_sum_one_third
end subroutine f90wrap_real_space_electrostatic_sum__get__one_third

! End of module real_space_electrostatic_sum defined in file real_space_electrostatic_sum.f90

