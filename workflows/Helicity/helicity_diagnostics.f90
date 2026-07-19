module helicity_diagnostics
    use helicity_numerics, only: dp, curl_vector, gradient, solve_poisson, volume_integral
    implicit none
    private
    public :: compute_error_diagnostics

contains

    subroutine compute_error_diagnostics(bx, by, bz, bp_x, bp_y, bp_z, &
                                         bj_x, bj_y, bj_z, ax, ay, az, &
                                         spacing, energy, divergence_error, curl_a_error)
        real(dp), intent(in) :: bx(:, :, :), by(:, :, :), bz(:, :, :)
        real(dp), intent(in) :: bp_x(:, :, :), bp_y(:, :, :), bp_z(:, :, :)
        real(dp), intent(in) :: bj_x(:, :, :), bj_y(:, :, :), bj_z(:, :, :)
        real(dp), intent(in) :: ax(:, :, :), ay(:, :, :), az(:, :, :)
        real(dp), intent(in) :: spacing, energy
        real(dp), intent(out) :: divergence_error, curl_a_error

        real(dp), parameter :: pi = 3.1415926_dp
        real(dp), allocatable :: work_x(:, :, :), work_y(:, :, :), work_z(:, :, :)
        real(dp), allocatable :: bpn_x(:, :, :), bpn_y(:, :, :), bpn_z(:, :, :)
        real(dp), allocatable :: bjn_x(:, :, :), bjn_y(:, :, :), bjn_z(:, :, :)
        real(dp), allocatable :: scalar(:, :, :)
        real(dp), allocatable :: zero_x(:, :), zero_y(:, :), zero_z(:, :)
        real(dp) :: divergence_energy
        integer :: nx, ny, nz

        nx = size(bx, 1)
        ny = size(bx, 2)
        nz = size(bx, 3)
        allocate(work_x(nx, ny, nz), work_y(nx, ny, nz), work_z(nx, ny, nz))
        allocate(bpn_x(nx, ny, nz), bpn_y(nx, ny, nz), bpn_z(nx, ny, nz))
        allocate(bjn_x(nx, ny, nz), bjn_y(nx, ny, nz), bjn_z(nx, ny, nz))
        allocate(scalar(nx, ny, nz))
        allocate(zero_x(ny, nz), zero_y(nx, nz), zero_z(nx, ny))
        zero_x = 0.0_dp
        zero_y = 0.0_dp
        zero_z = 0.0_dp

        call gradient(bp_x, 1, work_x)
        call gradient(bp_y, 2, work_y)
        call gradient(bp_z, 3, work_z)
        scalar = -(work_x + work_y + work_z)
        call solve_poisson('NNNNNN', nx-1, ny-1, nz-1, zero_x, zero_x, &
                           zero_y, zero_y, zero_z, zero_z, scalar)
        call gradient(scalar, 1, bpn_x)
        call gradient(scalar, 2, bpn_y)
        call gradient(scalar, 3, bpn_z)
        call zero_normal_components(bpn_x, bpn_y, bpn_z)

        call gradient(bj_x, 1, work_x)
        call gradient(bj_y, 2, work_y)
        call gradient(bj_z, 3, work_z)
        scalar = -(work_x + work_y + work_z)
        call solve_poisson('NNNNNN', nx-1, ny-1, nz-1, zero_x, zero_x, &
                           zero_y, zero_y, zero_z, zero_z, scalar)
        call gradient(scalar, 1, bjn_x)
        call gradient(scalar, 2, bjn_y)
        call gradient(scalar, 3, bjn_z)
        call zero_normal_components(bjn_x, bjn_y, bjn_z)

        ! Valori et al. energy decomposition, kept algebraically identical to
        ! the tested implementation.
        scalar = (bjn_x**2 + bjn_y**2 + bjn_z**2 + &
                  bpn_x**2 + bpn_y**2 + bpn_z**2) / (8.0_dp * pi) + &
                 ((bp_x-bpn_x)*bjn_x + (bp_y-bpn_y)*bjn_y + (bp_z-bpn_z)*bjn_z + &
                  (bp_x-bpn_x)*bpn_x + (bp_y-bpn_y)*bpn_y + (bp_z-bpn_z)*bpn_z + &
                  (bp_x-bpn_x)*(bj_x-bjn_x) + (bp_y-bpn_y)*(bj_y-bjn_y) + &
                  (bp_z-bpn_z)*(bj_z-bjn_z) + &
                  bpn_x*(bj_x-bjn_x) + bpn_y*(bj_y-bjn_y) + bpn_z*(bj_z-bjn_z) + &
                  bjn_x*(bj_x-bjn_x) + bjn_y*(bj_y-bjn_y) + bjn_z*(bj_z-bjn_z) + &
                  bjn_x*bpn_x + bjn_y*bpn_y + bjn_z*bpn_z) / (4.0_dp * pi)
        call volume_integral(scalar, spacing, divergence_energy)
        divergence_error = divergence_energy / energy

        call curl_vector(ax, ay, az, work_x, work_y, work_z, spacing)
        curl_a_error = sum(work_x**2 + work_y**2 + work_z**2) / &
                       sum(bx**2 + by**2 + bz**2)
    end subroutine compute_error_diagnostics


    subroutine zero_normal_components(field_x, field_y, field_z)
        real(dp), intent(inout) :: field_x(:, :, :), field_y(:, :, :), field_z(:, :, :)
        field_x(1, :, :) = 0.0_dp
        field_x(size(field_x, 1), :, :) = 0.0_dp
        field_y(:, 1, :) = 0.0_dp
        field_y(:, size(field_y, 2), :) = 0.0_dp
        field_z(:, :, 1) = 0.0_dp
        field_z(:, :, size(field_z, 3)) = 0.0_dp
    end subroutine zero_normal_components

end module helicity_diagnostics
