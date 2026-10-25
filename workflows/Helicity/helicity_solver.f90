module helicity_solver
    use helicity_numerics, only: dp, compute_vector_potential, gradient, solve_poisson, volume_integral
    use helicity_diagnostics, only: compute_error_diagnostics
    use helicity_io, only: magnetic_input
    implicit none
    private

    real(dp), parameter :: pi = 3.1415926_dp

    type, public :: helicity_result
        real(dp) :: energy = 0.0_dp
        real(dp) :: free_energy = 0.0_dp
        real(dp) :: h_v = 0.0_dp
        real(dp) :: h_j = 0.0_dp
        real(dp) :: magnetic_flux = 0.0_dp
        real(dp) :: divergence_error = 0.0_dp
        real(dp) :: curl_a_error = 0.0_dp
    end type helicity_result

    type, public :: helicity_workspace
        real(dp), allocatable :: scalar(:, :, :)
        real(dp), allocatable :: bp_x(:, :, :), bp_y(:, :, :), bp_z(:, :, :)
        real(dp), allocatable :: bj_x(:, :, :), bj_y(:, :, :), bj_z(:, :, :)
        real(dp), allocatable :: ax(:, :, :), ay(:, :, :), az(:, :, :)
        real(dp), allocatable :: ap_x(:, :, :), ap_y(:, :, :), ap_z(:, :, :)
    end type helicity_workspace

    public :: compute_helicity

contains

    subroutine compute_helicity(snapshot, workspace, diagnostics, result)
        type(magnetic_input), intent(in) :: snapshot
        type(helicity_workspace), intent(inout) :: workspace
        logical, intent(in) :: diagnostics
        type(helicity_result), intent(out) :: result

        call ensure_workspace(workspace, snapshot%nx, snapshot%ny, snapshot%nz)

        workspace%scalar = 0.0_dp
        call solve_poisson('NNNNNN', snapshot%nx-1, snapshot%ny-1, snapshot%nz-1, &
                           snapshot%bd_ax, snapshot%bd_bx, snapshot%bd_ay, snapshot%bd_by, &
                           snapshot%bd_az, snapshot%bd_bz, workspace%scalar)
        call gradient(workspace%scalar, 1, workspace%bp_x)
        call gradient(workspace%scalar, 2, workspace%bp_y)
        call gradient(workspace%scalar, 3, workspace%bp_z)

        ! The potential and input fields must have identical normal components.
        workspace%bp_x(1, :, :) = snapshot%bx(1, :, :)
        workspace%bp_x(snapshot%nx, :, :) = snapshot%bx(snapshot%nx, :, :)
        workspace%bp_y(:, 1, :) = snapshot%by(:, 1, :)
        workspace%bp_y(:, snapshot%ny, :) = snapshot%by(:, snapshot%ny, :)
        workspace%bp_z(:, :, 1) = snapshot%bz(:, :, 1)
        workspace%bp_z(:, :, snapshot%nz) = snapshot%bz(:, :, snapshot%nz)

        workspace%bj_x = snapshot%bx - workspace%bp_x
        workspace%bj_y = snapshot%by - workspace%bp_y
        workspace%bj_z = snapshot%bz - workspace%bp_z
        call compute_vector_potential(workspace%bp_x, workspace%bp_y, workspace%bp_z, &
                                      workspace%ap_x, workspace%ap_y, workspace%ap_z, snapshot%spacing)
        call compute_vector_potential(snapshot%bx, snapshot%by, snapshot%bz, &
                                      workspace%ax, workspace%ay, workspace%az, snapshot%spacing)

        workspace%scalar = (snapshot%bx**2 + snapshot%by**2 + snapshot%bz**2) / (8.0_dp * pi)
        call volume_integral(workspace%scalar, snapshot%spacing, result%energy)
        workspace%scalar = (workspace%bj_x**2 + workspace%bj_y**2 + workspace%bj_z**2) / (8.0_dp * pi)
        call volume_integral(workspace%scalar, snapshot%spacing, result%free_energy)
        workspace%scalar = (workspace%ax+workspace%ap_x)*workspace%bj_x + &
                           (workspace%ay+workspace%ap_y)*workspace%bj_y + &
                           (workspace%az+workspace%ap_z)*workspace%bj_z
        call volume_integral(workspace%scalar, snapshot%spacing, result%h_v)
        workspace%scalar = (workspace%ax-workspace%ap_x)*workspace%bj_x + &
                           (workspace%ay-workspace%ap_y)*workspace%bj_y + &
                           (workspace%az-workspace%ap_z)*workspace%bj_z
        call volume_integral(workspace%scalar, snapshot%spacing, result%h_j)
        result%magnetic_flux = sum(abs(snapshot%bz(:, :, 1))) * snapshot%spacing**2

        if (diagnostics) then
            call compute_error_diagnostics(snapshot%bx, snapshot%by, snapshot%bz, &
                                           workspace%bp_x, workspace%bp_y, workspace%bp_z, &
                                           workspace%bj_x, workspace%bj_y, workspace%bj_z, &
                                           workspace%ax, workspace%ay, workspace%az, snapshot%spacing, &
                                           result%energy, result%divergence_error, result%curl_a_error)
        end if
    end subroutine compute_helicity


    subroutine ensure_workspace(workspace, nx, ny, nz)
        type(helicity_workspace), intent(inout) :: workspace
        integer, intent(in) :: nx, ny, nz

        if (allocated(workspace%scalar)) then
            if (any(shape(workspace%scalar) /= [nx, ny, nz])) &
                error stop 'workspace dimensions do not match input dimensions'
            return
        end if
        allocate(workspace%scalar(nx, ny, nz))
        allocate(workspace%bp_x(nx, ny, nz), workspace%bp_y(nx, ny, nz), workspace%bp_z(nx, ny, nz))
        allocate(workspace%bj_x(nx, ny, nz), workspace%bj_y(nx, ny, nz), workspace%bj_z(nx, ny, nz))
        allocate(workspace%ax(nx, ny, nz), workspace%ay(nx, ny, nz), workspace%az(nx, ny, nz))
        allocate(workspace%ap_x(nx, ny, nz), workspace%ap_y(nx, ny, nz), workspace%ap_z(nx, ny, nz))
    end subroutine ensure_workspace

end module helicity_solver
