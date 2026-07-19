module helicity_numerics
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter, public :: dp = real64
    public :: compute_vector_potential, curl_vector, gradient
    public :: solve_poisson, volume_integral

contains

    subroutine gradient(field, axis, derivative, spacing)
        real(dp), intent(in) :: field(:, :, :)
        integer, intent(in) :: axis
        real(dp), intent(out) :: derivative(:, :, :)
        real(dp), intent(in), optional :: spacing
        real(dp) :: h
        integer :: nx, ny, nz

        nx = size(field, 1)
        ny = size(field, 2)
        nz = size(field, 3)
        h = 1.0_dp
        if (present(spacing)) h = spacing

        select case (axis)
        case (1)
            derivative(2:nx-1, :, :) = (field(3:nx, :, :) - field(1:nx-2, :, :)) * 0.5_dp / h
            derivative(1, :, :) = (-3.0_dp * field(1, :, :) + 4.0_dp * field(2, :, :) - field(3, :, :)) * 0.5_dp / h
            derivative(nx, :, :) = (3.0_dp * field(nx, :, :) - 4.0_dp * field(nx-1, :, :) + &
                                     field(nx-2, :, :)) * 0.5_dp / h
        case (2)
            derivative(:, 2:ny-1, :) = (field(:, 3:ny, :) - field(:, 1:ny-2, :)) * 0.5_dp / h
            derivative(:, 1, :) = (-3.0_dp * field(:, 1, :) + 4.0_dp * field(:, 2, :) - field(:, 3, :)) * 0.5_dp / h
            derivative(:, ny, :) = (3.0_dp * field(:, ny, :) - 4.0_dp * field(:, ny-1, :) + &
                                     field(:, ny-2, :)) * 0.5_dp / h
        case (3)
            derivative(:, :, 2:nz-1) = (field(:, :, 3:nz) - field(:, :, 1:nz-2)) * 0.5_dp / h
            derivative(:, :, 1) = (-3.0_dp * field(:, :, 1) + 4.0_dp * field(:, :, 2) - field(:, :, 3)) * 0.5_dp / h
            derivative(:, :, nz) = (3.0_dp * field(:, :, nz) - 4.0_dp * field(:, :, nz-1) + &
                                     field(:, :, nz-2)) * 0.5_dp / h
        case default
            error stop "gradient axis must be 1, 2, or 3"
        end select
    end subroutine gradient


    subroutine curl_vector(ax, ay, az, bx, by, bz, spacing)
        real(dp), intent(in) :: ax(:, :, :), ay(:, :, :), az(:, :, :)
        real(dp), intent(out) :: bx(:, :, :), by(:, :, :), bz(:, :, :)
        real(dp), intent(in) :: spacing
        real(dp), allocatable :: derivative_1(:, :, :), derivative_2(:, :, :)
        integer :: nx, ny, nz

        nx = size(ax, 1)
        ny = size(ax, 2)
        nz = size(ax, 3)
        if (any(shape(ay) /= shape(ax)) .or. any(shape(az) /= shape(ax))) &
            error stop "curl_vector input shapes do not match"
        if (any(shape(bx) /= shape(ax)) .or. any(shape(by) /= shape(ax)) .or. &
            any(shape(bz) /= shape(ax))) error stop "curl_vector output shapes do not match"

        ! The tested implementation used six full-size temporary arrays. Reusing two
        ! preserves the same finite differences while reducing scratch memory by 2/3.
        allocate(derivative_1(nx, ny, nz), derivative_2(nx, ny, nz))

        call gradient(az, 2, derivative_1)
        call gradient(ay, 3, derivative_2)
        bx = (derivative_1 - derivative_2) / spacing

        call gradient(ax, 3, derivative_1)
        call gradient(az, 1, derivative_2)
        by = (derivative_1 - derivative_2) / spacing

        call gradient(ay, 1, derivative_1)
        call gradient(ax, 2, derivative_2)
        bz = (derivative_1 - derivative_2) / spacing
    end subroutine curl_vector


    subroutine compute_vector_potential(bx, by, bz, ax, ay, az, spacing)
        real(dp), intent(in) :: bx(:, :, :), by(:, :, :), bz(:, :, :)
        real(dp), intent(out) :: ax(:, :, :), ay(:, :, :), az(:, :, :)
        real(dp), intent(in) :: spacing
        real(dp), allocatable :: boundary_x(:, :), boundary_y(:, :)
        integer :: i, j, k, nx, ny, nz

        nx = size(bx, 1)
        ny = size(bx, 2)
        nz = size(bx, 3)
        allocate(boundary_x(nx, ny), boundary_y(nx, ny))

        ! DeVore gauge: Az=0. The trapezoidal integration is kept identical to
        ! the calculation used by the tested guojh version.
        az = 0.0_dp
        do i = 1, nx
            do j = 1, ny
                do k = nz, 1, -1
                    if (k == nz) then
                        ax(i, j, k) = 0.0_dp
                        ay(i, j, k) = 0.0_dp
                    else if (k == nz - 1) then
                        ax(i, j, k) = (by(i, j, nz-1) + by(i, j, nz)) / 2.0_dp
                        ay(i, j, k) = (bx(i, j, nz-1) + bx(i, j, nz)) / 2.0_dp
                    else
                        ax(i, j, k) = ax(i, j, k+2) + 2.0_dp * by(i, j, k+1)
                        ay(i, j, k) = ay(i, j, k+2) + 2.0_dp * bx(i, j, k+1)
                    end if
                end do
            end do
        end do

        do i = 1, nx
            do j = 1, ny
                if (j == 1) then
                    boundary_x(i, j) = 0.0_dp
                else if (j == 2) then
                    boundary_x(i, j) = (bz(i, 1, nz) + bz(i, 2, nz)) / 2.0_dp
                else
                    boundary_x(i, j) = boundary_x(i, j-2) + 2.0_dp * bz(i, j-1, nz)
                end if

                if (i == 1) then
                    boundary_y(i, j) = 0.0_dp
                else if (i == 2) then
                    boundary_y(i, j) = (bz(1, j, nz) + bz(2, j, nz)) / 2.0_dp
                else
                    boundary_y(i, j) = boundary_y(i-2, j) + 2.0_dp * bz(i-1, j, nz)
                end if
            end do
        end do

        do k = 1, nz
            ax(:, :, k) = -ax(:, :, k) - 0.5_dp * boundary_x
            ay(:, :, k) =  ay(:, :, k) + 0.5_dp * boundary_y
        end do
        ax = spacing * ax
        ay = spacing * ay
    end subroutine compute_vector_potential


    subroutine solve_poisson(boundary_type, nx, ny, nz, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, potential)
        use mkl_poisson
        character(6), intent(in) :: boundary_type
        integer, intent(in) :: nx, ny, nz
        real(dp), intent(in) :: bd_ax(ny+1, nz+1), bd_bx(ny+1, nz+1)
        real(dp), intent(in) :: bd_ay(nx+1, nz+1), bd_by(nx+1, nz+1)
        real(dp), intent(in) :: bd_az(nx+1, ny+1), bd_bz(nx+1, ny+1)
        real(dp), intent(inout) :: potential(nx+1, ny+1, nz+1)
        integer :: status_code
        integer :: ipar(128)
        real(dp) :: dpar(13 * (nx + ny) / 2 + 9)
        real(dp) :: q
        type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle

        ipar = 0
        q = 0.0_dp
        call D_INIT_HELMHOLTZ_3D(0.0_dp, real(nx, dp), 0.0_dp, real(ny, dp), &
                                 0.0_dp, real(nz, dp), nx, ny, nz, boundary_type, &
                                 q, ipar, dpar, status_code)
        call check_mkl_status("D_INIT_HELMHOLTZ_3D", status_code)

        call D_COMMIT_HELMHOLTZ_3D(potential, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, &
                                   xhandle, yhandle, ipar, dpar, status_code)
        call check_mkl_status("D_COMMIT_HELMHOLTZ_3D", status_code)

        call D_HELMHOLTZ_3D(potential, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, &
                            xhandle, yhandle, ipar, dpar, status_code)
        call check_mkl_status("D_HELMHOLTZ_3D", status_code)

        call FREE_HELMHOLTZ_3D(xhandle, yhandle, ipar, status_code)
        call check_mkl_status("FREE_HELMHOLTZ_3D", status_code)
    end subroutine solve_poisson


    subroutine check_mkl_status(operation, status_code)
        character(*), intent(in) :: operation
        integer, intent(in) :: status_code
        if (status_code < 0) then
            write(*, '(A, A, A, I0)') "MKL Poisson error in ", trim(operation), ": ", status_code
            error stop "MKL Poisson solver failed"
        else if (status_code > 0) then
            write(*, '(A, A, A, I0)') "MKL Poisson warning in ", trim(operation), ": ", status_code
        end if
    end subroutine check_mkl_status


    subroutine volume_integral(array, spacing, integral)
        real(dp), intent(in) :: array(:, :, :)
        real(dp), intent(in) :: spacing
        real(dp), intent(out) :: integral
        real(dp) :: corners, edges, faces
        integer :: i, j, k, nx, ny, nz

        nx = size(array, 1)
        ny = size(array, 2)
        nz = size(array, 3)

        corners = 0.125_dp * (array(1, 1, 1) + array(1, 1, nz) + &
                  array(1, ny, 1) + array(1, ny, nz) + array(nx, 1, 1) + &
                  array(nx, 1, nz) + array(nx, ny, 1) + array(nx, ny, nz))

        edges = 0.0_dp
        do i = 1, nz - 2
            edges = edges + 0.25_dp * (array(1, 1, i+1) + array(1, ny, i+1) + &
                    array(nx, 1, i+1) + array(nx, ny, i+1))
        end do
        do j = 1, ny - 2
            edges = edges + 0.25_dp * (array(1, j+1, 1) + array(1, j+1, nz) + &
                    array(nx, j+1, 1) + array(nx, j+1, nz))
        end do
        do k = 1, nx - 2
            edges = edges + 0.25_dp * (array(k+1, 1, 1) + array(k+1, 1, nz) + &
                    array(k+1, ny, 1) + array(k+1, ny, nz))
        end do

        faces = 0.0_dp
        do j = 1, ny - 2
            do k = 1, nz - 2
                faces = faces + 0.5_dp * (array(1, j+1, k+1) + array(nx, j+1, k+1))
            end do
        end do
        do i = 1, nx - 2
            do k = 1, nz - 2
                faces = faces + 0.5_dp * (array(i+1, 1, k+1) + array(i+1, ny, k+1))
            end do
        end do
        do i = 1, nx - 2
            do j = 1, ny - 2
                faces = faces + 0.5_dp * (array(i+1, j+1, 1) + array(i+1, j+1, nz))
            end do
        end do

        integral = (corners + edges + faces + &
                   sum(array(2:nx-1, 2:ny-1, 2:nz-1))) * spacing**3
    end subroutine volume_integral

end module helicity_numerics
