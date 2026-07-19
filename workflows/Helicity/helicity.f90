program helicity
    use helicity_numerics, only: dp, compute_vector_potential, gradient, solve_poisson, volume_integral
    use helicity_diagnostics, only: compute_error_diagnostics
    implicit none

    real(dp), parameter :: pi = 3.1415926_dp
    character(len=1024) :: input_dir
    integer :: first_frame, last_frame, error_estimation
    integer :: frame, nx, ny, nz, expected_nx, expected_ny, expected_nz
    integer :: input_unit, output_unit
    character(len=4) :: frame_text
    character(len=2048) :: magnetic_file, boundary_file, output_file
    real(dp) :: spacing, energy, free_energy, h_v, h_j, magnetic_flux
    real(dp) :: divergence_error, curl_a_error, frame_start, frame_end
    real(dp), allocatable :: bx(:, :, :), by(:, :, :), bz(:, :, :), scalar(:, :, :)
    real(dp), allocatable :: bp_x(:, :, :), bp_y(:, :, :), bp_z(:, :, :)
    real(dp), allocatable :: bj_x(:, :, :), bj_y(:, :, :), bj_z(:, :, :)
    real(dp), allocatable :: ax(:, :, :), ay(:, :, :), az(:, :, :)
    real(dp), allocatable :: ap_x(:, :, :), ap_y(:, :, :), ap_z(:, :, :)
    real(dp), allocatable :: bd_ax(:, :), bd_bx(:, :), bd_ay(:, :), bd_by(:, :), bd_az(:, :), bd_bz(:, :)

    call parse_arguments(input_dir, first_frame, last_frame, error_estimation)
    if (error_estimation == 1) call ensure_diagnostic_directories(input_dir)

    output_file = trim(input_dir) // 'out.binary'
    open(newunit=output_unit, file=trim(output_file), form='unformatted', &
         access='stream', status='replace', action='write')

    expected_nx = -1
    expected_ny = -1
    expected_nz = -1

    do frame = first_frame, last_frame
        call cpu_time(frame_start)
        write(frame_text, '(I4.4)') frame
        magnetic_file = trim(input_dir) // 'nlfff/nlfff_' // frame_text // '.binary'
        boundary_file = trim(input_dir) // 'bd/bd_' // frame_text // '.binary'
        call require_file(magnetic_file)
        call require_file(boundary_file)

        open(newunit=input_unit, file=trim(magnetic_file), form='unformatted', &
             access='stream', status='old', action='read')
        read(input_unit) spacing
        read(input_unit) nx, ny, nz
        if (frame == first_frame) then
            call allocate_fields(nx, ny, nz)
            expected_nx = nx
            expected_ny = ny
            expected_nz = nz
        else if (nx /= expected_nx .or. ny /= expected_ny .or. nz /= expected_nz) then
            error stop 'all snapshots must use the same dimensions'
        end if
        read(input_unit) bx, by, bz
        close(input_unit)

        open(newunit=input_unit, file=trim(boundary_file), form='unformatted', &
             access='stream', status='old', action='read')
        read(input_unit) bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz
        close(input_unit)

        scalar = 0.0_dp
        call solve_poisson('NNNNNN', nx-1, ny-1, nz-1, bd_ax, bd_bx, &
                           bd_ay, bd_by, bd_az, bd_bz, scalar)
        call gradient(scalar, 1, bp_x)
        call gradient(scalar, 2, bp_y)
        call gradient(scalar, 3, bp_z)

        ! Enforce equality of the normal components on all six faces.
        bp_x(1, :, :) = bx(1, :, :)
        bp_x(nx, :, :) = bx(nx, :, :)
        bp_y(:, 1, :) = by(:, 1, :)
        bp_y(:, ny, :) = by(:, ny, :)
        bp_z(:, :, 1) = bz(:, :, 1)
        bp_z(:, :, nz) = bz(:, :, nz)

        bj_x = bx - bp_x
        bj_y = by - bp_y
        bj_z = bz - bp_z
        call compute_vector_potential(bp_x, bp_y, bp_z, ap_x, ap_y, ap_z, spacing)
        call compute_vector_potential(bx, by, bz, ax, ay, az, spacing)

        scalar = (bx**2 + by**2 + bz**2) / (8.0_dp * pi)
        call volume_integral(scalar, spacing, energy)
        scalar = (bj_x**2 + bj_y**2 + bj_z**2) / (8.0_dp * pi)
        call volume_integral(scalar, spacing, free_energy)
        scalar = (ax+ap_x)*(bx-bp_x) + (ay+ap_y)*(by-bp_y) + (az+ap_z)*(bz-bp_z)
        call volume_integral(scalar, spacing, h_v)
        scalar = (ax-ap_x)*(bx-bp_x) + (ay-ap_y)*(by-bp_y) + (az-ap_z)*(bz-bp_z)
        call volume_integral(scalar, spacing, h_j)
        magnetic_flux = sum(abs(bz(:, :, 1))) * spacing**2

        write(output_unit) energy, free_energy, h_v, h_j, magnetic_flux

        if (error_estimation == 1) then
            call compute_error_diagnostics(bx, by, bz, bp_x, bp_y, bp_z, &
                                           bj_x, bj_y, bj_z, ax, ay, az, spacing, &
                                           energy, divergence_error, curl_a_error)
            write(output_unit) divergence_error, curl_a_error
            call write_vector_field(trim(input_dir) // 'potential/potential_' // frame_text // '.binary', &
                                    bp_x, bp_y, bp_z)
            call write_vector_field(trim(input_dir) // 'A/A_' // frame_text // '.binary', ax, ay, az)
            call write_vector_field(trim(input_dir) // 'A/Ap_' // frame_text // '.binary', ap_x, ap_y, ap_z)
        end if

        call cpu_time(frame_end)
        write(*, '(A, I0, A, F10.3, A)') 'frame ', frame, ' finished in ', frame_end-frame_start, ' s'
    end do

    close(output_unit)
    write(*, '(A, A)') 'results written to ', trim(output_file)

contains

    subroutine parse_arguments(directory, start_index, end_index, diagnostics)
        character(len=*), intent(out) :: directory
        integer, intent(out) :: start_index, end_index, diagnostics
        character(len=1024) :: argument
        integer :: argument_count, io_status

        directory = './helicity_input/'
        start_index = 0
        end_index = 79
        diagnostics = 0
        argument_count = command_argument_count()

        if (argument_count > 0) then
            call get_command_argument(1, argument)
            if (trim(argument) == '-h' .or. trim(argument) == '--help') then
                call print_usage()
                stop
            end if
            if (len_trim(argument) > 0) directory = trim(argument)
        end if
        if (argument_count > 1) then
            call get_command_argument(2, argument)
            read(argument, *, iostat=io_status) start_index
            if (io_status /= 0) error stop 'START must be an integer'
        end if
        if (argument_count > 2) then
            call get_command_argument(3, argument)
            read(argument, *, iostat=io_status) end_index
            if (io_status /= 0) error stop 'END must be an integer'
        end if
        if (argument_count > 3) then
            call get_command_argument(4, argument)
            read(argument, *, iostat=io_status) diagnostics
            if (io_status /= 0) error stop 'ERROR_ESTIMATION must be 0 or 1'
        end if
        if (argument_count > 4) error stop 'too many arguments; use --help'
        if (start_index < 0 .or. end_index < start_index) &
            error stop 'frame range must satisfy 0 <= START <= END'
        if (diagnostics /= 0 .and. diagnostics /= 1) &
            error stop 'ERROR_ESTIMATION must be 0 or 1'
        if (directory(len_trim(directory):len_trim(directory)) /= '/') &
            directory = trim(directory) // '/'
    end subroutine parse_arguments


    subroutine print_usage()
        print '(A)', 'Usage: ./helicity [INPUT_DIR] [START] [END] [ERROR_ESTIMATION]'
        print '(A)', '  INPUT_DIR         contains nlfff/ and bd/ (default: ./helicity_input)'
        print '(A)', '  START, END        inclusive frame numbers (default: 0, 79)'
        print '(A)', '  ERROR_ESTIMATION  0 or 1 (default: 0)'
    end subroutine print_usage


    subroutine allocate_fields(size_x, size_y, size_z)
        integer, intent(in) :: size_x, size_y, size_z

        if (min(size_x, size_y, size_z) < 3) error stop 'each grid dimension must be at least 3'
        allocate(bx(size_x, size_y, size_z), by(size_x, size_y, size_z), bz(size_x, size_y, size_z))
        allocate(scalar(size_x, size_y, size_z))
        allocate(bp_x(size_x, size_y, size_z), bp_y(size_x, size_y, size_z), bp_z(size_x, size_y, size_z))
        allocate(bj_x(size_x, size_y, size_z), bj_y(size_x, size_y, size_z), bj_z(size_x, size_y, size_z))
        allocate(ax(size_x, size_y, size_z), ay(size_x, size_y, size_z), az(size_x, size_y, size_z))
        allocate(ap_x(size_x, size_y, size_z), ap_y(size_x, size_y, size_z), ap_z(size_x, size_y, size_z))
        allocate(bd_ax(size_y, size_z), bd_bx(size_y, size_z))
        allocate(bd_ay(size_x, size_z), bd_by(size_x, size_z))
        allocate(bd_az(size_x, size_y), bd_bz(size_x, size_y))
    end subroutine allocate_fields


    subroutine require_file(filename)
        character(len=*), intent(in) :: filename
        logical :: exists
        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            write(*, '(A, A)') 'missing input file: ', trim(filename)
            error stop 'required input file does not exist'
        end if
    end subroutine require_file


    subroutine ensure_diagnostic_directories(directory)
        character(len=*), intent(in) :: directory
        character(len=4096) :: command
        integer :: exit_status
        command = 'mkdir -p "' // trim(directory) // 'A" "' // trim(directory) // 'potential"'
        call execute_command_line(trim(command), exitstat=exit_status)
        if (exit_status /= 0) error stop 'could not create A/ and potential/ directories'
    end subroutine ensure_diagnostic_directories


    subroutine write_vector_field(filename, field_x, field_y, field_z)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: field_x(:, :, :), field_y(:, :, :), field_z(:, :, :)
        integer :: file_unit
        open(newunit=file_unit, file=trim(filename), form='unformatted', &
             access='stream', status='replace', action='write')
        write(file_unit) size(field_x, 1), size(field_x, 2), size(field_x, 3)
        write(file_unit) field_x, field_y, field_z
        close(file_unit)
    end subroutine write_vector_field

end program helicity
