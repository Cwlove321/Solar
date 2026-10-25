module helicity_io
    use iso_fortran_env, only: int32, int64
    use helicity_numerics, only: dp
    implicit none
    private

    type, public :: helicity_config
        character(len=1024) :: input_dir = './helicity_input/'
        integer :: first_frame = 0
        integer :: last_frame = 79
        logical :: diagnostics = .false.
    end type helicity_config

    type, public :: magnetic_input
        integer :: nx = 0, ny = 0, nz = 0
        real(dp) :: spacing = 0.0_dp
        real(dp), allocatable :: bx(:, :, :), by(:, :, :), bz(:, :, :)
        real(dp), allocatable :: bd_ax(:, :), bd_bx(:, :)
        real(dp), allocatable :: bd_ay(:, :), bd_by(:, :)
        real(dp), allocatable :: bd_az(:, :), bd_bz(:, :)
    end type magnetic_input

    public :: ensure_diagnostic_directories, parse_arguments
    public :: read_input_frame, write_vector_field

contains

    subroutine parse_arguments(config)
        type(helicity_config), intent(out) :: config
        character(len=1024) :: argument
        integer :: argument_count, diagnostics_flag, io_status

        config = helicity_config()
        argument_count = command_argument_count()
        if (argument_count > 4) error stop 'too many arguments; use --help'

        if (argument_count > 0) then
            call get_command_argument(1, argument)
            if (trim(argument) == '-h' .or. trim(argument) == '--help') then
                call print_usage()
                stop
            end if
            if (len_trim(argument) > 0) config%input_dir = trim(argument)
        end if
        if (argument_count > 1) then
            call get_command_argument(2, argument)
            read(argument, *, iostat=io_status) config%first_frame
            if (io_status /= 0) error stop 'START must be an integer'
        end if
        if (argument_count > 2) then
            call get_command_argument(3, argument)
            read(argument, *, iostat=io_status) config%last_frame
            if (io_status /= 0) error stop 'END must be an integer'
        end if
        diagnostics_flag = 0
        if (argument_count > 3) then
            call get_command_argument(4, argument)
            read(argument, *, iostat=io_status) diagnostics_flag
            if (io_status /= 0) error stop 'ERROR_ESTIMATION must be 0 or 1'
        end if
        if (config%first_frame < 0 .or. config%last_frame < config%first_frame) &
            error stop 'frame range must satisfy 0 <= START <= END'
        if (diagnostics_flag /= 0 .and. diagnostics_flag /= 1) &
            error stop 'ERROR_ESTIMATION must be 0 or 1'
        config%diagnostics = diagnostics_flag == 1
        if (config%input_dir(len_trim(config%input_dir):len_trim(config%input_dir)) /= '/') &
            config%input_dir = trim(config%input_dir) // '/'
    end subroutine parse_arguments


    subroutine print_usage()
        print '(A)', 'Usage: ./helicity [INPUT_DIR] [START] [END] [ERROR_ESTIMATION]'
        print '(A)', '  INPUT_DIR         contains nlfff/ and bd/ (default: ./helicity_input)'
        print '(A)', '  START, END        inclusive frame numbers (default: 0, 79)'
        print '(A)', '  ERROR_ESTIMATION  0 or 1 (default: 0)'
    end subroutine print_usage


    subroutine read_input_frame(directory, frame, snapshot)
        character(len=*), intent(in) :: directory
        integer, intent(in) :: frame
        type(magnetic_input), intent(inout) :: snapshot
        character(len=4) :: frame_text
        character(len=2048) :: magnetic_file, boundary_file
        integer :: input_unit, io_status
        integer(int32) :: dimensions(3)
        integer(int64) :: expected_bytes

        write(frame_text, '(I4.4)') frame
        magnetic_file = trim(directory) // 'nlfff/nlfff_' // frame_text // '.binary'
        boundary_file = trim(directory) // 'bd/bd_' // frame_text // '.binary'
        call require_file(magnetic_file)
        call require_file(boundary_file)

        open(newunit=input_unit, file=trim(magnetic_file), form='unformatted', &
             access='stream', status='old', action='read', iostat=io_status)
        call check_io(io_status, 'could not open magnetic input')
        read(input_unit, iostat=io_status) snapshot%spacing, dimensions
        call check_io(io_status, 'magnetic input header is incomplete')
        call ensure_input_shape(snapshot, int(dimensions(1)), int(dimensions(2)), int(dimensions(3)))
        expected_bytes = 20_int64 + 24_int64 * int(snapshot%nx, int64) * &
                         int(snapshot%ny, int64) * int(snapshot%nz, int64)
        call require_size(magnetic_file, expected_bytes)
        read(input_unit, iostat=io_status) snapshot%bx, snapshot%by, snapshot%bz
        call check_io(io_status, 'magnetic input data is incomplete')
        close(input_unit)

        expected_bytes = 16_int64 * (int(snapshot%ny, int64) * int(snapshot%nz, int64) + &
                         int(snapshot%nx, int64) * int(snapshot%nz, int64) + &
                         int(snapshot%nx, int64) * int(snapshot%ny, int64))
        call require_size(boundary_file, expected_bytes)
        open(newunit=input_unit, file=trim(boundary_file), form='unformatted', &
             access='stream', status='old', action='read', iostat=io_status)
        call check_io(io_status, 'could not open boundary input')
        read(input_unit, iostat=io_status) snapshot%bd_ax, snapshot%bd_bx, &
             snapshot%bd_ay, snapshot%bd_by, snapshot%bd_az, snapshot%bd_bz
        call check_io(io_status, 'boundary input data is incomplete')
        close(input_unit)
    end subroutine read_input_frame


    subroutine ensure_input_shape(snapshot, nx, ny, nz)
        type(magnetic_input), intent(inout) :: snapshot
        integer, intent(in) :: nx, ny, nz

        if (min(nx, ny, nz) < 3) error stop 'each grid dimension must be at least 3'
        if (allocated(snapshot%bx)) then
            if (nx /= snapshot%nx .or. ny /= snapshot%ny .or. nz /= snapshot%nz) &
                error stop 'all snapshots must use the same dimensions'
            return
        end if
        snapshot%nx = nx
        snapshot%ny = ny
        snapshot%nz = nz
        allocate(snapshot%bx(nx, ny, nz), snapshot%by(nx, ny, nz), snapshot%bz(nx, ny, nz))
        allocate(snapshot%bd_ax(ny, nz), snapshot%bd_bx(ny, nz))
        allocate(snapshot%bd_ay(nx, nz), snapshot%bd_by(nx, nz))
        allocate(snapshot%bd_az(nx, ny), snapshot%bd_bz(nx, ny))
    end subroutine ensure_input_shape


    subroutine require_file(filename)
        character(len=*), intent(in) :: filename
        logical :: exists

        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            write(*, '(A, A)') 'missing input file: ', trim(filename)
            error stop 'required input file does not exist'
        end if
    end subroutine require_file


    subroutine require_size(filename, expected_bytes)
        character(len=*), intent(in) :: filename
        integer(int64), intent(in) :: expected_bytes
        integer(int64) :: actual_bytes

        inquire(file=trim(filename), size=actual_bytes)
        if (actual_bytes /= expected_bytes) then
            write(*, '(A, A)') 'invalid input size: ', trim(filename)
            write(*, '(A, I0, A, I0)') 'expected bytes=', expected_bytes, ', actual bytes=', actual_bytes
            error stop 'input file size does not match its declared grid'
        end if
    end subroutine require_size


    subroutine check_io(io_status, message)
        integer, intent(in) :: io_status
        character(len=*), intent(in) :: message
        if (io_status /= 0) then
            write(*, '(A)') trim(message)
            error stop 'binary I/O failed'
        end if
    end subroutine check_io


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
        integer :: file_unit, io_status
        integer(int32) :: dimensions(3)

        dimensions = int(shape(field_x), int32)
        open(newunit=file_unit, file=trim(filename), form='unformatted', &
             access='stream', status='replace', action='write', iostat=io_status)
        call check_io(io_status, 'could not open diagnostic output')
        write(file_unit, iostat=io_status) dimensions, field_x, field_y, field_z
        call check_io(io_status, 'could not write diagnostic output')
        close(file_unit)
    end subroutine write_vector_field

end module helicity_io
