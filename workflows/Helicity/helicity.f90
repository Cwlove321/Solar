program helicity
    use helicity_io, only: helicity_config, magnetic_input, ensure_diagnostic_directories, &
                           parse_arguments, read_input_frame, write_vector_field
    use helicity_solver, only: helicity_result, helicity_workspace, compute_helicity
    implicit none

    type(helicity_config) :: config
    type(magnetic_input) :: snapshot
    type(helicity_workspace) :: workspace
    type(helicity_result) :: result
    integer :: frame, output_unit
    real :: frame_start, frame_end
    character(len=4) :: frame_text
    character(len=2048) :: output_file

    call parse_arguments(config)
    if (config%diagnostics) call ensure_diagnostic_directories(config%input_dir)

    output_file = trim(config%input_dir) // 'out.binary'
    open(newunit=output_unit, file=trim(output_file), form='unformatted', &
         access='stream', status='replace', action='write')

    do frame = config%first_frame, config%last_frame
        call cpu_time(frame_start)
        call read_input_frame(config%input_dir, frame, snapshot)
        call compute_helicity(snapshot, workspace, config%diagnostics, result)
        write(output_unit) result%energy, result%free_energy, result%h_v, &
                           result%h_j, result%magnetic_flux

        if (config%diagnostics) then
            write(output_unit) result%divergence_error, result%curl_a_error
            write(frame_text, '(I4.4)') frame
            call write_vector_field(trim(config%input_dir) // 'potential/potential_' // &
                                    frame_text // '.binary', workspace%bp_x, workspace%bp_y, workspace%bp_z)
            call write_vector_field(trim(config%input_dir) // 'A/A_' // frame_text // '.binary', &
                                    workspace%ax, workspace%ay, workspace%az)
            call write_vector_field(trim(config%input_dir) // 'A/Ap_' // frame_text // '.binary', &
                                    workspace%ap_x, workspace%ap_y, workspace%ap_z)
        end if

        call cpu_time(frame_end)
        write(*, '(A, I0, A, F10.3, A)') 'frame ', frame, ' finished in ', &
                                         frame_end-frame_start, ' s'
    end do

    close(output_unit)
    write(*, '(A, A)') 'results written to ', trim(output_file)
end program helicity
