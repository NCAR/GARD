!>------------------------------------------------------------
!!  Main Downscaling Code
!!
!!  Makes calls to read configuration options, setup data structures, run the model, and output the final data
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
program stat_down
    use data_structures
    use model_constants
    use config_mod,     only : read_config
    use atm_mod,        only : read_atm
    use init_mod,       only : model_init
    use obs_mod,        only : read_obs
    use geo,            only : geo_LUT
    use io_routines,    only : io_write
    use downscaling_mod,only : downscale
    use output_mod,     only : write_output

    implicit none

    type(config)        :: options
    type(atm)           :: training_atm, predictions
    type(obs)           :: training_obs
    type(results)       :: output
    character(len=MAXSTRINGLENGTH) :: name

    integer :: i

    options = read_config()

    call model_init(options)

    ! read in the atmospheric predictor data (e.g. GCM or GEFS)
    write(*,*) "=========================================="
    write(*,*) ""
    write(*,*) "Reading predictor"
    predictions  = read_atm(options%prediction)

    ! read in the training atmospheric data (e.g. reanalysis or GEFS)
    write(*,*) ""
    write(*,*) "Reading training"
    training_atm = read_atm(options%training)
    ! read in the training surface data (e.g. Maurer et al., Newman et al., Livneh et al., DAYMET )
    write(*,*) ""
    write(*,*) "Reading obs"
    training_obs = read_obs(options%obs)

    write(*,*) ""
    write(*,*) "=========================================="
    write(*,*) "options         ", trim(options%name)
    write(*,*) "obs             ", trim(training_obs%name)
    write(*,*) "atm:training    ", trim(training_atm%name)
    write(*,*) "atm:predictions ", trim(predictions%name)
    write(*,*) "=========================================="
    write(*,*) ""
    write(*,*) "Developing geographic interpolation..."

    call geo_LUT(training_obs, training_atm, options%prediction%interpolation_method)
    call geo_LUT(training_obs, predictions, options%training%interpolation_method)

    if (options%debug) then
        write(*,*) "=========================================="
        write(*,*) ""
        write(*,*) "Writing input data"
        do i=1,size(training_obs%variables)
            if (trim(options%obs%preloaded) == "") then
                options%obs%preloaded = "obs_preload_"
            endif
            call io_write(trim(options%obs%preloaded)//trim(training_obs%variables(i)%name)//".nc","data", &
                            training_obs%variables(i)%data)
        enddo
        do i=1,size(training_atm%variables)
            if (trim(options%training%preloaded) == "") then
                options%training%preloaded = "training_preload_"
            endif
            call io_write(trim(options%training%preloaded)//trim(training_atm%variables(i)%name)//".nc","data", &
                            training_atm%variables(i)%data)
            if (trim(options%prediction%preloaded) == "") then
                options%prediction%preloaded = "prediction_preload_"
            endif
            call io_write(trim(options%prediction%preloaded)//trim(predictions%variables(i)%name)//".nc","data", &
                            predictions%variables(i)%data)
        enddo
    endif

    write(*,*) ""
    write(*,*) "=========================================="
    write(*,*) ""
    write(*,*) "Running Downscaling Code"
    call downscale(training_atm, training_obs, predictions, output, options)

    write(*,*) "=========================================="
    write(*,*) ""
    write(*,*) "Writing output"
    call write_output(output, options)

end program stat_down
