!>------------------------------------------------------------
!!  Main Downscaling Code
!!
!!  Makes calls to read configuration options, setup data structures, run the model, and output the final data
!!  
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
program downscale
    use data_structures
    use model_constants
    use config_mod,     only : read_config
    use atm_mod,        only : read_atm
    use init_mod,       only : model_init
    ! use obs_mod,        only : read_obs
    ! use stats_mod,      only : sdown
    ! use output_mod,     only : write_output
    
    implicit none
    
    type(config)        :: options
    type(atm)           :: training_atm, predictions
    ! type(obs)          :: training_obs
    ! type(results)      :: output
    
    integer :: i
    
    options = read_config()
    
    call model_init(options)
    
    ! read in the training atmospheric data (e.g. reanalysis or GEFS)
    training_atm = read_atm(options%training)
    ! read in the training surface data (e.g. Maurer et al. )
    ! training_obs = read_obs(options%obs)
    
    ! read in the atmospheric predictor data (e.g. GCM or GEFS)
    ! predictions  = read_atm(options%predictions)
    
    if (options%debug) then
        ! print*, "options ",         trim(options%name)
        ! print*, "obs ",             trim(training_obs%name)
        print*, "atm:training ",    trim(training_atm%name)
        ! print*, "atm:predictions ", trim(predictions%name)
    endif
    
    ! output = sdown(training_atm, training_obs, predictions, options)
    
    ! call write_output(output, options)
    
end program downscale
