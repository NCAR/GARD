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
    use config_mod,     only : config
    use gcm_mod,        only : gcm
    use obs_mod,        only : observations
    use reanalysis_mod, only : reanalysis
    ! use results_mod,    only :: results_obj
    ! use ds_model_mod,   only :: run_model
    ! use output_mod,     only :: output_handler
    
    implicit none
    
    type(config)        :: options
    type(gcm)           :: gcm_data
    type(observations)  :: obs_data
    type(reanalysis)    :: current_data
    ! type(output_handler):: output
    ! type(results_obj)   :: results
    
    integer :: i
    
    options = config()
    
    call obs_data%init(options)
    gcm_data = gcm(options)
    call current_data%init(options)
    
    print*, "options ",      trim(options%as_string())
    print*, "obs_data ",     trim(obs_data%as_string())
    print*, "gcm_data ",     trim(gcm_data%as_string())
    print*, "current_data ", trim(current_data%as_string())
    ! output       = output_handler(options)

    ! results      = results_obj(options)
    
    ! do i = options%first_point, options%last_point
    !     
    !     call run_model(i, gcm_data, obs_data, current_data, results, options)
    !     call output%write(i, results)
    !     
    ! end do
    ! 
    ! call output%close()

end program downscale
