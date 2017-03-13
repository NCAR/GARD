module atm_mod
    use data_structures
    use model_constants

    use gcm_mod,         only : read_gcm
    use gefs_mod,         only : read_gefs
    ! use reanalysis_mod,  only : read_reanalysis

    implicit none
    ! private
    ! public read_atm

contains

    function read_atm(options) result(atm_data)
        implicit none
        class(atm_config) :: options
        type(atm) :: atm_data

        ! determine the type of data we are trying to read and call the relevant function
        select case (options%data_type)
            case (kGCM_TYPE)
                atm_data = read_gcm(options)
            ! case (kREANALYSIS_TYPE)
            !     atm_data = read_reanalysis(options)
            case (kGEFS_TYPE)
                atm_data = read_gefs(options)
        end select

    end function read_atm

end module atm_mod
