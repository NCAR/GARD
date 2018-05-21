!>------------------------------------------------
!! Handle all IO for Observational data
!!
!! Loops through Obs variables reading data and computing basic statistics
!! Reads time and lat/lon variables from the obs files as well.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module obs_mod

    use data_structures

    implicit none

interface
    !>------------------------------------------------
    !! Initialize the Obs module
    !!
    !!------------------------------------------------
    module subroutine init_obs_io(options)
        implicit none
        type(config), intent(in) :: options

    end subroutine init_obs_io

    !>------------------------------------------------
    !! Read in the Obs data
    !!
    !! Uses variable names and filenames defined in the options data structure
    !!
    !! Loops through all input variables, then reads lat, lon, and time variables
    !!
    !!------------------------------------------------
    module function read_obs(options) result(obs_data)
        implicit none
        class(obs_config), intent(in) :: options
        type(obs) :: obs_data

    end function read_obs

end interface
end module obs_mod
