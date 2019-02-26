module atm_mod
    use data_structures

    implicit none

interface
    module function read_atm(options) result(atm_data)
        implicit none
        class(atm_config) :: options
        type(atm) :: atm_data

    end function read_atm

end interface

end module atm_mod
