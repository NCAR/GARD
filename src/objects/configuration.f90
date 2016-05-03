module config_mod
    
    use data_structures
    
    implicit none
    private
    
    character(len=MAXFILELENGTH), parameter :: DEFAULT_OPTIONS_FILENAME = "downscale_options.txt"
    
    interface config
        module procedure new_config
    end interface config
    
    type, public, extends(options_type):: config
        private
        character(len=MAXFILELENGTH) :: options_filename
      contains
        procedure, public :: as_string  => as_string
    end type config
    
contains
    function new_config() result(options)
        implicit none
        type(config) :: options
        character(len=MAXFILELENGTH) :: filename
        
        CALL get_command_argument(1, filename)
        
        if (len_trim(filename) > 0) then
            options%options_filename = filename
        else
            options%options_filename = DEFAULT_OPTIONS_FILENAME
        endif
        
        options%gcm%name="gcm"
        options%obs%name="obs"
        options%training%name="train"
        
    end function new_config

    function as_string(this)
        implicit none
        class(config), intent(in) :: this
        character(len=MAXSTRINGLENGTH) :: as_string
        
        as_string = this%options_filename
        
    end function as_string
    

end module config_mod
