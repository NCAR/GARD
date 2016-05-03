module obs_mod

    use data_structures
    use geo
    
    implicit none
    
    interface observations
        module procedure :: new_obs
    end interface observations
    
    type, public :: observations
        private
        character(len=MAXFILELENGTH) :: name
      contains
        procedure, public  :: init        => init
        procedure, public  :: as_string   => as_string
        procedure, public  :: set_name    => set_name
    end type observations

contains
    function new_obs()
        implicit none
        type(observations) :: new_obs
        
        call new_obs%set_name( "default-obs-name" )
        
    end function new_obs
    
    subroutine init(this, options)
        implicit none
        class(observations), intent(inout) :: this
        class(options_type), intent(in)  :: options
        
        call this%set_name( options%obs%name )
        
    end subroutine init

    
    function as_string(this) result(output)
        implicit none
        class(observations), intent(in) :: this
        character(len=MAXVARLENGTH) :: output
        
        output = this%name
        
    end function as_string

    subroutine set_name(this, name)
        implicit none
        class(observations), intent(inout) :: this
        character(len=*) :: name
        
        this%name = name
        
    end subroutine set_name
    
end module obs_mod
