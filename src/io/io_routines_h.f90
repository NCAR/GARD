!>------------------------------------------------------------
!!  Basic file input/output routines
!!
!!  @details
!!  Primary use is io_read2d/3d
!!  io_write* routines are more used for debugging
!!  model output is performed in the output module
!!
!!  Generic interfaces are supplied for io_read and io_write
!!  but most code still uses the explicit read/write 2d/3d/etc.
!!  this keeps the code a little more obvious (at least 2D vs 3D)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module io_routines

    implicit none
    ! maximum number of dimensions for a netCDF file
    integer,parameter::io_maxDims=10

    !>------------------------------------------------------------
    !! Generic interface to the netcdf read routines
    !!------------------------------------------------------------
    interface io_read
        module procedure io_read2d, io_read3d, io_read6d, io_read2di, io_read1d, io_read4d, io_read1dd
    end interface

    !>------------------------------------------------------------
    !! Generic interface to the netcdf write routines
    !!------------------------------------------------------------
    interface io_write
        module procedure io_write6d, io_write4d, io_write3d, io_write2d, io_write1d, io_write1dd, io_write3di,io_write4di
    end interface

    !>------------------------------------------------------------
    !! Generic interface to the netcdf read_attribute_TYPE routines
    !!------------------------------------------------------------
    interface io_read_attribute
        module procedure io_read_attribute_r, io_read_attribute_i, io_read_attribute_c
    end interface
    ! to be added as necessary
    !, io_read_attribute_d

    !>------------------------------------------------------------
    !! Generic interface to the netcdf add_attribute_TYPE routines
    !!------------------------------------------------------------
    interface io_add_attribute
        module procedure io_add_attribute_r, io_add_attribute_i, io_add_attribute_c
    end interface
    ! to be added as necessary
    !, io_add_attribute_d

!   All routines are public
interface

    !>------------------------------------------------------------
    !! Tests to see if a file exists
    !!
    !! @param filename name of file to look for
    !! @retval logical true if file exists, false if it doesn't
    !!
    !!------------------------------------------------------------
    module logical function file_exists(filename)
        character(len=*), intent(in) :: filename

    end function file_exists


    !>------------------------------------------------------------
    !! Read the dimensions of a variable in a given netcdf file
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to find the dimensions of
    !! @param[out] dims     Allocated array to store output
    !! @retval dims(:) dims[1]=ndims, dims[i+1]=length of dimension i for a given variable
    !!
    !!------------------------------------------------------------
    module subroutine io_getdims(filename,varname,dims)
        implicit none
        character(len=*), intent(in) :: filename,varname
        integer,intent(out) :: dims(:)

    end subroutine io_getdims

    !>------------------------------------------------------------
    !! Reads in a variable from a netcdf file, allocating memory in data_in for it.
    !!
    !! if extradim is provided specifies this index for any extra dimensions (dims>6)
    !!   e.g. we may only want one time slice from a 6d variable
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to read
    !! @param[out] data_in     Allocatable 6-dimensional array to store output
    !! @param   extradim    OPTIONAL: specify the position to read for any extra (e.g. time) dimension
    !! @retval data_in     Allocated 6-dimensional array with the netCDF data
    !!
    !!------------------------------------------------------------
    module subroutine io_read6d(filename,varname,data_in,extradim)
        implicit none
        ! This is the name of the data_in file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        real,intent(out),allocatable :: data_in(:,:,:,:,:,:)
        integer, intent(in),optional :: extradim

    end subroutine io_read6d

    !>------------------------------------------------------------
    !! Same as io_read6d but for 4-dimensional data
    !!
    !! Reads in a variable from a netcdf file, allocating memory in data_in for it.
    !!
    !! if extradim is provided specifies this index for any extra dimensions (dims>3)
    !!   e.g. we may only want one time slice from a 3d variable
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to read
    !! @param[out] data_in     Allocatable 4-dimensional array to store output
    !! @param   extradim    OPTIONAL: specify the position to read for any extra (e.g. time) dimension
    !! @retval data_in     Allocated 3-dimensional array with the netCDF data
    !!
    !!------------------------------------------------------------
    module subroutine io_read4d(filename,varname,data_in,extradim)
        implicit none
        ! This is the name of the data_in file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        real,intent(out),allocatable :: data_in(:,:,:,:)
        integer, intent(in),optional :: extradim

    end subroutine io_read4d


    !>------------------------------------------------------------
    !! Same as io_read6d but for 3-dimensional data
    !!
    !! Reads in a variable from a netcdf file, allocating memory in data_in for it.
    !!
    !! if extradim is provided specifies this index for any extra dimensions (dims>3)
    !!   e.g. we may only want one time slice from a 3d variable
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to read
    !! @param[out] data_in     Allocatable 3-dimensional array to store output
    !! @param   extradim    OPTIONAL: specify the position to read for any extra (e.g. time) dimension
    !! @retval data_in     Allocated 3-dimensional array with the netCDF data
    !!
    !!------------------------------------------------------------
    module subroutine io_read3d(filename,varname,data_in,extradim)
        implicit none
        ! This is the name of the data_in file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        real,intent(out),allocatable :: data_in(:,:,:)
        integer, intent(in),optional :: extradim

    end subroutine io_read3d


    !>------------------------------------------------------------
    !! Same as io_read3d but for 2-dimensional data
    !!
    !! Reads in a variable from a netcdf file, allocating memory in data_in for it.
    !!
    !! if extradim is provided specifies this index for any extra dimensions (dims>2)
    !!   e.g. we may only want one time slice from a 2d variable
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to read
    !! @param[out] data_in     Allocatable 2-dimensional array to store output
    !! @param   extradim    OPTIONAL: specify the position to read for any extra (e.g. time) dimension
    !! @retval data_in     Allocated 2-dimensional array with the netCDF data
    !!
    !!------------------------------------------------------------
    module subroutine io_read2d(filename,varname,data_in,extradim)
        implicit none
        ! This is the name of the data_in file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        real,intent(out),allocatable :: data_in(:,:)
        integer, intent(in),optional :: extradim

    end subroutine io_read2d

    !>------------------------------------------------------------
    !! Same as io_read2d but for integer data
    !!
    !! Reads in a variable from a netcdf file, allocating memory in data_in for it.
    !!
    !! if extradim is provided specifies this index for any extra dimensions (dims>2)
    !!   e.g. we may only want one time slice from a 2d variable
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to read
    !! @param[out] data_in     Allocatable 2-dimensional array to store output
    !! @param   extradim    OPTIONAL: specify the position to read for any extra (e.g. time) dimension
    !! @retval data_in     Allocated 2-dimensional array with the netCDF data
    !!
    !!------------------------------------------------------------
    module subroutine io_read2di(filename,varname,data_in,extradim)
        implicit none
        ! This is the name of the data_in file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        integer,intent(out),allocatable :: data_in(:,:)
        integer, intent(in),optional :: extradim

    end subroutine io_read2di

    !>------------------------------------------------------------
    !! Same as io_read3d but for 1-dimensional data
    !!
    !! Reads in a variable from a netcdf file, allocating memory in data_in for it.
    !!
    !! if extradim is provided specifies this index for any extra dimensions (dims>1)
    !!   e.g. we may only want one time slice from a 1d variable
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to read
    !! @param[out] data_in     Allocatable 1-dimensional array to store output
    !! @param   extradim    OPTIONAL: specify the position to read for any extra (e.g. time) dimension
    !! @retval data_in     Allocated 1-dimensional array with the netCDF data
    !!
    !!------------------------------------------------------------
    module subroutine io_read1d(filename,varname,data_in,extradim)
        implicit none
        ! This is the name of the data_in file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        real,intent(out),allocatable :: data_in(:)
        integer, intent(in),optional :: extradim
    end subroutine io_read1d

    !>------------------------------------------------------------
    !! Same as io_read1d but for double precision data
    !!
    !! Reads in a variable from a netcdf file, allocating memory in data_in for it.
    !!
    !! if extradim is provided specifies this index for any extra dimensions (dims>1)
    !!   e.g. we may only want one time slice from a 1d variable
    !!
    !! @param   filename    Name of NetCDF file to look at
    !! @param   varname     Name of the NetCDF variable to read
    !! @param[out] data_in     Allocatable 1-dimensional array to store output
    !! @param   extradim    OPTIONAL: specify the position to read for any extra (e.g. time) dimension
    !! @retval data_in     Allocated 1-dimensional array with the netCDF data
    !!
    !!------------------------------------------------------------
    module subroutine io_read1dd(filename,varname,data_in,extradim)
        implicit none
        ! This is the name of the data_in file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        double precision,intent(out),allocatable :: data_in(:)
        integer, intent(in),optional :: extradim
    end subroutine io_read1dd



    !>------------------------------------------------------------
    !! Write a 6-dimensional variable to a netcdf file
    !!
    !! Create a netcdf file:filename with a variable:varname and write data_out to it
    !!
    !! @param   filename    Name of NetCDF file to write/create
    !! @param   varname     Name of the NetCDF variable to write
    !! @param   data_out    6-dimensional array to write to the file
    !!
    !!------------------------------------------------------------
    module subroutine io_write6d(filename,varname,data_out, dimnames)
        implicit none
        ! This is the name of the file and variable we will write.
        character(len=*), intent(in) :: filename, varname
        real,intent(in) :: data_out(:,:,:,:,:,:)
        character(len=*), optional, dimension(6), intent(in) :: dimnames

    end subroutine io_write6d

    !>------------------------------------------------------------
    !! Same as io_write6d but for 4-dimensional data
    !!
    !! Write a 4-dimensional variable to a netcdf file
    !!
    !! Create a netcdf file:filename with a variable:varname and write data_out to it
    !!
    !! @param   filename    Name of NetCDF file to write/create
    !! @param   varname     Name of the NetCDF variable to write
    !! @param   data_out    4-dimensional array to write to the file
    !!
    !!------------------------------------------------------------
    module subroutine io_write4d(filename,varname,data_out, dimnames)
        implicit none
        ! This is the name of the file and variable we will write.
        character(len=*), intent(in) :: filename, varname
        real,intent(in) :: data_out(:,:,:,:)
        character(len=*), optional, dimension(4), intent(in) :: dimnames

    end subroutine io_write4d

    !>------------------------------------------------------------
    !! Same as io_write4d but for integer data
    !!
    !! Write a 4-dimensional variable to a netcdf file
    !!
    !! Create a netcdf file:filename with a variable:varname and write data_out to it
    !!
    !! @param   filename    Name of NetCDF file to write/create
    !! @param   varname     Name of the NetCDF variable to write
    !! @param   data_out    4-dimensional array to write to the file
    !!
    !!------------------------------------------------------------
    module subroutine io_write4di(filename,varname,data_out)
        implicit none
        ! This is the name of the file and variable we will write.
        character(len=*), intent(in) :: filename, varname
        integer,intent(in) :: data_out(:,:,:,:)

    end subroutine io_write4di


    !>------------------------------------------------------------
    !! Same as io_write6d but for 3-dimensional data
    !!
    !! Write a 3-dimensional variable to a netcdf file
    !!
    !! Create a netcdf file:filename with a variable:varname and write data_out to it
    !!
    !! @param   filename    Name of NetCDF file to write/create
    !! @param   varname     Name of the NetCDF variable to write
    !! @param   data_out    3-dimensional array to write to the file
    !!
    !!------------------------------------------------------------
    module subroutine io_write3d(filename,varname,data_out, dimnames)
        implicit none
        ! This is the name of the file and variable we will write.
        character(len=*), intent(in) :: filename, varname
        real,intent(in) :: data_out(:,:,:)
        character(len=*), optional, dimension(3), intent(in) :: dimnames

    end subroutine io_write3d

    !>------------------------------------------------------------
    !! Same as io_write3d but for integer arrays
    !!
    !! Write a 3-dimensional variable to a netcdf file
    !!
    !! Create a netcdf file:filename with a variable:varname and write data_out to it
    !!
    !! @param   filename    Name of NetCDF file to write/create
    !! @param   varname     Name of the NetCDF variable to write
    !! @param   data_out    3-dimensional array to write to the file
    !!
    !!------------------------------------------------------------
    module subroutine io_write3di(filename,varname,data_out)
        implicit none
        ! This is the name of the data file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        integer,intent(in) :: data_out(:,:,:)

    end subroutine io_write3di

    !>------------------------------------------------------------
    !! Same as io_write3d but for 2-dimensional arrays
    !!
    !! Write a 2-dimensional variable to a netcdf file
    !!
    !! Create a netcdf file:filename with a variable:varname and write data_out to it
    !!
    !! @param   filename    Name of NetCDF file to write/create
    !! @param   varname     Name of the NetCDF variable to write
    !! @param   data_out    2-dimensional array to write to the file
    !!
    !!------------------------------------------------------------
    module subroutine io_write2d(filename,varname,data_out, dimnames)
        implicit none
        ! This is the name of the data file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        real,intent(in) :: data_out(:,:)
        character(len=*), optional, dimension(2), intent(in) :: dimnames
        
    end subroutine io_write2d

    module subroutine io_write1dd(filename,varname,data_out,dimname)
        implicit none
        ! This is the name of the data file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        double precision,intent(in) :: data_out(:)
        character(len=*), intent(in), optional :: dimname

    end subroutine io_write1dd

    module subroutine io_write1d(filename,varname,data_out)
        implicit none
        ! This is the name of the data file and variable we will read.
        character(len=*), intent(in) :: filename, varname
        real,intent(in) :: data_out(:)

    end subroutine io_write1d


    !>------------------------------------------------------------
    !! Read a real type attribute from a named file from an optional variable
    !!
    !! If a variable name is given reads the named attribute of that variable
    !! otherwise the named attribute is assumed to be a global attribute
    !!
    !! @param   filename    netcdf file to read the attribute from
    !! @param   att_name    name of attribute to read
    !! @param   att_value   output value to be returned (real*4)
    !! @param   var_name    OPTIONAL name of variable to read attribute from
    !!
    !!------------------------------------------------------------
    module subroutine io_read_attribute_r(filename, att_name, att_value, var_name, error)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: att_name
        real*4, intent(out) :: att_value
        character(len=*), intent(in), optional :: var_name
        integer,          intent(out),optional :: error

    end subroutine  io_read_attribute_r

    !>------------------------------------------------------------
    !! Read a integer type attribute from a named file from an optional variable
    !!
    !! If a variable name is given reads the named attribute of that variable
    !! otherwise the named attribute is assumed to be a global attribute
    !!
    !! @param   filename    netcdf file to read the attribute from
    !! @param   att_name    name of attribute to read
    !! @param   att_value   output value to be returned (integer)
    !! @param   var_name    OPTIONAL name of variable to read attribute from
    !!
    !!------------------------------------------------------------
    module subroutine io_read_attribute_i(filename, att_name, att_value, var_name, error)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: att_name
        integer, intent(out) :: att_value
        character(len=*), intent(in), optional :: var_name
        integer,          intent(out),optional :: error

    end subroutine  io_read_attribute_i

    !>------------------------------------------------------------
    !! Read a character type attribute from a named file from an optional variable
    !!
    !! If a variable name is given reads the named attribute of that variable
    !! otherwise the named attribute is assumed to be a global attribute
    !!
    !! @param   filename    netcdf file to read the attribute from
    !! @param   att_name    name of attribute to read
    !! @param   att_value   output value to be returned (character)
    !! @param   var_name    OPTIONAL name of variable to read attribute from
    !!
    !!------------------------------------------------------------
    module subroutine io_read_attribute_c(filename, att_name, att_value, var_name, error)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: att_name
        character(len=*), intent(out) :: att_value
        character(len=*), intent(in), optional :: var_name
        integer,          intent(out),optional :: error

    end subroutine  io_read_attribute_c


    !>------------------------------------------------------------
    !! Write a real type attribute to a named file for an optional variable
    !!
    !! If a variable name is given writes the named attribute to that variable
    !! otherwise the named attribute is assumed to be a global attribute
    !!
    !! @param   filename    netcdf file to write the attribute to
    !! @param   att_name    name of attribute to write
    !! @param   att_value   output value to be written (real*4)
    !! @param   var_name    OPTIONAL name of variable to write attribute to
    !!
    !!------------------------------------------------------------
    module subroutine io_add_attribute_r(filename, att_name, att_value, varname)
        implicit none
        character(len=*), intent(in)           :: filename
        character(len=*), intent(in)           :: att_name
        real*4,           intent(in)           :: att_value
        character(len=*), intent(in), optional :: varname

    end subroutine io_add_attribute_r


    !>------------------------------------------------------------
    !! Write an integer type attribute to a named file for an optional variable
    !!
    !! If a variable name is given writes the named attribute to that variable
    !! otherwise the named attribute is assumed to be a global attribute
    !!
    !! @param   filename    netcdf file to write the attribute to
    !! @param   att_name    name of attribute to write
    !! @param   att_value   output value to be written (integer)
    !! @param   var_name    OPTIONAL name of variable to write attribute to
    !!
    !!------------------------------------------------------------
    module subroutine io_add_attribute_i(filename, att_name, att_value, varname)
        implicit none
        character(len=*), intent(in)           :: filename
        character(len=*), intent(in)           :: att_name
        integer,          intent(in)           :: att_value
        character(len=*), intent(in), optional :: varname

    end subroutine io_add_attribute_i


    !>------------------------------------------------------------
    !! Write an character type attribute to a named file for an optional variable
    !!
    !! If a variable name is given writes the named attribute to that variable
    !! otherwise the named attribute is assumed to be a global attribute
    !!
    !! @param   filename    netcdf file to write the attribute to
    !! @param   att_name    name of attribute to write
    !! @param   att_value   output value to be written (character)
    !! @param   var_name    OPTIONAL name of variable to write attribute to
    !!
    !!------------------------------------------------------------
    module subroutine io_add_attribute_c(filename, att_name, att_value, varname)
        implicit none
        character(len=*), intent(in)           :: filename
        character(len=*), intent(in)           :: att_name
        character(len=*), intent(in)           :: att_value
        character(len=*), intent(in), optional :: varname

    end subroutine io_add_attribute_c


    !>------------------------------------------------------------
    !! Find an available file unit number.
    !!
    !! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
    !! The UNIT value is returned by the function, and also by the optional
    !! argument. This allows the function to be used directly in an OPEN
    !! statement, and optionally save the result in a local variable.
    !! If no units are available, -1 is returned.
    !! Newer versions of fortran can do this automatically, but this keeps one thing
    !! a little more backwards compatible
    !!
    !! @param[out]  unit    OPTIONAL integer to store the file logical unit number
    !! @retval      integer a file logical unit number
    !!
    !!------------------------------------------------------------
    module integer function io_newunit(unit)
        implicit none
        integer, intent(out), optional :: unit
    end function io_newunit

end interface
end module io_routines
