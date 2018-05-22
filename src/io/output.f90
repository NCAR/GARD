submodule(output_mod) output_implementation
    use io_routines, only : io_write, io_add_attribute
    use string, only : str

    implicit none

    interface shift_z_dim
        module procedure shift_z_dim_3d, shift_z_dim_4d
    end interface shift_z_dim
contains

    subroutine add_coordinates(filename, dimnames, lat, lon, times, year_zero)
        implicit none
        character(len=MAXFILELENGTH),       intent(in) :: filename
        character(len=MAXFILELENGTH),       intent(in) :: dimnames(3)
        real,             dimension(:,:),   intent(in) :: lat, lon
        double precision, dimension(:),     intent(in) :: times
        integer,                            intent(in) :: year_zero

        call io_write(filename, "time", times, dimnames(3))
        call io_add_attribute(filename, "standard_name", "time", "time")
        call io_add_attribute(filename, "units", "days since Jan. 1, "//trim(str(year_zero)), "time")
        call io_add_attribute(filename, "axis", "T", "time")

        call io_write(filename, "lat", lat, dimnames(1:2))
        call io_add_attribute(filename, "standard_name", "latitude", "lat")
        call io_add_attribute(filename, "units", "degrees_north", "lat")
        call io_add_attribute(filename, "axis", "Y", "lat")

        call io_write(filename, "lon", lon, dimnames(1:2))
        call io_add_attribute(filename, "standard_name", "longitude", "lon")
        call io_add_attribute(filename, "units", "degrees_east", "lon")
        call io_add_attribute(filename, "axis", "X", "lon")

    end subroutine

    module subroutine write_output(output, options)
        implicit none
        type(config),  intent(in)   :: options
        type(results), intent(in)   :: output

        character(len=MAXFILELENGTH)            :: filename
        character(len=MAXFILELENGTH)            :: dimnames(3)
        real, dimension(:,:,:),     allocatable :: output_data
        real, dimension(:,:,:,:),   allocatable :: output_data_4d
        integer :: nvars, i, nx, ny, nt, nv
        integer :: Mem_Error
        double precision, dimension(:), allocatable :: times

        nvars = size(output%variables)

        nt = size(output%variables(1)%data, 1)
        nx = size(output%variables(1)%data, 2)
        ny = size(output%variables(1)%data, 3)
        nv = size(output%variables(1)%predictors, 4)

        allocate(times(nt))
        do i=1,nt
            times(i) = output%times(i)%mjd()
        enddo

        allocate( output_data(nx,ny,nt), stat=Mem_Error )
        if (Mem_Error /= 0) call memory_error(Mem_Error, "output_data", [nx,ny,nt])

        if (options%debug) then
            allocate( output_data_4d(nx,ny,nt,nv), stat=Mem_Error )
            if (Mem_Error /= 0) call memory_error(Mem_Error, "output_data_4d", [nx,ny,nt,nv])
        endif

        dimnames = [character(len=4) :: "x", "y", "time"]

        do i=1,nvars

            filename = trim(options%output_file)//trim(output%variables(i)%name)//".nc"
            call shift_z_dim(output%variables(i)%data, output_data)
            call io_write(filename, trim(output%variables(i)%name), output_data, dimnames)
            call io_add_attribute(filename, "description", "predicted mean value", trim(output%variables(i)%name))
            call io_add_attribute(filename, "coordinates", "lon lat time", trim(output%variables(i)%name))
            call add_coordinates(filename, dimnames, output%lat, output%lon, times, output%times(1)%year_zero)

            filename = trim(options%output_file)//trim(output%variables(i)%name)//"_errors.nc"
            call shift_z_dim(output%variables(i)%errors, output_data)
            call io_write(filename, trim(output%variables(i)%name)//"_error", output_data, dimnames)
            call io_add_attribute(filename, "description", "predicted error term", trim(output%variables(i)%name)//"_error")
            call add_coordinates(filename, dimnames, output%lat, output%lon, times, output%times(1)%year_zero)

            if (options%logistic_threshold/=kFILL_VALUE) then
                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_logistic.nc"
                call shift_z_dim(output%variables(i)%logistic, output_data)
                call io_write(filename, trim(output%variables(i)%name)//"_exceedence_probability", output_data, dimnames)
                call io_add_attribute(filename, "description", "probability of threshold exceedence", trim(output%variables(i)%name)//"_exceedence_probability")
                call io_add_attribute(filename, "threshold", options%logistic_threshold, trim(output%variables(i)%name)//"_exceedence_probability")

                call add_coordinates(filename, dimnames, output%lat, output%lon, times, output%times(1)%year_zero)
            endif

            if (options%debug) then
                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_predictors.nc"
                call shift_z_dim(output%variables(i)%predictors, output_data_4d)
                call io_write(filename, "predictors", output_data_4d)

                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_obs.nc"
                call shift_z_dim(output%variables(i)%obs, output_data)
                call io_write(filename, "obs", output_data)

                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_training.nc"
                call io_write(filename, "training", output%variables(i)%training)
                ! this quickly takes too much memory, so we won't bother swapping dimensions
                ! call shift_z_dim(output%variables(i)%training, output_data_4d)
                ! call io_write(filename, "data", output_data_4d)

                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_coef.nc"
                ! note other 4d vars are nt, nx, ny, nv, but coeff is nv, nt, nx, ny
                ! call shift_z_dim(output%variables(i)%coefficients, output_data_4d)
                ! for now skip the shift and just output as is.
                call io_write(filename, "coefficients", output%variables(i)%coefficients)
            else if ((options%write_coefficients).and.(options%pure_regression)) then
                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_coef.nc"
                call io_write(filename, "coefficients", output%variables(i)%coefficients)
            endif

        enddo

    end subroutine write_output

    ! for whatever reason reshape(data, [nx, ny, nt], order=[2,3,1]) doesn't work!
    subroutine shift_z_dim_3d(input, output)
        implicit none
        real, intent(in), dimension(:,:,:) :: input
        real, intent(inout), dimension(:,:,:), allocatable :: output
        integer :: i,j,k
        integer :: nx,ny,nt
        integer :: Mem_Error

        nt = size(input, 1)
        nx = size(input, 2)
        ny = size(input, 3)

        if ((size(output,1) /= nx).or.(size(output,2) /= ny).or.(size(output,3) /= nt)) then
            deallocate(output)
            allocate(output(nx,ny,nt), stat=Mem_Error)
            if (Mem_Error /= 0) call memory_error(Mem_Error, "output (shift_z)", [nx,ny,nt])
        endif

        do j=1,ny
            do i=1,nx
                output(i,j,:) = input(:,i,j)
            enddo
        enddo

    end subroutine shift_z_dim_3d

    subroutine shift_z_dim_4d(input, output)
        implicit none
        real, intent(in), dimension(:,:,:,:) :: input
        real, intent(inout), dimension(:,:,:,:), allocatable :: output
        integer :: i,j,k
        integer :: nx,ny,nt, nv
        integer :: Mem_Error

        nt = size(input, 1)
        nx = size(input, 2)
        ny = size(input, 3)
        nv = size(input, 4)

        if ((size(output,1) /= nx).or.(size(output,2) /= ny).or.(size(output,3) /= nt).or.(size(output,4) /= nv)) then
            deallocate(output)
            allocate(output(nx,ny,nt,nv), stat=Mem_Error)
            if (Mem_Error /= 0) call memory_error(Mem_Error, "output_4d (shift_z)", [nx,ny,nt,nv])
        endif

        do k=1,nv
            do j=1,ny
                do i=1,nx
                    output(i,j,1:nt,k) = input(:,i,j,k)
                enddo
            enddo
        enddo

    end subroutine shift_z_dim_4d

    subroutine memory_error(error, variable_name, dims)
        implicit none
        integer,          intent(in)                :: error
        character(len=*), intent(in)                :: variable_name
        integer,          intent(in), dimension(:)  :: dims

        write(*,*) "Error allocating memory for variable: ", trim(variable_name)
        write(*,*) "  ERROR        = ", error
        write(*,*) "  Dimensions   = ", dims

        stop "MEMORY ALLOCATION ERROR"

    end subroutine memory_error


end submodule output_implementation
