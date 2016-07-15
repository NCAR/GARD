module output_mod
    use io_routines, only : io_write
    use data_structures
    use string, only : str
    
    implicit none

    interface shift_z_dim
        module procedure shift_z_dim_3d, shift_z_dim_4d
    end interface shift_z_dim
contains
    subroutine write_output(output, options)
        implicit none
        type(config),  intent(in)   :: options
        type(results), intent(in)   :: output
        character(len=MAXFILELENGTH):: filename
        real, dimension(:,:,:), allocatable :: output_data
        real, dimension(:,:,:,:), allocatable :: output_data_4d
        integer :: nvars, i, nx, ny, nt, nv
        
        nvars = size(output%variables)
        
        nt = size(output%variables(1)%data, 1)
        nx = size(output%variables(1)%data, 2)
        ny = size(output%variables(1)%data, 3)
        nv = size(output%variables(1)%predictors, 4)

        allocate( output_data(nx,ny,nt) )
        allocate( output_data_4d(nx,ny,nt,nv) )

        do i=1,nvars
            
            filename = trim(options%output_file)//trim(output%variables(i)%name)//".nc"
            call shift_z_dim(output%variables(i)%data, output_data)
            call io_write(filename, "data", output_data)
            
            filename = trim(options%output_file)//trim(output%variables(i)%name)//"_errors.nc"
            call shift_z_dim(output%variables(i)%errors, output_data)
            call io_write(filename, "data", output_data)

            if (options%logistic_threshold/=kFILL_VALUE) then
                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_logistic.nc"
                call shift_z_dim(output%variables(i)%logistic, output_data)
                call io_write(filename, "data", output_data)
            endif

            if (options%debug) then
                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_predictors.nc"
                call shift_z_dim(output%variables(i)%predictors, output_data_4d)
                call io_write(filename, "data", output_data_4d)

                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_obs.nc"
                call shift_z_dim(output%variables(i)%obs, output_data)
                call io_write(filename, "data", output_data)

                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_training.nc"
                call shift_z_dim(output%variables(i)%training, output_data_4d)
                call io_write(filename, "data", output_data_4d)
                
                filename = trim(options%output_file)//trim(output%variables(i)%name)//"_coef.nc"
                call shift_z_dim(output%variables(i)%coefficients, output_data_4d)
                call io_write(filename, "data", output_data_4d)
            endif
        enddo
        
    end subroutine write_output

    ! for whatever reason reshape(data, [nx, ny, nt], order=[2,3,1]) doesn't work!
    subroutine shift_z_dim_3d(input, output)
        implicit none
        real, intent(in), dimension(:,:,:) :: input
        real, intent(inout), dimension(:,:,:) :: output
        integer :: i,j,k
        integer :: nx,ny,nt

        nt = size(input, 1)
        nx = size(input, 2)
        ny = size(input, 3)

        do i=1,nx
            do j=1,ny
                output(i,j,:) = input(:,i,j)
            enddo
        enddo

    end subroutine shift_z_dim_3d

    subroutine shift_z_dim_4d(input, output)
        implicit none
        real, intent(in), dimension(:,:,:,:) :: input
        real, intent(inout), dimension(:,:,:,:) :: output
        integer :: i,j,k
        integer :: nx,ny,nt, nv

        nt = size(input, 1)
        nx = size(input, 2)
        ny = size(input, 3)
        nv = size(input, 4)

        do k=1,nv
            do i=1,nx
                do j=1,ny
                    output(i,j,:,k) = input(:,i,j,k)
                enddo
            enddo
        enddo

    end subroutine shift_z_dim_4d

end module output_mod
