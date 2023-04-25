module observer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module: observer
!
! This module is used to "observe" the state vector of a 
! system of equations by writing that state in some format to a file
! for later viewing or plotting by the user.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  private

  public :: observer_init, observer_write, observer_finalize, numFrames

  ! TODO Does this count as interface modification?
  integer :: ncid          ! netCDF file handle
  integer :: dimIDX        ! dimension ID for X
  integer :: dimIDY        ! dimension ID for Y
  integer :: dimIDState    ! dimension ID for S
  integer :: dimIDFrame    ! dimension ID for time
  integer :: varIDState    ! variable ID for height
  integer :: stat          ! netCDF status
  integer :: frame         ! frame count
  integer :: numFrames    ! number of frames


  contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Name: observer_init
  !
  ! Description: Initializes the observer module by, e.g., opening 
  !   files for later writing. This routine must be called before the 
  !   first call to observer_write().
  !
  ! Input: none
  !
  ! Output: none
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine observer_init()
#ifdef USE_MPI
    use mpi
#endif
    use netcdf
    use equations

    implicit none
   
    !
    ! Code...
    !
    frame = 1

#ifdef USE_MPI
    write( *, * ) "Creating Parellel IO file"
    stat = nf90_create( "states-mpi.nc", IOR(NF90_NETCDF4, NF90_MPIIO), ncid, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL )
#else

    write( *, * ) "Creating Single process file"
    stat = nf90_create( "states-single.nc", NF90_CLOBBER, ncid )

#endif

    if ( stat /= NF90_NOERR ) then
      ! We had an error
      write( 0, * )  __FILE__, ":", __LINE__, " : ", trim( nf90_strerror( stat ) )
      stop 1
    end if
    
    ! Define a dimension [nCells] with size X
    stat = nf90_def_dim( ncid, "frame", numFrames + 1, dimIDFrame )
    if ( stat /= NF90_NOERR ) then
      ! We had an error
      write( 0, * )  __FILE__, ":", __LINE__, " : ", trim( nf90_strerror( stat ) )
      stop 1
    end if

    ! Define a dimension [nCells] with size X
    stat = nf90_def_dim( ncid, "xdim", totalX, dimIDX )
    if ( stat /= NF90_NOERR ) then
      ! We had an error
      write( 0, * )  __FILE__, ":", __LINE__, " : ", trim( nf90_strerror( stat ) )
      stop 1
    end if
    
    stat = nf90_def_dim( ncid, "ydim", totalY, dimIDY )
    if ( stat /= NF90_NOERR ) then
      ! We had an error
      write( 0, * )  __FILE__, ":", __LINE__, " : ", trim( nf90_strerror( stat ) )
      stop 1
    end if

    stat = nf90_def_dim( ncid, "sdim", get_system_size(), dimIDState )
    if ( stat /= NF90_NOERR ) then
      ! We had an error
      write( 0, * )  __FILE__, ":", __LINE__, " : ", trim( nf90_strerror( stat ) )
      stop 1
    end if

    ! Define a real-valued variable [state] dimensioned by nCells
    stat = nf90_def_var( ncid, "state", NF90_FLOAT, [dimIDX, dimIDY,dimIDState,dimIDFrame], varIDState )
    if ( stat /= NF90_NOERR ) then
      ! We had an error
      write( 0, * )  __FILE__, ":", __LINE__, " : ", trim( nf90_strerror( stat ) )
      stop 1
    end if

    ! No more dimensions or variables to define so exit define mode
    stat = nf90_enddef( ncid )

  end subroutine observer_init


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Name: observer_write
  !
  ! Description: Formats and writes the contents of the state vector s
  !   to a file.
  !
  ! Input: s -- the state vector
  !
  ! Output: none
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine observer_write(s)
    use netcdf
    use equations

    implicit none
    real, dimension(:), intent(inout), target :: s

    real, dimension(:,:,:,:), pointer :: sRemap
    integer                           :: writeX, writeY, i, j

    if ( frame > ( numFrames + 1 ) ) then
      write( 0, * )  "No space left to write into CDF file, current frames : ", frame
      stop 1
    end if
    sRemap( 1:nx, 1:ny, 1:get_system_size(), 1:1 ) => s(:)
    writeX = 1
    writeY = 1

#ifdef USE_MPI

    ! Define where we shall write to based on composite rank
    writeX = aRank * nx + 1
    writeY = bRank * ny + 1

    ! write( *, "(a,i,a,i,a,i,a,i,a,i,a,i,a)" ) "Writing out ", aRank, ", ", bRank, " from [", writeX, ":", writeX + nx, ", ", writeY, ":", writeY + ny, "]"
    ! do i = 1,nx
    !   do j = 1,ny
    !     write( *, "(a,1f4.2)", advance="no" ) " ", sRemap(i,j,3,1)
    !   end do
    !   write( *, * ) " "
    ! end do

#endif

    ! do i = 1,nx
    !   do j = 1,ny
    !     write( *, "(a,1f4.2)", advance="no" ) " ", sRemap(i,j,3,1)
    !   end do
    !   write( *, * ) " "
    ! end do
    ! write( *, * ) "Writing into frame ", frame
    stat = nf90_put_var( ncid, varIDState, sRemap, start=[ writeX, writeY, 1, frame ] )
    if ( stat /= NF90_NOERR ) then
      ! We had an error
      write( 0, * )  __FILE__, ":", __LINE__, " : ", trim( nf90_strerror( stat ) )
      stop 1
    end if
        

    ! Move forward in frames
    frame = frame + 1
    
  end subroutine observer_write


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Name: observer_finalize
  !
  ! Description: Finalizes the observer module by, e.g., closing any
  !   files that were opened by the module. This routine must be called 
  !   only once after all calls to observer_write() have been made.
  !
  ! Input: none
  !
  ! Output: none
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine observer_finalize()
    use netcdf

    implicit none

    !
    ! Code...
    !
    stat = nf90_close( ncid )

  end subroutine observer_finalize

end module observer
