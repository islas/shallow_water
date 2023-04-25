! Shallow water driver
module integration
  private

  public :: forwardEuler, rk3

  contains

  subroutine forwardEuler( iter, timestep, state, tendency )

    use equations
    
    implicit none

    integer, intent(in)    :: iter
    real,    intent(in)    :: timestep
    real,    intent(inout), dimension(:)  :: state
    real,    intent(inout), dimension(:)  :: tendency

    real    :: time

    ! Propagate time, multiply to avoid FP round off errors, may need to seed with initial time
    time = timestep * iter
    ! Implications of array assignment being done like this?
    call f( time, state, tendency )

    ! Forward propagate
    state(:) = state(:) + tendency(:) * timestep

  end subroutine forwardEuler



  subroutine rk3( iter, timestep, state, k1, k2, k3 )
    use equations
#ifdef USE_MPI
    use mpi
#endif
    implicit none

    integer, intent(in)    :: iter
    real,    intent(in)    :: timestep
    
    real, intent(inout), dimension(:)  :: state
    real, intent(out),   dimension(:)  :: k1
    real, intent(out),   dimension(:)  :: k2
    real, intent(out),   dimension(:)  :: k3

    real    :: time
#ifdef USE_MPI
    integer :: mpiErr
#endif

    ! Propagate time, multiply to avoid FP round off errors, may need to seed with initial time
    time = timestep * iter

    ! MPI could be used here to have 3 processes calculating each one
    ! This is k1
    call f( time, state(:), k1 )

    ! This is k2
    call f( time + timestep * 0.5, state(:) + k1(:) * timestep * 0.5, k2 )

    ! This is k3
    call f( time + timestep, state(:) - timestep * k1(:) + 2 * timestep * k2(:), k3 )

    ! New state
    state(:) = state(:) + timestep / 6.0 * ( k1(:) + 4.0 * k2(:) + k3(:) )
      
  end subroutine rk3
end module integration


program solver
  ! Load modules
  use observer
  use equations
  use integration
#ifdef USE_MPI
  use mpi
#endif


  implicit none

  real    :: timestep
  real, dimension(:), pointer  :: stateVector    => null()
  ! k1 will be default tendency for Forward Euler
  real, dimension(:), allocatable  :: k1Scratch, k2Scratch, k3Scratch
  
  integer :: i, iterations, writeInterval
  
#ifdef USE_MPI
  integer             :: mpiErr, processRank
#endif

  
  timestep      = 0.05
  iterations    = 200
  writeInterval = 2
  numFrames     = iterations / writeInterval

#ifdef USE_MPI 
  write( *, * ) "Running MPI"
  call mpi_init( mpiErr )
  call mpi_comm_rank( MPI_COMM_WORLD, processRank, mpiErr )
#endif

  ! Init things - both call this to have the correct sizing
  call set_initial_state( stateVector )

#ifdef USE_MPI
  ! Sync all processes
  call mpi_barrier( MPI_COMM_WORLD, mpiErr )
#endif

  call observer_init()
  call observer_write( stateVector )

  ! Allocate interpolation scratch
  allocate( k1Scratch(size(stateVector) ) )
  allocate( k2Scratch(size(stateVector) ) )
  allocate( k3Scratch(size(stateVector) ) )

  do i = 1, iterations
    ! write( *, * ) "Iteration : ", i
#ifdef USE_MPI
    if ( processRank == 0 ) then
#endif
    if ( mod( float( i ) * timestep, 1.0 ) == 0.0 ) then
      write( *, * ) "Time: ", float( i ) * timestep
    endif
#ifdef USE_MPI
    endif
#endif

    call rk3( i, timestep, stateVector, k1Scratch, k2Scratch, k3Scratch )
    ! call forwardEuler( i, timestep, stateVector, k1Scratch )
    
    if ( mod( i, writeInterval ) == 0.0 ) then
      call observer_write( stateVector )
    endif
    
  end do

  call observer_finalize()
  call equations_finalize( stateVector )

  deallocate( k1Scratch )
  deallocate( k2Scratch )
  deallocate( k3Scratch )

#ifdef USE_MPI 
  call mpi_finalize( mpiErr )
#endif

end program solver