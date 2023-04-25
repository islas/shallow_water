module equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module: equations
!
! This module implements a system of ODEs by providing routines
! to query the size of the system's state vector, get initial 
! conditions, and return the time-derivative of the state vector given
! the current state.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private

  public :: get_system_size, set_initial_state, equations_finalize, f, nx, ny, totalX, totalY, aRank, bRank
  protected :: nx, ny, totalX, totalY, aRank, bRank

  ! Parameters
  real    :: g = 9.8065 ! m/s
  real    :: perturb = 0.5

  ! From init state
  integer :: totalX, totalY
  integer :: nx, ny
  real    :: deltaX, deltaY
  real    :: h0, fwhm
#ifdef USE_MPI
  integer :: processRank, aRank, bRank, sqrtCluster

  ! Boundary syncing
  real, dimension(:), allocatable, target :: iNextBoundary, iPrevBoundary
  real, dimension(:), allocatable, target :: jNextBoundary, jPrevBoundary
  real, dimension(:), allocatable         :: iNextBoundaryBuffer, iPrevBoundaryBuffer
  real, dimension(:), allocatable         :: jNextBoundaryBuffer, jPrevBoundaryBuffer
#endif

  contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Name: get_system_size
  !
  ! Description: Returns the size of the state vector used by the 
  !   system implemented in this module.
  !
  ! Input: none
  !
  ! Return value: the size of the system's state vector
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function get_system_size()

    implicit none

    ! Return value
    integer :: get_system_size


    !
    ! Code...
    !
    get_system_size = 3


  end function get_system_size


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Name: handleMPIErr
  !
  ! Description: Process any MPI error.
  !
  ! Input: mpiErr return from MPI subroutine call
  !
  ! Output: error message if err
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine handleMPIErr(mpiErr)
#ifdef USE_MPI
    use mpi

    integer, intent(in) :: mpiErr
    character( len=MPI_MAX_ERROR_STRING ) :: errorString
    integer                               :: mpiStringErr, stringLen

    if ( mpiErr /= MPI_SUCCESS ) then
      call mpi_error_string( mpiErr, errorString, stringLen, mpiStringErr )
      write( *, * ) "MPI Error (", mpiErr, ") : ", errorString
      call mpi_finalize( mpiStringErr )
      stop 1
    endif

#endif

  end subroutine handleMPIErr


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Name: set_initial_state
  !
  ! Description: Initializes a system state vector. Upon returning, 
  !   the elements of s have been set to contain the initial condition 
  !   of the system.
  !
  ! MPI Version: 
  !   Divide the # procs into a sqrt (for convenience, nothing about this math
  !   inherently requires a sqrt) and partition the grid into a,b indices of the 
  !   process rank into composite rank of a,b such that e.g. 9 procs results in
  !   
  !   0,0  0,1  0,2   where 0,0 is rank 0, 1,0 is rank 1, ... 1,2 is rank 7 and 2,2 is rank 8
  !   1,0  1,1  1,2   (col major to fit fortran)
  !   2,0  2,1  2,2
  !   
  !   Next, make sure that totalX and totalY are divisible by your a,b total indices
  !   and then make subblocks nx,ny such that:
  !   0,0
  !      `.
  !        1  2  .  .  .  ny
  !        2  .
  !        .     .
  !        .        .
  !        .            .
  !        nx            nx,ny
  !
  !   For each a,b block
  !
  ! Input: none
  !
  ! Output: s -- the initial condition state vector
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_initial_state(s)
#ifdef USE_MPI
    use mpi
#endif

    implicit none
    
    real, dimension(:), intent(out), pointer :: s
    real, dimension(:,:,:), pointer :: sRemap
      
    ! Iterating vars
    integer             :: i, j
    real                :: startX, startY

#ifdef USE_MPI
    integer             :: clusterSize, mpiErr
#endif

    ! Init conditions
    totalY = 1024
    totalX = 1024
        
    h0     =  4
    fwhm   =  4

    deltaX = 0.25
    deltaY = 0.25


#ifdef USE_MPI
    ! TODO maybe let this be set externally
    call mpi_comm_size( MPI_COMM_WORLD, clusterSize, mpiErr )
    call mpi_comm_rank( MPI_COMM_WORLD, processRank, mpiErr )
    sqrtCluster = int( sqrt( float( clusterSize ) ) )

    if ( mod( totalX, sqrtCluster ) /= 0 .or. &
         mod( totalY, sqrtCluster ) /= 0 ) then
      write( *, * ) "MPI # procs sqrt must be a multiple of total X and Y posts"
      call mpi_finalize( mpiErr )
      stop 1
    endif

    nx = totalX / sqrtCluster
    ny = totalY / sqrtCluster

    ! Calculate composite a,b rank
    bRank = int( processRank / sqrtCluster )
    aRank = mod( processRank,  sqrtCluster )

    ! Caclulate apparent start
    startX = ( -( totalX / 2.0 ) + aRank * nx ) * deltaX
    startY = ( -( totalY / 2.0 ) + bRank * ny ) * deltaY

    ! Gather up boundary conditions
    allocate( iNextBoundary( ny * get_system_size() ) )
    allocate( iPrevBoundary( ny * get_system_size() ) )
    allocate( jNextBoundary( nx * get_system_size() ) )
    allocate( jPrevBoundary( nx * get_system_size() ) )

    allocate( iNextBoundaryBuffer( ny * get_system_size() ) )
    allocate( iPrevBoundaryBuffer( ny * get_system_size() ) )
    allocate( jNextBoundaryBuffer( nx * get_system_size() ) )
    allocate( jPrevBoundaryBuffer( nx * get_system_size() ) )

#else
    nx = totalX
    ny = totalY

    startX = -( totalX / 2.0 ) * deltaX
    startY = -( totalY / 2.0 ) * deltaY

#endif

    ! Allocate s & meshNext

    allocate( s( nx*ny*get_system_size() ) )
  
    sRemap( 1:nx, 1:ny, 1:get_system_size() ) => s( 1:nx*ny*get_system_size() )



    ! Set initial perturbation, can parallelize this
    sRemap(:,:,1:2) = 0.0
#ifdef _OPENMP 
    write( *, * ) "Running OpenMP"
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( i, j )
#endif
    do i = 1,nx
      do j = 1,ny
        ! Keep this indented
        sRemap(i,j,3)   = h0 - perturb * exp( -( ( startX - 1 + float( i - 1 ) * deltaX )**2 + ( startY + float( j - 1 ) * deltaY )**2 ) / fwhm**2 )
        ! write( *, "(a,1f4.2)", advance="no" ) " ", sRemap(i,j,3)
      end do
      ! write( *, * ) " "
    end do

#ifdef _OPENMP 
    !$OMP END PARALLEL DO
#endif

    ! do i = 1,nx
    !   do j = 1,ny
    !     write( *, "(a,1f4.2)", advance="no" ) " ", sRemap(i,j,3)
    !   end do
    !   write( *, * ) " "
    ! end do


    ! s has been modified

  end subroutine set_initial_state


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Name: f
  !
  ! Description: This function returns a tendency vector, ds/dt, for 
  !   the system implemented in this module, given the current state of 
  !   the system.
  !
  ! Input: t -- the time at which the state, s, is valid
  !        s -- a state vector
  !        tendency -- a ds/dt for the current state s at time t
  !
  ! Return value: the time-derivative, ds/dt, for the system at s
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine f( t, s, tendency )
#ifdef USE_MPI
    use mpi
#endif

    implicit none
    
    real, intent(in) :: t
    ! tendency mesh grid of size x * y * 3 < u,v,h > for next state
    real, dimension(:), intent(out), target  :: tendency
    real, dimension(:), intent(in),  target  :: s
    
    ! Iterating vars
    integer          :: i, j
    integer          :: iNext, iPrev
    integer          :: jNext, jPrev

    ! Remaps
    real, dimension(:,:,:), pointer         :: sRemap
    real, dimension(:,:,:), pointer         :: meshRemap
    
    ! Boundary re-map to allow in-place MPI substitution of original indexing
    real, dimension(:,:,:), pointer         :: iNextPtr, iPrevPtr
    real, dimension(:,:,:), pointer         :: jNextPtr, jPrevPtr

#ifdef USE_MPI
    
    integer                    :: mpiErr, bufferPos, packedBytes
    integer                    :: iNextReq, iPrevReq
    integer                    :: jNextReq, jPrevReq
    integer                    :: iNextReqMPIStatus( MPI_STATUS_SIZE ), iPrevReqMPIStatus( MPI_STATUS_SIZE )
    integer                    :: jNextReqMPIStatus( MPI_STATUS_SIZE ), jPrevReqMPIStatus( MPI_STATUS_SIZE )

    integer                    :: iNextReqSendMPIStatus( MPI_STATUS_SIZE ), iPrevReqSendMPIStatus( MPI_STATUS_SIZE )
    integer                    :: jNextReqSendMPIStatus( MPI_STATUS_SIZE ), jPrevReqSendMPIStatus( MPI_STATUS_SIZE )

    integer                    :: iNextReqSend, iPrevReqSend
    integer                    :: jNextReqSend, jPrevReqSend
#endif


    !
    ! Code...
    !
    ! Parallelize here
    ! Remap all to rank-3 
    sRemap   ( 0:nx-1, 0:ny-1, 0:get_system_size()-1 ) =>        s( 1:nx*ny*get_system_size() )
    meshRemap( 0:nx-1, 0:ny-1, 0:get_system_size()-1 ) => tendency( 1:nx*ny*get_system_size() )

    ! Parallelize here
#ifdef USE_MPI

    iNextReqMPIStatus(:)      = 1234
    iPrevReqMPIStatus(:)      = 1234
    jNextReqMPIStatus(:)      = 1234
    jPrevReqMPIStatus(:)      = 1234
    iNextReqSendMPIStatus(:)  = 1234
    iPrevReqSendMPIStatus(:)  = 1234
    jNextReqSendMPIStatus(:)  = 1234
    jPrevReqSendMPIStatus(:)  = 1234

    ! Tags are :
    ! 1 : inext
    ! 2 : iprev
    ! 3 : jnext
    ! 4 : jprev
    ! i is max, go to aRank + 1 adjacent tile
    call mpi_irecv( iNextBoundary, ny * get_system_size(), MPI_REAL,                                                          &
                    mod( ( sqrtCluster + bRank ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank + 1, sqrtCluster ), &
                    1, MPI_COMM_WORLD, iNextReq, mpiErr )
    call handleMPIErr( mpiErr )

    ! i is min, go to aRank - 1 adjacent tile
    call mpi_irecv( iPrevBoundary, ny * get_system_size(), MPI_REAL,                                                          &
                    mod( ( sqrtCluster + bRank ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank - 1, sqrtCluster ), &
                    2, MPI_COMM_WORLD, iPrevReq, mpiErr )
    call handleMPIErr( mpiErr )

    ! j is max, go to bRank + 1 adjacent tile
    call mpi_irecv( jNextBoundary, nx * get_system_size(), MPI_REAL,                                                          &
                    mod( ( sqrtCluster + bRank + 1 ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank, sqrtCluster ), &
                    3, MPI_COMM_WORLD, jNextReq, mpiErr )
    call handleMPIErr( mpiErr )

    ! j is min, go to bRank - 1 adjacent tile
    call mpi_irecv( jPrevBoundary, nx * get_system_size(), MPI_REAL,                                                          &
                    mod( ( sqrtCluster + bRank - 1 ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank, sqrtCluster ), &
                    4, MPI_COMM_WORLD, jPrevReq, mpiErr )
    call handleMPIErr( mpiErr )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Let's be nice and send out our boundary conditions as well
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Tags are FOR RECEIVER:
    ! 1 : inext
    ! 2 : iprev
    ! 3 : jnext
    ! 4 : jprev
    ! i is max, go to aRank + 1 adjacent tile
    bufferPos = 0
    call mpi_pack_size( ny * get_system_size(), MPI_REAL, MPI_COMM_WORLD, packedBytes, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_pack( sRemap(0,:,:), ny * get_system_size(), MPI_REAL, iNextBoundaryBuffer, packedBytes,             &
                   bufferPos, MPI_COMM_WORLD, mpiErr )
    call handleMPIErr( mpiErr )

    call mpi_isend( iNextBoundaryBuffer, ny * get_system_size(), MPI_REAL,                                                          &
                    mod( ( sqrtCluster + bRank ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank - 1, sqrtCluster ), &
                    1, MPI_COMM_WORLD, iNextReqSend, mpiErr )

    call handleMPIErr( mpiErr )

    ! i is min, go to aRank - 1 adjacent tile
    bufferPos = 0
    call mpi_pack_size( ny * get_system_size(), MPI_REAL, MPI_COMM_WORLD, packedBytes, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_pack( sRemap(nx - 1,:,:), ny * get_system_size(), MPI_REAL, iPrevBoundaryBuffer, packedBytes,             &
                   bufferPos, MPI_COMM_WORLD, mpiErr )
    call handleMPIErr( mpiErr )

    call mpi_isend( iPrevBoundaryBuffer, ny * get_system_size(), MPI_REAL,                                                    &
                    mod( ( sqrtCluster + bRank ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank + 1, sqrtCluster ), &
                    2, MPI_COMM_WORLD, iPrevReqSend, mpiErr )


    call handleMPIErr( mpiErr )

    ! j is max, go to bRank + 1 adjacent tile
    bufferPos = 0
    call mpi_pack_size( ny * get_system_size(), MPI_REAL, MPI_COMM_WORLD, packedBytes, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_pack( sRemap(:,0,:), ny * get_system_size(), MPI_REAL, jNextBoundaryBuffer, packedBytes,             &
                   bufferPos, MPI_COMM_WORLD, mpiErr )
    call handleMPIErr( mpiErr )
    
    call mpi_isend( jNextBoundaryBuffer, nx * get_system_size(), MPI_REAL,                                                          &
                    mod( ( sqrtCluster + bRank - 1 ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank, sqrtCluster ), &
                    3, MPI_COMM_WORLD, jNextReqSend, mpiErr )


    call handleMPIErr( mpiErr )

    ! j is min, go to bRank - 1 adjacent tile
    bufferPos = 0
    call mpi_pack_size( ny * get_system_size(), MPI_REAL, MPI_COMM_WORLD, packedBytes, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_pack( sRemap(:,ny - 1,:), ny * get_system_size(), MPI_REAL, jPrevBoundaryBuffer, packedBytes,             &
                   bufferPos, MPI_COMM_WORLD, mpiErr )
    call handleMPIErr( mpiErr )
    
    call mpi_isend( jPrevBoundaryBuffer, nx * get_system_size(), MPI_REAL,                                                     &
                    mod( ( sqrtCluster + bRank + 1 ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank, sqrtCluster ), &
                    4, MPI_COMM_WORLD, jPrevReqSend, mpiErr )


    call handleMPIErr( mpiErr )


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now retrieve data OUTSIDE threaded region - TBD learn MPI threading
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call mpi_wait( iNextReq, iNextReqMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_wait( iPrevReq, iPrevReqMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )

    
    call mpi_wait( jNextReq, jNextReqMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_wait( jPrevReq, jPrevReqMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )

    ! ! Write out block
    ! if ( processRank == 2 ) then
    !   ! write( *, * ) "Process Rank: ", processRank, " [A]", aRank, " [B]", bRank
    !   ! do i = nx-1,nx-1
    !   !   do j = 0,ny-1
    !   !     write( *, "(a,1f8.6)", advance="no" ) " ", sRemap(i,j,2)
    !   !   end do
    !   !   write( *, * ) " "
    !   ! end do
    !   ! write( *, * ) " "
    !   flush( 0 )
    !   write( *, * ) "Process Rank: ", processRank, " sent to ", mod( ( sqrtCluster + bRank ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank - 1, sqrtCluster ), &
    !                 " iNextBoundary: ", sRemap(0,:,:)
    !   write( *, * ) "Process Rank: ", processRank, " sent to ", mod( ( sqrtCluster + bRank ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank + 1, sqrtCluster ), &
    !                 " iPrevBoundary: ", sRemap(nx - 1,:,:)
    !   write( *, * ) "Process Rank: ", processRank, " sent to ", mod( ( sqrtCluster + bRank - 1 ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank, sqrtCluster ), &
    !                 " jNextBoundary: ", sRemap(:,0,:)
    !   write( *, * ) "Process Rank: ", processRank, " sent to ", mod( ( sqrtCluster + bRank + 1 ), sqrtCluster )  * sqrtCluster + mod( sqrtCluster + aRank, sqrtCluster ), &
    !                 " jPrevBoundary: ", sRemap(:,ny - 1,:)
    !   flush( 0 )
      
    !   write( *, * ) " iNextReqMPIStatus MPI_STATUS  : ", iNextReqMPIStatus(:)
    !   write( *, * ) " iPrevReqMPIStatus MPI_STATUS  : ", iPrevReqMPIStatus(:)
    !   write( *, * ) " jNextReqMPIStatus MPI_STATUS  : ", jNextReqMPIStatus(:)
    !   write( *, * ) " jPrevReqMPIStatus MPI_STATUS  : ", jPrevReqMPIStatus(:)
    !   flush( 0 )

    ! endif

    ! ! Sync all processes
    ! call mpi_barrier( MPI_COMM_WORLD, mpiErr )

    ! if ( processRank == 0 ) then
    !   flush( 0 )
    !   write( *, * ) "Process Rank: ", processRank, " jNextBoundary : ", jNextBoundary
    !   flush( 0 )
    !   write( *, * ) "Process Rank: ", processRank, " jPrevBoundary : ", jPrevBoundary
    !   flush( 0 )
    !   write( *, * ) "jNextBoundary status : [MPI_SOURCE] : ", jNextReqMPIStatus(MPI_SOURCE), &
    !                                          " [MPI_TAG] : ", jNextReqMPIStatus(MPI_TAG),    &
    !                                        " [MPI_ERROR] : ", jNextReqMPIStatus(MPI_ERROR),  &
    !                                        " MPI_STATUS  : ", jNextReqMPIStatus(:)
    !   flush( 0 )
    !   write( *, * ) "jPrevBoundary status : [MPI_SOURCE] : ", jPrevReqMPIStatus(MPI_SOURCE), &
    !                                          " [MPI_TAG] : ", jPrevReqMPIStatus(MPI_TAG),    &
    !                                        " [MPI_ERROR] : ", jPrevReqMPIStatus(MPI_ERROR),  &
    !                                        " MPI_STATUS  : ", jPrevReqMPIStatus(:)
    !   flush( 0 )
    !   ! call handleMPIErr( jNextReqMPIStatus(MPI_ERROR) )
    !   ! call handleMPIErr( jPrevReqMPIStatus(MPI_ERROR) )

    ! endif
    
    ! ! Sync all processes
    ! call mpi_barrier( MPI_COMM_WORLD, mpiErr )

    ! if ( processRank == 3 ) then
    !   flush( 0 )
    !   write( *, * ) "Process Rank: ", processRank, " iNextBoundary : ", iNextBoundary
    !   flush( 0 )
    !   write( *, * ) "Process Rank: ", processRank, " iPrevBoundary : ", iPrevBoundary
    !   flush( 0 )
    !   write( *, * ) "iNextBoundary status : [MPI_SOURCE] : ", iNextReqMPIStatus(MPI_SOURCE), &
    !                                          " [MPI_TAG] : ", iNextReqMPIStatus(MPI_TAG),    &
    !                                        " [MPI_ERROR] : ", iNextReqMPIStatus(MPI_ERROR),  &
    !                                        " MPI_STATUS  : ", iNextReqMPIStatus(:)
    !   flush( 0 )
    !   write( *, * ) "iPrevBoundary status : [MPI_SOURCE] : ", iPrevReqMPIStatus(MPI_SOURCE), &
    !                                          " [MPI_TAG] : ", iPrevReqMPIStatus(MPI_TAG),    &
    !                                        " [MPI_ERROR] : ", iPrevReqMPIStatus(MPI_ERROR),  &
    !                                        " MPI_STATUS  : ", iPrevReqMPIStatus(:)
    !   flush( 0 )
    !   ! call handleMPIErr( iNextReqMPIStatus(MPI_ERROR) )
    !   ! call handleMPIErr( iPrevReqMPIStatus(MPI_ERROR) )
    ! endif


#endif


#ifdef _OPENMP 
    ! write( *, * ) "Running OpenMP"
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( i, j, iNext, iPrev, jNext, jPrev, iNextPtr, iPrevPtr, jNextPtr, jPrevPtr ) 
#endif
    do i = 0, nx-1
      do j = 0, ny-1
        iNext = mod( nx + i + 1, nx )
        iPrev = mod( nx + i - 1, nx )
        jNext = mod( ny + j + 1, ny )
        jPrev = mod( ny + j - 1, ny )

        ! Remap sRemap to basically each indexing array
        iNextPtr => sRemap
        iPrevPtr => sRemap
        jNextPtr => sRemap
        jPrevPtr => sRemap


#ifdef USE_MPI
        ! This could be moved to outer loop
        if ( iNext == 0 ) then
          ! Remap to what we received, custom to bounds iNext
          ! write( *, * ) "i is ", i, "... remapping iNextPtr"
          iNextPtr( iNext:iNext, 0:ny-1, 0:get_system_size()-1 ) => iNextBoundary(:)
          
        else if ( iPrev == ( nx - 1 ) ) then
          ! write( *, * ) "i is ", i, "... remapping iPrevPtr"
          ! Remap to what we received, custom to bounds iPrev
          iPrevPtr( iPrev:iPrev, 0:ny-1, 0:get_system_size()-1 ) => iPrevBoundary(:)

        endif

        if ( jNext == 0 ) then
          ! write( *, * ) "j is ", j, "... remapping jNextPtr"
          ! Remap to what we received, custom to bounds jNext
          jNextPtr( 0:nx-1, jNext:jNext, 0:get_system_size()-1 ) => jNextBoundary(:)

        else if ( jPrev == ( ny - 1 ) ) then
          ! write( *, * ) "j is ", j, "... remapping jPrevPtr"
          ! Remap to what we received, custom to bounds jPrev
          jPrevPtr( 0:nx-1, jPrev:jPrev, 0:get_system_size()-1 ) => jPrevBoundary(:)

        endif
#endif

        ! Calc each now
        ! Possible optimizations include storing sRemap( iPrev:iNext, jPrev:jNext, : ) in temp vars but that may be resolved
        ! via cache hits of them being collocated
        meshRemap(i,j,0) = -sRemap(i,j,0) * ( iNextPtr(iNext,j,0) - iPrevPtr(iPrev,j,0) ) / ( 2.0 * deltaX ) -              &
                            sRemap(i,j,1) * ( jNextPtr(i,jNext,0) - jPrevPtr(i,jPrev,0) ) / ( 2.0 * deltaY ) -              &
                            g * ( iNextPtr(iNext,j,2) - iPrevPtr(iPrev,j,2) ) / ( 2.0 * deltaX )
        
        meshRemap(i,j,1) = -sRemap(i,j,0) * ( iNextPtr(iNext,j,1) - iPrevPtr(iPrev,j,1) ) / ( 2.0 * deltaX ) -              &
                            sRemap(i,j,1) * ( jNextPtr(i,jNext,1) - jPrevPtr(i,jPrev,1) ) / ( 2.0 * deltaY ) -              &
                            g * ( jNextPtr(i,jNext,2) - jPrevPtr(i,jPrev,2) ) / ( 2.0 * deltaY )

        meshRemap(i,j,2) = -( iNextPtr(iNext,j,2) * iNextPtr(iNext,j,0) - iPrevPtr(iPrev,j,2) * iPrevPtr(iPrev,j,0) ) / ( 2.0 * deltaX ) - &
                            ( jNextPtr(i,jNext,2) * jNextPtr(i,jNext,1) - jPrevPtr(i,jPrev,2) * jPrevPtr(i,jPrev,1) ) / ( 2.0 * deltaY )
      end do
    end do

#ifdef _OPENMP 
    !$OMP END PARALLEL DO
#endif

#ifdef USE_MPI
    ! Wait for message to be received
    call mpi_wait( iNextReqSend, iNextReqSendMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_wait( iPrevReqSend, iPrevReqSendMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )

    
    call mpi_wait( jNextReqSend, jNextReqSendMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )
    call mpi_wait( jPrevReqSend, jPrevReqSendMPIStatus, mpiErr )
    call handleMPIErr( mpiErr )

    ! call mpi_barrier( MPI_COMM_WORLD, mpiErr )

    ! if ( processRank == 2 ) then
    !   flush( 0 )
      
    !   write( *, * ) " iNextReqSendMPIStatus MPI_STATUS  : ", iNextReqSendMPIStatus(:)
    !   write( *, * ) " iPrevReqSendMPIStatus MPI_STATUS  : ", iPrevReqSendMPIStatus(:)
    !   write( *, * ) " jNextReqSendMPIStatus MPI_STATUS  : ", jNextReqSendMPIStatus(:)
    !   write( *, * ) " jPrevReqSendMPIStatus MPI_STATUS  : ", jPrevReqSendMPIStatus(:)
    !   flush( 0 )
    ! endif

    ! call mpi_barrier( MPI_COMM_WORLD, mpiErr )
#endif

  end subroutine f

  subroutine equations_finalize( s )
    real, dimension(:), intent(inout), pointer :: s

#ifdef USE_MPI
    deallocate( iNextBoundary )
    deallocate( iPrevBoundary )
    deallocate( jNextBoundary )
    deallocate( jPrevBoundary )

    deallocate( iNextBoundaryBuffer )
    deallocate( iPrevBoundaryBuffer )
    deallocate( jNextBoundaryBuffer )
    deallocate( jPrevBoundaryBuffer )
#endif

    deallocate( s )
  end subroutine equations_finalize

end module equations
