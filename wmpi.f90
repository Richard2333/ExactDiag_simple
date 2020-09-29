  module prec
     !>> A module controls the precision. 
     !> when the nnzmax is larger than 2,147,483,647 then li=8,
     !> otherwise, li=4. 
     !> warning: li=4 was tested, li=8 is not tested yet.
     integer,parameter :: li=4 ! long integer
     integer,parameter :: Dp=kind(1.0d0) ! double precision  
  end module prec
  module wmpi
     use prec

#if defined (MPI)
     include 'mpif.h'
#endif

     integer :: cpuid  ! CPU id for mpi
     integer :: num_cpu  ! Number of processors for mpi

#if defined (MPI)
     integer, parameter :: mpi_in= mpi_integer
     integer, parameter :: mpi_dp= mpi_double_precision
     integer, parameter :: mpi_dc= mpi_double_complex
     integer, parameter :: mpi_cmw= mpi_comm_world
#endif 

     !> Define a structure containing information for doing communication
     type WTParCSRComm
  
        !> mpi communicator
        integer :: comm
  
        !> how many cpus that we need to send data on
        integer :: NumSends
  
        !> which cpus that we need to send data on
        integer, pointer :: SendCPUs(:)
  
        !> when before we send the vector data to other cpus, we need to get the 
        !> data which should be sent, then put these data to a array called 
        !> x_buf_data(:). The array SendMapElements(:) gives the position of the 
        !> data in vector that should be sent.
        integer, pointer :: SendMapStarts(:)
  
        !> with this array, we can select the vector data that should be sent
        integer(li), pointer :: SendMapElements(:)
  
        !> How many cpus that we need to recieve data from
        integer :: NumRecvs
  
        !> Which cpus that we need to recieve data from
        integer, pointer :: RecvCPUs(:)
  
        !> When recieved data from other cpus, we need to arrange those data into 
        !> an array. The length of this 
        integer, pointer :: RecvVecStarts(:)
     end type WTParCSRComm
  
     !> Define a structure containing information for doing communication
     type WTParVecComm
  
        !> mpi communicator
        integer :: comm
  
        !> how many cpus that we need to send data on
        integer :: NumSends
  
        !> which cpus that we need to send data on
        integer, pointer :: SendCPUs(:)
  
        !> when before we send the vector data to other cpus, we need to get the 
        !> data which should be sent, then put these data to a array called 
        !> x_buf_data(:). The array SendMapElements(:) gives the position of the 
        !> data in vector that should be sent.
        integer, pointer :: RecvMapStarts(:)
  
        !> with this array, we can select the vector data that should be sent
        integer(li), pointer :: RecvMapElements(:)
  
        !> How many cpus that we need to recieve data from
        integer :: NumRecvs
  
        !> Which cpus that we need to recieve data from
        integer, pointer :: RecvCPUs(:)
  
        !> When recieved data from other cpus, we need to arrange those data into 
        !> an array. The length of this 
        integer, pointer :: SendVecStarts(:)
  
        integer :: NumRowsDiag
        integer :: NumRowsOffd
  
        integer(li), pointer :: RowMapOffd(:)
        integer(li), pointer :: LocalIndexOffd(:)
        integer(li), pointer :: LocalIndexDiag(:)
     end type WTParVecComm
  
  
  
     !> define a handle for comm, that can be created and destroyed
     type WTCommHandle
  
        type(WTParCSRComm), pointer :: sendrecv
  
        integer :: numrequest
        integer, pointer :: mpirequest(:)
        
        complex(dp), pointer :: senddata(:)
        complex(dp), pointer :: recvdata(:)
  
     end type WTCommHandle
  
     integer               :: BasisStart
     integer               :: BasisEnd
  
     contains
  
     !> generate partition for any vector with a given length 
     subroutine WTGeneratePartition(length, nprocs, part)

        implicit none

        integer(li), intent(in) :: length
        integer, intent(in) :: nprocs
        integer(li), intent(out) :: part(nprocs+1)
  
        integer :: i
        integer(li) :: div
        integer(li) :: mod1
  
        mod1= mod(length, nprocs)
        if (mod1.eq.0) then !< each cpu has the same load balance
           div= length/nprocs
           do i=0, nprocs-1
              part(i+1)=1+ i*div
           enddo
        else
           div= length/nprocs+ 1
           do i=0, nprocs
              if (i.ge. (nprocs-mod1)) then
                 part(i+1)= 1+ i*div- (nprocs-mod1)
              else
                 part(i+1)= 1+ i*(div-1) !< the main cpu will get smaller data
              endif
           enddo
        endif
        part(nprocs+1)= length+ 1
  
        return
     end subroutine WTGeneratePartition
  
     !> generate local partition for any vector with a given length 
     subroutine WTGenerateLocalPartition(length, nprocs, icpu, first, last)
        implicit none
        integer, intent(in) :: nprocs
        integer, intent(in) :: icpu
        integer(li), intent(in) :: length
        integer(li), intent(out) :: first
        integer(li), intent(out) :: last
  
        integer(li) :: div
        integer(li) :: mod1
  
        mod1= mod(length, nprocs)
        if (mod1.eq.0) then !< each cpu has the same load balance
           div= length/nprocs
           first=1+ icpu*div
           last=(1+ icpu)*div
        else
           div= length/nprocs+ 1
           if (icpu.ge. (nprocs-mod1)) then
              first= 1+ icpu*div- (nprocs-mod1)
              last= (1+ icpu)*div- (nprocs-mod1)
           else
              first= 1+ icpu*(div-1)
              last= (1+ icpu)*(div-1)!< the main cpu will get smaller data
           endif
        endif
  
        return
     end subroutine WTGenerateLocalPartition
  
     !> given the send data, recieve data and sendrecv list, we can use
     !> this subroutine to send and recieve data from other cpus.
     !> when finished this subroutine calls, we need call WTCommHandleDestroy
     !> to check whether the send recv operation is finished
     subroutine WTCommHandleCreate(SendRecv, SendData, RecvData, &
             CommHandle)
  
        implicit none
  
        !* in variables
        type(WTParCSRComm), intent(in), pointer :: SendRecv
        complex(dp), pointer :: SendData(:)
        complex(dp), pointer :: RecvData(:)
  
        !* out variables
        type(WTCommHandle), pointer :: CommHandle
  
        integer(li) :: VecStart
        integer(li) :: VecLen
  
        !> sendrecv data
        integer, pointer :: SendCPUs(:)
        integer, pointer :: SendMapStarts(:)
        integer(li), pointer :: SendMapElements(:)
        integer, pointer :: RecvCPUs(:)
        integer, pointer :: RecvVecStarts(:)
  
        integer :: NumRecvs
        integer :: NumSends
        integer :: NumRequest
        integer, pointer :: MpiRequest(:)
  
        integer :: ierr
        integer :: comm
        integer :: NCPUs
        integer :: cpu_id
  
        integer :: i, j
        integer :: icpu
  
        !* initialize null pointers
        SendCPUs=> Null()
        SendMapStarts=> Null()
        SendMapElements=> Null()
        RecvCPUs=> Null()
        RecvVecStarts=> Null()
        MpiRequest=> Null()
  
#if defined (MPI)
        comm= SendRecv%Comm
        call mpi_comm_size(comm, NCPUS, ierr)
        call mpi_comm_rank(comm, cpu_id, ierr)
#endif  
        NumSends= SendRecv%NumSends
        NumRecvs= SendRecv%NumRecvs
        NumRequest= NumSends+ NumRecvs
  
        SendCPUs=> SendRecv%SendCPUs
        SendMapStarts=> SendRecv%SendMapStarts
        SendMapElements=> SendRecv%SendMapElements
        RecvCPUs=> SendRecv%RecvCPUs
        RecvVecStarts=> SendRecv%RecvVecStarts
  
        MpiRequest=> Null()
        if (.not.associated(CommHandle)) allocate(CommHandle)
        allocate(CommHandle%MpiRequest(NumRequest))
        MpiRequest=> CommHandle%MpiRequest
  
        !> recieve data from other cpus using non-block communication
        j=1
        do i=1, NumRecvs
           icpu= RecvCPUs(i)
           VecStart= RecvVecStarts(i)
           VecLen= RecvVecStarts(i+1)- RecvVecStarts(i)
#if defined (MPI)
           call mpi_irecv(RecvData(VecStart), VecLen, mpi_dc, icpu, &
                          0, comm, MpiRequest(j), ierr)
#endif
           j=j+1
        enddo
  
        !> send data to other cpus using non-block communication
        do i=1, NumSends
           icpu= SendCPUs(i)
           VecStart= SendMapStarts(i)
           VecLen= SendMapStarts(i+1)- SendMapStarts(i)
#if defined (MPI)
           call mpi_isend(SendData(VecStart), VecLen, mpi_dc, icpu, &
                          0, comm, MpiRequest(j), ierr)
#endif
           j=j+1
        enddo
  
        CommHandle%NumRequest= NumRequest
        CommHandle%SendData=> SendData
        CommHandle%RecvData=> RecvData
        CommHandle%SendRecv=> SendRecv
  
        return
  
     end subroutine WTCommHandleCreate
  
     subroutine WTCommHandleDestroy(CommHandle)
  
        implicit none
  
        type(WTCommHandle), pointer :: CommHandle
  
        integer :: ierr
        integer :: NumRequest
        integer, pointer :: MpiRequest(:)
        integer, pointer :: MpiStatus(:)
  
        NumRequest= CommHandle%NumRequest
        MpiRequest=> CommHandle%MpiRequest
  
        if (.not.associated(CommHandle)) return
  
        if (NumRequest>0) then
#if defined (MPI)
           allocate(MpiStatus(mpi_status_size*NumRequest))
           call mpi_waitall(NumRequest, MpiRequest, MpiStatus, ierr)
#endif
        endif
  
        return
     end subroutine WTCommHandleDestroy
  
     !> define my mpi_allreduce for complex array
     subroutine mp_allreduce_z(comm, ndim, vec, vec_mpi)
        implicit none
  
        integer, intent(in) :: comm
        integer, intent(in) :: ndim
        complex(dp), intent(in) :: vec(ndim)
        complex(dp), intent(inout) :: vec_mpi(ndim)
  
        integer :: ierr
  
        vec_mpi= 0d0
  
#if defined (MPI)
        call mpi_allreduce(vec, vec_mpi, ndim, &
                           mpi_dc, mpi_sum, comm, ierr)
#else
        vec_mpi= vec
#endif
  
        return
     end subroutine mp_allreduce_z

  end module wmpi
