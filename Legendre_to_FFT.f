!**********************************************************
! This program takes a uniform 3-D grid and performs the  *
! various transforms on the matrix for the different      *
! processes within the model.                             *
!                                                         *
! Grid-Space: Column models.  Divide the grid by latitude *
!   and longitude, but keep all Z for each lat/lon        *
!   together.                                             *
! FFT: Divide the grid by Z and latitude, but keep all    *
!   longitude (m) for each Z / lat pair on a single       *
!   processor                                             *
! Legendre: Divide the grid by Z and longitude (m), but   *
!   keep all latitudes (n) for each Z / lon pair on a     *
!   single processor                                      *
!**********************************************************

      Subroutine Legendre_to_FFT(Legendre_Array, FFT_Array)

!----------------------------------------------------------
! Purpose:
!    Perform transformation from Legendre space to 
!       FFT space

! Revision History:
!   Aug 2015 Matthew Turner: Created.  Need to adapt
!                            for thinned grids.
!----------------------------------------------------------

! Module files:
      Use Grid_Config

      Implicit None

! Include Files
      Include 'mpif.h'

! Passed variables:
      Real( 8 ), Intent(In) :: Legendre_Array(:,:,:)  ! Input array in Legendre space
      Real( 8 ), Intent(InOut) :: FFT_Array(:,:,:)  ! output array in FFT space

! Local variables:
      ! Allocatables
      Integer, Allocatable :: MPI_Ind(:)   ! Array of rank numbers
                                           !   for MPI communication
      Real( 8 ), Allocatable :: Buffer(:,:,:)

      Integer :: i, proc,ii ! counter
      Integer :: count     ! number of bytes to be sent / received
      Integer :: Astat     ! Allocation status
      Character(80) :: Msg   ! Error message
      Character(24) :: ProcessName = 'Legendre_to_FFT'

      ! Lat/lon/z send/recv temporary values
      Integer :: Start_lat_send ! Starting latitude of sending process
      Integer :: Start_lat_recv ! Starting latitude of receiving process
      Integer :: End_lat_send ! Ending latitude of sending process
      Integer :: End_lat_recv ! Ending latitude of receiving process
      Integer :: Start_lon_send ! Starting longitude of sending process
      Integer :: Start_lon_recv ! Starting longitude of receiving process
      Integer :: End_lon_send ! Ending longitude of sending process
      Integer :: End_lon_recv ! Ending longitude of receiving process
      Integer :: Start_Z_send ! Starting Z of sending process
      Integer :: Start_Z_recv ! Starting Z of receiving process
      Integer :: End_Z_send ! Ending Z of sending process
      Integer :: End_Z_recv ! Ending Z of receiving process

      Integer :: status, mpi_err

!-------------------------------- Code Begins -----------------------------

      Allocate(MPI_Ind(Num_Procs_Row), Stat = Astat)
      If ( Astat /= 0 ) Then
         Msg = 'Error allocating MPI_Ind in '
         Call alloc_status(Msg, ProcessName)
      End If

      ! Determine the processor Ranks that this processor will
      ! communicate with.
      Do i = 1, Num_Procs_Row
         MPI_Ind(i) = Mod((Rank),Num_Procs_Col) + Num_Procs_Col *
     &                (i - 1)
      End Do

      ! Loop through MPI_Id and send / receive
      Do ii = 1, Num_Procs_Row
         ! Define starting and ending lat/lon/z for the sending process
         Start_lat_send = Start_lat(MPI_Ind(ii))
         End_lat_send = End_lat(MPI_Ind(ii))
         Start_lon_send = Start_zonal(MPI_Ind(ii))
         End_lon_send = End_zonal(MPI_Ind(ii))
         Start_Z_send = Start_Z(MPI_Ind(ii))
         End_Z_send = End_Z(MPI_Ind(ii))

         ! Check if this process is the indexed process (MPI_Ind(ii)).
         ! If so, send a snippet of the grid-space array to the other
         ! processors within MPI_Ind.
         ! If not, receive a snippet of the grid-space array from the
         ! indexed process
         If ( Rank == MPI_Ind(ii) ) Then
            ! Loop through the receiving processors and send the data
            Do proc = 1, Num_Procs_Row
               ! If the receiving processor Rank is not the rank of this
               ! processor
               If ( Rank /= MPI_Ind(proc) ) Then
                  ! Calculate the starting and ending z of the
                  ! receiving process
                  Start_lat_recv = Start_lat(MPI_Ind(proc))
                  End_lat_recv = End_lat(MPI_Ind(proc))
   
                  ! Calculate the number of bytes to send (for MPI_SEND)
                  count =  (End_lat_recv - Start_lat_recv + 1) *
     &                     (End_lon_send - Start_lon_send + 1) *
     &                     (End_Z_send - Start_Z_send + 1)
                  ! Allocate the Buffer array for the data to be sent
                  Allocate(Buffer(End_lat_recv - Start_lat_recv + 1,
     &                   End_lon_send - Start_lon_send + 1,
     &                   End_Z_send - Start_Z_send + 1),
     &                   Stat = Astat)
                  If ( Astat /= 0 ) Then
                     Msg = 'Error allocating Buffer in '
                     Call alloc_status(Msg, ProcessName)
                  End If

                  ! Fill the Buffer with the data for the Z of the
                  ! receiving array
                  Buffer = Legendre_Array(Start_lat_recv:End_lat_recv,
     &                                    :,:)

                  ! Send the data
                  Call MPI_SEND(Buffer, count, MPI_DOUBLE,
     &                       MPI_Ind(proc), 0, MPI_COMM_WORLD, mpi_err)
                  If (mpi_err /= MPI_SUCCESS) Then
                     Write(*,*) 'MPI Error sending.  Terminating'
                     Call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_err)
                  End If

                  ! Deallocate the Buffer.  This is necessary as each
                  ! processor can possibly have a different sized Buffer
                  Deallocate(Buffer)
               End If ! Rank /= MPI_Ind(proc)
            End Do  ! Loop over proc

            ! Loop through the longitudes of this process and fill
            ! FFT_Array
            Do proc = Start_lon_send, End_lon_send
               FFT_Array(:,proc,:) =
     &             Legendre_Array(Start_lat_send:End_lat_send,
     &                        proc - Start_lon_send + 1,:)
            End Do
         Else  ! If this process is not the indexed process
            ! Calculate the starting and ending lat/lon/z of the
            ! receiving process
            Start_lat_recv = Start_lat(Rank)
            End_lat_recv = End_lat(Rank)
            Start_lon_recv = Start_zonal(Rank)
            End_lon_recv = End_zonal(Rank)
            Start_Z_recv = Start_Z(Rank)
            End_Z_recv = End_Z(Rank)

            ! Calculate the number of bytes to receive (use this
            ! processor's lat / Z, and sending processor's lon)
            count =  (End_lat_recv - Start_lat_recv + 1) * 
     &               (End_lon_send - Start_lon_send + 1)
     &             * (End_Z_send - Start_Z_send + 1)

            ! Allocate the Buffer array for the data to be received
            Allocate(Buffer(End_lat_recv - Start_lat_recv + 1,
     &               End_lon_send - Start_lon_send + 1,
     &               End_Z_send - Start_Z_send + 1),
     &               Stat = Astat)
            If ( Astat /= 0 ) Then
               Msg = 'Error allocating Buffer in '
               Call alloc_status(Msg, ProcessName)
            End If

            ! Receive data from indexed processor
            Call MPI_RECV(Buffer, count, MPI_DOUBLE,
     &                   MPI_Ind(ii), 0, MPI_COMM_WORLD, status, mpi_err)
            If (mpi_err /= MPI_SUCCESS) Then
               Write(*,*) 'MPI Error receiving.  Terminating'
               Call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_err)
            End If
 
            ! Loop through longitudes of the received data and fill
            ! FFT_Array
            Do proc = Start_lon_send, End_lon_send
               FFT_Array(:,proc,:) = Buffer(:,proc - 
     &                                  Start_lon_send + 1, :)
            End Do

            ! Deallocate Buffer array
            Deallocate(Buffer)
         End If
         ! This statement should not be needed, but the code seems to
         ! continue looping beyond the limit with optimization turned on
!         If ( ii == Num_Procs_Row ) Exit
      End Do  ! Loop over MPI_Ind

      ! Debug, mdt :: print the rank and FFT_Array for the given
      ! processor
      If ( rank == dRank  ) then
         write(*,*) 'Rank iFFT = ', Rank
         write(*,*) 'fft shape = ', shape(FFT_Array)
         write(*,*) 'max(fft) = ',maxval(FFT_Array(:,:,1))
         write(*,*) 'min(fft) = ',minval(FFT_Array(:,:,1))
!         write(*,*) 'iFFT = '
!         do i = Start_Z(Rank),End_Z(Rank)
!         write(*,*) FFT_Array
!         end do
      end if

      End Subroutine Legendre_to_FFT

