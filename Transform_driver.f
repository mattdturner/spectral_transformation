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

      Program Transform_Driver

!----------------------------------------------------------
! Purpose:
!    Program Transform_Driver initializes the full grid
!       and MPI processes.  It then runs the actual
!       transformation separately for each processor.

! Revision History:
!   Aug 2015 Matthew Turner: Created.  Need to adapt
!                            for thinned grids.
!----------------------------------------------------------

! Module files:
      Use Grid_Config    ! Grid configuration variables

      Implicit None

! Include Files
      Include 'mpif.h'

! Local variables:
      Logical, Save :: Firsttime = .True.
      Integer :: Lat_per_proc  ! Number of latitude cells per processor
      Integer :: Lon_per_proc  ! Number of longitude cells per processor
      Integer :: Z_per_proc    ! Number of Z cells per processor
      Integer :: Zonal_per_proc    ! Number of lon cells per processor for
                                 !   legendre space
      Integer :: Total_per_proc  ! Number of latitude cells per processor for spectral
      Integer :: Rem_lat       ! Number of leftover latitude cells after
                               !   integer division
      Integer :: Rem_lon       ! Number of leftover longitude cells after
                               !   integer division
      Integer :: Rem_Z         ! Number of leftover Z cells after
                               !   integer division
      Integer :: Rem_Zonal       ! Number of leftover longitude cells after
                               !   integer division for legendre space
      Integer :: Rem_total  ! Number of leftover latitude cells after
                               !   integer division for spectral 

      Real( 8 ), Allocatable :: Full_Grid(:,:,:)   ! Dummy variable for data
      Real( 8 ), Allocatable :: Grid_Array(:,:,:)   ! Dummy variable for data
      Real( 8 ), Allocatable :: Pre_FFT_Array(:,:,:)   ! Dummy variable for data
      Real( 8 ), Allocatable :: FFT_Array(:,:,:)   ! Dummy variable for data
      Real( 8 ), Allocatable :: Pre_Legendre_Array(:,:,:)   ! Dummy variable for data
      Real( 8 ), Allocatable :: Legendre_Array(:,:,:)   ! Dummy variable for data
      Real( 8 ), Allocatable :: Spectral_Array(:,:,:)   ! Dummy variable for data
      Real, Allocatable :: asig(:), bsig(:), alat(:)
      Integer, Allocatable :: lix(:)

      Integer :: i, j, k, proc_index, cell_counter, tmp_index   ! Counter
      Integer :: ierr ! MPI error value
      Integer :: kx

      Integer Astat   ! Array allocation error status
      Character(80) :: Msg
      Character(24) :: ProcessName = 'Transform_Driver'
      
      ! Interface
      Interface
         Subroutine Grid_to_FFT(Grid_Array, Pre_FFT_Array)
            Implicit None
            Real( 8 ), Intent(In)    :: Grid_Array(:,:,:)
            Real( 8 ), Intent(InOut) :: Pre_FFT_Array(:,:,:)
         End Subroutine Grid_to_FFT
         Subroutine FFT_to_Legendre(FFT_Array, Legendre_Array)
            Implicit None
            Real( 8 ), Intent(In)    :: FFT_Array(:,:,:)
            Real( 8 ), Intent(InOut) :: Legendre_Array(:,:,:)
         End Subroutine FFT_to_Legendre
         Subroutine Legendre_to_Spectral(Legendre_Array, Spectral_Array)
            Implicit None
            Real( 8 ), Intent(In)    :: Legendre_Array(:,:,:)
            Real( 8 ), Intent(InOut) :: Spectral_Array(:,:,:)
         End Subroutine Legendre_to_Spectral
         Subroutine Spectral_to_Legendre(Spectral_Array, Legendre_Array)
            Implicit None
            Real( 8 ), Intent(In)    :: Spectral_Array(:,:,:)
            Real( 8 ), Intent(InOut) :: Legendre_Array(:,:,:)
         End Subroutine Spectral_to_Legendre
         Subroutine Legendre_to_FFT(Legendre_Array, FFT_Array)
            Implicit None
            Real( 8 ), Intent(In)    :: Legendre_Array(:,:,:)
            Real( 8 ), Intent(InOut) :: FFT_Array(:,:,:)
         End Subroutine Legendre_to_FFT
         Subroutine FFT_to_Grid(FFT_Array, Grid_Array)
            Implicit None
            Real( 8 ), Intent(In)    :: FFT_Array(:,:,:)
            Real( 8 ), Intent(InOut) :: Grid_Array(:,:,:)
         End Subroutine FFT_to_Grid
         Subroutine alloc_status(Msg, ProcessName)
            Implicit None
            Character(80), Intent(In) :: Msg
            Character(24), Intent(In) :: ProcessName
         End Subroutine alloc_status
      End Interface
                               
!-------------------------------- Code Begins -----------------------------

      ! Read in domain file
      Open(Unit = 10, File = 
     &     '/home/mturner/Transform/noggeom_thinned.txt',
     &     Status = 'Old', Action = 'Read')
      read(10,'(a80)') ! Skip 1st line
      ! Read in latitude, max_lon, levels
      read(10,'(4i10)') lon,lat,Z !, num_points

      zonal_waves = lon / 2
      total_waves = ( lon / 3 ) - 1

      ! Allocate variables to be read from domain file
      Allocate (asig(Z+1), 
     &          bsig(Z+1), 
     &          alat(lat),
     &          lix(lat), Stat = Astat )
      If (Astat /= 0) Then
         Msg = 'Error allocating asig and bsig in '
         Call alloc_status(Msg, ProcessName)
      End If

      read(10,'(a80)')  ! Skip line
      ! Read in asig and bsig
      read(10,'(i5,2f20.15)') (kx,asig(k),bsig(k),k=1,Z+1)
      read(10,'(a80)')  ! Skip line
      ! Read in latitude centers and number of lon per lat
      read(10,'(i5,5x,f20.15,i6)') (kx,alat(j),lix(j),j=1,lat)

      ! Close domain file
      Close(10)

      ! Define the domain parameters (Number of columns and rows of
      ! processors)
      !   This will eventually be read from a run script
      Num_Procs_Col = 6
      Num_Procs_Row = 4

      ! Initialize MPI
      Call MPI_INIT(ierr)
      If (ierr /= MPI_SUCCESS) Then
         Write(*,*) 'Error starting MPI.  Terminating'
         Call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      End If

      ! Get the number of processors this job is using
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, Num_Procs, ierr)
      If (ierr /= MPI_SUCCESS) Then
         Write(*,*) 'Error getting number of MPI processes. Terminating'
         Call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      End If
      If ( Num_Procs /= Num_Procs_Col * Num_Procs_Row ) Then
         Write(*,*) 'num_proc = ', Num_Procs
         Write(*,*) 'num_multiply = ', Num_Procs_Col, Num_Procs_Row
         Write(*,*) 'Number of processors does not match the number of',
     &              ' rows x columns.  Terminating'
         Call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      End If

      ! Get the rank of the processor this thread is running on.
      Call MPI_COMM_RANK(MPI_COMM_WORLD, Rank, ierr)
      If (ierr /= MPI_SUCCESS) Then
         Write(*,*) 'Error getting MPI Rank.  Terminating'
         Call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      End If

      If ( Firsttime ) Then
         Firsttime = .False.
         ! Allocate arrays
         Allocate ( Full_Grid(lat,lon,Z), Stat = Astat )
         If ( Astat /= 0 ) Then
            Msg = 'Error allocating Full_Grid in process '
            Call alloc_status(Msg, ProcessName)
         End If
         If ( .not. Grid_Init() ) Then
            Msg = 'Error allocating grid variables in process '
            Call alloc_status(Msg, ProcessName)
         End If
      End If ! Firsttime

      ! Full the full grid array with a counter
      cell_counter = 0
      DO k = 1, Z
         Do i = 1, lat
            Do j = 1, lon
               If ( j <= lix(i) ) Then
                  cell_counter = cell_counter + 1
                  Full_Grid(i,j,k) = cell_counter
               Else
                  Full_Grid(i,j,k) = -999
               End If
            End Do
         End Do
      End Do

      ! Determine which cells go on which processor
      Lat_per_proc = lat / Num_Procs_Row
      Rem_lat = Mod(lat,Num_Procs_Row)
      Lon_per_proc = lon / Num_Procs_Col
      Rem_lon = Mod(lon,Num_Procs_Col)
      Z_per_proc = Z / Num_Procs_Col
      Rem_Z = Mod(Z,Num_Procs_Col)
      Zonal_per_proc = lon / (2 * Num_Procs_Row)
      Rem_Zonal = Mod(lon,2 * Num_Procs_Row)
      Total_per_proc = (lon / 3 - 1) / Num_Procs_Col
      Rem_total = Mod((lon / 3 - 1),Num_Procs_Col)

      Do i = 1, Num_Procs_Row
         If ( i <= Rem_lat ) Then
            Lat_on_proc(i) = Lat_per_proc + 1
         Else
            Lat_on_proc(i) = Lat_per_proc
         End If
      End Do

      Do j = 1, Num_Procs_Col
         If ( j <= Rem_lon ) Then
            Lon_on_proc(j) = Lon_per_proc + 1
         Else
            Lon_on_proc(j) = Lon_per_proc
         End If
      End Do

      Do k = 1, Num_Procs_Col
         If ( k <= Rem_Z ) Then
            Z_on_proc(k) = Z_per_proc + 1
         Else
            Z_on_proc(k) = Z_per_proc
         End If
      End Do

      Do j = 1, Num_Procs_Row
         If ( j <= Rem_Zonal ) Then
            Zonal_on_proc(j) = Zonal_per_proc + 1
         Else
            Zonal_on_proc(j) = Zonal_per_proc
         End If
      End Do

      Do j = 1, Num_Procs_Col
         If ( j <= Rem_total ) Then
            Total_on_proc(j) = Total_per_proc + 1
         Else
            Total_on_proc(j) = Total_per_proc
         End If
      End Do

      Do proc_index = 1, Num_Procs
         tmp_index = Mod(proc_index,Num_Procs_Col)
         If ( tmp_index == 0 ) Then
            tmp_index = Num_Procs_Col
         End If

         If ( tmp_index == 1 ) Then
            Start_lon(proc_index-1) = 1
            End_lon(proc_index-1) = Lon_on_proc(tmp_index)
            Start_Z(proc_index-1) = 1
            End_Z(proc_index-1) = Z_on_proc(tmp_index)
            Start_total(proc_index-1) = 1
            End_total(proc_index-1) =
     &                       Total_on_proc(tmp_index)
         Else
            Start_lon(proc_index-1) = End_lon(proc_index-2) + 1
            End_lon(proc_index-1) = End_lon(proc_index-2) + 
     &                            Lon_on_proc(tmp_index)
            Start_Z(proc_index-1) = End_Z(proc_index-2) + 1
            End_Z(proc_index-1) = End_Z(proc_index-2) +
     &                            Z_on_proc(tmp_index)
            Start_total(proc_index-1) =
     &                         End_total(proc_index-2) + 1
            End_total(proc_index-1) =
     &                         End_total(proc_index-2) + 
     &                         Total_on_proc(tmp_index)
         End If

         tmp_index = (proc_index - 1) / Num_Procs_Col + 1

         If ( proc_index <= Num_Procs_Col ) Then
            Start_lat(proc_index-1) = 1
            End_lat(proc_index-1) = Lat_on_proc(tmp_index)
            Start_zonal(proc_index-1) = 1
            End_zonal(proc_index-1) = Zonal_on_proc(tmp_index)
         Else
            Start_lat(proc_index-1) = End_lat(proc_index-Num_Procs_Col-1) 
     &                            + 1
            End_lat(proc_index-1)   = End_lat(proc_index-Num_Procs_Col-1) 
     &                            + Lat_on_proc(tmp_index)
            Start_zonal(proc_index-1) = End_zonal(proc_index-Num_Procs_Col-1) 
     &                            + 1
            End_zonal(proc_index-1)   = End_zonal(proc_index-Num_Procs_Col-1) 
     &                            + Zonal_on_proc(tmp_index)
         End If

      End Do

      ! Allocate Pre_FFT_Array (based on grid on this processor)
      !   FFT Array will have all lon and subset of lat / Z
      Allocate (Pre_FFT_Array(End_lat(Rank)-Start_lat(Rank)+1,lon,
     &                    End_Z(Rank)-Start_Z(Rank)+1),
     &           Stat = Astat )
      If ( Astat /= 0 ) Then
         Msg = 'Error allocating Pre-FFT_Array in '
         Call alloc_status(Msg, ProcessName)
      End If

      If ( rank == dRank  ) then
         write(*,*) 'Rank Processor = ', Rank
         write(*,*) 'grd shape = ',
     &               shape(Full_Grid(Start_lat(Rank):End_lat(Rank),
     &                               Start_lon(Rank):End_lon(Rank),:))
      end if

      ! Allocate Grid_Array (based on grid on this processor)
      !   Grid Array will have all Z and subset of lon / lat
      Allocate (Grid_Array(Start_lat(Rank):End_lat(Rank),
     &                           Start_lon(Rank):End_lon(Rank),Z),
     &           Stat = Astat )
      If ( Astat /= 0 ) Then
         Msg = 'Error allocating Grid_Array in '
         Call alloc_status(Msg, ProcessName)
      End If

      Grid_Array = Full_Grid(Start_lat(Rank):End_lat(Rank),
     &                           Start_lon(Rank):End_lon(Rank),:)
      if ( rank == drank ) then
         write(*,*) 'max(grd) = ', maxval(Grid_Array(:,:,1))
      end if

      If ( Rank == dRank ) then
         write(*,*) 'Full_grid(1,34:40,1) = ', Full_Grid(1,34:40,1)
      End if

      ! Transform the data from Grid space to FFT space
      Call Grid_to_FFT(Full_Grid(Start_lat(Rank):End_lat(Rank), 
     &                           Start_lon(Rank):End_lon(Rank),:),
     &                 Pre_FFT_Array)

      ! Data is now in a form where FFT can be performed.
      ! Call FFT(...)
      ! Data should now be in (lat,m,Z), i.e. lon converted to m

      ! Allocate FFT_Array (based on grid on this processor)
      !   FFT Array will have all lon and subset of lat / Z
      Allocate (FFT_Array(End_lat(Rank)-Start_lat(Rank)+1,lon/2,
     &                    End_Z(Rank)-Start_Z(Rank)+1),
     &           Stat = Astat )
      If ( Astat /= 0 ) Then
         Msg = 'Error allocating FFT_Array in '
         Call alloc_status(Msg, ProcessName)
      End If

! FFT to Legendre

      ! Allocate Pre_Legendre_Array (based on grid on this processor)
      !   Pre_Legendre Array will have all lat and subset of lon / Z
      Allocate (Pre_Legendre_Array(lat,End_zonal(Rank)-Start_zonal(Rank)+1,
     &                    End_Z(Rank)-Start_Z(Rank)+1),
     &           Stat = Astat )
      If ( Astat /= 0 ) Then
         Msg = 'Error allocating Pre_Legendre_Array in '
         Call alloc_status(Msg, ProcessName)
      End If

      ! Transform the data from FFT space to Legendre space
      Call FFT_to_Legendre(FFT_Array, Pre_Legendre_Array)

      ! Data is now in a form where Legendre Transform can be performed.
      ! Call Legendre(...)
      ! Data should now be in (n,m,Z), i.e. lat converted to n
      Allocate(Legendre_Array((lon / 3) - 1,
     &                        End_zonal(Rank)-Start_zonal(Rank)+1,
     &                        End_Z(Rank)-Start_Z(Rank)+1),
     &           Stat = Astat )
      If ( Astat /= 0 ) Then
         Msg = 'Error allocating Legendre_Array in '
         Call alloc_status(Msg, ProcessName)
      End If

! Legendre to Spectral

      ! Allocate Spectral_Array (based on grid on this processor)
      !   Spectral Array will have all Z and subset of lat / lon
      Allocate (Spectral_Array(End_total(Rank) - 
     &                    Start_total(Rank)+1,
     &                    End_zonal(Rank)-Start_zonal(Rank)+1,
     &                    Z),
     &           Stat = Astat )
      If ( Astat /= 0 ) Then
         Msg = 'Error allocating Spectral_Array in '
         Call alloc_status(Msg, ProcessName)
      End If

      ! Transform the data from Legendre space to Spectral space
      Call Legendre_to_Spectral(Legendre_Array,Spectral_Array)

! Do spectral calculations

      ! Transform the data from Spectral space back to Legendre space
      Call Spectral_to_Legendre(Spectral_Array,Legendre_Array)

! Perform inverse Legendre (from n, m, z to lat, m, z)

      ! Transform the data from Legendre space back to FFT space
      Call Legendre_to_FFT(Pre_Legendre_Array,FFT_Array)

! Perform inverse FFT (from lat, m, z to lat, lon, z)

      ! Transform the data from FFT space back to Grid space
      Call FFT_to_Grid(Pre_FFT_Array,Grid_Array)

! Perform grid processes

! Finalize MPI Processes

      ! Tell the MPI library to realease all resources it is using
      Call MPI_FINALIZE(ierr)
      If (ierr /= MPI_SUCCESS) Then
         Write(*,*) 'Error Finalizing MPI.  Terminating'
         Call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      End If

      End Program
