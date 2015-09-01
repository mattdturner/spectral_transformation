      Module Grid_Config

      Implicit None

      Integer, Allocatable :: Start_lat(:)   ! Start latitude for
                                             !   processor
      Integer, Allocatable :: Start_lon(:)   ! Start longitude for
                                             !   processor
      Integer, Allocatable :: End_lat(:)     ! End latitude for
                                             !   processor
      Integer, Allocatable :: End_lon(:)     ! End longitude for
                                             !   processor
      Integer, Allocatable :: Start_Z(:)     ! Start Z for processor
      Integer, Allocatable :: End_Z(:)       ! End Z for processor
      Integer, Allocatable :: Start_leg(:)   ! Start longitude for
                                             !   processor for legendre
      Integer, Allocatable :: End_leg(:)     ! End longitude for
                                             !   processor for legendre
      Integer, Allocatable :: Start_lat_spectral(:)   ! Start latitude for
                                             !   processor for spectral
      Integer, Allocatable :: End_lat_spectral(:)     ! End latitude for
                                             !   processor for spectral
      Integer, Allocatable :: Lat_on_proc(:) ! # of latitude cells on
                                             !   each processor
      Integer, Allocatable :: Lon_on_proc(:) ! # of longitude cells on
                                             !   each processor
      Integer, Allocatable :: Z_on_proc(:)   ! # of Z cells on each 
                                             !   processor
      Integer, Allocatable :: Leg_on_proc(:) ! # of longitude cells on
                                             !   each processor for legendre
      Integer, Allocatable :: Spectral_Lat_on_proc(:) ! # of latitude
                                             !   cells on each 
                                             !   processor for spectral

      Integer :: Num_Procs_Col ! Number of processor columns
      Integer :: Num_Procs_Row ! Number of processor rows
      Integer :: Rank  ! MPI Processor Rank
      Integer :: dRank = 0  ! MPI Processor Rank for debug statements
      Integer :: Num_Procs  ! Number of MPI Processors
      Integer :: lat, lon, Z   ! Grid configuration variables

      Contains

         Function Grid_Init() Result (Success)

         Implicit None

         Integer :: Astat
         Character(80) :: Msg
         Character(24) :: ProcessName = 'Grid_Init'
         Logical :: Success

         ! This function is called only once at startup, and allocates the
         ! grid configuration and transformation variables
         Allocate (Lat_on_proc(Num_Procs_Row), Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Lat_on_proc in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Lon_on_proc(Num_Procs_Col), Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Lon_on_proc in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Leg_on_proc(Num_Procs_Row), Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Leg_on_proc in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Z_on_proc(Num_Procs_Col), Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Z_on_proc in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Spectral_Lat_on_proc(Num_Procs_Col), Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Spectral_Lat_on_proc in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Start_lat(0:Num_Procs_Col*Num_Procs_Row-1),
     &              Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Start_lat in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Start_lat_spectral(0:Num_Procs_Col*Num_Procs_Row-1),
     &              Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Start_lat_spectral in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Start_lon(0:Num_Procs_Col*Num_Procs_Row-1),
     &              Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Start_lon in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Start_leg(0:Num_Procs_Col*Num_Procs_Row-1),
     &              Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Start_leg in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (Start_Z(0:Num_Procs_Col*Num_Procs_Row-1),
     &              Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating Start_Z in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (End_lat(0:Num_Procs_Col*Num_Procs_Row-1), 
     &             Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating End_lat in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (End_lat_spectral(0:Num_Procs_Col*Num_Procs_Row-1),
     &                               Stat = Astat)
         If (Astat /= 0) Then
            Msg = 'Error allocating End_lat_spectral in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (End_lon(0:Num_Procs_Col*Num_Procs_Row-1),
     &             Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating End_lon in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (End_leg(0:Num_Procs_Col*Num_Procs_Row-1),
     &              Stat = Astat )
         If (Astat /= 0) Then
            Msg = 'Error allocating End_leg in '
            Call alloc_status(Msg, ProcessName)
         End If
         Allocate (End_Z(0:Num_Procs_Col*Num_Procs_Row-1), Stat = Astat)
         If (Astat /= 0) Then
            Msg = 'Error allocating End_Z in '
            Call alloc_status(Msg, ProcessName)
         End If

         Success = .True.
 
         Return
         End Function Grid_init

      End Module Grid_Config
