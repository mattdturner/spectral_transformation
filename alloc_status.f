      Subroutine alloc_status(Msg, ProcessName)

!----------------------------------------------------------
! Purpose:
!    Checks the status of array allocations and kills the
!       simulation if the array failed to allocate.

! Revision History:
!   Aug 2015 Matthew Turner: Created.  Need to adapt
!                            for thinned grids.
!----------------------------------------------------------

! Module Files:
      ! Currently empty

      Implicit None

! Include Files:
      ! Currently Empty

! Passed variables:
      Character(80), Intent(In) :: Msg
      Character(24), Intent(In) :: ProcessName

!------------------------------ Code Begins ---------------------------

      Write(*,*) Msg, ProcessName
      Write(*,*) Msg, ProcessName
      Write(*,*) Msg, ProcessName
      Call Exit

      End Subroutine alloc_status
