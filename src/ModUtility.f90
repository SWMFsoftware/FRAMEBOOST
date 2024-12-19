!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModUtilities

  ! Simple methods which are used by CON and can be used
  ! by the components of the SWMF too.
  !
  ! F77 and C++ codes need an F90 interface to access these utilities.

  implicit none

  private ! except

  public:: CON_stop
  public:: linear_scalar
  public:: linear_vector

contains
  !============================================================================
  real function linear_scalar(a_I, iMin, iMax, x, x_I, DoExtrapolate, &
       iCell, Dist)

    ! Interface for default precision real array

    integer, intent(in) :: iMin, iMax
    real,    intent(in) :: a_I(iMin:iMax)
    real,    intent(in), optional:: x, x_I(iMin:), Dist
    logical, intent(in), optional :: DoExtrapolate
    integer, intent(in), optional :: iCell

    integer :: i1, i2
    real    :: Dx1, Dx2

    character(len=*), parameter:: NameSub = 'linear'
    !--------------------------------------------------------------------------
    if ( present(iCell) .and. present(Dist) ) then

       ! Calculate the remaining cell indices and interpolation weights

       i1 = iCell
       i2 = i1 + 1
       Dx1 = Dist
       Dx2 = 1.0 - Dx1

    else

       ! Locate the point Xyz_D on the grid
       call find_cell(iMin, iMax, x, i1, Dx1, x_I, DoExtrapolate, &
            "Called from "//NameSub)

       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1

    end if

    ! Perform interpolation (or extrapolation)
    linear_scalar = Dx2*a_I(i1) + Dx1*a_I(i2)

  end function linear_scalar
  !============================================================================
  function linear_vector(a_VI, nVar, iMin, iMax, x, x_I, x8, x8_I, &
       DoExtrapolate, iCell, Dist)

    ! interpolate default precision real a_VI vector array

    integer, intent(in):: nVar, iMin, iMax
    real,    intent(in):: a_VI(nVar,iMin:iMax)

    real,         intent(in), optional :: x, x_I(iMin:), Dist
    real,         intent(in), optional :: x8, x8_I(iMin:)
    logical,      intent(in), optional :: DoExtrapolate
    integer,      intent(in), optional :: iCell

    ! return value
    real                :: linear_vector(nVar)

    integer :: i1, i2
    real    :: Dx1, Dx2

    character(len=*), parameter:: NameSub = 'linear_vector'
    !--------------------------------------------------------------------------
    if ( present(iCell) .and. present(Dist) ) then

       ! Calculate the remaining cell indices and interpolation weights

       i1 = iCell
       i2 = i1 + 1
       Dx1 = Dist
       Dx2 = 1.0 - Dx1

    else
       ! Locate the point Xyz_D on the grid
       if(present(x8))then
          call find_cell(iMin, iMax, x8, i1, Dx1, x8_I, DoExtrapolate, &
               "Called from "//NameSub)
       else
          call find_cell(iMin, iMax, x, i1, Dx1, x_I, DoExtrapolate, &
               "Called from "//NameSub)
       end if

       ! Calculate the remaining cell indices and interpolation weights
       i2 = i1 + 1; Dx2 = 1.0 - Dx1

    end if

    ! Perform interpolation (or extrapolation) for multiple variables
    linear_vector = Dx2*a_VI(:,i1) + Dx1*a_VI(:,i2)

  end function linear_vector
  !============================================================================
  subroutine CON_stop(String)

    character(len=*), intent(in):: String

    ! Stop execution after the following actions:
    !
    ! Write out error message with processor rank and String.
    !--------------------------------------------------------------------------
    write(*,'(a)') 'ERROR: '//String
    stop
    
  end subroutine CON_stop
  !============================================================================
  subroutine find_cell(MinCoord, MaxCoord, Coord, iCoord, dCoord, &
       Coord_I, DoExtrapolate, StringError, IsInside)

    ! double precision version of the subroutine above

    integer,          intent(in)           :: MinCoord, MaxCoord
    real,             intent(in)           :: Coord
    integer,          intent(out)          :: iCoord
    real,             intent(out), optional:: dCoord
    real,             intent(in),  optional:: Coord_I(MinCoord:)
    logical,          intent(in),  optional:: DoExtrapolate
    character(len=*), intent(in),  optional:: StringError
    logical,          intent(out), optional:: IsInside

    integer:: i, Di, nIter, MaxIter
    real:: Tolerance

    logical:: IsUniform

    character(len=*), parameter:: NameSub = 'find_cell8'
    !--------------------------------------------------------------------------
    if(present(Coord_I)) then
       IsUniform = size(Coord_I) < 2
    else
       IsUniform = .true.
    endif

    if(IsUniform)then
       ! Uniform grid case with normalized coordinate

       iCoord = min(MaxCoord-1, max(MinCoord, floor(Coord)))
       dCoord = Coord - iCoord
       Tolerance = (MaxCoord - MinCoord)*1d-12

       if(Coord < MinCoord - Tolerance)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) NameSub, ': ', StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub// &
                  ': normalized coordinate is too small!')
          elseif(.not.DoExtrapolate)then
             ! Use lefttmost cell (first order accurate)
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
       elseif(Coord > MaxCoord + Tolerance)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub// &
                  ': normalized coordinate is too large!')
          elseif(.not.DoExtrapolate)then
             ! Use rightmost cell (first order accurate)
             dCoord = 1.0
          endif
          if(present(IsInside)) IsInside = .false.
       else
          if(present(IsInside)) IsInside = .true.
       end if

    elseif(Coord_I(MinCoord) < Coord_I(MaxCoord))then

       ! Monotone increasing coordinates
       Tolerance = (Coord_I(MaxCoord) - Coord_I(MinCoord))*1d-12

       if(present(IsInside)) IsInside = .true.

       if(Coord < Coord_I(MinCoord) - Tolerance)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord - Coord_I(iCoord)) &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       elseif(Coord <= Coord_I(MinCoord))then
          ! Very near the edge
          iCoord = MinCoord
          dCoord = 0.0
          RETURN
       end if

       if(Coord > Coord_I(MaxCoord) + Tolerance)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord - Coord_I(iCoord))  &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       elseif(Coord >= Coord_I(MaxCoord))then
          ! Very near the edge
          iCoord = MaxCoord - 1
          dCoord = 1.0
          RETURN
       end if

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       nIter   = 0
       MaxIter = MaxCoord - MinCoord ! Even a linear search would end
       do
          Di = (Di + 1)/2
          if(Coord < Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord > Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
          nIter = nIter + 1
          if (nIter > MaxIter) then
             write(*,*) NameSub,': ERROR in monotone increasing: '
             write(*,*) 'Tolerance=', Tolerance
             write(*,*) 'Coord_I(MinCoord), Coord_I(MaxCoord), Coord=', &
                  Coord_I(MinCoord), Coord_I(MaxCoord), Coord
             write(*,*) 'i, Di, Coord_I(i:i+1)=', i, Di, Coord_I(i:i+1)
             call CON_stop(NameSub//': maximum iteration exceeded!')
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord             - Coord_I(iCoord)) &
               /   (Coord_I(iCoord+1) - Coord_I(iCoord))
       end if
    else

       ! Monotone decreasing coordinates
       Tolerance = (Coord_I(MinCoord) - Coord_I(MaxCoord))*1d-12

       if(Coord < Coord_I(MaxCoord) - Tolerance)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord_I(iCoord) - Coord) &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MinCoord) + Tolerance)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord_I(iCoord) - Coord)  &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       nIter = 0
       MaxIter = MaxCoord - MinCoord ! Even a linear search ends by this
       do
          Di = (Di + 1)/2
          if(Coord > Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord < Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
          nIter = nIter + 1
          if (nIter > MaxIter) then
             write(*,*) NameSub,': ERROR in monotone decreasing: '
             write(*,*) 'Tolerance=', Tolerance
             write(*,*) 'Coord_I(MinCoord), Coord_I(MaxCoord), Coord=', &
                  Coord_I(MinCoord), Coord_I(MaxCoord), Coord
             write(*,*) 'i, Di, Coord_I(i:i+1)=', i, Di, Coord_I(i:i+1)
             call CON_stop(NameSub//': maximum iteration exceeded!')
          end if
       end do
       iCoord = i
       if(Coord_I(iCoord+1) == Coord_I(iCoord))then
          dCoord = 0.0
       else
          dCoord = (Coord_I(iCoord) - Coord  ) &
               /   (Coord_I(iCoord) - Coord_I(iCoord+1))
       end if
    end if

  end subroutine find_cell
  !============================================================================
end module ModUtilities
