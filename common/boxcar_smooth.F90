!-------------------------------------------------------------------------------
! Name: boxcar_smooth.F90
!
! Purpose:
! Subroutine which calculates a boxcar (i.e. moving window) average on a 2D
! array. Can be used to effectively reduce the resolution of satellite imagery,
! without degrading the number of pixels (useful for dealing with instruments
! where the collocation of different channels/views aren't collocated to the
! resolution of the imagery.
!
! Arguments:
! Name           Type    In/Out/Both Description
! width          integer In          Half width of the averaging window. Average
!                                    for index i will be calculated over the
!                                    range i-width,i+width. Same width is used
!                                    in both dimensions
! data           sreal   Both        Data array to be smoothed
!
! History:
! 2023/07/21, GT: Initial version
!
! Bugs:
! None known
!-------------------------------------------------------------------------------
  subroutine boxcar_smooth(width, nx, ny, data, missing, datasq)
    use common_constants_m

    implicit none
    
    integer(kind=lint), intent(in)     :: width
    integer,            intent(in)     :: nx, ny
    real(kind=sreal),   intent(inout)  :: data(nx,ny)
    real(kind=sreal),   intent(in)     :: missing
    real(kind=sreal),   intent(out)    :: datasq(nx,ny)
    
    integer(kind=lint)                 :: nsx, nsy
    integer(kind=lint)                 :: i, j, fw
    real(kind=sreal), allocatable      :: p1(:,:), p2(:,:)
    real(kind=sreal), allocatable      :: ps1(:,:), ps2(:,:)
    real(kind=sreal), allocatable      :: c1(:,:), c2(:,:)
    logical                            :: mask(2*width+1)
   
    ! Calculate full-width of the boxcar
    fw = 2*width + 1

    ! Define smooth output sizes
    nsx = nx - fw + 1
    nsy = ny - fw + 1
    
    allocate(p1(nx,nsy))
    allocate(p2(nsx,nsy))
    allocate(c1(nx,nsy))
    allocate(c2(nsx,nsy))
    allocate(ps1(nx,nsy))
    allocate(ps2(nsx,nsy))

    ! Pre calculate the data^2 array for efficiency
    datasq(:,:) = data(:,:)**2

    ! Do the partial sums (this is the most numerically efficient way of doing
    ! this I've found) (Order nsx*nsy operations)
    do i = 1, nx
       mask = data(i, 1:fw) .ne. missing
       p1(i, 1) = sum(data(i, 1:fw), MASK=mask)
       c1(i, 1) = count(mask, KIND=sreal)
       ps1(i, 1) = sum(datasq(i, 1:fw), MASK=mask)
       do j = 2, nsy
          p1(i, j) = p1(i, j-1)
          c1(i, j) = c1(i, j-1)
          ps1(i, j) = ps1(i, j-1)
          if (data(i, j-1) .ne. missing) then
             p1(i,j) = p1(i,j) - data(i, j-1)
             c1(i,j) = c1(i,j) - 1
             ps1(i,j) = ps1(i,j) - datasq(i, j-1)
          end if
          if (data(i, j+fw-1) .ne. missing) then
             p1(i,j) = p1(i,j) + data(i, j+fw-1)
             c1(i,j) = c1(i,j) + 1
             ps1(i,j) = ps1(i,j) + datasq(i, j+fw-1)
          end if
       end do
    end do
    do j = 1, nsy
       mask = p1(1:fw, j) .ne. missing
       p2(1, j) = sum(p1(1:fw, j), MASK=mask)
       !c2(1, j) = count(mask, KIND=sreal)
       c2(1, j) = sum(c1(1:fw, j))
       ps2(1, j) = sum(ps1(1:fw, j), MASK=mask)
       do i = 2, nsx
          p2(i, j) = p2(i-1, j)
          c2(i, j) = c2(i-1, j) - c1(i-1, j) + c1(i+fw-1, j)
          ps2(i, j) = ps2(i-1, j)
          if (p1(i-1, j) .ne. missing) then
             p2(i, j) = p2(i, j) - p1(i-1, j)
             c2(i, j) = c2(i, J) - c1(i-1, j)
             ps2(i, j) = ps2(i, j) - ps1(i-1, j)
          end if
          if (p1(i+fw-1, j) .ne. missing) then
             p2(i, j) = p2(i, j) + p1(i+fw-1, j)
             c2(i, j) = c2(i, j) + c1(i+fw-1, j)
             ps2(i, j) = ps2(i, j) + ps1(i+fw-1, j)
          end if
       end do
    end do
    c2(:,:) = c2(:,:) / 2.0
    ! Calculate standard deviation
    datasq(:,:) = sreal_fill_value
    ps2(:,:) = sqrt(c2(:,:)*ps2(:,:) - p2(:,:)**2) / c2(:,:)
    datasq(width+1:width+nsx, width+1:width+nsy) = ps2(:,:)
         
    ! Calculate the mean
    data(:,:) = sreal_fill_value
    p2(:,:) = p2(:,:) / c2(:,:)
    data(width+1:width+nsx, width+1:width+nsy) = p2(:,:)
    
    deallocate(p1)
    deallocate(p2)
    deallocate(c1)
    deallocate(c2)
    deallocate(ps1)
    deallocate(ps2)

  end subroutine  boxcar_smooth
