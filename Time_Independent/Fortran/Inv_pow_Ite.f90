!-------------------------------------------------!
! Script uses inverse power iteration with shift  !
! procedue to converge to the nearest eigen value !
!-------------------------------------------------!
!Input: mat.dat                                   !
!Output: dia_IPI.dat                              !
!-------------------------------------------------!
!pd_target: Targeted eigen value                  !
!shift_val: Shift or guess eigen value            !
!-------------------------------------------------!

program main
    implicit none (type, external)
    external :: sgesv, dnrm2
    real,allocatable  :: a1(:,:), b1(:), pivot1(:), b2(:), a11(:,:), la(:), Im(:,:), a12(:,:) 
    real     :: norm, pd, pd_target, pd_diff, shift_val !a(3,3)
    real     :: result, tarray(2)
    integer  :: rc, n, n1, i, j, k     


    open(10,file='mat.dat',status='old')
    open(101,file='Dia_IPI.dat')
    read(10,*)n1
    allocate (a1(n1,n1), a11(n1,n1), a12(n1,n1), b1(n1), pivot1(n1), la(n1), Im(n1,n1))
    do i=1,n1
       read(10,*)a1(:,i)
    enddo
    close(10)

    Im = 0.0d0

    do i = 1, n1
      Im(i,i) = 1.0d0
    enddo

    pd_target = 1.5d0

    shift_val = 1.3d0 !! shift value 

    Im = shift_val*Im

    !print *, Im
    !print *, 'Im_matrix'

    a1 = transpose(a1) !! Tranpose is done to convert the matrix to the correct format executable by sgesv 
    a11 = a1           !! However it does not matter here, as this is a symmetric matrix

    !print *, 'check'
    !print *, a11
    a12 = a1

    !print *, 'A reduced/shifted'
    a11 = a11 - Im ! The shifted matrix
    !print *, (a11)

    b1 = 1

    !print *, 'loop start'
    iloop: do i = 1, 100
       a1 = a11
    
       call sgesv(n1, 1, a1, n1, pivot1, b1, n1, rc)
       if (rc /= 0) then
           print '(a, i0)', 'Error: ', rc
           stop
       end if
   
       norm = 0.
       do j = 1, n1
       !    print *, b1(j)
           norm = norm + (b1(j)**2)
       enddo
       norm = sqrt(norm)
   
       !print *, norm
       !print *, 'normalised'
       b2 = b1 * (1/norm) !matmul(b1, norm)
       !print *, (b2)

       b1 = b2 

!! Here while finding the eigen value, original A-matrix (a1) is used not (A - l.I)(the shifted matrix), because the obtained eigen vector corresponds to the 
!! eigen value of shift, and will be obtained using the original A matrix

       la = matmul(a12, b1)
       pd = 0.00d0

       do k = 1, n1
          pd = pd + la(k) * b1(k)
       enddo

       !print *, 'see pd'
       !print *, i, ' ', (pd)
       write(101,10)i, pd
       10 format(I6,'   ',f14.8)

    enddo iloop

    print *, 'Converged eigenvalue at loop end'
    print *, pd

    deallocate(a1, a11, a12, b1, pivot1, la, Im)

    print *, 'Execution Time in Sec'
   
    call ETIME(tarray, result)
    print *, result  !! Run time since start in seconds.
    print *, tarray(1) !! User time in seconds.
    print *, tarray(2) !! System time in seconds.


end program main
