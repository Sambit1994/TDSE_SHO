!--------------------------------------------------------------------!
! Input from file 'mat.dat' (matrix to be diagonalized):             !
! line 1      : order of the symmetric matrix [M]                    !
! lines 2-n+1 : rows of the matrix                                   !
!--------------------------------------------------------------------!
! Output in file 'Diag_mat.dat':                                          !
! - eigenvalues                                                      !
! - eigenvectors (diagonalizing matrix [D])                          !
! - the original matrix [M] transformed by [D]; [1/D][M][D]          !
!--------------------------------------------------------------------!
! Matrix diagonalization uses LAPACK routine dsyev.f interfaced with !
! diasym.f90.                                                        !
!--------------------------------------------------------------------!


 program diagonalize
!!!!!!!!!!!!!!!!!!!!!
 implicit none

 integer :: i,n
 real(8), allocatable :: m0(:,:),m1(:,:),m2(:,:),eig(:),V1(:,:),V2(:,:),pdt(:,:)
 real                 :: result, tarray(2)

 open(10,file='mat.dat',status='old')
 read(10,*)n
 print *, n
 allocate (m0(n,n))
 allocate (m1(n,n))
 allocate (m2(n,n))
 allocate (eig(n))
 allocate (V1(n,n), V2(n,n), pdt(n,n)) !pdt(n))
 do i=1,n
    read(10,*)m0(:,i)
 enddo
 close(10)

 m1(:,:)=m0(:,:)
! print *, m1


 call diasym(m1,eig,n)
   
 open(10,file='Diag_mat.dat',status='replace')
 write(10,*)'Eigenvalues:'
 do i=1,n
    write(10,10)i,eig(i)
    10 format(I3,'   ',f14.8)
    !print *, i, eig(i)
 enddo
 write(10,*)
 write(10,*)'Eigenvectors:'
 do i=1,n
    write(10,20)i,m1(:,i)
    20 format(i3,'   ',10000f14.8)
 enddo

! print *, 'Check1'
! print *, V1
 V1 = m1
 V2 = transpose(V1)

! print *, 'Check2'
! print *, V2

 pdt = matmul(V1, V2) !! Checking the orthonormality

! print *, 'Check3'
! print *, pdt

 write(10,*)
 
 m2=matmul(transpose(m1),m0) !!! since the eigen vector matrix is orthogonal, it's inverse(of orthogonal matrix) is same as its' transpose
 m0=matmul(m2,m1)

 write(10,*)'Transformed matrix (check):'
 do i=1,n
    write(10,30)m0(:,i)
    30 format(10000f14.8)
 enddo
 write(10,*)

 write(10,*)'Orthonormality of Eigen vectors:'
 do i=1,n
    write(10,40)pdt(:,i)
    40 format(10000f14.8)
 enddo
 write(10,*)

 close(10)

 deallocate(m0); deallocate(m1); deallocate(m2);  deallocate(eig)
 deallocate(V1); deallocate(V2); deallocate(pdt) 

 print *, 'Execution Time in Sec'

 call ETIME(tarray, result)
 print *, result  !! Run time since start in seconds.
 print *, tarray(1) !! User time in seconds.
 print *, tarray(2) !! System time in seconds.

 end program diagonalize
!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
 subroutine diasym(a,eig,n)
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!!!!!!!!!!!!!!!!!!!!!!!!
