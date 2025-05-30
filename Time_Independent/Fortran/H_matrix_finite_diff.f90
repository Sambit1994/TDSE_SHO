!------------------------------------------------------------------------!
! Script uses                                                            !
!------------------------------------------------------------------------!
! span_i, span_f: Edges of box                                           !
! h: step size                                                           !
! The H matrix is generated using::                                      !
! Five-point central difference formula for second derivatives -- KE     !
! x**2/2 -- PE                                                           !
! The bump funxtions to potentails are commented out                     !
! The H matrix is written to mat.dat                                     !
!------------------------------------------------------------------------!

program SHO_matrix

implicit none

   integer :: i, j, N_start, N_end, k
   integer :: l, m, N_it 
   real :: numbers(5), V1, V2, V3, V4, V5, V_const 
   real :: matrix(3,3), h, x_val, V_pot, span_i, span_f
   real, dimension(5) :: arr
   double precision, allocatable  :: S(:,:)

   open(unit=101, file='mat.dat')

   span_i = -10.00d0
   span_f = 10.00d0

   h = 0.1d0 !1.0d0 !0.05d0 
   x_val = span_i - h 

   N_it = ((span_f-span_i)/h) + 1
   print *, N_it
   write(101, *) N_it

   N_start = 1
   N_end   = N_it

   V_const = 1.00d0/(24.00d0 * (h**2))

   V1 = 1.00d0 * V_const
   V2 = -16.00d0 * V_const
   V3 = 30.00d0 * V_const 
   V4 = -16.00d0 * V_const 
   V5 = 1.00d0 * V_const 

   allocate(S(1:N_end, 1:N_end))

   S = 0
!   print *, S

   do i = N_start, N_end
      x_val = x_val + h
      !print *, x_val

      !f0 = 2*(exp(-10*(x_val**2)))  
      !f1 = 3*(exp(-15*(x_val**2)))
      !f2 = 2*(exp(-5*(x_val**2)))
      !f3 = 2*(exp(-30*(x_val**2)))
      !f4 = 10*(exp(-30*(x_val**2)))
      !f5 = 10*(exp(-5*(x_val**2)))

      V_pot = ((x_val**2)/2) !+ (10*(exp(-30*(x_val**2))))
      !print *, V_pot

      j = i-2
      k = i-1
      l = i+1
      m = i+2

      if ((k .eq. (N_start-1)) .and. (j .eq. (N_start-2))) then
          S(i,i) = V3 + V_pot
          S(i,l) = V4 
          S(i,m) = V5
      elseif ((k .eq. (N_start)) .and. (j .eq. (N_start-1))) then
          S(i,k) = V2
          S(i,i) = V3 + V_pot 
          S(i,l) = V4
          S(i,m) = V5
      elseif ((l .eq. (N_end)) .and. (m .eq. (N_end+1))) then
          S(i,j) = V1
          S(i,k) = V2
          S(i,i) = V3 + V_pot 
          S(i,l) = V4
      elseif ((l .eq. (N_end+1)) .and. (m .eq. (N_end+2))) then
          S(i,j) = V1
          S(i,k) = V2
          S(i,i) = V3 + V_pot 
      else
          S(i,j) = V1
          S(i,k) = V2
          S(i,i) = V3 + V_pot 
          S(i,l) = V4
          S(i,m) = V5 
      endif

   enddo


   do i = 1, N_end
      !print *, (S(i, j), j = 1, N_end)
      write(101, *) (S(i, j), j = 1, N_end) 
      !print *, 'Check'
   end do

   deallocate(S)


   close(101)

end program
