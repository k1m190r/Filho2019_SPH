program flat_plate_diffusion

!-----------------------------------------------------------------
! these fortran source files have been used in heat diffusion in a
! homogeneous flat plate (sect. 4.2).
!-----------------------------------------------------------------
! the reader should follow the instructions below in order to
! create the program executable:
! 1. within the working directory, create and open a new file in
! the fortran editor˜
! 2. copy and paste this program and save as ''name.f90''
! 3. create a folder called output within the working directory
! 4. create a subfolder called temperature within the folder output
! 5. compile and run the program
! output: in the file called temperature_differences.dat are the
! points where the smallest and the largest temperature differences
! occur (comparing the sph results and solution provided by series).
!
!-----------------------------------------------------------------
! variables:
! aux_p - variable used in memory allocation
! dist - distance between the centres of mass of a fixed and a
!     neighbour particle
! dt - time step (in seconds)
! dx - horizontal distance between the centres of mass of
!     particles
! dy - vertical distance between the centres of mass of
!     particles
! eps - accuracy of simulation
! finish - variable used to record the simulation time
! hsml - smoothing lengths of particles
! int_function - interpolation function (kernel)
! itimestep - current iteration
! lx - length of the plate
! ly - width of the plate
! mass - mass of particles
! n_bound - number of boundary particles
! n_cells - number of cells in each cartesian direction
! neighbour - matrix containing the neighbours of each particle
!     of the domain
! n_global - sum of the number of particles at domain and
!     boundary particles
! np_side - number of particles per side of the plate
! ntotal - total number of particles at domain
! response - variable used to verify the achievement of accuracy
! rho - densities of particles
! start - variable used to record the simulation time
! step_out - step for recording results
! support_radius - radius of the domain of influence
! time - simulation time
! t_0 - initial temperatures of particles
! t - temperature provided by sph method at each numerical
!     iteration
! te - prescribed temperature at right side of the plate
! tn - prescribed temperature at top of the plate
! ts - prescribed temperature at bottom of the plate
! tw - prescribed temperature at left side of the plate
!     temp_laplacian - sph laplacian of temperature
! x - coordinates of particles
!-----------------------------------------------------------------

   implicit none
   integer, dimension (:,:), allocatable :: neighbour
   double precision, dimension (:), allocatable :: &
      rho, mass, hsml, t_0, t, temp_laplacian
   double precision, dimension (:,:), allocatable :: &
      x, dist, dx, dy
   integer int_function, aux_p, ntotal, itimestep, d, m, i, &
      hours, minutes, hours_seg, np_side, n_global, n_cells, n_bound, step_out
   double precision ts, tn, tw, te, lx, ly, xl, yl, factor
   double precision start, finish, time, seconds, dt, eps, support_radius
   character (len=15) response

   call cpu_time (start)

   open (01, file="output/domain_particles_positions.dat")
   open (02, file="output/boundary_particles_positions.dat")
   open (03, file="output/initial_temperature.dat")

!------------------------plate geometry---------------------------
   lx = 1.0
   ly = 1.0

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!---------------defining the temperatures at boundaries-----------
!------------------(dirichlet boundary conditions)----------------
!-----------------------------------------------------------------

   ts = 100.0 !temperature at the lower boundary
   tn = 0.0 !temperature at the upper side boundary
   tw = 0.0 !temperature at the left side boundary
   te = 0.0 !temperature at the right side boundary

!-------------------time step and accuracy------------------------
   dt = 1.e-5
   eps = 1.e-6

!-----------------------------------------------------------------
   print *, ' '
   print *, 'enter with the number of particles per side of the plate:'
   print *, ' '
   read  *, np_side
   print *, ' '
   print *, 'choose the interpolation function (kernel):'
   print *, ' '
   print *, '1 - lucy''s quartic kernel'
   print *, '2 - cubic spline kernel'
   print *, '3 - new quartic kernel'
   print *, '4 - quintic spline kernel'
   print *, ' '
   read  *, int_function
   print *, ' '
   print *, 'enter with the step for the output files:'
   print *, ' '
   read  *, step_out

   !distance between the centres of mass

   aux_p = np_side**2 + 1500 !total number of particles (for memory allocation)
   allocate ( &
      neighbour(aux_p,150), rho(aux_p), hsml(aux_p), t_0(aux_p), t(aux_p), &
      temp_laplacian(aux_p), x(aux_p,2), dist(aux_p,150), &
      dx(aux_p,150), dy(aux_p,150), mass(aux_p) )

   response = 'non_convergence'
   itimestep = 1

   print *, ' '
   print *,'-------------------------------------------------------'
   print *,'--flat plane diffusion simulation using the sph method-'
   print *,'-------------------------------------------------------'
   print *, ' '

   dx = lx/np_side
   dy = ly/np_side

   call temperature_series ( aux_p, np_side, ts, tn, tw, te, dx, dy, dt )

   call input ( &
      aux_p, x, mass, rho, np_side, ntotal, t_0, lx, ly, support_radius, n_cells )

   call boundary_particles ( &
      aux_p, itimestep, np_side, ts, tn, tw, te, &
      lx, ly, dx, dy, n_bound, ntotal, n_global, mass,x,rho,t_0 )

   call neighbour_search ( &
      aux_p, x, n_global, ntotal, n_bound, neighbour, dx, dy, support_radius, dist )

!-----------------------------------------------------------------
!--------------------until the accuracy is achieved---------------
!-----------------------------------------------------------------

   do while (response.eq.'non_convergence')
      call sph_temperature_laplacian ( &
         aux_p, int_function, factor, hsml, ntotal, dist, mass, rho, &
         support_radius, neighbour, t_0, temp_laplacian, itimestep )

      call temporal_integration ( &
         aux_p, itimestep, x, ntotal, dt, temp_laplacian, t_0, t, step_out )

      call convergence ( aux_p, itimestep,t_0,t,eps,response, ntotal )

      call temperature_differences ( itimestep, dt, ntotal, step_out )

      if (response.eq.'convergence') then
         print *, ' accuracy achieved in the ', itimestep, 'a. iteration'
         print *, ' end of simulation '
      end if

   end do

   call cpu_time(finish)

   time = finish-start
   hours_seg = 3600; !hours em seconds
   hours = (time/hours_seg) !resultado da hora
   minutes = (time -(hours_seg*hours))/60
   seconds = (time -(hours_seg*hours)-(minutes*60))

   print *,' '
   print *, ' -----------------------------------------------------'
   write (*,*) 'time processing (cpu time) '
   write (*,*) hours ,' hours', minutes,' minutes', seconds, ' seconds.'
   print *, ' -----------------------------------------------------'

   deallocate (neighbour, rho, hsml, t_0, t, temp_laplacian, x, dist, dx, dy, mass)

   open (776,file= 'output/simulation_parameters.dat')
   rewind (776)
   write (776,*) ' '
   if (int_function.eq.1) then !--------lucy's quartic kernel----
      write (776,*) 'interpolation function: lucy''s quartic kernel'
      write (776,*) ' '
   else
      if (int_function.eq.2) then !------cubic spline kernel-------
         write (776,*) 'interpolation function: cubic spline kernel'
         write (776,*) ' '
      else
         if (int_function.eq.3) then !----new quartic kernel--------
            write (776,*) 'interpolation function: new quartic kernel'
            write (776,*) ' '
         else
            if (int_function.eq.4) then !-----quintic spline kernel---
               write (776,*) 'interpolation function: quintic spline kernel'
               write (776,*) ' '
            end if
         end if
      end if
   end if

   open  (25,file='number_of_iterations.dat')
   write (25,*) itimestep
   close (25)

   write (776,*) 'number of particles per side of plate = ', np_side
   write (776,*) ' '
   write (776,*) 'time step (dt, em seconds) = ', dt
   write (776,*) ' '
   write (776,*) 'accuracy ( abs(t(m+1) - t(m)) ) = ', eps
   write (776,*) ' '
   write (776,*) 'time for the physical diffusion =', itimestep*dt
   write (776,*) ' '
   write (776,*) 'time processing (cpu time) = ', &
      hours ,' hours', minutes,' minutes', seconds, ' seconds.'

   close (01)
   close (02)
   close (03)
   close (776)

end program


!******************** solution by series *************************

subroutine temperature_series(aux_p, np_side, ts, tn, tw, te, dx, dy, dt)

!-----------------------------------------------------------------
! this subroutine defines the positions of the centres of mass of
! particles at domain and makes the calculation of the
! temperatures at steady-state using the solution by series
!
! x - coordinates of particles
! mass - mass of particles
! kmax - number of terms in series
!-----------------------------------------------------------------

   implicit none

   double precision, dimension (:,:), allocatable :: t , x1
   double precision, dimension (:), allocatable :: posx, posy, temp
   double precision lx, ly, dx, dy, aux, ts, tn, tw, te, x, y, sinx, txy, pi, dt
   integer aux_p, np_side, i, d, j, k, mp1, np1, kmax, cont, n_bound

   mp1 = np_side+1
   np1 = np_side+1

   allocate ( t(mp1,np1), posx(aux_p), posy(aux_p), temp(aux_p), x1(1500,3) )

   pi = acos(-1.)
   kmax = 90

   open (unit=10, file='output/temperatures.dat')
   open (unit=11, file='output/verifying_temp_series.dat')

!-----------------------------------------------------------------
!definition of positions occupied by centres of mass of particles

   cont = 0
   do i = 1, np_side
      x = (i-1)*dx + dx/2
      do j = 1, np_side
         cont = cont + 1
         y = (j-1)*dy + dy/2
         posx(cont) = x
         posy(cont) = y

         !calculation of temperatures

         txy=0.
         do k = 1, kmax, 2
            sinx = sin(k*pi*x)
            aux = ((sinh(k*pi*(y-1.))) /(sinh(k*pi)))
            txy = txy + ((-4.*ts)/(k*pi))* sinx * aux
         end do

         temp(cont)=txy

      end do

   end do

   do i = 1,cont
      write (11,23) i, posx(i), posy(i), temp(i)
   end do

23 format (1x,i10, 1x, d21.14, 1x, d21.14, 1x, d21.14)
   close (11)
   open  (23, file='n_part_series.dat')
   write (23, *) cont
   write (23, *) dt
   close (23)

   !defining a line of boundary particles and temperatures
   n_bound = 0
   j = 1
   do i = 1, mp1
      n_bound = n_bound+1
      x1(n_bound, 1) = (i-1)*dx
      x1(n_bound, 2) = 0.
      x1(n_bound, 3) = ts !lower boundary
   end do

   j = np1
   do i = 1,mp1
      n_bound = n_bound+1
      x1(n_bound, 1) = (i-1)*dx
      x1(n_bound, 2) = (np1-1)*dy
      x1(n_bound, 3)=tn !upper boundary
   end do

   i = 1
   do j=2, np_side
      n_bound = n_bound+1
      x1(n_bound, 1) = 0.
      x1(n_bound, 2) = (j-1)*dy
      x1(n_bound, 3) = tw !left boundary
   end do

   i = mp1
   do j=2, np_side
      n_bound=n_bound+1
      x1(n_bound, 1) = (mp1-1)*dx
      x1(n_bound, 2) = (j-1)*dy
      x1(n_bound, 3) = te !right boundary
   end do

   open (12, file='output/boundary_particles.dat')
   do i=1, n_bound
      write (12,23) i, (x1(i,d),d=1,3)
   end do
   close (12)

end subroutine temperature_series


!*************************** input *******************************

subroutine input(aux_p, x, mass, rho, np_side, ntotal, t_0, lx, ly, &
   support_radius, n_cells)

!-----------------------------------------------------------------
! this subroutine makes the distribution of particles at domain
! and the initial properties are set
!
! x - positions of particles
! mass - mass of particles
! rho - densities of particles
! hsml - smoothing lengths of particles
! ntotal - total number of particles
! ncells - total number of cells in each direction
! t_0 - initial temperature

   implicit none
   integer np_side, aux_p , n_cells, ntotal
   double precision lx, ly, xl, yl, dx, dy, support_radius
   integer i, k, d
   double precision mass(aux_p), rho(aux_p), t_0(aux_p)
   double precision x(aux_p,2), m(aux_p,4)

   n_cells = np_side/5
   yl = ly
   xl = lx
   dx = xl/np_side
   dy = yl/np_side

   open  (23, file='n_part_series.dat')
   read  (23, *) k
   close (23)

   open (unit=11,file='output/verifying_temp_series.dat') !particles at domain
   do i=1,k
      read(11,*) (m(i,d),d=1,4)
   end do

   do i=1,k
      x(i,1) = m(i,2)
      x(i,2) = m(i,3)
   end do

   close (11)

!---------------------defining the support radius ----------------
   support_radius = (xl/real(n_cells))/2.
!-----------------------------------------------------------------

!---------------------defining the properties of particles -------
   do i=1,aux_p
      rho (i) = 1.
      mass(i) = dx*dy*rho(i)
      t_0(i) = 0.
      write (03, 654) t_0(i)
   end do

654 format (1x, d21.14)

   do i=1, k
      write (01, 1013) (x(i,d), d = 1, 2)
   end do

1013 format (2(1x,d21.14))
   ntotal = k

   write (*,*)'---------------------------------------------------'
   write (*,*)' total number of particles at domain : '
   write (*,*) ntotal
   write (*,*)'---------------------------------------------------'
   write (*,*) ' '
   write (*,*) ' press any key to continue... '

   pause

end subroutine input


!************************ boundary particles**********************

subroutine boundary_particles(aux_p, itimestep, np_side, &
   ts, tn, tw, te, lx, ly, dx, dy, n_bound, ntotal, n_global, mass,x,rho,t_0)

!-----------------------------------------------------------------
! this subroutine defines a line of particles at the contour of
! the flat plate and the dirichlet boundary conditions
!
! itimestep - current time step
! ntotal - total number of particles at domain
! n_bound - number of boundary particles
! n_global - sum of the number of particles at domain and
! boundary particles
! hsml - smoothing length
! mass - masses of partic¸ces
! x - positions of particles
! rho - densities of particles
! t_0 - initial temperature
!-----------------------------------------------------------------

   implicit none

   integer aux_p, n_global, np_side, ntotal, n_bound, itimestep
   integer i, j, d, k
   double precision ts, tn, tw, te, lx, ly, dx,dy,max_x, max_y
   double precision :: rho(aux_p), mass(aux_p), hsml(aux_p), t_0(aux_p)
   double precision :: x(aux_p,2)

   n_bound = 0

!--defining a line of particles on the upper side and setting the physical properties
   do i = 1, 2*np_side -1
      n_bound = n_bound + 1
      x(ntotal + n_bound,1) = i*dx/2
      x(ntotal + n_bound,2) = ly
      rho (ntotal + n_bound) = 1.
      mass(ntotal + n_bound) = rho (ntotal + n_bound) * dx * dy
      t_0(ntotal + n_bound)= tn
   end do

!--defining a line of particles on the lower side and setting the physical properties
   do i = 1, 2*np_side -1
      n_bound = n_bound + 1
      x(ntotal + n_bound,1) = i*dx/2
      x(ntotal + n_bound,2) = 0.
      rho (ntotal + n_bound) = 1.
      mass(ntotal + n_bound) = rho (ntotal + n_bound) * dx * dy
      t_0(ntotal + n_bound)= ts
   end do

!--defining a line of particles on the left side and setting the physical properties
   do i = 1, 2*np_side + 1
      n_bound = n_bound + 1

      x(ntotal + n_bound,1) = 0.
      x(ntotal + n_bound,2) = (i-1)*dx/2
      rho (ntotal + n_bound) = 1.
      mass(ntotal + n_bound) = rho (ntotal + n_bound) * dx * dy
      t_0(ntotal + n_bound)= tw
   end do

!--defining a line of particles on the right side and setting the physical properties
   do i = 1, 2*np_side + 1
      n_bound = n_bound + 1
      x(ntotal + n_bound, 1) = lx
      x(ntotal + n_bound, 2) = (i-1)*dx/2
      rho (ntotal + n_bound) = 1.
      mass(ntotal + n_bound) = rho (ntotal + n_bound) * dx * dy
      t_0(ntotal + n_bound)= te
   end do

   max_y= 0.
   max_y =0.

!------------------------ output files----------------------------
   do i=ntotal+1, ntotal+ n_bound

      if (x(i,1).gt.max_x) max_x = x(i,1)
      if (x(i,2).gt.max_y) max_y = x(i,2)

      write (01, 1016) (x(i,d), d = 1, 2)
      n_global = ntotal + n_bound
      write (02, 1016) (x(i-ntotal,d), d = 1, 2)
      write (03, 654) t_0(i)
   end do

1016 format (2(1x, d21.14))
654 format (1x, d21.14)

   open  (17, file='output/geometry.dat')
   write (17, *) 0., max_x
   write (17, *) 0., max_y
   close (17)

end subroutine boundary_particles


!********************** neighbour search *************************

subroutine neighbour_search(aux_p, x, n_global, ntotal, n_bound, neighbour, &
   dx, dy, support_radius, dist)

!-----------------------------------------------------------------
! this subroutine finds the direct search for neighbour particles
! of each fixed particle
!
! ntotal - number of particles at domain
! n_bound - number of boundary particles
! n_global - sum of the number of particles at domain and
! boundary particles
! hsml - smoothing length
! x - positions of all particles
! n_neigh - number of neighbour particles
! dist - distance between a fixed particle and a neighbour particle
! neighbour - matrix of neighbour particles
!-----------------------------------------------------------------

   implicit none

   integer aux_p, n_global, ntotal, n_bound, n_neigh, i, j
   integer neighbour(aux_p,150)
   double precision r, dx_local, dy_local, support_radius
   double precision x(aux_p,2), dist(aux_p,150), dx(aux_p,150), dy(aux_p,150)

   do i=1,aux_p
      do j=1,150
         neighbour(i,j) = 0
         dist(i,j) = 0.
         dx(i,j) = 0.
         dy(i,j) = 0.
      end do
   end do

   do i=1,ntotal
      n_neigh = 1
      neighbour(i,n_neigh) = i

      do j = 1, n_global
         dx_local = x(i,1) - x(j,1)
         dy_local = x(i,2) - x(j,2)
         r = sqrt(dx_local**2 + dy_local**2)

         if (r.le.support_radius) then
            n_neigh = n_neigh +1
            neighbour(i,n_neigh) = j
            dx(i,n_neigh) = dx_local
            dy(i,n_neigh) = dy_local
            dist(i,n_neigh) = r
         end if
      end do

   end do

   open   (222,file='output/neighbouring_particles.dat')
   rewind (222)

   do i = 1, ntotal
      n_neigh = 1
      do while (n_neigh.ne.150)
         write (222, 150), neighbour(i,n_neigh)
         n_neigh = n_neigh +1
      end do
   end do

150 format (i10)
   close   (222)

end subroutine neighbour_search


!****************** laplacian of temperature *********************

subroutine sph_temperature_laplacian(aux_p, int_function, factor, hsml, ntotal, &
   dist, mass, rho, support_radius, neighbour, t_0, temp_laplacian, itimestep)

!-----------------------------------------------------------------
! this subroutine makes calculations of the divergent of the heat
! flux !using a kernel chosen by user
!
! dist - distance between a fixed particle and a neighbouring
! particle
! ntotal - number of particles at domain
! k_scale - scaling factor (depends on kernel)
! factor - normalization constant of kernel
! temp_laplacian - laplacian of temperature provided by sph
! method
! int_function - kernel interpolation
!
! kernel options:
! 1. lucy's quartic kernel
! 2. cubic spline kernel
! 3. new quartic kernel
! 4. quintic spline kernel
!-----------------------------------------------------------------

   implicit none

   integer aux_p, i, j, k, d, ntotal, n_global, n_neigh, int_function, itimestep
   integer neighbour(aux_p,150)
   double precision temp_laplacian (aux_p), aux, q, factor, r, support_radius, k_scale
   double precision rho(aux_p), mass(aux_p), hsml(aux_p), t_0(aux_p)
   double precision dist(aux_p,150), dx(aux_p,150), dy(aux_p,150), pi

   pi = acos(-1.)

   if (itimestep.eq.1) then

      do i=1,aux_p
         hsml(i) = 0.
      end do

      do i=1,ntotal

         if (int_function.eq.1) then !----------lucy's quartic kernel
            k_scale = 1.0
            hsml(i) = support_radius/k_scale
            factor = 5.0/(pi*hsml(i)*hsml(i))
         else
            if (int_function.eq.2) then !--------- cubic spline kernel
               k_scale = 2.0
               hsml(i) = support_radius/k_scale
               factor = 15./(7.*pi*hsml(i)*hsml(i))
            else
               if (int_function.eq.3) then !--------new quartic kernel
                  k_scale = 1.0
                  hsml(i) = support_radius/k_scale
                  factor = 15./ (7.*pi*hsml(i)*hsml(i))
               else
                  if (int_function.eq.4) then !-----quintic spline kernel
                     k_scale = 3.0
                     hsml(i) = support_radius/k_scale
                     factor = 7./(478.*pi*hsml(i)*hsml(i))
                  end if
               end if
            end if
         end if

      end do
   end if

!---------- sph temperature laplacian calculus(in the polar coordinates system)
   do i=1,ntotal
      n_neigh = 2
      temp_laplacian(i) = 0.

      do while (neighbour(i,n_neigh).ne.0)
         j = neighbour(i,n_neigh)
         r = dist(i,n_neigh)
         q = r/hsml(i)

!-----------------------------------------------------------------
!--------------------lucy's quartic kernel------------------------
!-----------------------------------------------------------------

         if (int_function.eq.1) then
            aux = 2*((t_0(i)-t_0(j))) * factor * ( (-12./hsml(i)**2.) + &
               (24.*r)/hsml(i)**3. - (12*r**2.)/hsml(i)**4. ) * (mass(j)/rho(j))

            temp_laplacian(i) = temp_laplacian(i) + aux
         end if

!-----------------------------------------------------------------
!---------------------cubic spline kernel-------------------------
!-----------------------------------------------------------------

         if (int_function.eq.2) then
            if ((q.ge.0.).and.(q.le.1.)) then
               aux = 2.*(( t_0(i)- t_0(j))) * factor * ( (-2.)/(hsml(i)**2.) + &
                  (3.*r)/(2.*hsml(i)**3.) ) * (mass(j)/rho(j))

               temp_laplacian(i) = temp_laplacian(i) + aux
            else
               if ((q.gt.1.).and.(q.le.2.)) then
                  aux = 2.*( (t_0(i)- t_0(j)) ) * factor * (-1./(2.*hsml(i))) * &
                     ((2. - (r/hsml(i) ) )**2 ) * (1./r) * (mass(j)/rho(j))

                  temp_laplacian(i) = temp_laplacian(i) + aux
               end if
            end if
         end if

!-----------------------------------------------------------------
!---------------------new quartic kernel--------------------------
!-----------------------------------------------------------------

         if (int_function.eq.3) then
            aux = 2*((t_0(i)-t_0(j))) * factor * ( -18./(8.*hsml(i)**2.) + &
               (57.*r)/(24.*hsml(i)**3) - (20.*r*r)/(32.*hsml(i)**4.) ) * &
               (mass(j)/rho(j))

            temp_laplacian(i) = temp_laplacian(i)+ aux
         end if

!-----------------------------------------------------------------
!---------------------quintic spline kernel-----------------------
!-----------------------------------------------------------------

         if (int_function.eq.4) then

            if ((q.ge.0.).and.(q.le.1.)) then
               aux = 2.*( (t_0(i)-t_0(j)) ) * factor * ( -5./hsml(i)) * &
                  ( 10.* (r**3)/hsml(i)**4 - 24.*(r**2)/hsml(i)**3 + 24./hsml(i) ) * &
                  (mass(j)/rho(j))

               temp_laplacian(i) = temp_laplacian(i) + aux
            else
               if ((q.gt.1.).and.(q.le.2.)) then
                  aux = 2.*( (t_0 (i)- t_0 (j)) ) * factor * ( -5./hsml(i) ) * &
                     ( - 5.* (r**3/hsml(i)**4) + 36.*(r**2/hsml(i)**3) - &
                     90.*(r/hsml(i)**2) + 84./hsml(i) - 15./r )* (mass(j)/rho(j))

                  temp_laplacian(i) = temp_laplacian(i) + aux
               else
                  if ((q.gt.2.).and.(q.le.3.)) then
                     aux = 2.*( (t_0 (i)- t_0 (j)) ) * factor * ( -5./hsml(i) ) * &
                        ((r**3)/hsml(i)**4 -(12.*r**2)/hsml(i)**3 + &
                        (54.*r)/hsml(i)**2 - 108./hsml(i) + 81./r ) * (mass(j)/rho(j))

                     temp_laplacian(i) = temp_laplacian(i) + aux
                  end if
               end if
            end if
         end if
         n_neigh = n_neigh +1
      end do
   end do

end subroutine sph_temperature_laplacian


!******************* temporal integration ************************

subroutine temporal_integration( &
   aux_p, itimestep, x, ntotal, dt, temp_laplacian, t_0, t, step_out)

!-----------------------------------------------------------------
! this subroutine performs the integration of the temperature of
! each particle at the domain
!
! ntotal - total number of particles at domain
! dt - timestep
! t_0 - temperature at current iteration
! t - temperature in the next iteration
!-----------------------------------------------------------------

   implicit none
   integer aux_p, ntotal, itimestep, d, step_out, i, j, aux, funit
   double precision temp_laplacian(aux_p), t_0(aux_p), t(aux_p), x(aux_p,2), dt
   character(len=20) aux2
   character(len=35) nome
   open (unit=66, file='output/verifying_temp_sph.dat')
   rewind (66)

!--------calculation of temperature of each particle at domain---

   do i=1,ntotal
      t(i) = t_0(i) + temp_laplacian(i)*dt
   end do

!-------------------------- output files -------------------------

   aux = itimestep

   if ((aux.eq.1).or.(mod(aux, step_out).eq.0)) then
      write (aux2,*) aux
      funit = aux+150

      if (aux.lt.10) then
         nome='output/temperature/000000'//trim(adjustl(aux2))//'.dat'
      end if
      if ((aux.ge.10).and.(aux.lt.100)) then
         nome='output/temperature/00000'//trim(adjustl(aux2))//'.dat'
      end if
      if ((aux.ge.100).and.(aux.lt.1000)) then
         nome='output/temperature/0000'//trim(adjustl(aux2))//'.dat'
      end if
      if ((aux.ge.1000).and.(aux.lt.10000)) then
         nome='output/temperature/000'//trim(adjustl(aux2))//'.dat'
      end if
      if ((aux.ge.10000).and.(aux.le.100000)) then
         nome='output/temperature/00'//trim(adjustl(aux2))//'.dat'
      end if
      if ((aux.ge.100000).and.(aux.lt.1000000)) then
         nome='output/temperature/0'//trim(adjustl(aux2))//'.dat'
      end if
      if (aux.eq.1000000) then
         nome='output/temperature'//trim(adjustl(aux2))//'.dat'
      end if

      open (funit, file=nome,status='unknown')

      do i = 1, ntotal
         write (funit,89) i, x(i,1), x(i,2), t(i)
      end do

      close (funit)

   end if

   do i=1,ntotal
      write (66,89) i, (x(i,d), d = 1, 2), t(i)
   end do

89 format (1x,i10, 1x,d21.14, 1x,d21.14, 1x,d21.14)
   close  (66)

   open  (unit=26, file='step_out.dat')
   write (26, *) step_out
   close (26)

end subroutine temporal_integration


!************************convergence *****************************

subroutine convergence(aux_p, itimestep,t_0,t,eps,response, ntotal)

!-----------------------------------------------------------------
! this subroutine verifies the convergence of the sph solution
! convergence criteria: abs(t_0(i) - t(i)) $<$= accuracy
! where i each particle at domain
!
! t_0 - temperature at current iteration
! t - temperature in the next iteration
! ntotal - total number of particles at domain
! dt - timestep
!-----------------------------------------------------------------

   implicit none
   integer aux_p, itimestep, i, ntotal
   double precision t_0(aux_p), t(aux_p), error(ntotal), eps, max_error
   character (len=15) response

   do i=1,ntotal
      error(i) = abs(t(i) - t_0(i))
   end do

   max_error = maxval(error)

   print *, 'maximum error = ', max_error, ' at ', itimestep, 'a. iteration'
   if (max_error.gt.eps) then
      response = 'non_convergence'
      itimestep = itimestep+1
   else
      response = 'convergence'
   end if

!------------- updating temperatures to the next iteration -------
   do i=1,ntotal
      t_0(i) = t(i)
   end do

end subroutine convergence


!************ point-to-point temperature difference **************

subroutine temperature_differences(itimestep, dt, ntotal, step_out)

!-----------------------------------------------------------------
! this subroutine calculates the differences between the
! temperatures obtained by the sph method and provided by series
!
!-----------------------------------------------------------------

   integer i, d, ntotal, itimestep, part_higher, part_smaller, step_out
   double precision series_solution(ntotal,4), sph_solution(ntotal,4), &
      temp_differences(ntotal), smaller_difference, higher_difference, dt

   open (unit=11, file='output/verifying_temp_series.dat')
   open (unit=66, file='output/verifying_temp_sph.dat')
   open (unit=73, file='output/temperature_differences.dat')
   open (unit=74, file='output/final_diferences.dat')
   open (unit=75, file='output/step_differences.dat')

   rewind (73)
   rewind (74)

   do i=1,ntotal
      read (11,*) (series_solution(i,d), d=1,4)
      read (66,*) (sph_solution(i,d), d=1,4)
   end do

   do i=1,ntotal
      temp_differences(i) = abs(series_solution(i,4) - sph_solution(i,4))
      write (73,55) (series_solution(i,d), d=2,3), temp_differences(i)
      write (74,55) (series_solution(i,d), d=2,3), temp_differences(i)

      if (mod(itimestep,step_out).eq.0) then
         write (75,55) (series_solution(i,d), d=2,3), temp_differences(i)
      end if

   end do

55 format (1x,d21.14, 1x,d21.14, 1x, d21.14)

   higher_difference = 0.
   smaller_difference = 1500.

   do i=1,ntotal
      if (temp_differences(i).gt.higher_difference) then
         part_higher = i
         higher_difference = temp_differences(i)
      end if

      if (temp_differences(i).lt.smaller_difference) then
         part_smaller = i
         smaller_difference = temp_differences(i)
      end if
   end do

   write (73,*) ''
   write (73,*) 'results at the ', itimestep , 'a. iteration'
   write (73,*) 'time = ' , itimestep*dt, ' seconds'
   write (73,*) ''
   write (73,*) '----------------------------------------'
   write (73,*) 'higher temperature difference'
   write (73,*) '----------------------------------------'
   write (73,*) ' position = ', (series_solution(part_higher,d),d=2,3)
   write (73,*) ' series temperature --- sph temperature --- difference'
   write (73,*) series_solution(part_higher,4), sph_solution(part_higher,4), &
      higher_difference
   write (73,*) '-----------------------------------'
   write (73,*) '----------------------------------------'
   write (73,*) 'smaller temperature difference'
   write (73,*) '----------------------------------------'
   write (73,*) ' position = ', (series_solution(part_smaller,d),d=2,3)
   write (73,*) ' series temperature --- sph temperature --- difference'
   write (73,*) series_solution(part_smaller,4), sph_solution(part_smaller,4), &
      smaller_difference
   write (73,*) '---------------------------------------'

   close (11)
   close (66)
   close (73)
   close (74)
   close (75)

end subroutine temperature_differences
