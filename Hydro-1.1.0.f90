PROGRAM HydroSim
  IMPLICIT NONE
  INTEGER :: n, i, j, ndims, nneighs, kernTot, qint
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IDs
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nears
  DOUBLE PRECISION :: t, G, pi, dist, temperature, dt, kb, dx, dy, dz, twrite, tcount, tlim, q, vdotr, viscalpha, viscbeta, muab, etasq, adiabrat, avsndspd
  DOUBLE PRECISION :: avdens
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ranPos, ranVel, m, kern0, kern1, dens, press, h, sndspd
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r, v, a, kern, kerndiv, dists, visc

  ! kern0 and kern1 contain values for the kernel and its first derivative for a given h at evenly spaced q-values
  ! kernTot dictates the spacing between the q-values used for kern0 and kern1. The lengths of kern0 and kern1 should be equal to kernTot
  ! dist is used as a temporary distance value in various places throughout
  ! dists contains distances between every pair of particles. I am keeping this as a 2D array with shape (n,n) just in case I end up needing to keep the distances stored somewhere so that they can be accessed again later on in the same timestep. If I find that this is unneeded, then I will change it to a 1D array with length n, which is all that is needed when initially calculating distances, nearest neighbours, kernels, etc
  ! qint is used when calculating approximate values for the kernel and its first derivative for a given particle. It is found by multiplying the calculated q value by kernTot/2, then adding 0.5 and casting q as an integer. Casting a value to an integer will truncate to the decimal point, the +0.5 term causes any values between x.5 and (x+1).0 to end up falling between (x+1).0 and (x+1).5, which truncates to x+1. This essentially rounds these values up to the nearest integer
  ! dx, dy and dz store the separations along the x, y and z axes respectively when calculating components of acceleration
  ! twrite is the ideal simulated time between successive file writes
  ! tcount keeps track of the progress towards the next file write
  ! tlim is the ideal total simulated time
  ! vdotr is the dot product of (v_i - v_j) and (r_i - r_j), used in viscosity calculations
  ! adiabrat is the adiabatic ratio

  OPEN(7, file='positions.txt', status='replace')
  OPEN(8, file='densities.txt', status='replace')
  OPEN(9, file='pressures.txt', status='replace')
  OPEN(10, file='time.txt', status='replace')
  
  n = 1000
  ndims = 3
  nneighs = 32
  kernTot = 1000
  
  ALLOCATE(IDs(1:n))
  ALLOCATE(nears(1:nneighs, 1:n))
  ALLOCATE(ranPos(1:ndims))
  ALLOCATE(ranVel(1:ndims))
  ALLOCATE(m(1:n))
  ALLOCATE(r(1:ndims, 1:n))
  ALLOCATE(v(1:ndims, 1:n))
  ALLOCATE(a(1:ndims, 1:n))
  ALLOCATE(kern(1:nneighs, 1:n))
  ALLOCATE(kerndiv(1:nneighs, 1:n))
  ALLOCATE(kern0(1:kernTot))
  ALLOCATE(kern1(1:kernTot))
  ALLOCATE(dists(1:n, 1:n))
  ALLOCATE(dens(1:n))
  ALLOCATE(press(1:n))
  ALLOCATE(visc(1:nneighs, 1:n))
  ALLOCATE(h(1:n)) ! Values of h are needed when calculating viscosity, which is only done after all densities are calculated as the average density of the particle pair is also needed
  ALLOCATE(sndspd(1:n)) ! Sound speed at both particles i and j is needed, so will store all sound speeds at any given time

  t = 0.
  tcount = 0.
  twrite = 100.
  tlim = 10000.
  dt = 0.001
  dx = 0.
  dy = 0.
  dz = 0.
  G = 6.67430e-11
  kb = 1.38065e-23
  pi = 4. * ATAN(1.)
  r = 0.
  v = 0.
  a = 0.
  kern = 0.
  kerndiv = 0.
  IDs = 0.
  nears = 0.
  h = 0.
  ranPos = 0.
  ranVel = 0.
  m = 1.
  dist = 0.
  dists = 0.
  kern0 = 0.
  kern1 = 0.
  dens = 0.
  press = 0.
  temperature = 300. ! Setting temperature to 300K for now because it is a convenient value that is close to room temperature
  q = 0.
  visc = 0.
  vdotr = 0.
  viscalpha = 1.
  viscbeta = 2. ! viscalpha and viscbeta should be near 1 and 2 respectively
  muab = 0.
  etasq = 0.
  adiabrat = 5./3. ! Value for a monatomic ideal gas with 3 translational degrees of freedom
  avsndspd = 0.
  avdens = 0.
  


  DO i=1,n ! Randomising initial conditions and setting up IDs array
     IDs(i) = i
10   CALL random_NUMBER(ranPos)
     r(:,i) = 2. * (ranPos(:) - 0.5) ! Offsets all coordinates so that some particles will have negative coordinates
     dist = SQRT((r(1,i)**2) + (r(2,i)**2) + (r(3,i)**2)) ! Finding distance from origin
     IF (dist > 1.) THEN
        GOTO 10 ! Re-randomises the initial position if the particle is further than 1 unit from the origin. This forces the particles into a sphere rather than a cube
     END IF
     CALL random_NUMBER(ranVel)
     v(:,i) = 2. * (ranVel(:) - 0.5) ! Offsets all velocities so that some particles will have negative components
  END DO

  DO
     ! Resetting some values that are calculated by summing over contributions from neighbouring particles
     dens = 0.
     a = 0.
     
     DO i=1,n ! Finding nearest neighbours, densities, pressures, sound speeds and kernels for each particle
        DO j=1,n ! Finding distances from particle i to each particle
           IF (i==j) CYCLE ! Finding the distance between a given particle and itself is completely redundant
           dists(j,i) = SQRT(((r(1,j) - r(1,i))**2) + ((r(2,j) - r(2,i))**2) + ((r(3,j) - r(3,i))**2))
        END DO ! At this point, the distances between particle i and all other particles has been calculated

        CALL heapsort(n, dists(:,i), IDs(:)) ! Sorts the IDs array in ascending order according to the distance of each particle from particle i. This means that position 1 in the IDs array should contain particle i's ID, and the nth nearest neighbour is in position (n+1)
        DO j=1,nneighs
           nears(j,i) = IDs(j+1) ! Note the use of j+1 as the relevant index for the IDs array due to IDs(1) being the ID of the current particle
        END DO ! At this point, the IDs of the nearest neighbours for particle i have been collected into the nears array

        h(i) = nears(32,i)/2. ! 2h is the distance to the 32nd nearest neighbour
        CALL kernel(kernTot, ndims, h(i), kern0(:), kern1(:)) ! Fills the kern0 and kern1 arrays with appropriate values

        DO j=1,nneighs ! This fills the kern and kerndiv arrays with approximate values, then works out the density at particle i by summing over contributions from the nearest neighbours
           dist = dists(nears(j,i),i) ! nears(j,i) should contain the ID for the jth nearest particle to particle i
           q = dist/h(i) ! Calculates q-value
           qint = INT((q*(kernTot/2.)) + 0.5) ! Adding 0.5 produces the intended result of rounding (q*(kernTot/2.)) to the nearest integer instead of truncating it
           IF (qint==0) THEN ! This is an edge-case where the approximate q-value corresponds to an index of zero, which is out of range for the kern0 and kern1 arrays. In this case, the kernel and its first derivative will be calculated explicitly with q=0, removing all of the q-terms from the equations
              IF (ndims==1) THEN
                 kern(j,i) = 2./(3. * h(i))
              ELSE IF (ndims==2) THEN
                 kern(j,i) = 10./(7. * pi * (h(i)**2))
              ELSE ! Assumes that the simulation will only ever be 1-, 2- or 3-dimensional
                 kern(j,i) = 1./(pi * (h(i)**3))
              END IF
              kerndiv(j,i) = 0. ! First derivative of the kernel is trivially zero for q=0     
           ELSE IF (qint >= 1 .AND. qint <= 1000) THEN
              kern(j,i) = kern0(qint)
              kerndiv(j,i) = kern1(qint)
           ELSE ! This bypasses a segmentation fault by telling the program to ignore that specific neighbour if the index will be out of range. This doesn't actually fix the segmentation fault (which likely arises in how qint is initially calculated), but it does prevent it from occuring during runtime
              kern(j,i) = 0.
              kerndiv(j,i) = 0.
           END IF

           dens(i) = dens(i) + (m(nears(j,i)) * kern(j,i)) ! This line is why density is universally set to zero at the start of each iteration of the main DO loop          
        END DO ! At this point, kern(:,i) and kerndiv(:,i) will both be filled with the appropriate values for particle i. The density at particle i has also been calculated

        press(i) = kb * dens(i) * temperature / m(i) ! From an ideal gas law, pressure (per unit volume) = number density x kb x temperature. Number density is the mass density divided by the mean particle mass, so P = rho*kb*T/m. Currently all particles have the same mass, so the mean particle mass is equivalent to the mass of a given particle i
        sndspd(i) = SQRT(adiabrat * press(i) / dens(i))
     END DO

     DO i=1,n ! Artificial Viscosity
        DO j=1,nneighs
           vdotr = dot_PRODUCT((v(:,i)-v(:,nears(j,i))),(r(:,i)-r(:,nears(j,i))))
           IF (vdotr > 0) THEN
              etasq = 0.01 * (h(i)**2)
              muab = (h(i) * vdotr) / ((dists(nears(j,i),i)**2) * etasq) ! dists(nears(j,i),i) is defined as (r_j - r_i). The viscosity formula uses (r_i - r_j)^2, which is defined the other way around, but this doesn't actually matter due to it being squared
              avsndspd = 0.5 * (sndspd(i) + sndspd(nears(j,i)))
              avdens = 0.5 * (dens(i) + dens(nears(j,i)))
              visc(j,i) = (((-1. * viscalpha * avsndspd * muab) + (viscbeta * (muab**2))) / avdens)
           ELSE
              visc(j,i) = 0.
           END IF
        END DO
     END DO
     
     DO i=1,n ! Acceleration needs to be calculated after everything else as the contribution to the acceleration of particle i from particle j depends on the pressures and densities of both i and j. Acceleration components can be found by multiplying the acceleration by dx/a, where a is the total separation and dx is the separation along a given axis
        DO j=1, nneighs
           dx = r(1,nears(j,i)) - r(1,i)
           dy = r(2,nears(j,i)) - r(2,i)
           dz = r(3,nears(j,i)) - r(3,i)
           dist = SQRT((dx**2)+(dy**2)+(dz**2))

           ! The following code assumes that the acceleration array has been set to zero at the start of the main DO loop
           ! This is probably quite slow. Could speed it up by introducing a variable that stores everything apart from the d/dist term so that this only needs to be calculated once per loop
           a(1,i) = a(1,i) - (m(nears(j,i)) * ((press(nears(j,i))/(dens(nears(j,i))**2)) + (press(i)/(dens(i)**2)) + visc(j,i)) * kerndiv(j,i) * (dx/dist))
           a(2,i) = a(2,i) - (m(nears(j,i)) * ((press(nears(j,i))/(dens(nears(j,i))**2)) + (press(i)/(dens(i)**2)) + visc(j,i)) * kerndiv(j,i) * (dy/dist))
           a(3,i) = a(3,i) - (m(nears(j,i)) * ((press(nears(j,i))/(dens(nears(j,i))**2)) + (press(i)/(dens(i)**2)) + visc(j,i)) * kerndiv(j,i) * (dz/dist))
        END DO
     END DO



     DO i=1,n ! Advancing the simulation
        r(:,i) = r(:,i) + (v(:,i) * dt)
        v(:,i) = v(:,i) + (a(:,i) * dt)
        t = t + dt
        tcount = tcount + dt
     END DO

     IF (t>tlim) THEN ! Writing to files
        WRITE(6,*) 'Writing final values to files'
        WRITE(7,*) r(1,:), r(2,:), r(3,:)
        WRITE(8,*) dens(:)
        WRITE(9,*) press(:)
        WRITE(10,*) t
        CLOSE(7)
        CLOSE(8)
        CLOSE(9)
        CLOSE(10)
        EXIT
     ELSE IF (tcount>=twrite) THEN
        WRITE(6,*) 'Writing to files', (t*100./tlim), '% complete'
        WRITE(7,*) r(1,:), r(2,:), r(3,:)
        WRITE(8,*) dens(:)
        WRITE(9,*) press(:)
        WRITE(10,*) t
        tcount = tcount - twrite
     END IF
     
  END DO

END PROGRAM HydroSim

! -----Subroutines-----
SUBROUTINE heapsort(psort,measureof,pwhichhas)
  ! does a heapsort (by APW)
  IMPLICIT NONE
  INTEGER, INTENT(IN)              :: psort              ! number of values to be sorted.
  DOUBLE PRECISION,    INTENT(IN)  :: measureof(1:psort) ! values to be sorted.
  INTEGER, INTENT(OUT)             :: pwhichhas(1:psort) ! identifier of value.
  INTEGER                          :: rank               ! rank of value.
  INTEGER                          :: ranknow            ! dummy rank.
  INTEGER                          :: ranktest           ! dummy rank.
  !
  DO rank=2,psort                ! THIS DO-LOOP BUILDS THE BINARY HEAP
     ranknow=rank
1    IF (ranknow==1) CYCLE
     ranktest=ranknow/2
     IF (measureof(pwhichhas(ranktest))>=measureof(pwhichhas(ranknow))) CYCLE
     CALL swapi(pwhichhas(ranknow),pwhichhas(ranktest))
     ranknow=ranktest
     GOTO 1
  END DO
  !
  DO rank=psort,2,-1             ! AND THIS DO-LOOP INVERTS THE BINARY HEAP
     CALL swapi(pwhichhas(rank),pwhichhas(1))
     ranknow=1
2    ranktest=2*ranknow
     IF (ranktest>=rank) CYCLE
     IF ((measureof(pwhichhas(ranktest+1))>measureof(pwhichhas(ranktest))) &
          & .AND.(ranktest+1<rank)) ranktest=ranktest+1
     IF (measureof(pwhichhas(ranktest))<=measureof(pwhichhas(ranknow))) CYCLE
     CALL swapi(pwhichhas(ranknow),pwhichhas(ranktest))
     ranknow=ranktest
     GOTO 2
  END DO
  !
  RETURN
END SUBROUTINE heapsort

!
! ===========================================================================
! ===========================================================================
!
SUBROUTINE swapi(item1,item2)
  !
  IMPLICIT NONE
  INTEGER :: item0,item1,item2
  !
  item0=item1; item1=item2; item2=item0
  !
  RETURN
END SUBROUTINE swapi

  
SUBROUTINE kernel(kernTot, dims, hval, kern0, kern1)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: kernTot, dims
  INTEGER :: i
  DOUBLE PRECISION, INTENT(IN) :: hval
  DOUBLE PRECISION, INTENT(OUT) :: kern0(1:kernTot), kern1(1:kernTot)
  DOUBLE PRECISION :: h, qval, pi, prefac1, prefac2
  
  qval = 1./500. ! Starts at qval = 1/500 due to how "qint" is calculated in the main program. qint will take on a value of 1 if 0.001<=q<0.003, so qval=0.002 will be the q-value for the first values in kern0 and kern1 to allow qint to correspond to the intended value in these arrays
  pi = 4. * ATAN(1.)
  IF (dims==1) THEN
     prefac1 = 2./(3. * hval)
  ELSE IF (dims==2) THEN
     prefac1 = 10./(7. * pi * (hval**2))
  ELSE ! Assumes that the simulation will only ever be 1-, 2- or 3-dimensional
     prefac1 = 1./(pi * (hval**3))
  END IF
  prefac2 = prefac1/4. ! When 1<qval<=2, there is an additional prefactor of 1/4
  
  DO i=1,kernTot ! Note that the following assumes that qval is between 0 and 2 (inclusive), which it should be within this subroutine
     IF (qval<=1.) THEN
        kern0(i) = prefac1 * (1 - (1.5 * (qval**2)) + (0.75 * (qval**3)))
        kern1(i) = prefac1 * ((-3. * qval) + (2.25 * (qval**2)))
     ELSE
        kern0(i) = prefac2 * ((2. - qval)**3)
        kern1(i) = -3. * prefac2 * ((2. - qval)**2)
     END IF
     qval = qval + (2./kernTot) ! Increments the value of qval for the next loop
  END DO

  RETURN
END SUBROUTINE kernel

        
     
  
     
     
  
  
