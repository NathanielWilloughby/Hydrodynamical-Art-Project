PROGRAM Hydro
  IMPLICIT NONE
  INTEGER :: n, i, j, k, ndims, nneighs, kernTot, qint, counter, kdum
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IDs
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nears
  DOUBLE PRECISION :: tempDist, q, pi, kb, adiabrat, etasq, vdotr, muab, avsndspd, avdens, tempAcc, t, tcount, dt, tlim, twrite, viscalpha, viscbeta, ran3
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m, ranPos, dens, h, kern0, kern1, press, temps, dr, sndspd
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r, v, a, dists, kern, kerndiv, visc

  OPEN(7, file='positions.txt', status='replace')
  OPEN(8, file='velocities.txt', status='replace')
  OPEN(9, file='densities.txt', status='replace')
  OPEN(10, file='pressures.txt', status='replace')
  OPEN(11, file='temperatures.txt', status='replace')
  OPEN(12, file='time.txt', status='replace')

  n = 1000
  ndims = 3
  nneighs = 32 ! Number of nearest neighbours to consider
  kernTot = 1000 ! Number of values in the smoothed kernel arrays

  ALLOCATE(IDs(1:n))
  ALLOCATE(m(1:n))
  ALLOCATE(r(1:ndims, 1:n))
  ALLOCATE(v(1:ndims, 1:n))
  ALLOCATE(a(1:ndims, 1:n))
  ALLOCATE(ranPos(1:ndims))
  ALLOCATE(dens(1:n))
  ALLOCATE(dists(1:n, 1:n))
  ALLOCATE(nears(1:nneighs, 1:n))
  ALLOCATE(h(1:n))
  ALLOCATE(kern0(1:kernTot))
  ALLOCATE(kern1(1:kernTot))
  ALLOCATE(kern(1:nneighs, 1:n))
  ALLOCATE(kerndiv(1:nneighs, 1:n))
  ALLOCATE(press(1:n))
  ALLOCATE(temps(1:n)) ! Temperatures stored in an array just in case the system is set up to be non-isothermal, or if particle temperatures are allowed to differ from particle to particle
  ALLOCATE(visc(1:nneighs, 1:n))
  ALLOCATE(dr(1:ndims))
  ALLOCATE(sndspd(1:n))

  t = 0.
  tcount = 0.
  dt = 0.001
  tlim = 2000.
  twrite = 1.
  counter = 0 ! Stores the current timestep number, which is useful for debugging

  m = 1. ! Universally setting masses to 1kg for now
  r = 0.
  v = 0.
  a = 0.
  ranPos = 0.
  tempDist = 0. ! Used throughout to temporarily store distances
  dens = 0.
  dists = 0.
  nears = 0 ! IDs of the nearest neighbours to each particle
  h = 0. ! 2h is the distance of the 32nd nearest neighbour from a given particle
  kern0 = 0. ! Smoothed kernel array
  kern1 = 0. ! Smoothed kernel derivative array
  q = 0.
  pi = 4. * ATAN(1.)
  kern = 0.
  kerndiv = 0.
  press = 0.
  temps = 300. ! Currently setting every particle to be just above room temperature
  kb = 1.38065e-23
  adiabrat = 5./3. ! This can be changed depending on the system
  visc = 0.
  etasq = 0.
  muab = 0.
  vdotr = 0.
  avsndspd = 0.
  avdens = 0.
  tempAcc = 0.
  dr = 0.
  viscalpha = 1. ! The exact values of these parameters are not too important as long as they are close to 1 and 2 respectively
  viscbeta = 2.
  sndspd = 0.
  qint = 0
  kdum = 762 ! Arbitrary random seed. This is used for testing when the same initial conditions are desired
  

  DO i=1,n ! Placing particles and setting up IDs array
     IDs(i) = i ! Fills up IDs array. This array will be constantly re-ordered throughout, but the array needs to be set up such that the IDs initially correspond to the particle at the same index in the position array. After IDs is re-ordered, the IDs will be sorted in ascending distance from a given particle. This is then used to fill the nears array with the nneighs nearest IDs to particle i

10   CALL random_NUMBER(ranPos)
     r(:,i) = 2. * (ranPos(:) - 0.5) ! Offsets all coordinates so that some particles will be at negative coordinates
     tempDist = 0. ! Setting to zero because I am calculating it by summing over contributions in each dimension then taking a square root, which keeps it valid for any integer number of dimensions
     DO j=1,ndims
        tempDist = tempDist + (r(j,i)**2)
     END DO
     tempDist = SQRT(tempDist)
     IF (tempDist > 1.) THEN
        GOTO 10 ! Re-randomises initial position if particle is outside of a sphere of radius 1 unit centred on the origin
     END IF
  END DO

  DO
     ! Resetting some values that are calculated by summing over contributions
     dens = 0.
     a = 0.
     
     counter = counter + 1 ! This goes at the start of the loop because I want counter to equal 1 during the first timestep

     DO i=1,n ! Nearest neighbours, densities, pressures, sound speeds and kernels for each particle
        DO j=1,n ! Distances between each particle
           tempDist = 0.
           IF (i==j) CYCLE ! Distance between a particle and itself is redundant
           DO k=1,ndims ! Adding the square of each component together
              tempDist = tempDist + ((r(k,j) - r(k,i))**2)
           END DO
           dists(j,i) = SQRT(tempDist)
           
        END DO ! Distances from particle i should now have been calculated

        CALL heapsort(n, dists(:,i), IDs(:)) ! Sorts the IDs array in ascending order by distance from particle i

        DO j=1,nneighs ! Fills nears(:,i) with the nearest neighbours to particle i
           nears(j,i) = IDs(j+1) ! IDs(1) is the ID of particle i, so need to start from IDs(2)
        END DO ! Nearest neighbours should now have been calculated for particle i

        h(i) = dists(nears(32,i),i)/2.
        

        CALL kernel(kernTot, ndims, h(i), kern0(:), kern1(:)) ! Fills up the smoothed kernel arrays for the current particle

        DO j=1,nneighs ! Filling kern and kerndiv arrays with appropriate values and working out the contribution to the density of particle i by particle j
           tempDist = dists(nears(j,i),i)
           q = tempDist/h(i)
           
           qint = INT((q*(kernTot/2.)) + 0.5) ! Casting as an integer truncates the value, so adding 0.5 before truncating essentially simulates rounding to the nearest integer
           

           IF (qint==0) THEN ! Edge-case as the first value in the smoothed arrays has a non-zero q

              IF (ndims==1) THEN
                 kern(j,i) = 2./(3. * h(i))
              ELSE IF (ndims==2) THEN
                 kern(j,i) = 10./(7. * pi * (h(i)**2))
              ELSE ! Assumes simulation will only ever be 1-, 2- or 3-dimensional
                 kern(j,i) = 1./(pi * (h(i)**3))
              END IF

              kerndiv(j,i) = 0. ! First derivative of the kernel is trivially zero for q=0
           ELSE IF ((qint<0) .OR. (qint>kernTot)) THEN ! Out of range, will just tell program to ignore this particle's contribution by setting the kernel and its derivative to zero
              kern(j,i) = 0.
              kerndiv(j,i) = 0.
           ELSE
              kern(j,i) = kern0(qint)
              kerndiv(j,i) = kern1(qint)
           END IF

           dens(i) = dens(i) + (m(nears(j,i)) * kern(j,i))
        END DO ! The density, kernel and derivative of the kernel should now have been fully computed for particle i

        press(i) = kb * dens(i) * temps(i) / m(i) ! Ideal gas law. This can be changed if different equations of state or relations between density and pressure are desired

        sndspd(i) = SQRT(adiabrat * press(i) / dens(i))
     END DO ! Nearest neighbours, kernels, first kernel derivatives, densities, pressures and sound speeds for all particles should now have been found

     DO i=1,n ! Artificial Viscosity

        etasq = 0.01 * (h(i)**2)
        DO j=1,nneighs
           vdotr = dot_PRODUCT((v(:,i)-v(:,nears(j,i))),(r(:,i)-r(:,nears(j,i))))
           IF (vdotr > 0.) THEN
              muab = (h(i) * vdotr) / ((dists(nears(j,i),i)**2) * etasq) ! dists(nears(j,i),i) is defined as (r_j - r_i) when using hard boundaries. The viscosity formula uses (r_i - r_j), which is the other way around. However, this doesn't matter as the value is squared anyway
              avsndspd = 0.5 * (sndspd(i) + sndspd(nears(j,i)))
              avdens = 0.5 * (dens(i) + dens(nears(j,i)))
              visc(j,i) = (((-1. * viscalpha * avsndspd * muab) + (viscbeta * (muab**2))) / avdens)
           ELSE
              visc(j,i) = 0.
           END IF
        END DO
     END DO ! Artificial viscosity from the nearest neighbours for each particle should now have been computed

     DO i=1,n ! Accceleration needs to be calculated after everything else as the contribution to the acceleration of particle i from neighbour j depends on the pressures and densities at both i and j. Components of acceleration can be found by multiplying the acceleration by dr(k)/a, where a is the total separation and dr(k) is the separation along a given axis k

        DO j=1, nneighs
           tempDist = 0.
           tempAcc = m(nears(j,i)) * ((press(nears(j,i))/(dens(nears(j,i))**2)) + (press(i)/(dens(i)**2)) + visc(j,i)) * kerndiv(j,i)
           DO k=1, ndims ! Technically can just do this entire DO loop in one line when using hard boundaries, but I am keeping the components separate to ease the conversion to looping boundaries when desired
              dr(k) = r(k,nears(j,i)) - r(k,i)
              tempDist = tempDist + (dr(k)**2)
           END DO
           tempDist = SQRT(tempDist)

           ! The following assumes that the acceleration array has been set to zero at the start of the main DO loop
           a(:,i) = a(:,i) - (tempAcc * (dr(:)/tempDist))
        END DO
     END DO ! Accelerations should now have been computed for each particle

     DO i=1,n ! Advancing the simulation
        DO k=1,ndims ! Again, all components of position and velocity could be calculated in one line each for hard boundaries, but I am separating the components so that it is easier to convert to periodic boundaries if desired
           r(k,i) = r(k,i) + (v(k,i) * dt)          
           v(k,i) = v(k,i) + (a(k,i) * dt)           
        END DO
     END DO

     t = t + dt
     tcount = tcount + dt

     
     IF (t>=tlim) THEN ! Writing to files
        WRITE(6,*) 'Writing final values to files'
        WRITE(7,*) r(:,:)
        WRITE(8,*) v(:,:)
        WRITE(9,*) dens(:)
        WRITE(10,*) press(:)
        WRITE(11,*) temps(:)
        WRITE(12,*) t
99      CLOSE(7)
        CLOSE(8)
        CLOSE(9)
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        EXIT
     END IF
     
     IF (tcount>=twrite) THEN
        WRITE(6,*) 'Writing to files', (t*100./tlim), '% complete'
        WRITE(7,*) r(:,:)
        WRITE(8,*) v(:,:)
        WRITE(9,*) dens(:)
        WRITE(10,*) press(:)
        WRITE(11,*) temps(:)
        WRITE(12,*) t
        tcount = tcount - twrite
     END IF

  END DO

END PROGRAM Hydro

        
        
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

     
SUBROUTINE heapsort(psort,measureof,pwhichhas)
  ! does a heapsort (by APW)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: psort              ! number of values to be sorted.
  DOUBLE PRECISION,    INTENT(IN)  :: measureof(1:psort) ! values to be sorted.
  INTEGER, INTENT(OUT) :: pwhichhas(1:psort) ! identifier of value.
  INTEGER              :: rank               ! rank of value.
  INTEGER              :: ranknow            ! dummy rank.
  INTEGER              :: ranktest           ! dummy rank.
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

!
! ===========================================================================
      DOUBLE PRECISION FUNCTION ran3(idum)
      INTEGER :: idum
      INTEGER :: MBIG,MSEED,MZ
      DOUBLE PRECISION :: FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER :: i,iff,ii,inext,inextp,k
      INTEGER :: mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      IF(idum.LT.0.OR.iff.EQ.0)THEN
        iff=1
        mj=MSEED-iabs(idum)
        mj=MOD(mj,MBIG)
        ma(55)=mj
        mk=1
        DO 11 i=1,54
          ii=MOD(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          IF(mk.LT.MZ)mk=mk+MBIG
          mj=ma(ii)
11      CONTINUE
        DO 13 k=1,4
          DO 12 i=1,55
            ma(i)=ma(i)-ma(1+MOD(i+30,55))
            IF(ma(i).LT.MZ)ma(i)=ma(i)+MBIG
12        CONTINUE
13      CONTINUE
        inext=0
        inextp=31
        idum=1
      ENDIF
      inext=inext+1
      IF(inext.EQ.56)inext=1
      inextp=inextp+1
      IF(inextp.EQ.56)inextp=1
      mj=ma(inext)-ma(inextp)
      IF(mj.LT.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      RETURN
      END FUNCTION ran3
!


















