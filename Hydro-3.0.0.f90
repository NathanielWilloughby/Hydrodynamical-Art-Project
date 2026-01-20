SUBROUTINE HydroSim(n, ndims, rin, vin, tempin, tlim, dt, cin, press, dens, temps, r, v, counter)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, ndims, cin
  INTEGER, INTENT(OUT) :: counter
  INTEGER :: i, j, k, nneighs, kernTot, qint, recomp, recomplim, nbneighs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IDs, bneighsIDs, tempIDs
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nearIDs
  DOUBLE PRECISION, INTENT(IN) :: rin(1:ndims,1:n), vin(1:ndims,1:n), tempin(1:n)
  DOUBLE PRECISION, INTENT(OUT) :: press(1:n), dens(1:n), temps(1:n), r(1:ndims,1:n), v(1:ndims,1:n)
  DOUBLE PRECISION :: tempDist, q, pi, kb, adiabrat, etasq, vdotr, muab, avsndspd, avdens, tempAcc, t, tcount, dt, tlim, twrite, viscalpha, viscbeta
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m, ranPos, h, kern0, kern1, dr, sndspd, allDists, nearDists
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: a, kern, kerndiv, visc, dists

  nneighs = 32 ! Number of nearest neighbours to consider
  nbneighs = 100 ! Number of nearest neighbours to identify every few timesteps
  kernTot = 1000 ! Number of values in the smoothed kernel arrays

  ALLOCATE(IDs(1:n))
  ALLOCATE(bneighsIDs(1:nbneighs))
  ALLOCATE(m(1:n))
  ALLOCATE(a(1:ndims, 1:n))
  ALLOCATE(ranPos(1:ndims))
  ALLOCATE(allDists(1:n))
  ALLOCATE(h(1:n))
  ALLOCATE(kern0(1:kernTot))
  ALLOCATE(kern1(1:kernTot))
  ALLOCATE(kern(1:nneighs, 1:n))
  ALLOCATE(kerndiv(1:nneighs, 1:n))
  ALLOCATE(visc(1:nneighs, 1:n))
  ALLOCATE(dr(1:ndims))
  ALLOCATE(sndspd(1:n))
  ALLOCATE(nearIDs(1:nbneighs, 1:n))
  ALLOCATE(nearDists(1:nbneighs))
  ALLOCATE(tempIDs(1:nbneighs))
  ALLOCATE(dists(1:nbneighs, 1:n))

  t = 0.
  dt = 0.001
  tlim = 1.
  recomplim = 10 ! Number of timesteps between each time the nbneighs nearest neighbours should be identified
  recomp = recomplim ! The nbneighs nearest neighbours will be identified in the first timestep

  m = 1.
  a = 0.
  ranPos = 0.
  tempDist = 0. ! Used throughout to temporarily store distances
  dists = 0.
  h = 0. ! 2h is the distance of the 32nd nearest neighbour from a given particle
  kern0 = 0. ! Smoothed kernel array
  kern1 = 0. ! Smoothed kernel derivative array
  q = 0.
  pi = 4. * ATAN(1.)
  kern = 0.
  kerndiv = 0.
  kb = 1.38065e-23
  adiabrat = 1.4 ! This can be changed depending on the system
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
  nearIDs = 0 ! This stores the actual IDs of the nbneighs nearest neighbours, and is not touched by any sorting algorithms. Values in this array are assigned explicitly when needed
  bneighsIDs = 0
  nearDists = 0.
  tempIDs = 0
  allDists = 0.
  r(:,:) = rin(:,:)
  v(:,:) = vin(:,:)
  temps(:) = tempin(:)
  counter = cin

  DO i=1,nbneighs ! Setting up an array of the first 100 natural numbers. These correspond to indices for values in the nearIDs array, and will be constantly re-ordered throughout according to distance from particle i
     bneighsIDs(i) = i
  END DO
       

  DO i=1,n ! Setting up IDs array
     IDs(i) = i ! Fills up IDs array. This array will be constantly re-ordered throughout, but the array needs to be set up such that the IDs initially correspond to the particle at the same index in the position array
  END DO

  DO
     ! Resetting some values that are calculated by summing over contributions
     dens = 0.
     a = 0.
     
     counter = counter + 1 ! This goes at the start of the loop because I want counter to equal 1 during the first timestep
    
     IF (recomp==recomplim) THEN ! This IF-ELSE statement finds the IDs of the nearest neighbours to particle i, as well as their distances to particle i
        recomp = 1
        DO i=1,n
           DO j=1,n
              tempDist = 0.
              DO k=1,ndims
                 tempDist = tempDist + ((r(k,j) - r(k,i))**2)
              END DO
              allDists(j) = SQRT(tempDist)
           END DO
           

           CALL heapsort(n, allDists(:), IDs(:))
           DO j=1,nbneighs
              nearIDs(j,i) = IDs(j+1)
              dists(j,i) = allDists(nearIDs(j,i))
           END DO ! nearIDs(:,i) will now contain the IDs of the nbneighs nearest particles to particle i. dists(:,i) will now contain the distances to these particles
        END DO
     ELSE
        recomp = recomp + 1
        DO i=1,n
           DO j=1,nbneighs
              tempDist = 0.
              DO k=1,ndims
                 tempDist = tempDist + ((r(k,nearIDs(j,i)) - r(k,i))**2)
              END DO
              nearDists(j) = SQRT(tempDist)
           END DO ! nearDists(j) now contains the distances between particle i and the nbneighs nearest neighbours that were most recently identified, sorted by the distance from particle i in the previous timestep

           CALL heapsort(nbneighs, nearDists(:), bneighsIDs(:))
           DO j=1,nbneighs
              tempIDs(j) = nearIDs(bneighsIDs(j),i)
              dists(j,i) = nearDists(bneighsIDs(j))
           END DO

           nearIDs(:,i) = tempIDs(:) ! Updates nearIDs, needed so that the rest of the timestep can occur in the same way regardless of which part of the current IF-ELSE statement was executed
        END DO
     END IF ! dists now contains the distances of the nbneighs nearest neighbours to each particle in ascending order per particle. nearIDs contains the IDs of the nbneighs nearest neighbours to each particle, also sorted in ascending order by distance

     

     DO i=1,n ! Densities, pressures, sound speeds and kernels for each particle
        
        h(i) = dists(nneighs,i)/2.
        
        CALL kernel(kernTot, ndims, h(i), kern0(:), kern1(:)) ! Fills up the smoothed kernel arrays for the current particle

        DO j=1,nneighs ! Filling kern and kerndiv arrays with appropriate values and working out the contribution to the density of particle i by particle j
           q = dists(j,i)/h(i)
           
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

           dens(i) = dens(i) + (m(nearIDs(j,i)) * kern(j,i))           
        END DO ! The density, kernel and derivative of the kernel should now have been fully computed for particle i

        press(i) = kb * dens(i) * temps(i) / m(i)

        sndspd(i) = SQRT(adiabrat * press(i) / dens(i))       
     END DO ! Kernels, first kernel derivatives, densities, pressures and sound speeds for all particles should now have been found
     
        
     DO i=1,n ! Artificial Viscosity

        etasq = 0.01 * (h(i)**2)
        DO j=1,nneighs
           vdotr = 0.
           DO k=1,ndims
              vdotr = vdotr + ((v(k,i) - v(k,nearIDs(j,i))) * (r(k,i) - r(k,nearIDs(j,i))))
           END DO
           IF (vdotr < 0.) THEN              
              muab = (h(i) * vdotr) / ((dists(j,i)**2) * etasq) ! dists(j,i) can be used even though the viscosity formula uses a distance vector that is defined the other way around. This is because the dot product of the viscosity distance vector with itself is taken, which is independent of the direction in which the distance is defined as the power of 2 removes the difference in sign
              avsndspd = 0.5 * (sndspd(i) + sndspd(nearIDs(j,i)))
              avdens = 0.5 * (dens(i) + dens(nearIDs(j,i)))
              visc(j,i) = (((-1. * viscalpha * avsndspd * muab) + (viscbeta * (muab**2))) / avdens)
      
           ELSE
              visc(j,i) = 0.
           END IF
        END DO
     END DO ! Artificial viscosity from the nearest neighbours for each particle should now have been computed     

     DO i=1,n ! Accceleration needs to be calculated after everything else as the contribution to the acceleration of particle i from neighbour j depends on the pressures and densities at both i and j. Components of acceleration can be found by multiplying the acceleration by dr(k)/a, where a is the total separation and dr(k) is the separation along a given axis k

        DO j=1, nneighs
           tempDist = 0.
           tempAcc = m(nearIDs(j,i)) * ((press(nearIDs(j,i))/(dens(nearIDs(j,i))**2)) + (press(i)/(dens(i)**2))) * kerndiv(j,i)
           DO k=1,ndims
              dr(k) = r(k,i) - r(k,nearIDs(j,i))
              tempDist = tempDist + (dr(k)**2)     
           END DO
           tempDist = SQRT(tempDist)
           
           ! The following assumes that the acceleration array has been set to zero at the start of the main DO loop
           DO k=1,ndims
              a(k,i) = a(k,i) - (tempAcc * (dr(k)/tempDist))
           END DO
           
        END DO

        
     END DO ! Accelerations should now have been computed for each particle
     
     DO i=1,n ! Advancing the simulation       
        DO k=1,ndims
           r(k,i) = r(k,i) + (v(k,i) * dt)           
           v(k,i) = v(k,i) + (a(k,i) * dt)
        END DO
     END DO


     t = t + dt
     
     IF (t>=tlim) THEN
        RETURN
        EXIT
     END IF
  END DO

END SUBROUTINE HydroSim

SUBROUTINE sprSetup(n, ndims, rad, r)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, ndims
  INTEGER :: i, k
  DOUBLE PRECISION, INTENT(IN) :: rad
  DOUBLE PRECISION, INTENT(OUT) :: r(1:ndims,1:n)
  DOUBLE PRECISION :: tempDist, ranPos(1:ndims)

  tempDist = 0.
  ranPos = 0.

  DO i=1,n
10   CONTINUE
     DO k=1,ndims
        CALL random_NUMBER(ranPos(k))
     END DO

     r(:,i) = 2. * (ranPos(:) - 0.5) ! Modifies all positions so that some particles will have negative coordinates
     r(:,i) = r(:,i) * rad ! Scales all of the positions to account for the possibility of a starting sphere with a radius that isn't 1m
     tempDist = 0.
     DO k=1,ndims
        tempDist = tempDist + (r(k,i)**2)
     END DO
     tempDist = SQRT(tempDist)
     IF (tempDist > rad) THEN
        GOTO 10 ! Re-randomised initial position if particle is outside the sphere of radius rad centred on the origin
     END IF
  END DO

END SUBROUTINE sprSetup

SUBROUTINE cubSetup(n, ndims, lim, r)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, ndims
  INTEGER :: i, k
  DOUBLE PRECISION, INTENT(IN) :: lim(1:ndims)
  DOUBLE PRECISION, INTENT(OUT) :: r(1:ndims,1:n)
  DOUBLE PRECISION :: ranPos(1:ndims)

  ranPos = 0.

  DO i=1,n
     DO k=1,ndims
        CALL random_NUMBER(ranPos(k))
        r(k,i) = ranPos(k) * lim(k) ! Scales positions based on the size of the initial cuboid
     END DO
  END DO

  RETURN
END SUBROUTINE cubSetup

        
        
        
        
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


















