
C#==============================================================================
C  Diffusion in a 1-dimensional finite difference grid 
C  all inputs are vectors
C#==============================================================================

      SUBROUTINE diff1D (N, C, Cup, Cdown, fluxup, fluxdown, aup, adown,        &
     &                   BcUp, BcDown, D, VF, A, dx, dxaux,                     &
     &                   Flux, dC) 
      IMPLICIT NONE
      INTEGER N                  ! length of C
C input
      DOUBLE PRECISION C(N)

C Boundary concentrations (used if Bc..=2,4), fluxes (used if Bc= 1) 
C and convection coeff (used if Bc=4)
      DOUBLE PRECISION Cup, Cdown, fluxup, fluxdown, aup, adown

C Diffusion, volume fraction, surface Area
      DOUBLE PRECISION D(N+1), VF(N+1), A(N+1)

C grid size, distance from mid to mid
      DOUBLE PRECISION dx(N), dxaux(N+1)

C boundary concitions (1= flux, 2=conc, 3 = 0-grad, 4=convect)
      INTEGER BcUp, BcDown   

C output: fluxes and rate of change
      DOUBLE PRECISION Flux(N+1), dC(N)

C locals 
      INTEGER I
      DOUBLE PRECISION AVF

C -------------------------------------------------------------------------------

C Flux - first internal cells

      DO I = 2,N
        Flux(I) = -VF(I)*D(I) * (C(I)-C(I-1)) /dxaux(I)
      ENDDO

C Then the outer cells 
C upstream boundary
      IF (BcUp .EQ. 1) THEN
        Flux(1) = fluxup

      ELSE IF (BcUp .EQ. 2) THEN
        Flux(1) = -VF(1)*D(1) * (C(1)-Cup) /dxaux(1)

      ELSE IF (BcUp .EQ. 3) THEN
        Flux(1) = 0.D0

      ELSE 
        Flux(1) = aup * (Cup - C(1))

      ENDIF
       

C downstream boundary
      IF (BcDown .EQ. 1) THEN
        Flux(N+1) = fluxdown

      ELSE IF (BcDown .EQ. 2) THEN
        Flux(N+1) = -VF(N+1)*D(N+1) * (Cdown-C(N)) /dxaux(N+1)

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) =0.D0

      ELSE
        Flux(N+1) = -adown * (Cdown-C(N))

      ENDIF


C Rate of change = negative flux gradient
      DO I = 1,N
        AVF   = 0.25 * (A(I)+A(I+1)) * (VF(I)+VF(I+1))
        dC(I) = -(A(I+1)*Flux(I+1) - A(I)*Flux(I)) / AVF / dx(I)
      ENDDO
    
      RETURN
      END SUBROUTINE diff1D




C#==============================================================================
C  Diffusion in a one-dimensional finite difference grid 
C  all inputs except C are scalars
C#==============================================================================

      SUBROUTINE diff1Dscal (N, C, Cup, Cdown, fluxup, fluxdown,                &
     &                       aup, adown,BcUp, BcDown, D, VF, dx,                &
     &                       Flux, dC) 
      IMPLICIT NONE
      INTEGER N                  ! length of C
C input
      DOUBLE PRECISION C(N), Cup, Cdown, fluxup, fluxdown, aup, adown
      DOUBLE PRECISION D, VF, dx

C boundary conditions (1= flux, 2 = conc, 3 = 0 grad, 4= convective)
      INTEGER BcUp, BcDown   

C output
      DOUBLE PRECISION Flux(N+1), dC(N)

C locals
      INTEGER I

C -------------------------------------------------------------------------------

C Flux - first internal cells
      DO I = 2,N
        Flux(I) = -VF*D * (C(I)-C(I-1)) /dx
      ENDDO

C Then the outer cells
C upstream
      IF (BcUp .EQ. 1) THEN
        Flux(1) = fluxup

      ELSE IF (BcUp .EQ. 2) THEN
        Flux(1) = -VF*D * (C(1)-Cup) /dx

      ELSE IF (BcUp .EQ. 3) THEN
        Flux(1) = 0.D0

      ELSE
        Flux(1) = -aup * (C(1)-Cup)

      ENDIF

C downstream       
      IF (BcDown .EQ. 1) THEN
        Flux(N+1) = fluxdown

      ELSE IF (BcDown .EQ. 2) THEN
        Flux(N+1) = -VF*D * (Cdown-C(N)) /dx

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) = 0.D0

      ELSE
        Flux(N+1) = -adown * (Cdown-C(N))

      ENDIF

C rate of change = negative flux gradient
      DO I = 1,N
        dC(I) = -(Flux(I+1) - Flux(I)) /VF /dx
      ENDDO
    
      RETURN
      END
