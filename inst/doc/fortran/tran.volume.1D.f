C#==============================================================================
C transport in a 1-dimensional finite volume grid (no lateral input)
C all inputs are vectors, upstream weighted advection; 
C Flow is one value, no lateral input
C compile        system("R CMD SHLIB tran.1D.vol.f")
C load           dyn.load("tran.1D.vol.dll")
C#==============================================================================

      SUBROUTINE tran1dvol (N, C, Cup, Cdown, Fup, Fdown,                        &
     &                   BcUp, BcDown, D, Flow, Volume,                          &
     &                   Flux, dC, Fin, Fout) 
      IMPLICIT NONE
      INTEGER N                  ! length of C
C input concentration
      DOUBLE PRECISION C(N)

C Boundary concentrations (used if Bc = 2), fluxes (used if Bc = 1) 
      DOUBLE PRECISION Cup, Cdown, Fup, Fdown 

C Diffusion, advection weighing, Volume
      DOUBLE PRECISION D(N+1),  Volume(N)

C Flows
      DOUBLE PRECISION Flow 

C boundary concitions (1= flux, 2=conc, 3 = 0-gradient)
      INTEGER BcUp, BcDown   

C output: flux vector, rate of change, and flux in and out (mass/sec)
      DOUBLE PRECISION Flux(N+1), dC(N), Fin, Fout

C locals 
      INTEGER I

C -------------------------------------------------------------------------------

C Flux - first internal cells

      DO I = 2,N
        Flux(I) = -D(I) * (C(I)-C(I-1))  

        IF (Flow  > 0) THEN
          Flux (I) = Flux(I) + Flow *C(I-1)
        ELSE IF (Flow  < 0) THEN
          Flux (I) = Flux(I) - Flow *C(I)
        ENDIF 

      ENDDO

C Then the outer cells 
C upstream boundary
      IF (BcUp .EQ. 1) THEN
        Flux(1) = Fup

      ELSE IF (BcUp .EQ. 2) THEN
        Flux(1) = -D(1) * (C(1)-Cup) 

	      IF (Flow  > 0) THEN
          Flux (1) = Flux(1) + Flow * Cup
	      ELSE IF (Flow  < 0) THEN
          Flux (1) = Flux(1) + Flow * C(1)
        ENDIF 

      ELSE IF (BcUp .EQ. 3) THEN
        Flux(1) = 0.D0

      ELSE 
        Flux(1) = 0.D0

      ENDIF

C downstream boundary
      IF (BcDown .EQ. 1) THEN
        Flux(N+1) = Fdown

      ELSE IF (BcDown .EQ. 2) THEN
        Flux(N+1) = -D(N+1) * (Cdown-C(N)) 

	      IF (Flow  > 0) THEN
          Flux (N+1) = Flux(N+1) + Flow*C(N)
    	  ELSE IF (Flow < 0) THEN
          Flux (N+1) = Flux(N+1) + Flow*Cdown
        ENDIF 

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) =0.D0

      ELSE
        Flux(N+1) =0.D0

      ENDIF


C Rate of change = negative flux gradient
      DO I = 1,N
        dC(I) = -(Flux(I+1) - Flux(I)) / Volume(I)
      ENDDO

      Fin  = Flux(1) 
      Fout = Flux(N+1) 
	    
      RETURN
      END SUBROUTINE tran1Dvol


C#==============================================================================
C transport in a 1-dimensional finite volume grid (no lateral input)
C all inputs are vectors, upstream weighted advection; Flow is a vector
C#==============================================================================

      SUBROUTINE tran1dvolb (N, C, Cup, Cdown, Fup, Fdown,                        &
     &                   BcUp, BcDown, D, Flow, Clat, Volume,                     &
     &                   Flux, dC, Fin, Fout) 
      IMPLICIT NONE
      INTEGER N                  ! length of C
C input concentration
      DOUBLE PRECISION C(N)

C Boundary concentrations (used if Bc.= 2), fluxes (used if Bc = 1) 
      DOUBLE PRECISION Cup, Cdown, Fup, Fdown 

C Diffusion, advection weighing, Volume
      DOUBLE PRECISION D(N+1), AFDW(N+1), Volume(N)

C Flows at upstream interface, lateral and lateral concentrations
      DOUBLE PRECISION Flow(N+1), Clat(N)

C boundary concitions (1= flux, 2=conc, 3 = 0-gradient)
      INTEGER BcUp, BcDown   

C output: flux vector, rate of change, and flux in and out (mass/sec)
      DOUBLE PRECISION Flux(N+1), dC(N), Flat(N), Fin, Fout

C locals 
      INTEGER I

C -------------------------------------------------------------------------------

C Flux - first internal cells

      DO I = 2,N
        Flux(I) = -D(I) * (C(I)-C(I-1))  

        IF (Flow(I) > 0) THEN
          Flux (I) = Flux(I) + Flow(I) *C(I-1)
        ELSE IF (Flow(I) < 0) THEN
          Flux (I) = Flux(I) - Flow(I) *C(I)
        ENDIF 

      ENDDO

C Then the outer cells 
C upstream boundary
      IF (BcUp .EQ. 1) THEN
        Flux(1) = Fup

      ELSE IF (BcUp .EQ. 2) THEN
        Flux(1) = -D(1) * (C(1)-Cup) 

        IF (Flow(1) > 0) THEN
          Flux (1) = Flux(1) + Flow(1)*Cup
        ELSE IF (Flow(1) < 0) THEN
          Flux (1) = Flux(1) + Flow(1)*C(1)
        ENDIF 

      ELSE IF (BcUp .EQ. 3) THEN
        Flux(1) = 0.D0

      ELSE 
        Flux(1) = 0.D0

      ENDIF
       

C downstream boundary
      IF (BcDown .EQ. 1) THEN
        Flux(N+1) = Fdown

      ELSE IF (BcDown .EQ. 2) THEN
        Flux(N+1) = -D(N+1) * (Cdown-C(N)) 

        IF (Flow(N+1) > 0) THEN
          Flux (N+1) = Flux(N+1) + Flow(N+1)*C(N)
        ELSE IF (Flow(N+1) < 0) THEN
          Flux (N+1) = Flux(N+1) + Flow(N+1)*Cdown
        ENDIF 

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) =0.D0

      ELSE
        Flux(N+1) =0.D0

      ENDIF


C Rate of change = negative flux gradient
      DO I = 1,N
        FLat(I) = (Flow(I) - Flow(I+1))*Clat(I)
        dC(I) = -(Flux(I+1) - Flux(I) + Flat(I)) / Volume(I) 
      ENDDO
    
      Fin  = Flux(1) 
      Fout = Flux(N+1) 
	    
      RETURN
      END SUBROUTINE tran1Dvolb

