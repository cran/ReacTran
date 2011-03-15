C#==============================================================================
C  Diffusion in a 2-dimensional finite difference grid 
C  all inputs are vectors/matrices
C#==============================================================================

      SUBROUTINE diff2D (Nx, Ny, C, Cxup, Cxdown, Cyup, Cydown,                 &
     &                   Fxup, Fxdown, Fyup, Fydown,                            &
     &                   axup, axdown, ayup, aydown,                            &
     &                   BcXup, BcXdown, BcYup, BcYdown,                        &
     &                   D_x, D_y, VF_x, VF_y, A_x, A_y,                        &
     &                   dx, dxaux, dy, dyaux, FluxX, FluxY, dC) 
      IMPLICIT NONE
      INTEGER Nx,Ny                  ! dimension of C
C input
      DOUBLE PRECISION C(Nx,Ny)

C when Bc.. = 2; boundary concentrations should be a vector then, else a scalar
      DOUBLE PRECISION Cxup(*),Cxdown(*),Cyup(*),Cydown(*)

C when Bc.. = 1; boundary fluxes should be a vector, else a scalar
      DOUBLE PRECISION Fxup(*), Fxdown(*), Fyup(*), Fydown(*)

C when Bc.. = 4; should be a vector, else a scalar
      DOUBLE PRECISION axup(*), axdown(*), ayup(*), aydown(*)

C diffusion coefficients, volume fractions, surface Areas
      DOUBLE PRECISION D_x(Nx+1,Ny), D_y(Nx,Ny+1), VF_x(Nx+1,Ny),               &
     &                 VF_y(Nx,Ny+1), A_x(Nx+1,Ny),A_y(Nx,Ny+1)

C Grid sizes and auxillary distances (from mid to mid)
      DOUBLE PRECISION dx(Nx,Ny),dxaux(Nx+1,Ny),dy(Nx,Ny),dyaux(Nx,Ny+1)

C Boundaries: 1= flux, 2 = conc, 3 = 0 grad, 4 = convection
      INTEGER  BcXup, BcXdown, BcYup, BcYdown  

C Output: Fluxes in X- and Y-direction, rate of change
      DOUBLE PRECISION FluxX(Nx+1, Ny),  FluxY(Nx, Ny+1), dC(Nx,Ny)

C locals
      DOUBLE PRECISION AVx, AVy 
      INTEGER I, J 
C -------------------------------------------------------------------------------

C ## transport in X-direction ##
C First diffusion internal cells

      DO I = 2, Nx
        DO J = 1, Ny
          FluxX(I,J) = -VF_x(I,J)*D_x(I,J)*(C(I,J)-C(I-1,J)) /dxaux(I,J)   
        ENDDO 
      ENDDO
      FluxX(1,:)    = 0.D0
      FluxX(Nx+1,:) = 0.D0

C Then the outer cells - depending on boundary
C upstream X-boundary

      IF (BcXup .EQ. 1) THEN
        DO J = 1, Ny 
          FluxX(1,J) = Fxup(J)
        ENDDO

      ELSE IF (BcXup .EQ. 2) THEN     
        DO J = 1, Ny 
          FluxX(1,J) = -VF_x(1,J)*D_x(1,J)*(C(1,J)-Cxup(J))/dxaux(1,J)
        ENDDO

      ELSE IF (BcXup .EQ. 3) THEN
        FluxX(1,:) = 0.D0

      ELSE IF (BcXup .EQ. 4) THEN     
        DO J = 1, Ny 
          FluxX(1,J) = -axup(J)*(C(1,J)-Cxup(J))
        ENDDO

      ENDIF
       
C downstream X-boundary

      IF (BcXdown .EQ. 1) THEN
        DO J = 1, Ny 
          FluxX(Nx+1,J) = Fxdown(J)
        ENDDO

      ELSE IF (BcXdown .EQ. 2) THEN     
        DO J = 1, Ny 
          FluxX(Nx+1,J) = -VF_x(Nx+1,J)*D_x(Nx+1,J)*(Cxdown(J)-C(Nx,J))          &
     &                                           /dxaux(Nx+1,J)
        ENDDO

      ELSE IF (BcXdown .EQ. 3) THEN
        FluxX(Nx+1,:) = 0.D0

      ELSE IF (BcXdown .EQ. 4) THEN     
        DO J = 1, Ny 
          FluxX(Nx+1,J) = -axdown(J)*(Cxdown(J)-C(Nx,J))
        ENDDO

      ENDIF

C Rate of change in X-direction

      DO I = 1,Nx
        DO J = 1,Ny
          AVx=0.25*(A_x(I,J)+A_x(I+1,J))*(VF_x(I,J)+VF_x(I+1,J))*dx(I,J)
          dC(I,J) = -(A_x(I+1,J)*FluxX(I+1,J)-A_x(I,J)*FluxX(I,J))/AVx            
        ENDDO
      ENDDO

C ## transport in Y-direction ##
C First diffusion internal cells

      DO I = 1, Nx
        DO J = 2, Ny
          FluxY(I,J) = -VF_y(I,J)*D_y(I,J)*(C(I,J)-C(I,J-1)) /dyaux(I,J) 
        ENDDO 
      ENDDO
      FluxY(:,1)    = 0.D0
      FluxY(:,Ny+1) = 0.D0
       
C upstream Y-boundary

      IF (BcYup.EQ. 1) THEN
        DO I = 1, Nx
          FluxY(I,1) = Fyup(I)
        ENDDO

      ELSE IF (BcYup .EQ. 2) THEN
        DO I = 1, Nx
          FluxY(I,1) = -VF_y(I,1)*D_y(I,1)*(C(I,1)-Cyup(I))/dyaux(I,1)
        ENDDO

      ELSE IF (BcYup .EQ. 3) THEN
        FluxY(:,1) = 0.D0

      ELSE IF (BcYup .EQ. 4) THEN
        DO I = 1, Nx
          FluxY(I,1) = -ayup(I)*(C(I,1)-Cyup(I))
        ENDDO

      ENDIF

C downstream Y-boundary

      IF (BcYdown.EQ. 1) THEN
        DO I = 1, Nx
          FluxY(I,Ny+1) = Fydown(I)
        ENDDO

      ELSE IF (BcYdown .EQ. 2) THEN
        DO I = 1, Nx
          FluxY(I,Ny+1) = -VF_y(I,Ny+1)*D_y(I,Ny+1)*(Cydown(I)-C(I,Ny))          &
     &                                           /dyaux(I,Ny+1)
        ENDDO

      ELSE IF (BcYdown .EQ. 3) THEN
        FluxY(:,Ny+1) = 0.D0

      ELSE IF (BcYdown .EQ. 4) THEN
        DO I = 1, Nx
          FluxY(I,Ny+1) = -aydown(I)*(Cydown(I)-C(I,Ny)) 
        ENDDO

      ENDIF

C rate of change, negative flux gradient

      DO I = 1,Nx
        DO J = 1,Ny
          AVy=0.25*(A_y(I,J)+A_y(I,J+1))*(VF_y(I,J)+VF_y(I,J+1))*dy(I,J)

          dC(I,J) = dC(I,J)                                                     &
     &              -(A_y(I,J+1)*FluxY(I,J+1)-A_y(I,J)*FluxY(I,J))/AVy 
        ENDDO
      ENDDO
    
      RETURN
      END SUBROUTINE diff2D

C#==============================================================================
C  Diffusion in a 2-dimensional finite difference grid - inputs are scalars
C#==============================================================================

      SUBROUTINE diff2Dscal (Nx, Ny, C, Cxup, Cxdown, Cyup, Cydown,             &
     &                   Fxup, Fxdown, Fyup, Fydown,                            &
     &                   axup, axdown, ayup, aydown,                            &
     &                   BcXup, BcXdown, BcYup, BcYdown,                        &
     &                   D, VF, dx, dy, FluxX, FluxY, dC) 
      IMPLICIT NONE
      INTEGER Nx, Ny                  ! dimension of C
C input
      DOUBLE PRECISION C(Nx,Ny),Cxup, Cxdown, Cyup, Cydown 
      DOUBLE PRECISION Fxup, Fxdown, Fyup, Fydown
      DOUBLE PRECISION axup, axdown, ayup, aydown

      DOUBLE PRECISION D, VF, dx, dy

C boundaries: 1= flux, 2 = conc, 3 = 0 grad, 4=convection
      INTEGER BcXup, BcXdown, BcYup, BcYdown  

C output
      DOUBLE PRECISION FluxX(Nx+1, Ny), FluxY(Nx, Ny+1), dC(Nx,Ny)

C locals
      INTEGER I, J

C -------------------------------------------------------------------------------

C ## transport in X-direction ##
C First diffusion internal cells

      DO I = 2, Nx
        DO J = 1, Ny
          FluxX(I,J) = -VF*D*( (C(I,J)-C(I-1,J)) /dx)
        ENDDO 
      ENDDO

      FluxX(1,:)    = 0.D0
      FluxX(Nx+1,:) = 0.D0

C Then the outer cells - depending on boundary
C upstream X-boundary

      IF (BcXup .EQ. 1) THEN
        DO J = 1, Ny 
          FluxX(1,J) = Fxup
        ENDDO

      ELSE IF (BcXup .EQ. 2) THEN     
        DO J = 1, Ny 
          FluxX(1,J) = -VF*D*(C(1,J)-Cxup)/dx
        ENDDO

      ELSE IF (BcXup .EQ. 3) THEN
        FluxX(1,:) = 0.D0

      ELSE IF (BcXup .EQ. 4) THEN     
        DO J = 1, Ny 
          FluxX(1,J) = -axup*(C(1,J)-Cxup)
        ENDDO

      ENDIF
       
C downstream X-boundary

      IF (BcXdown .EQ. 1) THEN
        DO J = 1, Ny 
          FluxX(Nx+1,J) = Fxdown
        ENDDO

      ELSE IF (BcXdown .EQ. 2) THEN     
        DO J = 1, Ny 
          FluxX(Nx+1,J) = -VF*D*(Cxdown-C(Nx,J)) /dx
        ENDDO

      ELSE IF (BcXdown .EQ. 3) THEN
        FluxX(Nx+1,:) = 0.D0

      ELSE IF (BcXdown .EQ. 4) THEN     
        DO J = 1, Ny 
          FluxX(Nx+1,J) = -axdown*(Cxdown-C(Nx,J))
        ENDDO

      ENDIF

C Rate of change = negative flux gradient

      DO I = 1,Nx
        DO J = 1,Ny
          dC(I,J) = -(FluxX(I+1,J)-FluxX(I,J))/VF/dx   
        ENDDO
      ENDDO

C ## transport in Y-direction ##
C first internal cells

      DO I = 1, Nx
        DO J = 2, Ny
          FluxY(I,J) = -VF*D*( (C(I,J)-C(I,J-1)) /dy)
        ENDDO 
      ENDDO

      FluxY(:,1)    = 0.D0
      FluxY(:,Ny+1) = 0.D0
       
C upstream Y-boundary

      IF (BcYup.EQ. 1) THEN
        DO I = 1, Nx
          FluxY(I,1) = Fyup 
        ENDDO

      ELSE IF (BcYup .EQ. 2) THEN
        DO I = 1, Nx
          FluxY(I,1) = -VF*D*(C(I,1)-Cyup)/dy
        ENDDO
 
      ELSE IF (BcYup .EQ. 3) THEN
        FluxY(:,1) = 0.D0
 
      ELSE IF (BcYup .EQ. 4) THEN
        DO I = 1, Nx
          FluxY(I,1) = -ayup*(C(I,1)-Cyup)
        ENDDO
 
      ENDIF

C downstream Y-boundary

      IF (BcYdown.EQ. 1) THEN
        DO I = 1, Nx
          FluxY(I,Ny+1) = Fydown
        ENDDO

      ELSE IF (BcYdown .EQ. 2) THEN
        DO I = 1, Nx
          FluxY(I,Ny+1) = -VF*D*(Cydown-C(I,Ny))/dy
        ENDDO

      ELSE IF (BcYdown .EQ. 3) THEN
        FluxY(:,Ny+1) = 0.D0

      ELSE IF (BcYdown .EQ. 4) THEN
        DO I = 1, Nx
          FluxY(I,Ny+1) = -aydown*(Cydown-C(I,Ny)) 
        ENDDO

      ENDIF


C Rate of change = negative flux gradient

      DO I = 1,Nx
        DO J = 1,Ny
          dC(I,J) = dC(I,J) -(FluxY(I,J+1)-FluxY(I,J))/VF/dy
        ENDDO
      ENDDO
    
      RETURN
      END SUBROUTINE diff2Dscal


C#==============================================================================
C  Diffusion in a 2-dimensional polar coordinate system - all inputs are vectors
C  C is matrix
C#==============================================================================

      SUBROUTINE diffpolar (Nr, Ntet, C, Crup, Crdown, Ctetup, Ctetdown,        &
     &                   Frup, Frdown, Ftetup, Ftetdown,                        &
     &                   Bcrup, Bcrdown, Bctetup, Bctetdown,                    &
     &                   D_r, D_tet, r, tet, Fluxr, Fluxtet, dC) 
      IMPLICIT NONE
C dimension of 
      INTEGER Nr, Ntet                   
C concentration
      DOUBLE PRECISION C(Nr,Ntet)

C used when Bc.. = 2, should be a vector then, else one number
      DOUBLE PRECISION Crup(*), Crdown(*), Ctetup(*), Ctetdown(*)

C used when Bc.. = 1, should be a vectorthen , else one number
      DOUBLE PRECISION Frup(*), Frdown(*), Ftetup(*), Ftetdown(*)

      DOUBLE PRECISION D_r(Nr+1), D_tet(Ntet+1),                                &
     &                 r(Nr+1), tet(Ntet+1), rc(Nr), tetc(Ntet)
      DOUBLE PRECISION dr(Nr),draux(Nr+1),dtet(Ntet),dtetaux(Ntet+1)

c     boundaries: 1= flux, 2 = conc, 3 = zero- grad, 5 = cyclic boundary
      INTEGER  BcRup, BcRdown, BcTetup, BcTetdown  

C output
      DOUBLE PRECISION Fluxr(Nr+1,Ntet), Fluxtet(Nr,Ntet+1), dC(Nr,Ntet)

C locals
      INTEGER I, J 

C -------------------------------------------------------------------------------

C The grid

      DO I = 1, Nr
        rc(I) = 0.5*(r(I)+r(I+1))
        dr(I) = r(I+1) - r(I)
      ENDDO

      DO I = 1, Nr-1
        draux(I+1) = (rc(I+1) - rc(I))
      ENDDO
      draux(1)    = rc(1)-r(1)
      draux(Nr+1) = r(Nr+1) - rc(Nr)
      

      DO I = 1, Ntet
        tetc(I) = 0.5*(tet(I)+tet(I+1))
        dtet(I) = tet(I+1) - tet(I)
      ENDDO

      DO I = 1, Ntet-1
        dtetaux(I+1) = (tetc(I+1) - tetc(I))
      ENDDO
      dtetaux(1)      = tetc(1)-tet(1)
      dtetaux(Ntet+1) = tet(Ntet+1) - tetc(Ntet)
      
C ## transport in r-direction ##
C First diffusion internal cells
      
      DO I = 2, Nr
        DO J = 1, Ntet
          Fluxr(I,J) = -D_r(I)*(C(I,J)-C(I-1,J)) /draux(I)   
        ENDDO 
      ENDDO

      Fluxr(1,:)    = 0.D0
      Fluxr(Nr+1,:) = 0.D0

C Then the outer cells - depending on boundary type
C upstream r-boundary
      
      IF (Bcrup .EQ. 1) THEN
        DO J = 1, Ntet 
          Fluxr(1,J) = Frup(J)
        ENDDO

      ELSE IF (Bcrup .EQ. 2) THEN     
        DO J = 1, Ntet 
          Fluxr(1,J) = -D_r(1) * (C(1,J)-Crup(J)) /draux(1)
        ENDDO

      ELSE IF (Bcrup .EQ. 3) THEN
        Fluxr(1,:) = 0.D0

      ELSE IF (Bcrup .EQ. 5) THEN
        DO J = 1, Ntet
          Crup(J) = (C(1,J)*draux(Nr+1) + C(Nr,J)*draux(1))  /                  &
     &              (draux(1)+draux(Nr+1))
          Crdown(J) = Crup(J)
          Fluxr(1,J) = -D_r(1)*(C(1,J)-Crup(J))/draux(1)
        ENDDO

      ENDIF
       
C downstream r-boundary

      IF (Bcrdown .EQ. 1) THEN
        DO J = 1, Ntet 
          Fluxr(Nr+1,J) = Frdown(J)
        ENDDO

      ELSE IF (Bcrdown .EQ. 2 .OR. Bcrdown .EQ. 5) THEN       ! if 5: Crdown already estimated
        DO J = 1, Ntet 
          Fluxr(Nr+1,J) = -D_r(Nr+1) * (Crdown(J)-C(Nr,J)) /draux(Nr+1)
        ENDDO

      ELSE IF (Bcrdown .EQ. 3) THEN
        Fluxr(Nr+1,:) = 0.D0

      ENDIF

C Rate of change in r-direction

      DO I = 1,Nr
        DO J = 1,Ntet
          dC(I,J) = -(Fluxr(I+1,J)*r(I+1)-Fluxr(I,J)*r(I))/dr(I)/rc(I)            
        ENDDO
      ENDDO

C ## transport in teta -direction ##
C First diffusion internal cells

      DO I = 1, Nr
        DO J = 2, Ntet
          Fluxtet(I,J) = -D_tet(J)*(C(I,J)-C(I,J-1)) /dtetaux(J)/rc(I) 
        ENDDO 
      ENDDO

      Fluxtet(:,1)    = 0.D0
      Fluxtet(:,Ntet+1) = 0.D0
       
C upstream teta-boundary

      IF (BcTetup.EQ. 1) THEN
        DO I = 1, Nr
          Fluxtet(I,1) = Ftetup(I)
        ENDDO

      ELSE IF (BcTetup .EQ. 2) THEN
        DO I = 1, Nr
          Fluxtet(I,1) = -D_tet(1)*(C(I,1)-Ctetup(I))/dtetaux(1)/rc(I) 
        ENDDO

      ELSE IF (BcTetup .EQ. 3) THEN
        Fluxtet(:,1) = 0.D0

      ELSE IF (BcTetup .EQ. 5) THEN
        DO I = 1, Nr
          Ctetup(I) = (C(I,1)*dtetaux(Ntet+1) + C(I,Ntet)*dtetaux(1)) /         &
     &              (dtetaux(1)+dtetaux(Ntet+1))
          Ctetdown(I) = Ctetup(I)
          Fluxtet(I,1) = -D_tet(1)*(C(I,1)-Ctetup(I))/dtetaux(1)/rc(I) 
        ENDDO

      ENDIF

C downstream teta-boundary

      IF (BcTetdown.EQ. 1) THEN
        DO I = 1, Nr
          Fluxtet(I,Ntet+1) = Ftetdown(I)
        ENDDO

      ELSE IF (BcTetdown .EQ. 2 .OR. BcTetdown .EQ. 5) THEN
        DO I = 1, Nr
          Fluxtet(I,Ntet+1) = -D_tet(Ntet+1)*(Ctetdown(I)-C(I,Ntet))            & 
     &                                           /dtetaux(Ntet+1)/rc(I) 
        ENDDO

      ELSE IF (BcTetdown .EQ. 3) THEN
        Fluxtet(:,Ntet+1) = 0.D0

      ENDIF

C rate of change, negative flux gradient

      DO I = 1,Nr
        DO J = 1,Ntet

          dC(I,J) = dC(I,J)                                                     &
     &              -(Fluxtet(I,J+1)-Fluxtet(I,J))/dtet(J)/rc(I) 
        ENDDO
      ENDDO
    
      RETURN
      END


