
C#==============================================================================
C Diffusion in a 3-dimensional spherical coordinate system (r, theta, phi) 
C all inputs are vectors/scalars
C#==============================================================================

      SUBROUTINE diffspher (Nr, Ntet, Nphi, C,                                   &
     &             Crup, Crdown, Ctetup, Ctetdown, Cphiup, Cphidown,             &
     &             Frup, Frdown, Ftetup, Ftetdown, Fphiup, Fphidown,             &
     &             Bcrup, Bcrdown, Bctetup, Bctetdown,                           &
     &             Bcphiup, Bcphidown, D_r, D_tet, D_phi,                        &
     &             r, tet, phi, FluxrUp, FluxrDown,                              &
     &             FluxtetUp, FluxtetDown, FluxphiUp, FluxPhiDown, dC) 

      IMPLICIT NONE
C dimension of 
      INTEGER Nr, Ntet, Nphi                   
C concentration
      DOUBLE PRECISION C(Nr, Ntet, Nphi)

C used when Bc.. = 2
      DOUBLE PRECISION Crup, Crdown, Ctetup, Ctetdown, Cphiup, Cphidown

C used when Bc.. = 1
      DOUBLE PRECISION Frup, Frdown, Ftetup, Ftetdown, Fphiup, Fphidown

      DOUBLE PRECISION D_r(Nr+1), D_tet(Ntet+1), D_phi(Nphi+1),                  &
     &                 r(Nr+1), tet(Ntet+1), phi(Nphi+1),                        &
     &                 rc(Nr), tetc(Ntet), phic(Nphi)
      DOUBLE PRECISION dr(Nr),draux(Nr+1),dtet(Ntet),dtetaux(Ntet+1),            &
     &                 dphi(Nphi), dphiaux(Nphi+1)

c     boundaries: 1= flux, 2 = conc, 3 = zero- grad, 5 = cyclic boundary
      INTEGER  BcRup, BcRdown, BcTetup, BcTetdown, BcPhiup, BcPhidown

C output
      DOUBLE PRECISION FluxrUp(Ntet,Nphi), FluxrDown(Ntet,Nphi)
      DOUBLE PRECISION Fluxtetup(Nr,Nphi), FluxtetDown(Nr,Nphi)
      DOUBLE PRECISION Fluxphiup(Nr,Ntet), FluxphiDown(Nr,Ntet)
      DOUBLE PRECISION dC(Nr,Ntet,Nphi)

C locals
      INTEGER I, J, K 
      DOUBLE PRECISION Flux(Nr+1, Ntet+1, Nphi+1), Cbound

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
      

      DO I = 1, Nphi
        phic(I) = 0.5*(phi(I)+phi(I+1))
        dphi(I) = phi(I+1) - phi(I)
      ENDDO

      DO I = 1, Nphi-1
        dphiaux(I+1) = (phic(I+1) - phic(I))
      ENDDO
      dphiaux(1)      = phic(1)-phi(1)
      dphiaux(Nphi+1) = phi(Nphi+1) - phic(Nphi)

C ## transport in r-direction ##
C First diffusion internal cells
      
      DO I = 2, Nr
        DO J = 1, Ntet
          DO K = 1, Nphi
            Flux(I,J,K) = -D_r(I)*(C(I,J,K)-C(I-1,J,K)) /draux(I)   
          ENDDO
        ENDDO 
      ENDDO

      Flux(1,:,:)    = 0.D0
      Flux(Nr+1,:,:) = 0.D0

C Then the outer cells - depending on boundary type
C upstream r-boundary
      
      IF (Bcrup .EQ. 1) THEN
        DO J = 1, Ntet 
          DO K = 1, Nphi
            Flux(1,J,K) = Frup
          ENDDO
        ENDDO

      ELSE IF (Bcrup .EQ. 2) THEN     
        DO J = 1, Ntet 
          DO K = 1, Nphi
            Flux(1,J,K) = -D_r(1) * (C(1,J,K)-Crup) /draux(1)
          ENDDO  
        ENDDO

      ELSE IF (Bcrup .EQ. 3) THEN
        Flux(1,:,:) = 0.D0

      ELSE IF (Bcrup .EQ. 5) THEN
        DO J = 1, Ntet
          DO K = 1, Nphi
          Cbound = (C(1,J,K)*D_r(1)*draux(Nr+1) +                                &
     &              C(Nr,J,K)*D_r(Nr+1)*draux(1))  /                             &
     &              (draux(1)*D_r(Nr+1)+draux(Nr+1)*D_r(1))
          Flux(1,J,K)    = -D_r(1)    * (C(1,J,K)-Cbound ) /draux(1)
          Flux(Nr+1,J,K) = -D_r(Nr+1) * (Cbound-C(Nr,J,K)) /draux(Nr+1)
          ENDDO
        ENDDO

      ENDIF

      DO J = 1, Ntet
        DO K = 1, Nphi
          Fluxrup(J,K) = Flux(1,J,K)
        ENDDO
      ENDDO
       
C downstream r-boundary

      IF (Bcrdown .EQ. 1) THEN
        DO J = 1, Ntet 
          DO K = 1, Nphi
            Flux(Nr+1,J,K) = Frdown
          ENDDO
        ENDDO

      ELSE IF (Bcrdown .EQ. 2) THEN       ! if 5: Crdown already estimated
        DO J = 1, Ntet 
          DO K = 1, Nphi
            Flux(Nr+1,J,K) =-D_r(Nr+1) * (Crdown-C(Nr,J,K)) /draux(Nr+1)
          ENDDO
        ENDDO

      ELSE IF (Bcrdown .EQ. 3) THEN
        Flux(Nr+1,:,:) = 0.D0

      ENDIF

      DO J = 1, Ntet
        DO K = 1, Nphi
          Fluxrdown(J,K) = Flux(Nr+1,J,K)
        ENDDO
      ENDDO

C Rate of change in r-direction

      DO I = 1,Nr
        DO J = 1,Ntet
          DO K = 1, Nphi
            IF(rc(I) .NE. 0.D0) THEN
              dC(I,J,K) =                                                        &
     &          -(Flux(I+1,J,K)*(r(I+1)**2.)-Flux(I,J,K)*(r(I)**2.))             &
     &                  /dr(I)/(rc(I)**2.)
            ELSE
 	        dC(I,J,K) = 0.D0
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C ## transport in teta -direction ##
C First diffusion internal cells

      DO I = 1, Nr
        DO J = 2, Ntet
          DO K = 1, Nphi
            IF(rc(I) .NE. 0.D0) THEN                                             
              Flux(I,J,K) = -D_tet(J)*(C(I,J,K)-C(I,J-1,K)) /dtetaux(J)            &
     &                    /rc(I) 
	      ELSE
              Flux(I,J,K) = 0.D0
	      ENDIF
          ENDDO
        ENDDO 
      ENDDO

      Flux(:,1,:)      = 0.D0
      Flux(:,Ntet+1,:) = 0.D0
       
C upstream teta-boundary

      IF (BcTetup .EQ. 1) THEN
        DO I = 1, Nr
          DO K = 1, Nphi
            Flux(I,1,K) = Ftetup
          ENDDO
        ENDDO

      ELSE IF (BcTetup .EQ. 2) THEN
        DO I = 1, Nr
          DO K = 1, Nphi
            IF(rc(I) .NE. 0.D0) THEN                                             
              Flux(I,1,K) = -D_tet(1)*(C(I,1,K)-Ctetup)/dtetaux(1)/rc(I)
            ENDIF 
          ENDDO  
        ENDDO

      ELSE IF (BcTetup .EQ. 3) THEN
        Flux(:,1,:) = 0.D0

      ELSE IF (BcTetup .EQ. 5) THEN
        DO I = 1, Nr
          DO K = 1, Nphi
            
		  IF(rc(I) .NE. 0.D0) THEN

              Cbound = (C(I,1,K)   *D_tet(1)     *dtetaux(Ntet+1) +              &
     &               C(I,Ntet,K)*D_tet(Ntet+1)*dtetaux(1)        )/              &
     &              (D_tet(Ntet+1)*dtetaux(1)+D_tet(1)*dtetaux(Ntet+1))
              Flux(I,1,K)      = -D_tet(1)     *(C(I,1,K)-Cbound)                &
     &                                         /dtetaux(1)/rc(I) 
              Flux(I,Ntet+1,K) = -D_tet(Ntet+1)*(Cbound-C(I,Ntet,K))             & 
     &                                       /dtetaux(Ntet+1)/rc(I) 
	      ENDIF

          ENDDO
        ENDDO

      ENDIF

      DO I = 1, Nr
        DO K = 1, Nphi
          Fluxtetup(I,K) = Flux(I,1,K)
        ENDDO
      ENDDO

C downstream teta-boundary

      IF (BcTetdown .EQ. 1) THEN
        DO I = 1, Nr
          DO K = 1, Nphi
            Flux(I,Ntet+1,K) = Ftetdown
          ENDDO
        ENDDO

      ELSE IF (BcTetdown .EQ. 2) THEN
        DO I = 1, Nr
          DO K = 1, Nphi
		  IF(rc(I) .NE. 0.D0) THEN
            Flux(I,Ntet+1,K) = -D_tet(Ntet+1) * (Ctetdown-C(I,Ntet,K))           & 
     &                                           /dtetaux(Ntet+1)/rc(I) 
            ENDIF 
          ENDDO 
        ENDDO

      ELSE IF (BcTetdown .EQ. 3) THEN
        Flux(:,Ntet+1,:) = 0.D0

      ENDIF

      DO I = 1, Nr
        DO K = 1, Nphi
          Fluxtetdown(I,K) = Flux(I,Ntet+1,K)
        ENDDO
      ENDDO

C rate of change, negative flux gradient

      DO I = 1,Nr
        DO J = 1,Ntet
          DO K = 1, Nphi
            if (rc(I) * sin(tetc(J)) .NE. 0.D0)                                  &
     &      dC(I,J,K) = dC(I,J,K)                                                &
     &        -(Flux(I,J+1,K)*sin(tet(J+1))-Flux(I,J,K)*sin(tet(J)))             &
     &           /dtet(J)/rc(I)/sin(tetc(J))
          ENDDO 
        ENDDO
      ENDDO
    

C ## transport in phi-direction ##
C First diffusion internal cells

      DO I = 1, Nr
        DO J = 1, Ntet
          DO K = 2, Nphi
		  IF(rc(I) * sin(tet(J)) .NE. 0.D0) THEN
              Flux(I,J,K) = -D_phi(J)*(C(I,J,K)-C(I,J,K-1))                      &
     &                           /dphiaux(K)  /rc(I)/sin(tet(j))  
            ELSE
	        Flux(I,J,K) = 0.D0
            ENDIF
          ENDDO
        ENDDO 
      ENDDO

      Flux(:,:,1)      = 0.D0
      Flux(:,:,Nphi+1) = 0.D0
       
C upstream phi-boundary

      IF (BcPhiup .EQ. 1) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            Flux(I,J,1) = Fphiup
          ENDDO
        ENDDO

      ELSE IF (BcPhiup .EQ. 2) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            IF (rc(I) * sin(tet(J)) .NE. 0.D0)                                   &
     &        Flux(I,J,1) = -D_phi(1)*(C(I,J,1)-Cphiup)                          &
     &                            /dphiaux(1)/rc(I) /sin(tet(j))
          ENDDO  
        ENDDO

      ELSE IF (BcPhiup .EQ. 3) THEN
        Flux(:,:,1) = 0.D0

      ELSE IF (BcPhiup .EQ. 5) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            IF(rc(I) * sin(tet(j)) .NE. 0.D0) THEN                                             

          Cbound = (C(I,J,1)   *D_phi(1)     *dphiaux(Nphi+1)                    &
     &            + C(I,J,Nphi)*D_phi(Nphi+1)*dphiaux(1))/                       &
     &              (D_phi(Nphi+1)*dphiaux(1)+D_phi(1)*dphiaux(Nphi+1))
          Flux(I,J,1)      = -D_phi(1)*(C(I,J,1)-Cbound)                         &
     &                             /dphiaux(1)     /rc(I) /sin(tet(j))
          Flux(I,J,Nphi+1) = -D_phi(Nphi+1)*(Cbound-C(I,J,Nphi))                 & 
     &                             /dphiaux(Nphi+1)/rc(I) /sin(tet(j))               
            ENDIF
          ENDDO                                      
        ENDDO

      ENDIF

      DO I = 1, Nr
        DO J = 1, Ntet
          Fluxphiup(I,J) = Flux(I,J,1)
        ENDDO
      ENDDO

C downstream phi-boundary

      IF (Bcphidown .EQ. 1) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            Flux(I,J,Nphi+1) = Fphidown
          ENDDO
        ENDDO

      ELSE IF (Bcphidown .EQ. 2) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            if (rc(I) * sin(tet(J)) .NE. 0.D0)                                   &
     &      Flux(I,J,Nphi+1) = -D_phi(Nphi+1) * (Cphidown-C(I,J,Nphi))           & 
     &                             /dphiaux(Nphi+1)/rc(I) /sin(tet(j))               
          ENDDO 
        ENDDO

      ELSE IF (Bcphidown .EQ. 3) THEN
        Flux(:,:,Nphi+1) = 0.D0

      ENDIF

      DO I = 1, Nr
        DO J = 1, Ntet
          Fluxphidown(I,J) = Flux(I,J,Nphi+1)
        ENDDO
      ENDDO

C rate of change, negative flux gradient

      DO I = 1, Nr
        DO J = 1, Ntet
          DO K = 1, Nphi
            if (rc(I) .NE. 0.D0 .AND. sin(tetc(J)) .NE. 0.D0)                    &
     &      dC(I,J,K) = dC(I,J,K)                                                &
     &        -(Flux(I,J,K+1)-Flux(I,J,K))/dphi(K)/rc(I)/sin(tetc(J))

          ENDDO 
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE diffspher 


C#==============================================================================
C Diffusion in a 3-dimensional cylindrical coordinate system (r, theta, z)
C all inputs are vectors/scalars
C CHECK THIS!!!
C#==============================================================================

      SUBROUTINE diffcyl   (Nr, Ntet, Nz, C,                                     &
     &             Crup, Crdown, Ctetup, Ctetdown, Czup, Czdown,                 &
     &             Frup, Frdown, Ftetup, Ftetdown, Fzup, Fzdown,                 &
     &             Bcrup, Bcrdown, Bctetup, Bctetdown,                           &
     &             BCzUp, BCzDown, D_r, D_tet, D_z,                              &
     &             r, tet, z, FluxrUp, FluxrDown,                                &
     &             FluxtetUp, FluxtetDown, FluxzUp, FluxzDown, dC) 

      IMPLICIT NONE
C dimension of 
      INTEGER Nr, Ntet, Nz                   
C concentration
      DOUBLE PRECISION C(Nr, Ntet, Nz)

C used when Bc.. = 2
      DOUBLE PRECISION Crup, Crdown, Ctetup, Ctetdown, Czup, Czdown

C used when Bc.. = 1
      DOUBLE PRECISION Frup, Frdown, Ftetup, Ftetdown, Fzup, Fzdown

      DOUBLE PRECISION D_r(Nr+1), D_tet(Ntet+1), D_z(Nz+1),                      &
     &                 r(Nr+1), tet(Ntet+1), z(Nz+1),                            &
     &                 rc(Nr), tetc(Ntet), zc(Nz)
      DOUBLE PRECISION dr(Nr), draux(Nr+1), dtet(Ntet),dtetaux(Ntet+1),          &
     &                 dz(Nz), dzaux(Nz+1)

c     boundaries: 1= flux, 2 = conc, 3 = zero- grad, 5 = cyclic boundary
      INTEGER  BcRup, BcRdown, BcTetup, BcTetdown, BCzUp, BCzDown

C output
      DOUBLE PRECISION FluxrUp(Ntet,Nz), FluxrDown(Ntet,Nz)
      DOUBLE PRECISION Fluxtetup(Nr,Nz), FluxtetDown(Nr,Nz)
      DOUBLE PRECISION FluxZup(Nr,Ntet), FluxZDown(Nr,Ntet)
      DOUBLE PRECISION dC(Nr,Ntet,Nz)

C locals
      INTEGER I, J, K 
      DOUBLE PRECISION Flux(Nr+1, Ntet+1, Nz+1), Cbound

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
      

      DO I = 1, Nz
        zc(I) = 0.5*(z(I)+z(I+1))
        dz(I) = z(I+1) - z(I)
      ENDDO

      DO I = 1, Nz-1
        dzaux(I+1) = (zc(I+1) - zc(I))
      ENDDO
      dzaux(1)    = zc(1)-z(1)
      dzaux(Nz+1) = z(Nz+1) - zc(Nz)

C ## transport in r-direction ##
C First diffusion internal cells
      
      DO I = 2, Nr
        DO J = 1, Ntet
          DO K = 1, Nz
            Flux(I,J,K) = -D_r(I)*(C(I,J,K)-C(I-1,J,K)) /draux(I)   
          ENDDO
        ENDDO 
      ENDDO

      Flux(1,:,:)    = 0.D0
      Flux(Nr+1,:,:) = 0.D0

C Then the outer cells - depending on boundary type
C upstream r-boundary
      
      IF (Bcrup .EQ. 1) THEN
        DO J = 1, Ntet 
          DO K = 1, Nz
            Flux(1,J,K) = Frup
          ENDDO
        ENDDO

      ELSE IF (Bcrup .EQ. 2) THEN     
        DO J = 1, Ntet 
          DO K = 1, Nz
            Flux(1,J,K) = -D_r(1) * (C(1,J,K)-Crup) /draux(1)
          ENDDO  
        ENDDO

      ELSE IF (Bcrup .EQ. 3) THEN
        Flux(1,:,:) = 0.D0

      ELSE IF (Bcrup .EQ. 5) THEN
        DO J = 1, Ntet
          DO K = 1, Nz

          Cbound = (C(1,J,K)*D_r(1)*draux(Nr+1) +                                &
     &              C(Nr,J,K)*D_r(Nr+1)*draux(1))  /                             &
     &              (draux(1)*D_r(Nr+1)+draux(Nr+1)*D_r(1))

          Flux(1,J,K)    = -D_r(1)    * (C(1,J,K)-Cbound ) /draux(1)
          Flux(Nr+1,J,K) = -D_r(Nr+1) * (Cbound-C(Nr,J,K)) /draux(Nr+1)
          ENDDO
        ENDDO

      ENDIF

      DO J = 1, Ntet
        DO K = 1, Nz
          Fluxrup(J,K) = Flux(1,J,K)
        ENDDO
      ENDDO
       
C downstream r-boundary

      IF (Bcrdown .EQ. 1) THEN
        DO J = 1, Ntet 
          DO K = 1, Nz
            Flux(Nr+1,J,K) = Frdown
          ENDDO
        ENDDO

      ELSE IF (Bcrdown .EQ. 2) THEN       ! if 5: Crdown already estimated
        DO J = 1, Ntet 
          DO K = 1, Nz
            Flux(Nr+1,J,K) =-D_r(Nr+1) * (Crdown-C(Nr,J,K)) /draux(Nr+1)
          ENDDO
        ENDDO

      ELSE IF (Bcrdown .EQ. 3) THEN
        Flux(Nr+1,:,:) = 0.D0

      ENDIF

      DO J = 1, Ntet
        DO K = 1, Nz
          Fluxrdown(J,K) = Flux(Nr+1,J,K)
        ENDDO
      ENDDO

C Rate of change in r-direction

      DO I = 1,Nr
        DO J = 1,Ntet
           DO K = 1, Nz
            IF (rc(I) .NE. 0.D0)                                                 &
     &      dC(I,J,K) =                                                          &
     &          -(Flux(I+1,J,K)*r(I+1) - Flux(I,J,K)*r(I))/dr(I)/(rc(I))
          ENDDO
        ENDDO
      ENDDO

C ## transport in teta -direction ##
C First diffusion internal cells

      DO I = 1, Nr
        DO J = 2, Ntet
          DO K = 1, Nz
            Flux(I,J,K) = -D_tet(J)*(C(I,J,K)-C(I,J-1,K)) /dtetaux(J)            &
     &                    /rc(I) 
          ENDDO
        ENDDO 
      ENDDO

      Flux(:,1,:)      = 0.D0
      Flux(:,Ntet+1,:) = 0.D0
       
C upstream teta-boundary

      IF (BcTetup .EQ. 1) THEN
        DO I = 1, Nr
          DO K = 1, Nz
            Flux(I,1,K) = Ftetup
          ENDDO
        ENDDO

      ELSE IF (BcTetup .EQ. 2) THEN
        DO I = 1, Nr
          DO K = 1, Nz
           Flux(I,1,K) = -D_tet(1)*(C(I,1,K)-Ctetup)/dtetaux(1)/rc(I)
          ENDDO  
        ENDDO

      ELSE IF (BcTetup .EQ. 3) THEN
        Flux(:,1,:) = 0.D0

      ELSE IF (BcTetup .EQ. 5) THEN
        DO I = 1, Nr
          DO K = 1, Nz

          Cbound = (C(I,1,K)   *D_tet(1)     *dtetaux(Ntet+1) +                  &
     &              C(I,Ntet,K)*D_tet(Ntet+1)*dtetaux(1))/                       &
     &              (D_tet(Ntet+1)*dtetaux(1)+D_tet(1)*dtetaux(Ntet+1))

          Flux(I,1,K)      = -D_tet(1)     *(C(I,1,K)-Cbound)                    &
     &                                       /dtetaux(1)/rc(I) 
          Flux(I,Ntet+1,K) = -D_tet(Ntet+1)*(Cbound-C(I,Ntet,K))                 & 
     &                                       /dtetaux(Ntet+1)/rc(I) 
          ENDDO
        ENDDO

      ENDIF

      DO I = 1, Nr
        DO K = 1, Nz
          Fluxtetup(I,K) = Flux(I,1,K)
        ENDDO
      ENDDO

C downstream teta-boundary

      IF (BcTetdown .EQ. 1) THEN
        DO I = 1, Nr
          DO K = 1, Nz
            Flux(I,Ntet+1,K) = Ftetdown
          ENDDO
        ENDDO

      ELSE IF (BcTetdown .EQ. 2) THEN
        DO I = 1, Nr
          DO K = 1, Nz
            Flux(I,Ntet+1,K) = -D_tet(Ntet+1) * (Ctetdown-C(I,Ntet,K))           & 
     &                                           /dtetaux(Ntet+1)/rc(I) 
          ENDDO 
        ENDDO

      ELSE IF (BcTetdown .EQ. 3) THEN
        Flux(:,Ntet+1,:) = 0.D0

      ENDIF

      DO I = 1, Nr
        DO K = 1, Nz
          Fluxtetdown(I,K) = Flux(I,Ntet+1,K)
        ENDDO
      ENDDO

C rate of change, negative flux gradient

      DO I = 1,Nr
        DO J = 1,Ntet
          DO K = 1, Nz
            IF(rc(I) .NE. 0.D0 )                                                 &
     &      dC(I,J,K) = dC(I,J,K)                                                &
     &        -(Flux(I,J+1,K)-Flux(I,J,K)) /dtet(J)/rc(I)
          ENDDO 
        ENDDO
      ENDDO
    

C ## transport in z-direction ##
C First diffusion internal cells

      DO I = 1, Nr
        DO J = 1, Ntet
          DO K = 2, Nz
            Flux(I,J,K) = -D_z(J)*(C(I,J,K)-C(I,J,K-1)) /dzaux(K) 
          ENDDO
        ENDDO 
      ENDDO

      Flux(:,:,1)      = 0.D0
      Flux(:,:,Nz+1) = 0.D0
       
C upstream z-boundary

      IF (BCzUp .EQ. 1) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            Flux(I,J,1) = Fzup
          ENDDO
        ENDDO

      ELSE IF (BCzUp .EQ. 2) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
           Flux(I,J,1) = -D_z(1)*(C(I,J,1)-Czup)  /dzaux(1) 
          ENDDO  
        ENDDO

      ELSE IF (BCzUp .EQ. 3) THEN
        Flux(:,:,1) = 0.D0

      ELSE IF (BCzUp .EQ. 5) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
          Cbound = (C(I,J,1 )*D_z(1)   *dzaux(Nz+1) +                            &
     &              C(I,J,Nz)*D_z(Nz+1)*dzaux(1)) /                              &
     &              (D_z(Nz+1)*dzaux(1)+D_z(1)*dzaux(Nz+1))
          Flux(I,J,1)    = -D_z(1)   *(C(I,J,1)-Cbound)  /dzaux(1)    
          Flux(I,J,Nz+1) = -D_z(Nz+1)*(Cbound-C(I,J,Nz)) /dzaux(Nz+1)              
          ENDDO                                      
        ENDDO

      ENDIF

      DO I = 1, Nr
        DO J = 1, Ntet
          FluxZup(I,J) = Flux(I,J,1)
        ENDDO
      ENDDO

C downstream z-boundary

      IF (BCzDown .EQ. 1) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            Flux(I,J,Nz+1) = Fzdown
          ENDDO
        ENDDO

      ELSE IF (BCzDown .EQ. 2) THEN
        DO I = 1, Nr
          DO J = 1, Ntet
            Flux(I,J,Nz+1) = -D_Z(Nz+1) * (CZdown-C(I,J,Nz))                     & 
     &                             /dZaux(NZ+1)               
          ENDDO 
        ENDDO

      ELSE IF (BCzDown .EQ. 3) THEN
        Flux(:,:,Nz+1) = 0.D0

      ENDIF

      DO I = 1, Nr
        DO J = 1, Ntet
          FluxZdown(I,J) = Flux(I,J,Nz+1)
        ENDDO
      ENDDO

C rate of change, negative flux gradient

      DO I = 1, Nr
        DO J = 1, Ntet
          DO K = 1, Nz

            dC(I,J,K) = dC(I,J,K)                                                &
     &        -(Flux(I,J,K+1)-Flux(I,J,K))/dz(K) 
          ENDDO 
        ENDDO
      ENDDO

      RETURN
      END

