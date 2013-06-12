C#==============================================================================
C transport in a 3-dimensional finite volume grid 
C compile        system("R CMD SHLIB tran.volume.3D.f")
C load           dyn.load("tran.volume.3D.dll")
C#==============================================================================

      SUBROUTINE tran3dvol (Nx, Ny, Nz, C ,                                     &
     &                 V_xup, V_xdown, V_yup, V_ydown, V_zup, V_zdown,          &
     &                 BcxUp, BcxDown, BcyUp, BcyDown, BczUp, BczDown,          &
C NOT      &                 D_x, D_y, D_z,                                           &
     &                 Flow_x, Flow_y, Flow_z, Volume, dC,                      &
     &                 Fxup, Fxdown, Fyup, Fydown, Fzup, Fzdown) 
       IMPLICIT NONE
       INTEGER Nx,Ny,Nz                  ! dim of C

C input concentration
       DOUBLE PRECISION C(Nx, Ny, Nz), dC(Nx, Ny, Nz)

C Boundary concentrations (used if Bc = 2) or fluxes (used if Bc = 1) 
       DOUBLE PRECISION V_xup(Ny, Nz), V_xdown(Ny, Nz)
       DOUBLE PRECISION V_yup(Nx, Nz), V_ydown(Nx, Nz)
       DOUBLE PRECISION V_zup(Nx, Ny). V_zdown(Nx, Ny) 

       DOUBLE PRECISION Volume(Nx, Ny, Nz)

C Flows at upstream interface, lateral and lateral concentrations
       DOUBLE PRECISION Flow_x(Nx+1, Ny, Nz)
       DOUBLE PRECISION Flow_y(Nx, Ny+1, Nz), Flow_Z(Nx, Ny) !Z just bottom flux
C       DOUBLE PRECISION D_x(Nx+1, Ny, Nz), D_Y(Nx, Ny+1, Nz), D_Z(Nx, Ny, Nz+1)
       DOUBLE PRECISION Flow_z, Flow_z_prev, Flux

C boundary concitions (1= flux, 2=conc, 3 = 0-gradient)
      INTEGER BcxUp, BcxDown, BcyUp, BcyDown, BczUp, BczDown

C locals 
       INTEGER I, J, K
       CHARACTER (LEN = 100) msg
C -------------------------------------------------------------------------------
       dC = 0.d0

C transport in x-direction
       DO I = 1, Nx+1
         DO J = 1, Ny
           DO K = 1, Nz

	           IF (Flow_x(I, J, K) > 0) THEN
	             IF (I > 1) THEN 
                 Flux  = Flow_x(I, J, K)  * C(I-1, J, K)
               ENDIF
             ELSE   
               IF (I < Nx +1) THEN
                 Flux  = Flow_x(I, J, K) * C(I, J, K)
               ENDIF
	           ENDIF 

             IF ( I .EQ. 1) THEN   
               IF(BcxUp .EQ. 3) THEN 
                 Flux  = Flow_x(I, J, K)  * C(I, J, K)
               ELSE IF(BcxUp .EQ. 1) THEN 
                 Flux  = V_xup(J, K)
               ELSE
                 Flux =  Flow_x(I, J, K) * V_xup(J, K)  
               ENDIF       
             ENDIF  

             IF (I .EQ. Nx + 1) THEN   
               IF(BcxDown .EQ. 3) THEN 
                 Flux  = Flow_x(I, J, K) * C(I-1, J, K)
               ELSE IF(BcxDown .EQ. 1) THEN 
                 Flux  = V_xdown(J, K)
               ELSE
                 Flux  = Flow_x(I, J, K) * V_xdown(J, K)
               ENDIF
	           ENDIF 

             IF (I < Nx + 1) THEN
               dC(I, J, K)   = dC(I, J, K)   + Flux/Volume(I, J, K)
             ENDIF
             IF (I > 1) THEN
               dC(I-1, J, K) = dC(I-1, J, K) - Flux/volume(I-1, J, K)
             ENDIF  	    
           
           ENDDO
         ENDDO
       ENDDO

C transport in y-direction
       DO I = 1, Nx
         DO J = 1, Ny+1
           DO K = 1, Nz
             IF (Flow_y(I, J, K) > 0) THEN
               IF (J > 1) THEN
                 Flux  = Flow_y(I, J, K)  * C(I, J-1, K)
               ENDIF
             ELSE   
               IF (J < Ny + 1) THEN 
                 Flux  = Flow_y(I, J, K) * C(I, J, K)
               ENDIF 
	           ENDIF 

             IF (J .EQ. 1) THEN
               IF(BcyUp .EQ. 3) THEN 
                 Flux  = Flow_y(I, J, K)  * C(I, J, K)
               ELSE IF(BcyUp .EQ. 1) THEN 
                 Flux  = V_yup(I, K)
               ELSE
                 Flux =  Flow_y(I, J, K) * V_yup(I, K)  
               ENDIF 
	           ENDIF
	           IF (J .EQ. Ny+1) THEN
               IF(BcyDown .EQ. 3) THEN   
                 Flux  = Flow_y(I, J, K) * C(I, J-1, K)
               ELSE IF(BcyDown .EQ. 1) THEN 
                 Flux  = V_ydown(I, K)
               ELSE
                 Flux =  Flow_y(I, J, K) * V_ydown(I, K)  
               ENDIF
	           ENDIF
             IF (J < Ny + 1) THEN 
               dC(I, J, K)   = dC(I, J, K) + Flux/Volume(I, J, K)
             ENDIF 
             IF (J > 1) THEN
               dC(I, J-1, K) = dC(I, J-1, K) -Flux/volume(I, J-1, K)	    
             ENDIF
           
           ENDDO
         ENDDO
       ENDDO

C transport in z-direction       !!! flux_z upstream = given
       DO I = 1, Nx
         DO J = 1, Ny
           Flow_z_prev = Flow_z(I,J)
           DO K = 1, Nz+1
             IF (K .EQ. 1) 
               Flowz = Flow_z(I,J)
             ELSE  
               Flowz = Flow_z_prev -(Flow_x(I+1,J,K-1) - Flow_x(I,J,K-1) 
     &                             +Flow_y(I,J+1,K-1) - Flow_y(I,J,K-1))              
             ENDIF
             Flow_z_prev = Flowz

	           IF (Flowz > 0) THEN
	             IF(K > 1) THEN 
                 Flux  = Flowz  * C(I, J, K-1)
               ENDIF       
             ELSE   
               IF (K < Nz + 1) THEN
                 Flux = Flowz * C(I, J, K)
               ENDIF   
	           ENDIF 
             IF (K .EQ. 1)
               IF(BczUp .EQ. 1) THEN 
                 Flux  = V_zup(I, J)
               ELSE IF(BczUp .EQ. 3) THEN 
                 Flux  = Flowz  * C(I, J, K)
               ELSE
                 Flux =  Flowz * V_zup(I, J)  
               ENDIF  
             ENDIF
             IF (K == Nz + 1) THEN
               IF(BczDown .EQ. 3) THEN 
                 Flux = Flowz * C(I, J, K-1)
               ELSE IF(BczDown .EQ. 1) THEN 
                 Flux  = V_zdown(I, J)
               ELSE
                 Flux =  Flowz * V_zdown(I, J)  
               ENDIF
             ENDIF
             
             IF (K < Nz + 1) THEN
               dC(I, J, K) = dC(I, J, K) + Flux/Volume(I, J, K)
             ENDIF
             dC(I, J, K-1) = dC(I, J, K-1) -Flux/volume(I, J, K-1)	    
            
           ENDDO

         ENDDO
       ENDDO

       RETURN
       END SUBROUTINE tran3Dvol



C#==============================================================================
C Horizontal fluxes in a 3-dimensional finite volume grid
C#==============================================================================

      SUBROUTINE flux3DHor (Nx, Ny, Nz, C, Flow_x, Flow_y,                    &
     &                        Flux_x, Flux_y) 
       IMPLICIT NONE
       INTEGER Nx,Ny,Nz                  ! length of C

C input concentration
       DOUBLE PRECISION C(Nx, Ny, Nz) 
       DOUBLE PRECISION Volume(Nx, Ny, Nz)

C Flows at upstream interface, lateral and lateral concentrations
       DOUBLE PRECISION Flow_x(Nx+1, Ny, Nz), Flow_y(Nx, Ny+1, Nz)
       DOUBLE PRECISION Flow_z, Flow_z_prev, Flux
       DOUBLE PRECISION Flux_x(Nx+1, Ny), Flux_y(Nx, Ny+1)

C locals 
       INTEGER I, J, K
       CHARACTER (LEN = 100) msg
C -------------------------------------------------------------------------------

C transport in x-direction
       DO I = 1, Nx+1
         DO J = 1, Ny
           Flux_x(I,J) = 0.d0
           DO K = 1, Nz
	           IF (Flow_x(I, J, K) > 0) THEN
	             IF (I > 1) THEN 
                 Flux  = Flow_x(I, J, K)  * C(I-1, J, K)
               ENDIF
             ELSE   
               IF (I < Nx +1) THEN
                 Flux  = Flow_x(I, J, K) * C(I, J, K)
               ENDIF
	           ENDIF 

             IF ( I .EQ. 1) THEN   
               IF(BcxUp .EQ. 3) THEN 
                 Flux  = Flow_x(I, J, K)  * C(I, J, K)
               ELSE IF(BcxUp .EQ. 1) THEN 
                 Flux  = V_xup(J, K)
               ELSE
                 Flux =  Flow_x(I, J, K) * V_xup(J, K)  
               ENDIF       
             ENDIF  

             IF (I .EQ. Nx + 1) THEN   
               IF(BcxDown .EQ. 3) THEN 
                 Flux  = Flow_x(I, J, K) * C(I-1, J, K)
               ELSE IF(BcxDown .EQ. 1) THEN 
                 Flux  = V_xdown(J, K)
               ELSE
                 Flux  = Flow_x(I, J, K) * V_xdown(J, K)
               ENDIF
	           ENDIF 

             Flux_x(I, J) = Flux_x(I, J) + Flux
           
           ENDDO
         ENDDO
       ENDDO

C transport in y-direction
       DO I = 1, Nx
         DO J = 1, Ny+1
           Flux_y(I,J) = 0.D0
           DO K = 1, Nz
             IF (Flow_y(I, J, K) > 0) THEN
               IF (J > 1) THEN
                 Flux  = Flow_y(I, J, K)  * C(I, J-1, K)
               ENDIF
             ELSE   
               IF (J < Ny + 1) THEN 
                 Flux  = Flow_y(I, J, K) * C(I, J, K)
               ENDIF 
	           ENDIF 

             IF (J .EQ. 1) THEN
               IF(BcyUp .EQ. 3) THEN 
                 Flux  = Flow_y(I, J, K)  * C(I, J, K)
               ELSE IF(BcyUp .EQ. 1) THEN 
                 Flux  = V_yup(I, K)
               ELSE
                 Flux =  Flow_y(I, J, K) * V_yup(I, K)  
               ENDIF 
	           ENDIF
	           IF (J .EQ. Ny+1) THEN
               IF(BcyDown .EQ. 3) THEN   
                 Flux  = Flow_y(I, J, K) * C(I, J-1, K)
               ELSE IF(BcyDown .EQ. 1) THEN 
                 Flux  = V_ydown(I, K)
               ELSE
                 Flux =  Flow_y(I, J, K) * V_ydown(I, K)  
               ENDIF
	           ENDIF
	          
             Flux_y(I, J) = Flux_y(I, J) + Flux

           
           ENDDO
         ENDDO
       ENDDO

       RETURN
       END SUBROUTINE flux3DHor




C#==============================================================================
C Internal Fluxes across an internal boundary in 3-D domain
C#==============================================================================

      SUBROUTINE Flux3DBound (Nx, Ny, Nz, C, Flow_x, Flow_y, NBnd, Bnd,         &
     &                      Flux_x, Flux_y, Flux) 
       IMPLICIT NONE
       INTEGER Nx, Ny, Nz, NBnd                  ! dimension of C
       INTEGER Bnd(NBnd, 2)                      ! Position of boundary

C input concentration
       DOUBLE PRECISION C(Nx, Ny, Nz)

C Flows at upstream interface
       DOUBLE PRECISION Flow_x(Nx+1, Ny, Nz), Flow_y(Nx, Ny+1, Nz)

C RETURN ELEMENTS
       DOUBLE PRECISION Flux_x, Flux_y, Flux

C locals 
       INTEGER I, J, K, N, prevI, prevJ
       CHARACTER (LEN = 100) msg
C -------------------------------------------------------------------------------
C transport in x-direction
       prevI = 0
       prevJ = 0
       Flux_x = 0.d0 
       Flux_y = 0.d0

C Fluxes in North-South (y-direction)       
 
       DO N = 1, Nbnd
         I = BND(N, 1)
         J = BND(N, 2)
        
           IF (J .EQ. prevJ) CYCLE
           prevJ = J        
           
           DO K = 1, Nz

	           IF (Flow_x(I, J, K) > 0) THEN
	             IF (I > 1) THEN 
                 Flux_x  = Flux_x + Flow_x(I, J, K)  * C(I-1, J, K)
               ENDIF       
             ELSE   
               IF (I < Nx +1) THEN
                 Flux_x  = Flux_x + Flow_x(I, J, K) * C(I, J, K)
               ENDIF
	           ENDIF 
           ENDDO
       ENDDO

C transport in x-direction
       DO N = 1, Nbnd
         I = BND(N, 1)
         J = BND(N, 2)
        
           IF (I .EQ. prevI) CYCLE
           prevI = I        
           DO K = 1, Nz
	          
             IF (Flow_y(I, J, K) > 0) THEN
               IF (J > 1) THEN
                 Flux_y  = Flux_y + Flow_y(I, J, K)  * C(I, J-1, K)
               ENDIF 
             ELSE   
               IF (J < Ny + 1) THEN 
                 Flux_y  = Flux_y + Flow_y(I, J, K) * C(I, J, K)
               ENDIF
	           ENDIF 
           ENDDO
       ENDDO
       Flux = +Flux_x +Flux_y
       RETURN
       END SUBROUTINE Flux3DBound
