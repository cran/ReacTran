C#==============================================================================
C transport in a 2-dimensional finite volume grid 
C compile        system("R CMD SHLIB tran.volume.2D.f")
C load           dyn.load("tran.volume.2D.dll")
C#==============================================================================

      SUBROUTINE tran2dvol (Nx, Ny, C,                                          &
     &                 V_xup, V_xdown, V_yup, V_ydown, V_z,                     &
     &                 BcxUp, BcxDown, BcyUp, BcyDown, Bcz,                     &
C NOT      &                 D_x, D_y, D_z,                                           &
     &                 Flow_x, Flow_y, Volume, dC,                              &
     &                 Fxup, Fxdown, Fyup, Fydown, Fz) 
       IMPLICIT NONE
       INTEGER Nx,Ny                  ! dim of C

C input concentration
       DOUBLE PRECISION C(Nx, Ny), dC(Nx, Ny)

C Boundary concentrations (used if Bc = 2) or fluxes (used if Bc = 1) 
       DOUBLE PRECISION V_xup(Ny), V_xdown(Ny)
       DOUBLE PRECISION V_yup(Nx), V_ydown(Nx, Nz)
       DOUBLE PRECISION V_z(Nx, Ny)  

       DOUBLE PRECISION Volume(Nx, Ny)

C Flows at upstream interface, lateral and lateral concentrations
       DOUBLE PRECISION Flow_x(Nx+1, Ny)
       DOUBLE PRECISION Flow_y(Nx, Ny+1)
C       DOUBLE PRECISION D_x(Nx+1, Ny), D_Y(Nx, Ny+1), D_Z(Nx, Ny)
       DOUBLE PRECISION Flow_z, Flux

C boundary concitions (1= flux, 2=conc, 3 = 0-gradient)
      INTEGER BcxUp, BcxDown, BcyUp, BcyDown, Bcz

C locals 
       INTEGER I, J, K
       CHARACTER (LEN = 100) msg
C -------------------------------------------------------------------------------
       dC = 0.d0

C transport in x-direction
       DO I = 1, Nx+1
         DO J = 1, Ny
	           IF (Flow_x(I, J) > 0) THEN
	             IF (I > 1) THEN 
                 Flux  = Flow_x(I, J)  * C(I-1, J)
               ENDIF
             ELSE   
               IF (I < Nx +1) THEN
                 Flux  = Flow_x(I, J) * C(I, J)
               ENDIF
	           ENDIF 

             IF ( I .EQ. 1) THEN   
               IF(BcxUp .EQ. 3) THEN 
                 Flux  = Flow_x(I, J)  * C(I, J)
               ELSE IF(BcxUp .EQ. 1) THEN 
                 Flux  = V_xup(J)
               ELSE
                 Flux =  Flow_x(I, J) * V_xup(J)  
               ENDIF       
             ENDIF  

             IF (I .EQ. Nx + 1) THEN   
               IF(BcxDown .EQ. 3) THEN 
                 Flux  = Flow_x(I, J) * C(I-1, J)
               ELSE IF(BcxDown .EQ. 1) THEN 
                 Flux  = V_xdown(J)
               ELSE
                 Flux  = Flow_x(I, J) * V_xdown(J)
               ENDIF
	           ENDIF 

             IF (I < Nx + 1) THEN
               dC(I, J)   = dC(I, J)   + Flux/Volume(I, J)
             ENDIF
             IF (I > 1) THEN
               dC(I-1, J) = dC(I-1, J) - Flux/volume(I-1, J)
             ENDIF  	    
         ENDDO
       ENDDO

C transport in y-direction
       DO I = 1, Nx
         DO J = 1, Ny+1
             IF (Flow_y(I, J) > 0) THEN
               IF (J > 1) THEN
                 Flux  = Flow_y(I, J)  * C(I, J-1)
               ENDIF
             ELSE   
               IF (J < Ny + 1) THEN 
                 Flux  = Flow_y(I, J) * C(I, J)
               ENDIF 
	           ENDIF 

             IF (J .EQ. 1) THEN
               IF(BcyUp .EQ. 3) THEN 
                 Flux  = Flow_y(I, J)  * C(I, J)
               ELSE IF(BcyUp .EQ. 1) THEN 
                 Flux  = V_yup(I)
               ELSE
                 Flux =  Flow_y(I, J) * V_yup(I)  
               ENDIF 
	           ENDIF
	           IF (J .EQ. Ny+1) THEN
               IF(BcyDown .EQ. 3) THEN   
                 Flux  = Flow_y(I, J) * C(I, J-1)
               ELSE IF(BcyDown .EQ. 1) THEN 
                 Flux  = V_ydown(I)
               ELSE
                 Flux =  Flow_y(I, J) * V_ydown(I)  
               ENDIF
	           ENDIF
             IF (J < Ny + 1) THEN 
               dC(I, J)   = dC(I, J) + Flux/Volume(I, J)
             ENDIF 
             IF (J > 1) THEN
               dC(I, J-1) = dC(I, J-1) -Flux/volume(I, J-1)	    
             ENDIF
         ENDDO
       ENDDO

C transport in z-direction       !!! flux_z upstream = given
       DO I = 1, Nx
         DO J = 1, Ny
            Flowz =  -(Flow_x(I+1,J) - Flow_x(I,J) 
     &                   +Flow_y(I,J+1) - Flow_y(I,J))              

               IF(Bcz .EQ. 1) THEN 
                 Flux  = V_z(I, J)
               ELSE IF(Bcz .EQ. 3) THEN 
                 Flux  = Flowz  * C(I, J)
               ELSE
                 Flux =  Flowz * V_zup(I, J)  
               ENDIF  

               dC(I, J) = dC(I, J) + Flux/Volume(I, J)
         ENDDO
       ENDDO

       RETURN
       END SUBROUTINE tran2Dvol



C#==============================================================================
C Horizontal fluxes in a 2-dimensional finite volume grid
C#==============================================================================

      SUBROUTINE flux2DHor (Nx, Ny, C, Flow_x, Flow_y,                          &
     &                        Flux_x, Flux_y) 
       IMPLICIT NONE
       INTEGER Nx, Ny                   ! length of C

C input concentration
       DOUBLE PRECISION C(Nx, Ny) 
       DOUBLE PRECISION Volume(Nx, Ny)

C Flows at upstream interface, lateral and lateral concentrations
       DOUBLE PRECISION Flow_x(Nx+1, Ny), Flow_y(Nx, Ny+1)
       DOUBLE PRECISION Flow_z, Flux
       DOUBLE PRECISION Flux_x(Nx+1, Ny), Flux_y(Nx, Ny+1)

C locals 
       INTEGER I, J, K
       CHARACTER (LEN = 100) msg
C -------------------------------------------------------------------------------

C transport in x-direction
       DO I = 1, Nx+1
         DO J = 1, Ny
	           IF (Flow_x(I, J, K) > 0) THEN
	           IF (Flow_x(I, J) > 0) THEN
	             IF (I > 1) THEN 
                 Flux  = Flow_x(I, J)  * C(I-1, J)
               ENDIF
             ELSE   
               IF (I < Nx +1) THEN
                 Flux  = Flow_x(I, J) * C(I, J)
               ENDIF
	           ENDIF 

             IF ( I .EQ. 1) THEN   
               IF(BcxUp .EQ. 3) THEN 
                 Flux  = Flow_x(I, J)  * C(I, J)
               ELSE IF(BcxUp .EQ. 1) THEN 
                 Flux  = V_xup(J)
               ELSE
                 Flux =  Flow_x(I, J) * V_xup(J)  
               ENDIF       
             ENDIF  
            Flux_x(i,j) = Flux
         ENDDO
       ENDDO

C transport in y-direction
       DO I = 1, Nx
         DO J = 1, Ny+1
             IF (Flow_y(I, J) > 0) THEN
               IF (J > 1) THEN
                 Flux  = Flow_y(I, J)  * C(I, J-1)
               ENDIF
             ELSE   
               IF (J < Ny + 1) THEN 
                 Flux  = Flow_y(I, J) * C(I, J)
               ENDIF 
	           ENDIF 

             IF (J .EQ. 1) THEN
               IF(BcyUp .EQ. 3) THEN 
                 Flux  = Flow_y(I, J)  * C(I, J)
               ELSE IF(BcyUp .EQ. 1) THEN 
                 Flux  = V_yup(I)
               ELSE
                 Flux =  Flow_y(I, J) * V_yup(I)  
               ENDIF 
	           ENDIF
	           IF (J .EQ. Ny+1) THEN
               IF(BcyDown .EQ. 3) THEN   
                 Flux  = Flow_y(I, J) * C(I, J-1)
               ELSE IF(BcyDown .EQ. 1) THEN 
                 Flux  = V_ydown(I)
               ELSE
                 Flux =  Flow_y(I, J) * V_ydown(I)  
               ENDIF
	           ENDIF
           
            Flux_y(i,j) = Flux
         ENDDO
       ENDDO

       RETURN
       END SUBROUTINE fluxHor




C#==============================================================================
C Internal Fluxes across an internal boundary in 3-D domain
C#==============================================================================

C      SUBROUTINE Flux2DBound (Nx, Ny, C, Flow_x, Flow_y, NBnd, Bnd,           &
C     &                      Flux_x, Flux_y, Flux) 
C       END SUBROUTINE FluxBound
