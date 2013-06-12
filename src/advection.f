c karline: changes to take into account zero-grad and conc boundary conditions

      subroutine advection (N, Y, dt, h, hint, v, Bcup, Bcdown,              &
     &       Yup, Ydown, VFint, VF, Aint, A, method, mode, split,            &
     &       dY, cu, it)

c-----------------------------------------------------------------------------------------
c based on the advection scheme in the GOTM model, code from 2006-11-06
c fluxes are defined on the interfaces, in an upstream-biased way.
c slope-delimeters are applied to obtain monotonic and positive schemes
c also in the presence of large gradients. 
c there are 5 different slope delimeters; first-order upstream, 
c 3rd order upstream-biased polynomial, 3rd order with superbee delimiter, 
c 3rd order with MUSCL limiter
c as described in Pietrzak 98

c Karline: made changes to make it work for negative ww...
c          added volume fraction, surface area; these two properties will generally = 1.
cc-----------------------------------------------------------------------------------------

      IMPLICIT NONE
c  number of vertical layers, time step
      INTEGER                  :: N
      DOUBLE PRECISION         :: dt

c  layer thickness (m), distance from mid to mid of each layer
      DOUBLE PRECISION         :: h(N), hint(0:N)

c  advection speed in the direction of the axis
      DOUBLE PRECISION         :: v(0:N), ww(0:N)  

c  volume fraction and surface at interface and in middle of layers
      DOUBLE PRECISION         :: VFint(0:N), Aint(0:N), VF(N), A(N)  

c  type of upper and lower Boundary Condition (only 1 and 2 used in R)
      INTEGER                  :: Bcdown, Bcup
      integer, parameter       :: Flux     = 1
      integer, parameter       :: Value    = 2
      integer, parameter       :: ZeroGrad = 3
      integer, parameter       :: zeroDivergence = 4

c  value of upper and lower bnd conc
      DOUBLE PRECISION         :: Ydown, Yup

c  type of advection scheme, slope delimeters 
      INTEGER                  :: method
      integer, parameter       :: UPSTREAM =5  
      integer, parameter       :: P2       =4   
      integer, parameter       :: P2_PDM   =3
      integer, parameter       :: Superbee =2
      integer, parameter       :: MUSCL    =1

c  advection mode: 0= non-conservative (e.g.water flow), 1= conservative, e.g. sinking
      INTEGER                  :: mode

c  splitting mode: >0 split, <= 0: do not split, one step and return)
      INTEGER                  :: split

c  concentration to be transported
      DOUBLE PRECISION         :: Y(N)

c  rate of change due to advection
      DOUBLE PRECISION         :: dY(N)

      DOUBLE PRECISION         :: one6th=1.0d0/6.0d0
c maximal number of iterations
      INTEGER, parameter       :: itmax=10000

c  LOCAL VARIABLES:
      integer                  :: i,k,it, istart, istop
      DOUBLE PRECISION         :: x,r,Phi,limit=0.d0
      DOUBLE PRECISION         :: Yu,Yc,Yd
      DOUBLE PRECISION         :: c,cmax
      DOUBLE PRECISION         :: cu(0:N)
c-----------------------------------------------------------------------------------------
c  copy of current value of state variables 
      do k =1,N
         dy(k) = y(k)
      enddo

      istart = 1
      istop = N-1
      IF (Bcdown == Value) THEN
        istop = N
      ELSEIF (Bcdown == ZeroGrad) THEN
        istop = N
        Ydown = y(N)
      ENDIF
      IF (Bcup == Value) THEN
        istart = 0
      ELSEIF (Bcup == ZeroGrad) THEN
        istart = 0
        Yup = y(1)
      ENDIF

c  initialize maximum Courant number
      cmax = 0.d0

c convert to per fraction (in case VFint != 1)
       do k=0,N
         ww(k) = v(k)*VFint(k)  
       enddo

      IF (split .GT. 0) THEN

c  compute maximum Courant number; estimate nr of iterations 
       do k=0,N
         c = dabs(ww(k))*dt/hint(k)
         if (c.gt.cmax) cmax=c
       enddo

       if (cmax . GT. 1) then 
        it = min(itmax, int(cmax)+1)   
        it = max(1,it)
       else 
        it = 1
       endif   
     
      ELSE
        it = 1
      ENDIF
            
c  (time) splitting loop
      do i=1,it

c      initialize upstream interface fluxes with zero
        cu   = 0.d0


c     spatial loop - karline : changed into 1:N-1
        do k = istart, istop

c        positive speed 
           if (ww(k) .gt. 0.d0) then
              if(k > 0) THEN
                c=ww(k)/dble(it)*dt/hint(k-1)          ! courant number
              else
                c=ww(k)/dble(it)*dt/hint(k)            ! courant number
              endif
              
              if (k .gt. 1) then
                 Yu=Y(k-1)                             ! upstream value
              else
                 Yu=Yup
              end if
              
              if (k .gt. 0) then
                 Yc = Y(k)                             ! central value
              else
                 Yc = Yup
              end if
              
              if (k .lt. N) then
                Yd = Y(k+1)                           ! downstream value
              else
                Yd = Ydown
              endif
                
c        negative speed
           else
             
             c=-ww(k)/dble(it)*dt/hint(k)              ! courant number

              if (k .gt. 0) then
                Yd = Y(k)                                 ! downstream value
              else
                 Yd = Yup
              end if
                  
              if (k .lt. N) then
                Yc = Y(k+1)                               ! central value
              else
                Yc = Ydown
              endif
                    
             if (k .lt. N-1) then
               Yu = Y(k+2)                             ! upstream value
             else
               Yu = Yup                               
             endif  

           end if
                 
           if (abs(Yd - Yc) .gt. 1e-10) then           ! slope ratio
                r=(Yc - Yu)/(Yd - Yc)
           else
                r=(Yc - Yu)*1.e10
           end if

c        limit the flux according to different suggestions
           select case (method)
             case (UPSTREAM)
               limit = 0.d0
             case ((P2),(P2_PDM))
c        the flux-factor phi
               x    =  one6th*(1.-2.0*c)
               Phi  =  (0.5+x)+(0.5-x)*r
               if (method.eq.P2) then
                  limit=Phi
               else
                  limit=max(0.d0,min(Phi,2./(1.d0-c),2.*r/(c+1.e-10)))
               end if
             case (Superbee)
               limit=max(0.d0, min(1.d0, 2.0*r), min(r,2.*1.d0) )
             case (MUSCL)
               limit=max(0.d0,min(2.*1.d0,2.0*r,0.5*(1.d0+r)))
             case default
c               call rerror( 'unkown advection method')  ! should not happen
           end select

c        compute the limited flux  
           cu(k)  = ww(k) *(Yc+0.5d0*limit*(1.-c)*(Yd-Yc))

        end do

c       downstream boundary conditions
        select case (Bcdown)
          case (flux)
            cu(N) = Ydown                ! flux OUT of the domain is positive
          case (zeroDivergence)
            cu(N) = cu(N-1)
c          case default
c            call rwarn('unkown downstream boundary condition type')
        end select


c       upstream boundary conditions
        select case (Bcup)
          case (flux)
            cu(0) =   Yup                 ! flux into the domain is positive
          case (zeroDivergence)
           cu(0) = cu(1)
c          case default
c             call rwarn('unkown upstream boundary condition type')
        end select

c     the advection step

        if (mode.eq.0) then ! conservative - KARLINE CHECK...
           do k=1,N
             Y(k)=Y(k)-1.d0/dble(it)*dt*((Aint(k)*cu(k)-                       &
     &                       Aint(k-1)*cu(k-1))/ h(k)/A(k)/VF(k)          &
     &          -Y(k)*(ww(k)-ww(k-1))/h(k)/A(k)/VF(k))
           enddo
        else                ! non-conservative
           do k=1,N
             Y(k)=Y(k)-1.d0/dble(it)*dt*((Aint(k)*cu(k)-                       &
     &                       Aint(k-1)*cu(k-1))/h(k)/A(k)/VF(k))
           enddo
        end if

      end do ! end of the iteration loop

c     rate of change due to advection 
      do k =1,N
         dy(k) = (y(k)-dy(k))/dt
      enddo

c Still to do: integrate fluxes in time (now cu = cu of last step )
c            flux = 0. ! at start
c            do k=1,N
c             flux(k)=flux(k)+1.d0/dble(it)*dt*cu(k)
 
      return
      end subroutine advection

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
! ... extensively modified by Karline Soetaert
!-----------------------------------------------------------------------






c karline: original routine, differs with respect to boundary conditions 


      subroutine advection_ori(N, Y, dt, h, hint, v, Bcup, Bcdown,              &
     &       Yup, Ydown, VFint, VF, Aint, A, method, mode, split,           &
     &       dY, cu, it)

c-----------------------------------------------------------------------------------------
c based on the advection scheme in the GOTM model, code from 2006-11-06
c fluxes are defined on the interfaces, in an upstream-biased way.
c slope-delimeters are applied to obtain monotonic and positive schemes
c also in the presence of large gradients. 
c there are 5 different slope delimeters; first-order upstream, 
c 3rd order upstream-biased polynomial, 3rd order with superbee delimiter, 
c 3rd order with MUSCL limiter
c as described in Pietrzak 98

c Karline: made changes to make it work for negative ww...
c          added volume fraction, surface area; these two properties will generally = 1.
cc-----------------------------------------------------------------------------------------

      IMPLICIT NONE
c  number of vertical layers, time step
      INTEGER                  :: N
      DOUBLE PRECISION         :: dt

c  layer thickness (m), distance from mid to mid of each layer
      DOUBLE PRECISION         :: h(N), hint(0:N)

c  advection speed in the direction of the axis
      DOUBLE PRECISION         :: v(0:N), ww(0:N)  

c  volume fraction and surface at interface and in middle of layers
      DOUBLE PRECISION         :: VFint(0:N), Aint(0:N), VF(N), A(N)  

c  type of upper and lower Boundary Condition (only 1 and 2 used in R)
      INTEGER                  :: Bcdown, Bcup
      integer, parameter       :: Flux     = 1
      integer, parameter       :: Value    = 2
      integer, parameter       :: ZeroGrad = 3
      integer, parameter       :: zeroDivergence = 4

c  value of upper and lower bnd conc
      DOUBLE PRECISION         :: Ydown, Yup

c  type of advection scheme, slope delimeters 
      INTEGER                  :: method
      integer, parameter       :: UPSTREAM =5  
      integer, parameter       :: P2       =4   
      integer, parameter       :: P2_PDM   =3
      integer, parameter       :: Superbee =2
      integer, parameter       :: MUSCL    =1

c  advection mode: 0= non-conservative (e.g.water flow), 1= conservative, e.g. sinking
      INTEGER                  :: mode

c  splitting mode: >0 split, <= 0: do not split, one step and return)
      INTEGER                  :: split

c  concentration to be transported
      DOUBLE PRECISION         :: Y(N)

c  rate of change due to advection
      DOUBLE PRECISION         :: dY(N)

      DOUBLE PRECISION         :: one6th=1.0d0/6.0d0
c maximal number of iterations
      INTEGER, parameter       :: itmax=10000

c  LOCAL VARIABLES:
      integer                  :: i,k,it
      DOUBLE PRECISION         :: x,r,Phi,limit=0.d0
      DOUBLE PRECISION         :: Yu,Yc,Yd
      DOUBLE PRECISION         :: c,cmax
      DOUBLE PRECISION         :: cu(0:N)
c-----------------------------------------------------------------------------------------
c  copy of current value of state variables 
      do k =1,N
         dy(k) = y(k)
      enddo

c  initialize maximum Courant number
      cmax = 0.d0

c convert to per fraction (in case VFint != 1)
       do k=0,N
         ww(k) = v(k)*VFint(k)  
       enddo

      IF (split .GT. 0) THEN

c  compute maximum Courant number; estimate nr of iterations 
       do k=0,N
         c = dabs(ww(k))*dt/hint(k)
         if (c.gt.cmax) cmax=c
       enddo

       if (cmax . GT. 1) then 
        it = min(itmax,int(cmax)+1)   ! WAS: min(itmax,int(cmax)+1)
        it = max(1,it)
       else 
        it = 1
       endif   
     
      ELSE
        it = 1
      ENDIF
            
c  (time) splitting loop
      do i=1,it

c      initialize upstream interface fluxes with zero
        cu   = 0.d0


c     spatial loop - karline : changed into 1:N-1
        do k = 1, N-1

c        positive speed 
           if (ww(k) .gt. 0.d0) then

              c=ww(k)/dble(it)*dt/hint(k-1)            ! courant number

              if (k .gt. 1) then
                 Yu=Y(k-1)                             ! upstream value
              else
                 Yu=Y(k)
              end if
              
              Yc=Y(k  )                                ! central value
              
              Yd=Y(k+1)                              ! downstream value

c        negative speed
           else

             c=-ww(k)/dble(it)*dt/hint(k)              ! courant number

             Yd = Y(k)                                ! downstream value
                  
             Yc = Y(k+1)                                   ! central value

             if (k .lt. N-1) then
               Yu = Y(k+2) 
             else
               Yu = Y(N)                               
             endif  

           end if
                 
           if (abs(Yd-Yc) .gt. 1e-10) then          ! slope ratio
                r=(Yc-Yu)/(Yd-Yc)
           else
                r=(Yc-Yu)*1.e10
           end if

c        limit the flux according to different suggestions
           select case (method)
             case (UPSTREAM)
               limit=0.d0
             case ((P2),(P2_PDM))
c        the flux-factor phi
               x    =  one6th*(1.-2.0*c)
               Phi  =  (0.5+x)+(0.5-x)*r
               if (method.eq.P2) then
                  limit=Phi
               else
                  limit=max(0.d0,min(Phi,2./(1.d0-c),2.*r/(c+1.e-10)))
               end if
             case (Superbee)
               limit=max(0.d0, min(1.d0, 2.0*r), min(r,2.*1.d0) )
             case (MUSCL)
               limit=max(0.d0,min(2.*1.d0,2.0*r,0.5*(1.d0+r)))
             case default
c               call rerror( 'unkown advection method')  ! should not happen
           end select

c        compute the limited flux  
           cu(k)  = ww(k) *(Yc+0.5d0*limit*(1.-c)*(Yd-Yc))

        end do

c       downstream boundary conditions
        select case (Bcdown)
          case (flux)
            cu(N) = Ydown                ! flux OUT of the domain is positive
          case (value)
            if (ww(N).lt. 0.d0) then
              cu(N) =  ww(N)*Ydown
            else
              cu(N) =  ww(N)*Y(N)
            end if
          case (zerograd)
            cu(N) =  ww(N)*Y(N)
          case (zeroDivergence)
            cu(N) = cu(N-1)
          case default
            call rwarn('unkown downstream boundary condition type')
        end select


c       upstream boundary conditions
        select case (Bcup)
          case (flux)
            cu(0) =   Yup                 ! flux into the domain is positive
          case (value)
            if(ww(0) .gt. 0.d0) then       
              cu(0) =  ww(0)*Yup
            else
              cu(0) =  ww(0)*Y(1)
            end if
          case (ZeroGrad)
              cu(0) =  ww(0)*Y(1)
          case (zeroDivergence)
           cu(0) = cu(1)
          case default
             call rwarn('unkown upstream boundary condition type')
        end select

c     the advection step

        if (mode.eq.0) then ! conservative - KARLINE CHECK...
           do k=1,N
             Y(k)=Y(k)-1.d0/dble(it)*dt*((Aint(k)*cu(k)-                       &
     &                       Aint(k-1)*cu(k-1))/ h(k)/A(k)/VF(k)          &
     &          -Y(k)*(ww(k)-ww(k-1))/h(k)/A(k)/VF(k))
           enddo
        else                ! non-conservative
           do k=1,N
             Y(k)=Y(k)-1.d0/dble(it)*dt*((Aint(k)*cu(k)-                       &
     &                       Aint(k-1)*cu(k-1))/h(k)/A(k)/VF(k))
           enddo
        end if

      end do ! end of the iteration loop

c     rate of change due to advection 
      do k =1,N
         dy(k) = (y(k)-dy(k))/dt
      enddo

c Still to do: integrate fluxes in time (now cu = cu of last step )
c            flux = 0. ! at start
c            do k=1,N
c             flux(k)=flux(k)+1.d0/dble(it)*dt*cu(k)
 
      return
      end subroutine advection_ori

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
! ... extensively modified by Karline Soetaert
!-----------------------------------------------------------------------






      subroutine advectvol(N, Y, dt, V, Vint, flow, Bcup, Bcdown,           &
     &                     Yup, Ydown,  method,mode,dY, cu, it)

c-----------------------------------------------------------------------------------------
c Similar as above, but for volumetric transport:
c     use flow   = v*A rather than v
c         volume = h*A rather than h
c         hint  -> c(volume(1),volume) ... for now ... not so important
c-----------------------------------------------------------------------------------------

      IMPLICIT NONE
c  number of vertical layers, time step
      INTEGER                  :: N
      DOUBLE PRECISION         :: dt

c  layer thickness (m), distance from mid to mid
      DOUBLE PRECISION         :: V(N), Vint(0:N)

c  vertical advection speed 
      DOUBLE PRECISION         :: flow(0:N)

c  type of upper and lower Boundary Condition  (1 and 2 used)
      INTEGER                  :: Bcdown, Bcup
      integer, parameter       :: Flux     =1
      integer, parameter       :: Value    =2
      integer, parameter       :: ZeroGrad =3         ! not used
      integer, parameter       :: zeroDivergence =4   ! not used

c  value of upper and lower bnd conc
      DOUBLE PRECISION         :: Ydown, Yup

c  type of advection scheme, slope delimeters 
      INTEGER                  :: method
      integer, parameter       :: UPSTREAM =5  
      integer, parameter       :: P2       =4   
      integer, parameter       :: P2_PDM   =3
      integer, parameter       :: Superbee =2
      integer, parameter       :: MUSCL    =1

c  advection mode 0: non-conservative (e.g. water flow), 1: conservative, eg.g. sinking
      INTEGER                  :: mode

c  concentration to be transported
      DOUBLE PRECISION         :: Y(N)

c  rate of change due to advection
      DOUBLE PRECISION         :: dY(N)

      DOUBLE PRECISION         :: one6th=1.0d0/6.0d0
      INTEGER, parameter       :: itmax=100

c  LOCAL VARIABLES:
      integer                  :: i,k,it
      DOUBLE PRECISION         :: x,r,Phi,limit=0.d0
      DOUBLE PRECISION         :: Yu,Yc,Yd
      DOUBLE PRECISION         :: c,cmax
      DOUBLE PRECISION         :: cu(0:N)
c-----------------------------------------------------------------------------------------

c  initialize maximum Courant number
      cmax = 0.d0

c  copy of current value of state variables 
      do k =1,N
          dy(k) = y(k)
      enddo

c  compute maximum Courant number; estimate nr of iterations 
      do k=0,N
        c=dabs(flow(k))*dt/Vint(k)
        if (c.gt.cmax) cmax=c
      enddo

      it=min(itmax,int(cmax)+1)

c  (time) splitting loop
      do i=1,it

c  initialize interface fluxes with zero
      cu   = 0.d0

c     spatial loop  
        do k=1,N-1

c        positive speed 
           if (flow(k) .gt. 0.d0) then

              c=flow(k)/dble(it)*dt/Vint(k-1)          ! courant number

              if (k .gt. 1) then
                 Yu=Y(k-1)                             ! upstream value
              else
                 Yu=Y(k)
              end if
              
              Yc=Y(k  )                                ! central value
              
              Yd=Y(k+1)                              ! downstream value
                 
c        negative speed
           else

             c=-flow(k)/dble(it)*dt/Vint(k)            ! courant number

             if (k .lt. N-1) then
               Yu=Y(k+2)                               ! upstream value
             else
               Yu=Y(N)
             end if
             
             Yc=Y(k+1)                                   ! central value
             
             Yd=Y(k)

           end if
           if (abs(Yd-Yc) .gt. 1e-10) then          ! slope ratio
              r=(Yc-Yu)/(Yd-Yc)
           else
              r=(Yc-Yu)*1.e10
           end if

c        limit the flux according to different suggestions, phi = flux-factor
           select case (method)
             case (UPSTREAM)
               limit=0.d0
             case ((P2),(P2_PDM))
c         - for quickest
               x    =  one6th*(1.d0-2.d0*c)
               Phi  =  (0.5d0+x)+(0.5d0-x)*r

               if (method.eq.P2) then
                  limit=Phi
               else
                  limit=max(0.d0,min(Phi,2./(1.-c),2.*r/(c+1.e-10)))
               end if
             case (Superbee)
               limit=max(0.d0, min(1.d0, 2.0*r), min(r,2.*1.d0) )
             case (MUSCL)
               limit=max(0.d0,min(2.*1.d0,2.0*r,0.5*(1.0+r)))
             case default
c               call rerror( 'unkown advection method')
           end select

c        compute the limited flux  
           cu(k)  = flow(k) *(Yc+0.5d0*limit*(1.-c)*(Yd-Yc))

        end do

c       downstream boundary conditions
        select case (Bcdown)
          case (flux)
            cu(N) = Ydown                  ! flux OUT of the domain is positive
          case (value)
            if (flow(N).lt. 0.d0) then
              cu(N) =  flow(N)*Ydown
            else
              cu(N) =  flow(N)*Y(N)
            end if
          case (ZeroGrad)
              cu(N) =  flow(N)*Y(N)
          case (zeroDivergence)
            cu(N) = cu(N-1)
          case default
            call rwarn('unkown downstream boundary condition type')
        end select


c       upstream boundary conditions
        select case (Bcup)
          case (flux)
            cu(0) =   Yup                   ! flux into the domain is positive
          case (value)
            if(flow(0) .gt. 0.d0) then      ! Karline: CHECK!
              cu(0) =  flow(0)*Yup
            else
              cu(0) =  flow(0)*Y(1)
            end if
          case (ZeroGrad)
              cu(0) =  flow(0)*Y(1)
          case (zeroDivergence)
           cu(0) = cu(1)
          case default
             call rwarn('unkown upstream boundary condition type')
        end select

c     the advection step

        if (mode.eq.0) then ! non-conservative - KARLINE CHECK...
           do k=1,N
             Y(k)=Y(k)-1.d0/dble(it)*dt*((cu(k)- cu(k-1))/ V(k)                &
     &          -Y(k)*(flow(k)-flow(k-1))/V(k))
           enddo
        else                ! conservative - this is actually used
           do k=1,N
             Y(k)=Y(k)-1.d0/dble(it)*dt*((cu(k)- cu(k-1))/V(k))
           enddo
        end if

      end do ! end of the iteration loop

c     rate of change due to advection 
      do k =1,N
         dy(k) = (y(k)-dy(k))/dt
      enddo

c Still to do: integrate fluxes in time (now cu = cu of last step )
c            flux = 0. ! at start
c            do k=1,N
c             flux(k)=flux(k)+1.d0/real(it)*dt*cu(k)
 
      return
      end subroutine advectvol
